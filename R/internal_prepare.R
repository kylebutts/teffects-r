#' Mark a treatment variable in a treatment-effects formula
#'
#' `treat()` identifies the treatment assignment and its reference level in a
#' `teffects` formula. The estimators parse the call directly, so `treat()` does
#' not accept any arguments beyond the treatment variable and its reference
#' level.
#'
#' @param treatment A treatment variable.
#' @param ref The reference (control) treatment value. When `NULL`, estimators
#'   use `0` or `FALSE` when present and otherwise use the first observed value.
#'
#' @return The result of [fixest::i()] with the requested reference level. This
#'   is only used when `treat()` is evaluated directly; within `teffects`
#'   formulas the call is parsed as a treatment marker before model-matrix
#'   construction.
#' @export
treat <- function(treatment, ref = NULL) {
  fixest::i(treatment, ref = ref)
}


# `outcome_design` and `treatment_design` say how much of each model's design
# to build: "reduced" drops linearly dependent columns, "full" keeps the model
# matrix as specified, and "none" skips the design entirely.
.teffects_prepare <- function(
  fml,
  data,
  treatment_fml = NULL,
  outcome_design = c("reduced", "full", "none"),
  treatment_design = c("reduced", "full", "none"),
  additional_complete = NULL,
  additional_formulas = list(),
  weights = NULL
) {
  fml <- fixest::xpd(fml)
  weights <- .teffects_observation_weights(weights, data)
  outcome_design <- match.arg(outcome_design)
  treatment_design <- match.arg(treatment_design)

  ## Parse for `treat()` call
  formula_terms <- stats::terms(fml)
  term_labels <- attr(formula_terms, "term.labels")
  term_calls <- lapply(term_labels, str2lang)
  is_treat_call <- vapply(
    term_calls,
    function(term_call) {
      if (!is.call(term_call)) {
        return(FALSE)
      }
      call_head <- term_call[[1L]]
      identical(call_head, quote(treat)) ||
        (is.call(call_head) &&
          identical(call_head[[1L]], quote(`::`)) &&
          identical(call_head[[2L]], quote(teffects)) &&
          identical(call_head[[3L]], quote(treat)))
    },
    logical(1)
  )
  treat_position <- seq_along(is_treat_call)[is_treat_call]
  if (length(treat_position) != 1L) {
    stop("`fml` must contain exactly one `treat()` term.", call. = FALSE)
  }
  treat_call <- term_calls[[treat_position]]
  treat_call[[1L]] <- quote(treat)
  # `treat()` only accepts a treatment variable and an optional reference, so
  # reject anything else before `match.call()` can emit a less specific error.
  treat_arguments <- as.list(treat_call)[-1L]
  treat_names <- names(treat_arguments)
  if (is.null(treat_names)) {
    treat_names <- rep("", length(treat_arguments))
  }
  if (
    length(treat_arguments) > 2L ||
      any(!treat_names %in% c("", "treatment", "ref"))
  ) {
    stop("`treat()` must only be a single variable.", call. = FALSE)
  }
  matched_treat <- match.call(treat, treat_call, expand.dots = FALSE)
  treatment_fml_from_treat <- stats::as.formula(
    call("~", matched_treat$treatment),
    env = environment(fml)
  )
  treatment_frame <- stats::model.frame(
    treatment_fml_from_treat,
    data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )
  treatment <- treatment_frame[[1L]]
  if (length(treatment) != nrow(data)) {
    stop("The treatment must have one value per row of `data`.", call. = FALSE)
  }
  reference <- if (!is.null(matched_treat$ref)) {
    eval(matched_treat$ref, data, environment(fml))
  } else {
    NULL
  }
  if (!is.null(reference)) {
    dreamerr::check_value(
      reference,
      "scalar",
      .arg_name = "ref",
      .message = "`ref` must identify one nonmissing treatment level."
    )
  }
  outcome_formula <- stats::reformulate(
    term_labels[-treat_position],
    response = paste(deparse(fml[[2L]]), collapse = ""),
    intercept = attr(formula_terms, "intercept")
  )
  environment(outcome_formula) <- environment(fml)
  outcome_terms <- stats::terms(outcome_formula)
  shared_design <- is.null(treatment_fml)
  # Every term in `fml` defines the common estimation sample, including outcome
  # covariates an estimator never consumes.
  outcome_frame <- stats::model.frame(
    outcome_terms,
    data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )
  outcome <- stats::model.response(outcome_frame)
  if (!is.numeric(outcome) || length(outcome) != nrow(data)) {
    stop(
      "The outcome must be one numeric value per row of `data`.",
      call. = FALSE
    )
  }
  outcome <- as.numeric(outcome)

  if (shared_design) {
    treatment_terms <- stats::delete.response(outcome_terms)
    treatment_design_frame <- outcome_frame
  } else {
    treatment_terms <- stats::terms(treatment_fml)
    treatment_design_frame <- stats::model.frame(
      treatment_terms,
      data,
      na.action = stats::na.pass,
      drop.unused.levels = FALSE
    )
  }

  additional_terms <- lapply(additional_formulas, stats::terms)
  additional_frames <- lapply(
    additional_terms,
    stats::model.frame,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )

  outcome_candidate <- if (identical(outcome_design, "none")) {
    NULL
  } else {
    .teffects_design_matrix(
      stats::delete.response(outcome_terms),
      outcome_frame,
      outcome_design
    )
  }
  treatment_candidate <- if (identical(treatment_design, "none")) {
    NULL
  } else {
    .teffects_design_matrix(
      treatment_terms,
      treatment_design_frame,
      treatment_design
    )
  }

  complete_parts <- c(
    list(outcome, treatment, outcome_frame, weights),
    if (!shared_design) list(treatment_design_frame),
    if (!is.null(outcome_candidate)) list(outcome_candidate),
    if (!shared_design && !is.null(treatment_candidate)) {
      list(treatment_candidate)
    },
    additional_frames
  )
  complete <- do.call(.teffects_complete_cases, complete_parts)
  if (!is.null(additional_complete)) {
    if (length(additional_complete) != nrow(data)) {
      stop(
        "`additional_complete` must have one value per row of `data`.",
        call. = FALSE
      )
    }
    complete <- complete & !is.na(additional_complete) & additional_complete
  }
  if (!any(complete)) {
    stop(
      "No complete observations remain after evaluating the models.",
      call. = FALSE
    )
  }

  additional_frames <- lapply(
    additional_frames,
    function(frame) frame[complete, , drop = FALSE]
  )

  outcome_matrix <- if (is.null(outcome_candidate)) {
    NULL
  } else {
    outcome_candidate[complete, , drop = FALSE]
  }
  if (identical(outcome_design, "reduced")) {
    outcome_matrix <- .teffects_independent_columns(outcome_matrix)
  }
  treatment_matrix <- if (is.null(treatment_candidate)) {
    NULL
  } else {
    treatment_candidate[complete, , drop = FALSE]
  }
  if (identical(treatment_design, "reduced")) {
    treatment_matrix <- .teffects_independent_columns(treatment_matrix)
  }

  compact_treatment <- treatment[complete]
  if (is.factor(compact_treatment)) {
    compact_treatment <- as.character(compact_treatment)
    if (!is.null(reference)) {
      reference <- as.character(reference)
    }
  }
  observed_levels <- unique(compact_treatment)
  control_level <- reference
  if (is.null(control_level)) {
    default <- if (is.logical(compact_treatment)) {
      FALSE
    } else if (is.numeric(compact_treatment)) {
      0
    } else {
      observed_levels[[1L]]
    }
    control_level <- if (default %in% observed_levels) {
      default
    } else {
      observed_levels[[1L]]
    }
  }
  if (!(control_level %in% observed_levels)) {
    stop(
      "The treatment reference level must be observed in the analysis sample.",
      call. = FALSE
    )
  }
  observed_levels <- c(
    control_level,
    sort(observed_levels[observed_levels != control_level])
  )
  treated_levels <- observed_levels[observed_levels != control_level]
  list(
    outcome = outcome[complete],
    treatment = compact_treatment,
    outcome_design = outcome_matrix,
    treatment_design = treatment_matrix,
    complete = complete,
    n = sum(complete),
    n_original = nrow(data),
    weights = weights[complete],
    treated_levels = treated_levels,
    control_level = control_level,
    additional_frames = additional_frames,
    additional_terms = additional_terms
  )
}

.teffects_complete_cases <- function(...) {
  inputs <- list(...)
  if (!length(inputs)) {
    stop("Supply at least one input.", call. = FALSE)
  }

  n <- vapply(
    inputs,
    function(x) if (is.null(dim(x))) length(x) else nrow(x),
    integer(1)
  )
  if (length(unique(n)) != 1L) {
    stop("All inputs must have the same number of rows.", call. = FALSE)
  }

  parts <- unlist(
    lapply(inputs, function(x) {
      if (is.data.frame(x)) unname(as.list(x)) else list(x)
    }),
    recursive = FALSE
  )
  if (!length(parts)) {
    return(rep.int(TRUE, n[[1L]]))
  }

  complete <- lapply(parts, function(x) {
    if (is.null(dim(x))) {
      if (is.numeric(x) || is.complex(x)) {
        return(is.finite(x))
      }
      return(!is.na(x))
    }
    if (inherits(x, "Matrix") || is.numeric(x) || is.complex(x)) {
      return(Matrix::rowSums(is.finite(x)) == ncol(x))
    }
    stats::complete.cases(x)
  })
  Reduce(`&`, complete)
}

.teffects_observation_weights <- function(weights, data) {
  if (is.null(weights)) {
    return(rep(1, nrow(data)))
  }
  if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(data)) {
      stop("The observation-weight column is not in `data`.", call. = FALSE)
    }
    weights <- data[[weights]]
  }
  if (!is.numeric(weights) || length(weights) != nrow(data)) {
    stop(
      "`weights` must be a numeric vector with one value per row of `data`.",
      call. = FALSE
    )
  }
  if (any(!is.finite(weights)) || any(weights < 0) || !any(weights > 0)) {
    stop(
      "`weights` must be finite, nonnegative, and positive for at least one observation.",
      call. = FALSE
    )
  }
  as.numeric(weights)
}

# Build the full candidate matrix from an evaluated model frame. Reduced designs
# use R's standard contrasts; full designs retain one indicator per factor level.
.teffects_design_matrix <- function(terms, data, mode) {
  contrasts_arg <- if (identical(mode, "full")) {
    factor_columns <- vapply(data, is.factor, logical(1))
    lapply(data[factor_columns], stats::contrasts, contrasts = FALSE)
  }
  if (!length(contrasts_arg)) {
    contrasts_arg <- NULL
  }

  design <- stats::model.matrix(
    terms,
    data,
    contrasts.arg = contrasts_arg
  )
  design <- Matrix::Matrix(design, sparse = TRUE)
  rownames(design) <- NULL
  design
}
