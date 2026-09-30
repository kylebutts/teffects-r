## Main ----

#' Treatment-effects estimation using lasso
#'
#' Estimates average treatment effects and potential-outcome means using
#' augmented inverse-probability weighting after plugin-lasso selection of
#' high-dimensional controls. Outcome controls are selected separately within
#' each treatment arm and all selected models are refitted without a penalty.
#'
#' @inheritParams teffects_estimators
#' @param selection Penalty-selection method. The rigorous plugin method is
#'   currently implemented.
#' @param xfolds Number of folds used for cross-fitting. One, the default,
#'   estimates the nuisance functions on the full estimation sample.
#' @param resamples Number of independently generated cross-fitting partitions
#'   to average. Values above one require `xfolds > 1`.
#' @param seed Optional integer seed used to generate cross-fitting partitions.
#' @param outcome_ainclude,treatment_ainclude Optional character vectors naming
#'   expanded model-matrix columns that must be retained in the corresponding
#'   post-lasso model. The intercept, when present, is always retained.
#' @return A `teffects_telasso` object.
#' @export
telasso <- function(
  fml,
  data,
  treatment_fml = NULL,
  omodel = "linear",
  tmodel = "logit",
  stat = c("ate", "att", "atc", "pomeans"),
  selection = "plugin",
  xfolds = 1L,
  resamples = 1L,
  seed = NULL,
  outcome_ainclude = NULL,
  treatment_ainclude = NULL,
  weights = NULL,
  vce = "robust",
  pstolerance = 1e-5,
  ...
) {
  dreamerr::check_arg(fml, "ts formula mbt")
  dreamerr::check_arg(data, "data.frame mbt")
  dreamerr::check_arg(treatment_fml, "NULL os formula")
  stat <- tolower(stat)
  dreamerr::check_set_arg(stat, "strict match")
  dreamerr::check_arg(
    outcome_ainclude,
    treatment_ainclude,
    "NULL character vector no na"
  )
  dreamerr::check_set_arg(
    omodel,
    "strict match(linear)",
    .message = "Only `omodel = \"linear\"` is implemented for `telasso()` yet."
  )
  dreamerr::check_set_arg(
    tmodel,
    "strict match(logit)",
    .message = "Only `tmodel = \"logit\"` is implemented for `telasso()` yet."
  )
  dreamerr::check_set_arg(
    selection,
    "strict match(plugin)",
    .message = "Only `selection = \"plugin\"` is implemented for `telasso()` yet."
  )
  dreamerr::check_set_arg(
    vce,
    "strict match(robust)",
    .message = "Only `vce = \"robust\"` is implemented for `telasso()` yet."
  )
  dreamerr::check_value(xfolds, "integer scalar GE{0}")
  dreamerr::check_value(resamples, "integer scalar GE{0}")
  if (resamples > 1L && xfolds == 1L) {
    stop("`resamples > 1` requires `xfolds > 1`.", call. = FALSE)
  }
  dreamerr::check_arg(
    seed,
    "NULL numeric scalar GE{0} LT{(.Machine$integer.max + 1)}",
    .message = "`seed` must be `NULL` or one nonnegative integer."
  )
  if (!is.null(seed)) {
    dreamerr::check_arg(
      seed,
      "class(numeric, integer)",
      .message = "`seed` must be `NULL` or one nonnegative integer."
    )
    dreamerr::check_arg(
      seed,
      "strict integer scalar",
      .message = "`seed` must be `NULL` or one nonnegative integer."
    )
    seed <- as.integer(seed)
    if (seed > .Machine$integer.max - resamples + 1L) {
      stop("`seed + resamples - 1` must fit in an integer.", call. = FALSE)
    }
  }
  dreamerr::check_value(
    pstolerance,
    "numeric scalar GE{0} LT{0.5}",
    .arg_name = "pstolerance",
  )

  .teffects_telasso(
    fml = fml,
    data = data,
    treatment_fml = treatment_fml,
    stat = stat,
    xfolds = xfolds,
    resamples = resamples,
    seed = seed,
    outcome_ainclude = outcome_ainclude,
    treatment_ainclude = treatment_ainclude,
    weights = weights,
    pstolerance = pstolerance,
    call = match.call()
  )
}


## Estimator ----

.teffects_telasso <- function(
  fml,
  data,
  treatment_fml,
  stat,
  xfolds,
  resamples,
  seed,
  outcome_ainclude,
  treatment_ainclude,
  weights,
  pstolerance,
  call
) {
  prepared <- .teffects_prepare(
    fml = fml,
    data = data,
    treatment_fml = treatment_fml,
    outcome_design = "full",
    treatment_design = "full",
    weights = weights
  )
  y <- prepared$outcome
  outcome_X <- prepared$outcome_design
  treatment <- prepared$treatment
  treatment_X <- prepared$treatment_design
  keep <- prepared$complete

  treated_levels <- prepared$treated_levels
  treatment_levels <- c(prepared$control_level, treated_levels)
  if (length(treatment_levels) != 2L) {
    stop("`telasso()` requires a binary treatment.", call. = FALSE)
  }
  control <- prepared$control_level
  treated <- treated_levels[[1L]]
  z <- as.numeric(treatment == treated)
  group_sizes <- vapply(
    treatment_levels,
    function(level) sum(treatment == level),
    numeric(1)
  )
  if (xfolds > min(group_sizes)) {
    stop(
      "`xfolds` cannot exceed the size of either treatment arm.",
      call. = FALSE
    )
  }

  outcome_ainclude <- .telasso_validate_ainclude(
    outcome_ainclude,
    colnames(outcome_X),
    "outcome_ainclude"
  )
  treatment_ainclude <- .telasso_validate_ainclude(
    treatment_ainclude,
    colnames(treatment_X),
    "treatment_ainclude"
  )

  estimates <- vector("list", resamples)
  influences <- vector("list", resamples)
  nuisance <- vector("list", resamples)
  for (r in seq_len(resamples)) {
    folds <- .telasso_folds(
      z,
      xfolds,
      if (is.null(seed)) NA_integer_ else seed + (r - 1L)
    )
    fitted <- .telasso_crossfit(
      y = y,
      z = z,
      outcome_X = outcome_X,
      treatment_X = treatment_X,
      folds = folds,
      outcome_ainclude = outcome_ainclude,
      treatment_ainclude = treatment_ainclude
    )
    .teffects_overlap_check(fitted$propensity, pstolerance)
    score <- .telasso_score(
      y = y,
      z = z,
      mu_0 = fitted$outcome_mean[, 1L],
      mu_1 = fitted$outcome_mean[, 2L],
      propensity = fitted$propensity,
      stat = stat,
      control = control,
      treated = treated
    )
    estimates[[r]] <- score$estimate
    influences[[r]] <- score$influence
    nuisance[[r]] <- c(
      fitted,
      list(
        folds = folds,
        pseudo_outcome = score$pseudo_outcome
      )
    )
  }

  estimate <- .telasso_average(estimates)
  influence <- .telasso_average(influences)
  names(estimate) <- colnames(influence) <- names(estimates[[1L]])
  propensity <- .telasso_average(lapply(nuisance, `[[`, "propensity"))
  outcome_mean <- .telasso_average(lapply(nuisance, `[[`, "outcome_mean"))
  pseudo_outcome <- .telasso_average(lapply(nuisance, `[[`, "pseudo_outcome"))

  .teffects_result(
    estimate,
    influence,
    "telasso",
    list(
      stat = stat,
      outcome_model_name = "linear",
      treatment_model_name = "logit",
      selection = "plugin",
      control = control,
      treated = treated,
      target = names(estimate),
      sample = which(keep),
      osample = .teffects_osample(prepared$n_original, keep),
      formula = fml,
      treatment_formula = treatment_fml,
      call = call,
      xfolds = xfolds,
      resamples = resamples,
      seed = seed,
      propensity_score = propensity,
      outcome_mean = outcome_mean,
      pseudo_outcome = pseudo_outcome,
      nuisance = nuisance
    )
  )
}


## Nuisance estimation ----

.telasso_average <- function(values) {
  if (length(values) == 1L) {
    return(values[[1L]])
  }
  Reduce(`+`, values) / length(values)
}

.telasso_validate_ainclude <- function(value, columns, argument) {
  if (is.null(value)) {
    return(character())
  }
  unknown <- setdiff(value, columns)
  if (length(unknown)) {
    stop(
      sprintf(
        "Unknown column%s in `%s`: %s.",
        if (length(unknown) == 1L) "" else "s",
        argument,
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  unique(value)
}

.telasso_folds <- function(z, xfolds, seed) {
  if (xfolds == 1L) {
    return(rep.int(1L, length(z)))
  }
  if (!is.na(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
    }
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    })
    set.seed(seed)
  }
  folds <- integer(length(z))
  for (arm in 0:1) {
    in_arm <- z == arm
    folds[in_arm] <- sample(rep(seq_len(xfolds), length.out = sum(in_arm)))
  }
  folds
}

.telasso_selected_columns <- function(x, y, family, intercept, always) {
  candidate <- setdiff(colnames(x), c("(Intercept)", always))
  if (length(candidate)) {
    varying <- vapply(
      candidate,
      function(column) {
        values <- as.numeric(x[, column])
        diff(range(values)) > sqrt(.Machine$double.eps)
      },
      logical(1)
    )
    candidate <- candidate[varying]
  }
  selected <- character()
  selector <- NULL
  if (length(candidate)) {
    dense_x <- as.matrix(x[, candidate, drop = FALSE])
    selector <- if (identical(family, "gaussian")) {
      hdm::rlasso(dense_x, y, post = FALSE, intercept = intercept)
    } else if (identical(family, "binomial")) {
      hdm::rlassologit(dense_x, y, post = FALSE, intercept = intercept)
    }
    selected <- candidate[selector$index]
  }
  columns <- unique(c(
    intersect("(Intercept)", colnames(x)),
    always,
    selected
  ))
  if (!length(columns)) {
    stop(
      "Plugin lasso selected no variables for a no-intercept model.",
      call. = FALSE
    )
  }
  list(columns = columns, selector = selector)
}

.telasso_post_design <- function(x, columns, always, model) {
  design <- .teffects_independent_columns(x[, columns, drop = FALSE])
  omitted_always <- setdiff(always, colnames(design))
  if (length(omitted_always)) {
    stop(
      sprintf(
        "Always-included column%s %s not estimable in the %s model.",
        if (length(omitted_always) == 1L) "" else "s",
        paste(omitted_always, collapse = ", "),
        model
      ),
      call. = FALSE
    )
  }
  design
}

.telasso_crossfit <- function(
  y,
  z,
  outcome_X,
  treatment_X,
  folds,
  outcome_ainclude,
  treatment_ainclude
) {
  fold_levels <- sort(unique(folds))
  n <- length(y)
  outcome_mean <- matrix(NA_real_, n, 2L)
  propensity <- rep(NA_real_, n)
  fits <- vector("list", length(fold_levels))
  outcome_intercept <- "(Intercept)" %in% colnames(outcome_X)
  treatment_intercept <- "(Intercept)" %in% colnames(treatment_X)

  for (k in seq_along(fold_levels)) {
    test <- folds == fold_levels[[k]]
    train <- if (length(fold_levels) == 1L) test else !test
    arm_fits <- vector("list", 2L)

    for (arm in 0:1) {
      arm_train <- train & z == arm
      arm_X <- outcome_X[arm_train, , drop = FALSE]
      arm_y <- y[arm_train]
      selected <- .telasso_selected_columns(
        arm_X,
        arm_y,
        family = "gaussian",
        intercept = outcome_intercept,
        always = outcome_ainclude
      )
      post_train_X <- .telasso_post_design(
        arm_X,
        selected$columns,
        outcome_ainclude,
        sprintf("outcome model for treatment arm %d", arm)
      )
      model <- tryCatch(
        .teffects_lm(post_train_X, arm_y),
        error = function(error) {
          stop(
            sprintf(
              "The post-lasso outcome model for treatment arm %d could not be fit: %s",
              arm,
              conditionMessage(error)
            ),
            call. = FALSE
          )
        }
      )
      outcome_mean[test, arm + 1L] <- as.numeric(
        outcome_X[test, colnames(post_train_X), drop = FALSE] %*%
          model$coefficients
      )
      arm_fits[[arm + 1L]] <- list(
        selected = colnames(post_train_X),
        selector = selected$selector,
        model = model
      )
    }

    treatment_train_X <- treatment_X[train, , drop = FALSE]
    treatment_train <- z[train]
    selected_treatment <- .telasso_selected_columns(
      treatment_train_X,
      treatment_train,
      family = "binomial",
      intercept = treatment_intercept,
      always = treatment_ainclude
    )
    post_treatment_train_X <- .telasso_post_design(
      treatment_train_X,
      selected_treatment$columns,
      treatment_ainclude,
      "treatment"
    )
    treatment_model <- tryCatch(
      .teffects_glm(
        post_treatment_train_X,
        treatment_train,
        "logit",
        compute_vcov = FALSE
      ),
      error = function(error) {
        stop(
          sprintf(
            "The post-lasso treatment model could not be fit: %s",
            conditionMessage(error)
          ),
          call. = FALSE
        )
      }
    )
    propensity[test] <- stats::plogis(as.numeric(
      treatment_X[test, colnames(post_treatment_train_X), drop = FALSE] %*%
        treatment_model$coefficients
    ))
    fits[[k]] <- list(
      fold = fold_levels[[k]],
      outcome = arm_fits,
      treatment = list(
        selected = colnames(post_treatment_train_X),
        selector = selected_treatment$selector,
        model = treatment_model
      )
    )
  }

  list(
    outcome_mean = outcome_mean,
    propensity = propensity,
    fits = fits
  )
}


## Orthogonal scores ----

.telasso_score <- function(
  y,
  z,
  mu_0,
  mu_1,
  propensity,
  stat,
  control,
  treated
) {
  if (stat %in% c("ate", "pomeans")) {
    # t = star: pseudo_k(O_i) = mu_k(X_i) + Z_{ki} / pi_k(X_i)
    #                             {Y_i - mu_k(X_i)}.
    pseudo0 <- mu_0 + (1 - z) * (y - mu_0) / (1 - propensity)
    pseudo1 <- mu_1 + z * (y - mu_1) / propensity
    theta <- c(mean(pseudo0), mean(pseudo1))
    influence <- cbind(pseudo0 - theta[[1L]], pseudo1 - theta[[2L]])
  } else if (identical(stat, "att")) {
    # t = 1: a_{1i} = z and q_t = mean(z). The control mean transports the
    # control residual with weight pi(X_i) / {1 - pi(X_i)}.
    p1 <- mean(z)
    pseudo0 <- (z *
      mu_0 +
      propensity * (1 - z) * (y - mu_0) / (1 - propensity)) /
      p1
    pseudo1 <- z * y / p1
    theta <- c(mean(pseudo0), mean(pseudo1))
    influence <- cbind(
      pseudo0 - z * theta[[1L]] / p1,
      pseudo1 - z * theta[[2L]] / p1
    )
  } else if (identical(stat, "atc")) {
    # t = 0: a_{0i} = 1 - z and q_t = mean(1 - z). The treated mean
    # transports the treated residual with weight {1 - pi(X_i)} / pi(X_i).
    p0 <- mean(1 - z)
    pseudo0 <- (1 - z) * y / p0
    pseudo1 <- ((1 - z) *
      mu_1 +
      (1 - propensity) * z * (y - mu_1) / propensity) /
      p0
    theta <- c(mean(pseudo0), mean(pseudo1))
    influence <- cbind(
      pseudo0 - (1 - z) * theta[[1L]] / p0,
      pseudo1 - (1 - z) * theta[[2L]] / p0
    )
  }

  level_names <- as.character(c(control, treated))
  if (identical(stat, "pomeans")) {
    estimate <- theta
    influence <- influence
    coefficient_names <- sprintf("POmean[%s]", level_names)
  } else if (stat %in% c("ate", "att", "atc")) {
    estimate <- c(theta[[2L]] - theta[[1L]], theta[[1L]])
    influence <- cbind(
      influence[, 2L] - influence[, 1L],
      influence[, 1L]
    )
    coefficient_names <- c(
      sprintf(
        "%s[%s vs %s]",
        toupper(stat),
        level_names[[2L]],
        level_names[[1L]]
      ),
      if (identical(stat, "att")) {
        sprintf("POmean[%s|%s]", level_names[[1L]], level_names[[2L]])
      } else if (stat %in% c("ate", "atc")) {
        sprintf("POmean[%s]", level_names[[1L]])
      }
    )
  }
  names(estimate) <- colnames(influence) <- coefficient_names
  list(
    estimate = estimate,
    influence = influence,
    pseudo_outcome = cbind(control = pseudo0, treated = pseudo1)
  )
}


## Methods ----

#' @export
print.teffects_telasso <- function(
  x,
  digits = max(3L, getOption("digits") - 2L),
  ...
) {
  cat("Treatment-effects lasso estimation\n")
  cat("Outcome model:", x$teffects$outcome_model_name, "\n")
  cat("Treatment model:", x$teffects$treatment_model_name, "\n")
  cat("Selection:", x$teffects$selection, "lasso\n")
  cat("Cross-fit folds:", x$teffects$xfolds, "\n")
  cat("Resamples:", x$teffects$resamples, "\n")
  cat("Statistic:", toupper(x$teffects$stat), "\n\n")
  stats::printCoefmat(
    x$coeftable,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )
  invisible(x)
}
