## Distance indexes ----

.teffects_make_quadratic_distance <- function(design, scaling, caliper = NULL) {
  design <- as.matrix(design)
  if (!ncol(design)) {
    coordinates <- matrix(numeric(), nrow(design), 0L)
  } else {
    decomposition <- eigen(scaling, symmetric = TRUE)
    root <- decomposition$vectors %*%
      diag(sqrt(pmax(decomposition$values, 0)), nrow = ncol(design)) %*%
      t(decomposition$vectors)
    coordinates <- design %*% root
  }
  structure(
    list(
      coordinates = coordinates,
      caliper = caliper
    ),
    class = "teffects_quadratic_distance"
  )
}

.teffects_matching_index <- function(
  distance,
  treatment,
  exact_group = NULL
) {
  n <- length(treatment)
  if (length(unique(treatment)) != 2L) {
    stop("Matching requires a binary treatment.", call. = FALSE)
  }
  exact_id <- if (is.null(exact_group)) {
    rep.int(1L, n)
  } else {
    as.integer(interaction(exact_group, drop = TRUE))
  }
  pool_id <- 2L * (exact_id - 1L) + treatment + 1L
  number_of_pools <- 2L * max(exact_id)
  pools <- split(
    seq_len(n),
    factor(pool_id, levels = seq_len(number_of_pools)),
    drop = FALSE
  )

  if (inherits(distance, "teffects_scalar_distance")) {
    values <- as.numeric(distance$values)
    pools <- lapply(pools, function(rows) {
      rows[order(values[rows], rows)]
    })
    pool_values <- lapply(pools, function(rows) values[rows])
    return(structure(
      list(
        type = "scalar",
        values = values,
        caliper = distance$caliper,
        treatment = treatment,
        exact_id = exact_id,
        pool_id = pool_id,
        pools = pools,
        pool_values = pool_values,
        n = n
      ),
      class = "teffects_matching_index"
    ))
  }

  if (!inherits(distance, "teffects_quadratic_distance")) {
    stop("Unsupported matching distance.", call. = FALSE)
  }
  if (is.null(distance$coordinates)) {
    distance <- .teffects_make_quadratic_distance(
      distance$design,
      distance$scaling,
      distance$caliper
    )
  }
  structure(
    list(
      type = "quadratic",
      coordinates = distance$coordinates,
      caliper = distance$caliper,
      treatment = treatment,
      exact_id = exact_id,
      pool_id = pool_id,
      pools = pools,
      n = n
    ),
    class = "teffects_matching_index"
  )
}

.teffects_query_pool_id <- function(index, rows, same) {
  treatment <- if (same) {
    index$treatment[rows]
  } else {
    1L - index$treatment[rows]
  }
  2L * (index$exact_id[rows] - 1L) + treatment + 1L
}


## Neighbor search ----

.teffects_select_neighbors <- function(
  candidate,
  candidate_distance,
  number,
  caliper,
  tolerance
) {
  keep <- is.finite(candidate_distance)
  if (!is.null(caliper)) {
    keep <- keep & candidate_distance <= caliper
  }
  candidate <- candidate[keep]
  candidate_distance <- candidate_distance[keep]
  if (!length(candidate)) {
    return(integer())
  }

  take <- min(number, length(candidate))
  cutoff <- sort.int(candidate_distance, partial = take)[[take]]
  expanded_cutoff <- cutoff + tolerance * max(1, abs(cutoff))
  if (!is.null(caliper)) {
    expanded_cutoff <- min(expanded_cutoff, caliper)
  }
  keep <- candidate_distance <= expanded_cutoff
  candidate <- candidate[keep]
  candidate_distance <- candidate_distance[keep]
  candidate[order(candidate_distance, candidate)]
}

.teffects_query_scalar_neighbors <- function(
  index,
  number,
  same,
  include_self,
  rows,
  tolerance
) {
  values <- index$values
  query_pool <- .teffects_query_pool_id(index, rows, same)
  exclude_self <- same && !include_self
  result <- vector("list", length(rows))

  for (position in seq_along(rows)) {
    i <- rows[[position]]
    pool <- index$pools[[query_pool[[position]]]]
    if (!length(pool)) {
      result[[position]] <- integer()
      next
    }
    scores <- index$pool_values[[query_pool[[position]]]]
    query <- values[[i]]
    if (is.null(index$caliper)) {
      available_start <- 1L
      available_end <- length(pool)
    } else {
      boundary_slop <- 8 *
        .Machine$double.eps *
        max(1, abs(query), index$caliper)
      available_start <- findInterval(
        query - index$caliper - boundary_slop,
        scores,
        left.open = TRUE
      ) +
        1L
      available_end <- findInterval(
        query + index$caliper + boundary_slop,
        scores
      )
      while (
        available_start <= available_end &&
          abs(scores[[available_start]] - query) > index$caliper
      ) {
        available_start <- available_start + 1L
      }
      while (
        available_end >= available_start &&
          abs(scores[[available_end]] - query) > index$caliper
      ) {
        available_end <- available_end - 1L
      }
    }
    if (available_start > available_end) {
      result[[position]] <- integer()
      next
    }

    available <- available_end - available_start + 1L
    if (exclude_self) {
      available <- available - 1L
    }
    if (!available) {
      result[[position]] <- integer()
      next
    }

    take <- min(number, available)
    insertion <- findInterval(query, scores)
    window <- seq.int(
      max(available_start, insertion - take - 1L),
      min(available_end, insertion + take + 1L)
    )
    if (exclude_self) {
      window <- window[pool[window] != i]
    }
    window_distance <- abs(scores[window] - query)
    cutoff <- sort.int(window_distance, partial = take)[[take]]
    expanded_cutoff <- cutoff + tolerance * max(1, abs(cutoff))
    if (!is.null(index$caliper)) {
      expanded_cutoff <- min(expanded_cutoff, index$caliper)
    }
    boundary_slop <- 8 *
      .Machine$double.eps *
      max(1, abs(query), expanded_cutoff)
    selected_start <- findInterval(
      query - expanded_cutoff - boundary_slop,
      scores,
      left.open = TRUE
    ) +
      1L
    selected_end <- findInterval(
      query + expanded_cutoff + boundary_slop,
      scores
    )
    selected_position <- seq.int(selected_start, selected_end)
    selected_distance <- abs(scores[selected_position] - query)
    within_cutoff <- selected_distance <= expanded_cutoff
    selected <- pool[selected_position[within_cutoff]]
    selected_distance <- selected_distance[within_cutoff]
    if (exclude_self) {
      keep <- selected != i
      selected <- selected[keep]
      selected_distance <- selected_distance[keep]
    }
    result[[position]] <- selected[order(selected_distance, selected)]
  }
  result
}

.teffects_quadratic_distances <- function(index, query, candidate) {
  if (!ncol(index$coordinates)) {
    return(numeric(length(candidate)))
  }
  difference <- sweep(
    index$coordinates[candidate, , drop = FALSE],
    2L,
    index$coordinates[query, ],
    "-"
  )
  sqrt(rowSums(difference^2))
}

.teffects_query_quadratic_exhaustive <- function(
  index,
  number,
  same,
  include_self,
  rows,
  tolerance
) {
  query_pool <- .teffects_query_pool_id(index, rows, same)
  result <- vector("list", length(rows))
  for (position in seq_along(rows)) {
    i <- rows[[position]]
    candidate <- index$pools[[query_pool[[position]]]]
    if (same && !include_self) {
      candidate <- candidate[candidate != i]
    }
    result[[position]] <- .teffects_select_neighbors(
      candidate = candidate,
      candidate_distance = .teffects_quadratic_distances(
        index,
        i,
        candidate
      ),
      number = number,
      caliper = index$caliper,
      tolerance = tolerance
    )
  }
  result
}

.teffects_query_quadratic_neighbors <- function(
  index,
  number,
  same,
  include_self,
  rows,
  tolerance
) {
  if (!ncol(index$coordinates)) {
    return(.teffects_query_quadratic_exhaustive(
      index,
      number,
      same,
      include_self,
      rows,
      tolerance
    ))
  }

  query_pool <- .teffects_query_pool_id(index, rows, same)
  row_groups <- split(seq_along(rows), query_pool)
  result <- vector("list", length(rows))

  for (positions in row_groups) {
    focal <- rows[positions]
    pool <- index$pools[[query_pool[[positions[[1L]]]]]]
    if (!length(pool)) {
      result[positions] <- rep(list(integer()), length(positions))
      next
    }

    extra_self <- as.integer(same && !include_self)
    search_k <- min(length(pool), number + 1L + extra_self)
    nearest <- FNN::get.knnx(
      data = index$coordinates[pool, , drop = FALSE],
      query = index$coordinates[focal, , drop = FALSE],
      k = search_k
    )
    nearest_index <- matrix(
      nearest$nn.index,
      nrow = length(focal),
      ncol = search_k
    )
    nearest_distance <- matrix(
      nearest$nn.dist,
      nrow = length(focal),
      ncol = search_k
    )

    for (j in seq_along(focal)) {
      i <- focal[[j]]
      candidate <- pool[nearest_index[j, ]]
      candidate_distance <- nearest_distance[j, ]
      if (same && !include_self) {
        keep <- candidate != i
        candidate <- candidate[keep]
        candidate_distance <- candidate_distance[keep]
      }
      keep <- is.finite(candidate_distance)
      if (!is.null(index$caliper)) {
        keep <- keep & candidate_distance <= index$caliper
      }
      candidate <- candidate[keep]
      candidate_distance <- candidate_distance[keep]
      if (!length(candidate)) {
        result[[positions[[j]]]] <- integer()
        next
      }

      candidate_order <- order(candidate_distance, candidate)
      candidate <- candidate[candidate_order]
      candidate_distance <- candidate_distance[candidate_order]
      take <- min(number, length(candidate))
      boundary_tie <-
        length(candidate_distance) > number &&
        candidate_distance[[number + 1L]] <=
          candidate_distance[[number]] +
            tolerance * max(1, abs(candidate_distance[[number]]))
      if (take < number || !boundary_tie) {
        result[[positions[[j]]]] <- candidate[seq_len(take)]
        next
      }

      all_candidate <- pool
      if (same && !include_self) {
        all_candidate <- all_candidate[all_candidate != i]
      }
      result[[positions[[j]]]] <- .teffects_select_neighbors(
        candidate = all_candidate,
        candidate_distance = .teffects_quadratic_distances(
          index,
          i,
          all_candidate
        ),
        number = number,
        caliper = index$caliper,
        tolerance = tolerance
      )
    }
  }
  result
}

.teffects_query_neighbors <- function(
  index,
  number,
  same = TRUE,
  include_self = TRUE,
  rows = seq_len(index$n),
  tolerance = sqrt(.Machine$double.eps),
  use_kd = TRUE
) {
  if (number < 1L) {
    stop("The number of neighbors must be positive.", call. = FALSE)
  }
  if (identical(index$type, "scalar")) {
    return(.teffects_query_scalar_neighbors(
      index,
      number,
      same,
      include_self,
      rows,
      tolerance
    ))
  }
  if (!use_kd) {
    return(.teffects_query_quadratic_exhaustive(
      index,
      number,
      same,
      include_self,
      rows,
      tolerance
    ))
  }
  .teffects_query_quadratic_neighbors(
    index,
    number,
    same,
    include_self,
    rows,
    tolerance
  )
}

# Retained as a compact test and diagnostic entry point.
.teffects_generic_neighbor_sets <- function(
  distance,
  treatment,
  number,
  same,
  include_self,
  exact_group,
  rows,
  tolerance
) {
  scalar <- inherits(distance, "teffects_scalar_distance")
  quadratic <- inherits(distance, "teffects_quadratic_distance")
  if (quadratic && is.null(distance$coordinates)) {
    distance <- .teffects_make_quadratic_distance(
      distance$design,
      distance$scaling,
      distance$caliper
    )
  }
  result <- vector("list", length(rows))
  for (position in seq_along(rows)) {
    i <- rows[[position]]
    eligible <- if (same) {
      treatment == treatment[[i]]
    } else {
      treatment != treatment[[i]]
    }
    if (!include_self) {
      eligible[[i]] <- FALSE
    }
    if (!is.null(exact_group)) {
      eligible <- eligible & exact_group == exact_group[[i]]
    }
    candidate <- which(eligible)
    candidate_distance <- if (scalar) {
      abs(distance$values[candidate] - distance$values[[i]])
    } else if (quadratic) {
      difference <- sweep(
        distance$coordinates[candidate, , drop = FALSE],
        2L,
        distance$coordinates[i, ],
        "-"
      )
      sqrt(rowSums(difference^2))
    } else {
      distance[i, candidate]
    }
    result[[position]] <- .teffects_select_neighbors(
      candidate,
      candidate_distance,
      number,
      if (scalar || quadratic) distance$caliper else NULL,
      tolerance
    )
  }
  result
}

.teffects_neighbor_sets <- function(
  distance,
  treatment,
  number,
  same = TRUE,
  include_self = TRUE,
  exact_group = NULL,
  rows = seq_along(treatment),
  tolerance = sqrt(.Machine$double.eps),
  use_kd = TRUE
) {
  supported_distance <- inherits(
    distance,
    c("teffects_scalar_distance", "teffects_quadratic_distance")
  )
  if (!supported_distance || length(unique(treatment)) != 2L) {
    return(.teffects_generic_neighbor_sets(
      distance,
      treatment,
      number,
      same,
      include_self,
      exact_group,
      rows,
      tolerance
    ))
  }
  index <- .teffects_matching_index(distance, treatment, exact_group)
  .teffects_query_neighbors(
    index,
    number,
    same,
    include_self,
    rows,
    tolerance,
    use_kd
  )
}


## Match plans ----

.teffects_sets_operator <- function(sets, focal, n) {
  degree <- lengths(sets)
  edge_focal <- rep.int(focal, degree)
  neighbor <- unlist(sets, use.names = FALSE)
  edge_weight <- rep.int(1 / degree, degree)
  operator <- Matrix::sparseMatrix(
    i = edge_focal,
    j = neighbor,
    x = edge_weight,
    dims = c(n, n),
    repr = "C"
  )
  list(
    operator = operator,
    focal = focal,
    degree = degree,
    edge_count = length(neighbor)
  )
}

.teffects_generated_matches <- function(
  sets,
  focal,
  treatment,
  row_map,
  prefix
) {
  directions <- list(
    att = seq_along(focal)[treatment[focal] == 1L],
    atc = seq_along(focal)[treatment[focal] == 0L]
  )
  directions <- directions[lengths(directions) > 0L]
  lapply(directions, function(positions) {
    direction_sets <- sets[positions]
    width <- max(lengths(direction_sets))
    result <- matrix(
      NA_character_,
      nrow = length(positions),
      ncol = width,
      dimnames = list(as.character(row_map[focal[positions]]), NULL)
    )
    for (i in seq_along(direction_sets)) {
      neighbors <- direction_sets[[i]]
      result[i, seq_along(neighbors)] <- as.character(row_map[neighbors])
    }
    colnames(result) <- paste0(prefix, seq_len(width))
    result
  })
}

.teffects_matching_target <- function(treatment, stat) {
  switch(
    stat,
    att = treatment == 1L,
    atc = treatment == 0L,
    ate = rep.int(TRUE, length(treatment))
  )
}

.teffects_build_match_plan <- function(
  index,
  treatment,
  stat,
  nneighbor,
  tolerance = sqrt(.Machine$double.eps),
  generate = NULL,
  row_map = seq_along(treatment)
) {
  target <- .teffects_matching_target(treatment, stat)
  focal <- which(target)
  sets <- .teffects_query_neighbors(
    index,
    nneighbor,
    same = FALSE,
    include_self = FALSE,
    rows = focal,
    tolerance = tolerance
  )
  if (any(lengths(sets) < nneighbor)) {
    stop("At least one observation could not be matched.", call. = FALSE)
  }

  result <- .teffects_sets_operator(sets, focal, length(treatment))
  names(result$degree) <- as.character(row_map[focal])
  reuse_k <- as.numeric(Matrix::colSums(result$operator))
  reuse_k_prime <- as.numeric(Matrix::colSums(
    result$operator * result$operator
  ))
  target_weight <- numeric(length(treatment))
  target_weight[focal] <- 1
  result$target <- target
  result$reuse <- list(k = reuse_k, k_prime = reuse_k_prime)
  result$weights <- target_weight + reuse_k
  result$generated_matches <- if (is.null(generate)) {
    NULL
  } else {
    .teffects_generated_matches(
      sets,
      focal,
      treatment,
      row_map,
      generate
    )
  }
  result
}

.teffects_build_neighborhood_plan <- function(
  index,
  number,
  same,
  include_self,
  rows = seq_len(index$n),
  tolerance = sqrt(.Machine$double.eps)
) {
  sets <- .teffects_query_neighbors(
    index,
    number,
    same,
    include_self,
    rows,
    tolerance
  )
  .teffects_sets_operator(sets, rows, index$n)
}
