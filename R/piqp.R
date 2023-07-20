# This file is part of PIQP-R. It is based on osqp-r
# (https://github.com/osqp/osqp-r) which is licensed
# under Apache License 2.0
#
# Copyright (c) 2023 EPFL
# Copyright (c) 2019 Paul Goulart, Bartolomeo Stellato
#
# This source code is licensed under the BSD 2-Clause License found in the
# LICENSE file in the root directory of this source tree.

#' PIQP Solver object
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom methods as is
#' @importFrom R6 R6Class
#' @param P dense or sparse matrix of class dgCMatrix or coercible into such, must be positive semidefinite
#' @param c numeric vector
#' @param A dense or sparse matrix of class dgCMatrix or coercible into such
#' @param b numeric vector
#' @param G dense or sparse matrix of class dgCMatrix or coercible into such
#' @param h numeric vector
#' @param x_lb numeric vector, with possibly infinite elements
#' @param x_ub numeric vector, with possibly infinite elements
#' @param settings list with optimization parameters, conveniently set with the function \code{\link{piqp_settings}}.
#' @param backend which backend to use, if auto and P, A or G are sparse then sparse backend is used (`"auto"`, `"sparse"` or `"dense"`) (`"auto"`)
#' @return An R6-object of class "piqp_model" with methods defined which can be further
#'   used to solve the problem with updated settings / parameters.
#' @seealso \code{\link{solve_piqp}}
#' @section Usage:
#' \preformatted{model = piqp(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL, x_lb = NULL, x_ub = NULL, settings = piqp_settings(), backend = c("auto", "sparse", "dense"))
#'
#' model$solve()
#' model$update(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL, x_lb = NULL, x_ub = NULL)
#' model$get_settings()
#' model$get_dims()
#' model$update_settings(new_settings = piqp_settings())
#'
#' print(model)
#' }
#' @details
#' Allows one to solve a parametric
#' problem with for example warm starts between updates of the parameter, c.f. the examples.
#' The object returned by \code{piqp} contains several methods which can be used to either update/get details of the
#' problem, modify the optimization settings or attempt to solve the problem.
#' @examples
#' ## example, adapted from PIQP documentation
#' library(piqp)
#' library(Matrix)
#'
#' P <- Matrix(c(6., 0.,
#'               0., 4.), 2, 2, sparse = TRUE)
#' c <- c(-1., -4.)
#' A <- Matrix(c(1., -2.), 1, 2, sparse = TRUE)
#' b <- c(1.)
#' G <- Matrix(c(1., 2., -1., 0.), 2, 2, sparse = TRUE)
#' h <- c(0.2, -1.)
#' x_lb <- c(-1., -Inf)
#' x_ub <- c(1., Inf)
#'
#' settings <- piqp_settings(verbose = TRUE)
#'
#' model <- piqp(P, c, A, b, G, h, x_lb, x_ub, settings)
#'
#' # Solve
#' res <- model$solve()
#' res$x
#'
#' # Define new data
#' A_new <- Matrix(c(1., -3.), 1, 2, sparse = TRUE)
#' h_new <- c(2., 1.)
#'
#' # Update model and solve again
#' model$update(A = A_new, h = h_new)
#' res <- model$solve()
#' res$x
#'
#' @export
piqp <- function(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL, x_lb = NULL, x_ub = NULL, settings = piqp_settings(), backend = c("auto", "sparse", "dense")) {

  # match possible options
  backend <- match.arg(backend)

  if (backend == "auto") {
    if ((!is.null(P) && is(P, 'sparseMatrix')) || !is.null(A) && is(A, 'sparseMatrix') || !is.null(G) && is(G, 'sparseMatrix')) {
      backend <- "sparse"
    } else {
      backend <- "dense"
    }
  }

  if (is.null(P) && is.null(c)) {
    stop("At least one of P and q must be supplied")
  }

  if (is.null(P)) {
    n <- length(q)
  } else {
    n <- dim(P)[1]
  }

  if (backend == "dense") {

    if (is.null(P)) {
      P <- matrix(0, n, n)
    } else {
      P <- as.matrix(P)
    }

    if (is.null(A)) {
      p <- 0
      A <- matrix(0, p, n)
      b <- numeric()
    } else {
      if (is.null(b))
        stop("b must be supplied if A is supplied")

      A <- as.matrix(A)
      p <- nrow(A)
      b <- as.numeric(b)
    }

    if (is.null(G)) {
      m <- 0
      G <- matrix(0, m, n)
      h <- numeric()
    } else {
      if (is.null(h))
        stop("h must be supplied if G is supplied")

      G <- as.matrix(G)
      m <- nrow(G)
      h <- as.numeric(h)
    }

  } else {

    if (is.null(P)) {
      P <- sparseMatrix(integer(), integer(), x = numeric(), dims = c(n, n))
    } else {
      P <- as(P, "dgCMatrix")
    }

    if (is.null(A)) {
      p <- 0
      A <- sparseMatrix(integer(), integer(), x = numeric(), dims = c(p, n))
      b <- numeric()
    } else {
      if (is.null(b))
        stop("b must be supplied if A is supplied")

      A <- as(A, "dgCMatrix")
      p <- nrow(A)
      b <- as.numeric(b)
    }

    if (is.null(G)) {
      m <- 0
      G <- sparseMatrix(integer(), integer(), x = numeric(), dims = c(m, n))
      h <- numeric()
    } else {
      if (is.null(h))
        stop("h must be supplied if G is supplied")

      G <- as(G, "dgCMatrix")
      m <- nrow(G)
      h <- as.numeric(h)
    }

  }

  if (is.null(c)) {
    c <- numeric(n)
  } else {
    c <- as.numeric(c)
  }

  if (!is.null(x_lb)) {
    x_lb <- as.numeric(x_lb)
  }

  if (!is.null(x_ub)) {
    x_ub <- as.numeric(x_ub)
  }

  stopifnot(dim(P) == c(n, n),
            length(c) == n,
            dim(A) == c(p, n),
            length(b) == p,
            dim(G) == c(m, n),
            length(h) == m,
            length(x_lb) == n,
            length(x_ub) == n)

  R6Class("piqp_model",
          public =
            list(
              initialize = function(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL, x_lb = NULL, x_ub = NULL, settings = list(), backend = NULL) {
                private$.backend <- backend
                if (private$.backend == "dense") {
                  private$.work <- piqp_setup_dense(P, c, A, b, G, h, x_lb, x_ub, settings)
                } else {
                  private$.work <- piqp_setup_sparse(P, c, A, b, G, h, x_lb, x_ub, settings)
                }
              },
              solve = function() {
                if (private$.backend == "dense") {
                  piqp_solve_dense(private$.work)
                } else {
                  piqp_solve_sparse(private$.work)
                }
              },
              update = function(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL, h = NULL, x_lb = NULL, x_ub = NULL) {
                if (private$.backend == "dense") {
                  dims <- piqp_get_dims_dense(private$.work)
                  stopifnot(length(c) %in% c(0, dims[[1]]),
                            length(b) %in% c(0, dims[[2]]),
                            length(h) %in% c(0, dims[[3]]),
                            length(x_lb) %in% c(0, dims[[1]]),
                            length(x_ub) %in% c(0, dims[[1]])
                            )
                  if (!is.null(P)) {
                    P <- as.matrix(P)
                  }
                  if (!is.null(A)) {
                    A <- as.matrix(A)
                  }
                  if (!is.null(G)) {
                    G <- as.matrix(G)
                  }
                  piqp_update_dense(private$.work, P, c, A, b, G, h, x_lb, x_ub)
                } else {
                  dims <- piqp_get_dims_sparse(private$.work)
                  stopifnot(length(c) %in% c(0, dims[[1]]),
                            length(b) %in% c(0, dims[[2]]),
                            length(h) %in% c(0, dims[[3]]),
                            length(x_lb) %in% c(0, dims[[1]]),
                            length(x_ub) %in% c(0, dims[[1]])
                            )
                  if (!is.null(P)) {
                    P <- as(P, "dgCMatrix")
                  }
                  if (!is.null(A)) {
                    A <- as(A, "dgCMatrix")
                  }
                  if (!is.null(G)) {
                    G <- as(G, "dgCMatrix")
                  }
                  piqp_update_sparse(private$.work, P, c, A, b, G, h, x_lb, x_ub)
                }
              },
              get_settings = function() {
                if (private$.backend == "dense") {
                  piqp_get_settings_dense(private$.work)
                } else {
                  piqp_get_settings_sparse(private$.work)
                }
              },
              get_dims = function() {
                if (private$.backend == "dense") {
                  piqp_get_dims_dense(private$.work)
                } else {
                  piqp_get_dims_sparse(private$.work)
                }
              },
              update_settings = function(new_settings = piqp_settings()) {
                stopifnot(is.list(new_settings))
                if (private$.backend == "dense") {
                  piqp_update_settings_dense(private$.work, new_settings)
                } else {
                  piqp_update_settings_sparse(private$.work, new_settings)
                }
              }
            ),
          private = list(.work = NULL, .backend = NULL)
  )$new(P, c, A, b, G, h, x_lb, x_ub, settings, backend)
}

#' @export
format.piqp_model <- function(x, ...) {
  dims <- x$GetDims()
  sprintf("PIQP-modelobject\n\nNumber of variables: %i\nNumber of constraints: %i", dims[[1]], dims[[2]])
}

#' @export
print.piqp_model <- function(x, ...)
  cat(format(x))


private <- NULL # to suppress cran note