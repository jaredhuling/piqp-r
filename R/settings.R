# This file is part of PIQP-R. It is based on osqp-r
# (https://github.com/osqp/osqp-r) which is licensed
# under Apache License 2.0
#
# Copyright (c) 2023 EPFL
# Copyright (c) 2019 Paul Goulart, Bartolomeo Stellato
#
# This source code is licensed under the BSD 2-Clause License found in the
# LICENSE file in the root directory of this source tree.
  
#' Settings for PIQP
#'
#' For further details please consult the PIQP documentation:
#' \url{https://predict-epfl.github.io/piqp/}
#' @param rho_init initial value for the primal proximal penalty parameter rho (`1e-6`)
#' @param delta_init initial value for the augmented lagrangian penalty parameter delta (`1e-4`)
#' @param eps_abs absolute tolerance (`1e-8`)
#' @param eps_rel relative tolerance (`1e-9`)
#' @param check_duality_gap check terminal criterion on duality gap (`TRUE`)
#' @param eps_duality_gap_abs absolute tolerance on duality gap (`1e-8`)
#' @param eps_duality_gap_rel relative tolerance on duality gap (`1e-9`)
#' @param reg_lower_limit lower limit for regularization (`1e-10`)
#' @param reg_finetune_lower_limit fine tune lower limit regularization (`1e-13`)
#' @param reg_finetune_primal_update_threshold threshold of number of no primal updates to transition to fine tune mode (`7L`)
#' @param reg_finetune_dual_update_threshold threshold of number of no dual updates to transition to fine tune mode (`5L`)
#' @param max_iter maximum number of iterations (`250L`)
#' @param max_factor_retires maximum number of factorization retires before failure (`10L`)
#' @param preconditioner_scale_cost scale cost in Ruiz preconditioner (`FALSE`)
#' @param preconditioner_iter maximum of preconditioner iterations (`10L`)
#' @param tau maximum interior point step length (`0.99`)
#' @param iterative_refinement_always_enabled always run iterative refinement and not only on factorization failure (`FALSE`)
#' @param iterative_refinement_eps_abs iterative refinement absolute tolerance (`1e-12`)
#' @param iterative_refinement_eps_rel iterative refinement relative tolerance (`1e-12`)
#' @param iterative_refinement_max_iter maximum number of iterations for iterative refinement (`10L`)
#' @param iterative_refinement_min_improvement_rate minimum improvement rate for iterative refinement (`5.0`)
#' @param iterative_refinement_static_regularization_eps static regularization for KKT system for iterative refinement (`1e-7`)
#' @param iterative_refinement_static_regularization_rel static regularization w.r.t. the maximum abs diagonal term of KKT system (`.Machine$double.eps * .Machine$double.eps`)
#' @param verbose verbose printing (`FALSE`)
#' @param compute_timings measure timing information internally (`FALSE`)
#' @export
piqp_settings <- function(rho_init = 1e-6,
                          delta_init = 1e-4,
                          eps_abs = 1e-8,
                          eps_rel = 1e-9,
                          check_duality_gap = TRUE,
                          eps_duality_gap_abs = 1e-8,
                          eps_duality_gap_rel = 1e-9,
                          reg_lower_limit = 1e-10,
                          reg_finetune_lower_limit = 1e-13,
                          reg_finetune_primal_update_threshold = 7L,
                          reg_finetune_dual_update_threshold = 5L,
                          max_iter = 250L,
                          max_factor_retires = 10L,
                          preconditioner_scale_cost = FALSE,
                          preconditioner_iter = 10L,
                          tau = 0.99,
                          iterative_refinement_always_enabled = FALSE,
                          iterative_refinement_eps_abs = 1e-12,
                          iterative_refinement_eps_rel = 1e-12,
                          iterative_refinement_max_iter = 10L,
                          iterative_refinement_min_improvement_rate = 5.0,
                          iterative_refinement_static_regularization_eps = 1e-7,
                          iterative_refinement_static_regularization_rel = .Machine$double.eps * .Machine$double.eps,
                          verbose = FALSE,
                          compute_timings = FALSE
                        ) {
  inpars = as.list(match.call())[-1]
  pars = sapply(simplify = FALSE, USE.NAMES = TRUE, names(inpars), function(nm) {
    checkpar(inpars[[nm]], default_piqp_settings[[nm]])
  })
  pars
}


default_piqp_settings <- list(rho_init = 1e-6,
                              delta_init = 1e-4,
                              eps_abs = 1e-8,
                              eps_rel = 1e-9,
                              check_duality_gap = TRUE,
                              eps_duality_gap_abs = 1e-8,
                              eps_duality_gap_rel = 1e-9,
                              reg_lower_limit = 1e-10,
                              reg_finetune_lower_limit = 1e-13,
                              reg_finetune_primal_update_threshold = 7L,
                              reg_finetune_dual_update_threshold = 5L,
                              max_iter = 250L,
                              max_factor_retires = 10L,
                              preconditioner_scale_cost = FALSE,
                              preconditioner_iter = 10L,
                              tau = 0.99,
                              iterative_refinement_always_enabled = FALSE,
                              iterative_refinement_eps_abs = 1e-12,
                              iterative_refinement_eps_rel = 1e-12,
                              iterative_refinement_max_iter = 10L,
                              iterative_refinement_min_improvement_rate = 5.0,
                              iterative_refinement_static_regularization_eps = 1e-7,
                              iterative_refinement_static_regularization_rel = .Machine$double.eps * .Machine$double.eps,
                              verbose = FALSE,
                              compute_timings = FALSE)


checkpar <- function(l, r) {
  l <- switch(typeof(r),
             integer = as.integer(l),
             double = as.numeric(l),
             logical = as.logical(l))
  if (length(l) != 1 || is.na(l))
    return(r)
  l
}
