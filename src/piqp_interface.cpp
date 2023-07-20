// This file is part of PIQP-R.
//
// Copyright (c) 2023 EPFL
//
// This source code is licensed under the BSD 2-Clause License found in the
// LICENSE file in the root directory of this source tree.

#include <RcppEigen.h>
#include "piqp_types.h"
#include "piqp/piqp.hpp"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

template<typename T>
piqp::optional<T> nullable2optional(Rcpp::Nullable<T> data)
{
    if (data.isNotNull())
    {
        return piqp::optional<T>(Rcpp::as<T>(data.get()));
    }
    return piqp::optional<T>();
}

void update_settings(piqp::Settings<double>& piqp_settings, const List &r_settings)
{
    CharacterVector names(r_settings.names());
    for (int i = 0; i < r_settings.size(); i++)
    {
        if (Rf_isNull(names[i]))
            continue;
        std::string n = as<std::string>(names[i]);

        if (n == "rho_init")
            piqp_settings.rho_init = as<double>(r_settings[i]);
        else if (n == "delta_init")
            piqp_settings.delta_init = as<double>(r_settings[i]);
        else if (n == "eps_abs")
            piqp_settings.eps_abs = as<double>(r_settings[i]);
        else if (n == "eps_rel")
            piqp_settings.eps_rel = as<double>(r_settings[i]);
        else if (n == "check_duality_gap")
            piqp_settings.check_duality_gap = as<bool>(r_settings[i]);
        else if (n == "eps_duality_gap_abs")
            piqp_settings.eps_duality_gap_abs = as<double>(r_settings[i]);
        else if (n == "eps_duality_gap_rel")
            piqp_settings.eps_duality_gap_rel = as<double>(r_settings[i]);
        else if (n == "reg_lower_limit")
            piqp_settings.reg_lower_limit = as<double>(r_settings[i]);
        else if (n == "reg_finetune_lower_limit")
            piqp_settings.reg_finetune_lower_limit = as<double>(r_settings[i]);
        else if (n == "reg_finetune_primal_update_threshold")
            piqp_settings.reg_finetune_primal_update_threshold = as<piqp::isize>(r_settings[i]);
        else if (n == "reg_finetune_dual_update_threshold")
            piqp_settings.reg_finetune_dual_update_threshold = as<piqp::isize>(r_settings[i]);
        else if (n == "max_iter")
            piqp_settings.max_iter = as<piqp::isize>(r_settings[i]);
        else if (n == "max_factor_retires")
            piqp_settings.max_factor_retires = as<piqp::isize>(r_settings[i]);
        else if (n == "preconditioner_scale_cost")
            piqp_settings.preconditioner_scale_cost = as<bool>(r_settings[i]);
        else if (n == "preconditioner_iter")
            piqp_settings.preconditioner_iter = as<piqp::isize>(r_settings[i]);
        else if (n == "tau")
            piqp_settings.tau = as<double>(r_settings[i]);
        else if (n == "iterative_refinement_always_enabled")
            piqp_settings.iterative_refinement_always_enabled = as<bool>(r_settings[i]);
        else if (n == "iterative_refinement_eps_abs")
            piqp_settings.iterative_refinement_eps_abs = as<double>(r_settings[i]);
        else if (n == "iterative_refinement_eps_rel")
            piqp_settings.iterative_refinement_eps_rel = as<double>(r_settings[i]);
        else if (n == "iterative_refinement_max_iter")
            piqp_settings.iterative_refinement_max_iter = as<piqp::isize>(r_settings[i]);
        else if (n == "iterative_refinement_min_improvement_rate")
            piqp_settings.iterative_refinement_min_improvement_rate = as<double>(r_settings[i]);
        else if (n == "iterative_refinement_static_regularization_eps")
            piqp_settings.iterative_refinement_static_regularization_eps = as<double>(r_settings[i]);
        else if (n == "iterative_refinement_static_regularization_rel")
            piqp_settings.iterative_refinement_static_regularization_rel = as<double>(r_settings[i]);
        else if (n == "verbose")
            piqp_settings.verbose = as<bool>(r_settings[i]);
        else if (n == "compute_timings")
            piqp_settings.compute_timings = as<bool>(r_settings[i]);
    }
}

// [[Rcpp::export]]
SEXP piqp_setup_dense(Eigen::Map<Mat> P, 
                      Eigen::Map<Vec> c,
                      Eigen::Map<Mat> A,
                      Eigen::Map<Vec> b,
                      Eigen::Map<Mat> G,
                      Eigen::Map<Vec> h,
                      Rcpp::Nullable<Eigen::Map<Vec>> x_lb,
                      Rcpp::Nullable<Eigen::Map<Vec>> x_ub,
                      const List& settings)
{
    piqp::DenseSolver<double>* solver = new piqp::DenseSolver<double>();

    if (settings.size()) {
        update_settings(solver->settings(), settings);
    }

    solver->setup(P, c, A, b, G, h, nullable2optional(x_lb), nullable2optional(x_ub));

    XPtr<piqp::DenseSolver<double>> solver_p(solver);
    return solver_p;
}

// [[Rcpp::export]]
SEXP piqp_setup_sparse(Eigen::Map<SparseMat> P, 
                       Eigen::Map<Vec> c,
                       Eigen::Map<SparseMat> A,
                       Eigen::Map<Vec> b,
                       Eigen::Map<SparseMat> G,
                       Eigen::Map<Vec> h,
                       Nullable<Eigen::Map<Vec>> x_lb,
                       Nullable<Eigen::Map<Vec>> x_ub,
                       const List& settings)
{
    piqp::SparseSolver<double>* solver = new piqp::SparseSolver<double>();

    if (settings.size()) {
        update_settings(solver->settings(), settings);
    }

    solver->setup(P, c, A, b, G, h, nullable2optional(x_lb), nullable2optional(x_ub));

    Rcpp::XPtr<piqp::SparseSolver<double>> solver_p(solver);
    return solver_p;
}

template<typename Solver>
List piqp_solve(Solver& solver)
{
    solver.solve();

    std::string status = piqp::status_to_string(solver.result().info.status);
    List info(24);
    CharacterVector info_names(24);
    info_names[0] = "status"; info[0] = status;
    info_names[1] = "iter"; info[1] = solver.result().info.iter;
    info_names[2] = "rho"; info[2] = solver.result().info.rho;
    info_names[3] = "delta"; info[3] = solver.result().info.delta;
    info_names[4] = "mu"; info[4] = solver.result().info.mu;
    info_names[5] = "sigma"; info[5] = solver.result().info.sigma;
    info_names[6] = "primal_step"; info[6] = solver.result().info.primal_step;
    info_names[7] = "dual_step"; info[7] = solver.result().info.dual_step;
    info_names[8] = "primal_inf"; info[8] = solver.result().info.primal_inf;
    info_names[9] = "primal_rel_inf"; info[9] = solver.result().info.primal_rel_inf;
    info_names[10] = "dual_inf"; info[10] = solver.result().info.dual_inf;
    info_names[11] = "dual_rel_inf"; info[11] = solver.result().info.dual_rel_inf;
    info_names[12] = "primal_obj"; info[12] = solver.result().info.primal_obj;
    info_names[13] = "dual_obj"; info[13] = solver.result().info.dual_obj;
    info_names[14] = "duality_gap"; info[14] = solver.result().info.duality_gap;
    info_names[15] = "duality_gap_rel"; info[15] = solver.result().info.duality_gap_rel;
    info_names[16] = "factor_retires"; info[16] = solver.result().info.factor_retires;
    info_names[17] = "reg_limit"; info[17] = solver.result().info.reg_limit;
    info_names[18] = "no_primal_update"; info[18] = solver.result().info.no_primal_update;
    info_names[19] = "no_dual_update"; info[19] = solver.result().info.no_dual_update;
    info_names[20] = "setup_time"; info[20] = solver.result().info.setup_time;
    info_names[21] = "update_time"; info[21] = solver.result().info.update_time;
    info_names[22] = "solve_time"; info[22] = solver.result().info.solve_time;
    info_names[23] = "run_time"; info[23] = solver.result().info.run_time;
    info.names() = info_names;

    List result = List::create(
        _("x") = solver.result().x,
        _("y") = solver.result().y,
        _("z") = solver.result().z,
        _("z_lb") = solver.result().z_lb,
        _("z_ub") = solver.result().z_ub,
        _("s") = solver.result().s,
        _("s_lb") = solver.result().s_lb,
        _("s_ub") = solver.result().s_ub,
        _("zeta") = solver.result().zeta,
        _("lambda") = solver.result().lambda,
        _("nu") = solver.result().nu,
        _("nu_lb") = solver.result().nu_lb,
        _("nu_ub") = solver.result().nu_ub,
        _("info") = info
    );

    return result;
}

// [[Rcpp::export]]
List piqp_solve_dense(SEXP solver_p)
{
    auto solver = as<XPtr<piqp::DenseSolver<double>>>(solver_p);
    return piqp_solve(*solver);
}

// [[Rcpp::export]]
List piqp_solve_sparse(SEXP solver_p)
{
    auto solver = as<XPtr<piqp::SparseSolver<double>>>(solver_p);
    return piqp_solve(*solver);
}

template<typename Solver>
IntegerVector piqp_get_dims(const Solver& solver)
{
    auto res = IntegerVector::create(_("n") = solver.result().x.rows(),
                                     _("p") = solver.result().y.rows(),
                                     _("m") = solver.result().z.rows());
    return res;
}

// [[Rcpp::export]]
IntegerVector piqp_get_dims_dense(SEXP solver_p)
{
    auto solver = as<XPtr<piqp::DenseSolver<double>>>(solver_p);
    return piqp_get_dims(*solver);
}

// [[Rcpp::export]]
IntegerVector piqp_get_dims_sparse(SEXP solver_p)
{
    auto solver = as<XPtr<piqp::SparseSolver<double>>>(solver_p);
    return piqp_get_dims(*solver);
}

// [[Rcpp::export]]
void piqp_update_dense(SEXP solver_p,
                       Rcpp::Nullable<Eigen::Map<Mat>> P, 
                       Rcpp::Nullable<Eigen::Map<Vec>> c,
                       Rcpp::Nullable<Eigen::Map<Mat>> A,
                       Rcpp::Nullable<Eigen::Map<Vec>> b,
                       Rcpp::Nullable<Eigen::Map<Mat>> G,
                       Rcpp::Nullable<Eigen::Map<Vec>> h,
                       Rcpp::Nullable<Eigen::Map<Vec>> x_lb,
                       Rcpp::Nullable<Eigen::Map<Vec>> x_ub)
{
    auto solver = as<XPtr<piqp::DenseSolver<double>>>(solver_p);
    solver->update(nullable2optional(P),
                   nullable2optional(c),
                   nullable2optional(A),
                   nullable2optional(b),
                   nullable2optional(G),
                   nullable2optional(h),
                   nullable2optional(x_lb),
                   nullable2optional(x_ub));
}

// [[Rcpp::export]]
void piqp_update_sparse(SEXP solver_p,
                        Rcpp::Nullable<Eigen::Map<SparseMat>> P, 
                        Rcpp::Nullable<Eigen::Map<Vec>> c,
                        Rcpp::Nullable<Eigen::Map<SparseMat>> A,
                        Rcpp::Nullable<Eigen::Map<Vec>> b,
                        Rcpp::Nullable<Eigen::Map<SparseMat>> G,
                        Rcpp::Nullable<Eigen::Map<Vec>> h,
                        Nullable<Eigen::Map<Vec>> x_lb,
                        Nullable<Eigen::Map<Vec>> x_ub)
{
    auto solver = as<XPtr<piqp::SparseSolver<double>>>(solver_p);
    solver->update(nullable2optional(P),
                   nullable2optional(c),
                   nullable2optional(A),
                   nullable2optional(b),
                   nullable2optional(G),
                   nullable2optional(h),
                   nullable2optional(x_lb),
                   nullable2optional(x_ub));
}

template<typename Solver>
List piqp_get_settings(Solver& solver)
{
    List settings(25);
    CharacterVector settings_names(25);
    settings_names[0] = "rho_init"; settings[0] = solver.settings().rho_init;
    settings_names[1] = "delta_init"; settings[1] = solver.settings().delta_init;
    settings_names[2] = "eps_abs"; settings[2] = solver.settings().eps_abs;
    settings_names[3] = "eps_rel"; settings[3] = solver.settings().eps_rel;
    settings_names[4] = "check_duality_gap"; settings[4] = solver.settings().check_duality_gap;
    settings_names[5] = "eps_duality_gap_abs"; settings[5] = solver.settings().eps_duality_gap_abs;
    settings_names[6] = "eps_duality_gap_rel"; settings[6] = solver.settings().eps_duality_gap_rel;
    settings_names[7] = "reg_lower_limit"; settings[7] = solver.settings().reg_lower_limit;
    settings_names[8] = "reg_finetune_lower_limit"; settings[8] = solver.settings().reg_finetune_lower_limit;
    settings_names[9] = "reg_finetune_primal_update_threshold"; settings[9] = solver.settings().reg_finetune_primal_update_threshold;
    settings_names[10] = "reg_finetune_dual_update_threshold"; settings[10] = solver.settings().reg_finetune_dual_update_threshold;
    settings_names[11] = "max_iter"; settings[11] = solver.settings().max_iter;
    settings_names[12] = "max_factor_retires"; settings[12] = solver.settings().max_factor_retires;
    settings_names[13] = "preconditioner_scale_cost"; settings[13] = solver.settings().preconditioner_scale_cost;
    settings_names[14] = "preconditioner_iter"; settings[14] = solver.settings().preconditioner_iter;
    settings_names[15] = "tau"; settings[15] = solver.settings().tau;
    settings_names[16] = "iterative_refinement_always_enabled"; settings[16] = solver.settings().iterative_refinement_always_enabled;
    settings_names[17] = "iterative_refinement_eps_abs"; settings[17] = solver.settings().iterative_refinement_eps_abs;
    settings_names[18] = "iterative_refinement_eps_rel"; settings[18] = solver.settings().iterative_refinement_eps_rel;
    settings_names[19] = "iterative_refinement_max_iter"; settings[19] = solver.settings().iterative_refinement_max_iter;
    settings_names[20] = "iterative_refinement_min_improvement_rate"; settings[20] = solver.settings().iterative_refinement_min_improvement_rate;
    settings_names[21] = "iterative_refinement_static_regularization_eps"; settings[21] = solver.settings().iterative_refinement_static_regularization_eps;
    settings_names[22] = "iterative_refinement_static_regularization_rel"; settings[22] = solver.settings().iterative_refinement_static_regularization_rel;
    settings_names[23] = "verbose"; settings[23] = solver.settings().verbose;
    settings_names[24] = "compute_timings"; settings[24] = solver.settings().compute_timings;
    settings.names() = settings_names;

    return settings;
}

// [[Rcpp::export]]
List piqp_get_settings_dense(SEXP solver_p)
{
    auto solver = as<XPtr<piqp::DenseSolver<double>>>(solver_p);
    return piqp_get_settings(*solver);
}

// [[Rcpp::export]]
List piqp_get_settings_sparse(SEXP solver_p)
{
    auto solver = as<XPtr<piqp::SparseSolver<double>>>(solver_p);
    return piqp_get_settings(*solver);
}

// [[Rcpp::export]]
void piqp_update_settings_dense(SEXP solver_p, const List& settings)
{
    auto solver = as<XPtr<piqp::DenseSolver<double>>>(solver_p);
    
    if (settings.size()) {
        update_settings(solver->settings(), settings);
    }
}

// [[Rcpp::export]]
void piqp_update_settings_sparse(SEXP solver_p, const List& settings)
{
    auto solver = as<XPtr<piqp::SparseSolver<double>>>(solver_p);
    
    if (settings.size()) {
        update_settings(solver->settings(), settings);
    }
}
