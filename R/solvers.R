.solver_extract_pars <- function(sol, solver = "solnp")
{
    switch(solver,
           "solnp" = sol$pars,
           "nloptr" = sol$solution,
           "optim" = sol$par)

}

.default_options <- function(solver, control = list(trace = 1))
{
    control <- switch(solver,
                      "solnp" = .solnp_defaults(control),
                      "nloptr" = .nloptr_defaults(control),
                      "optim" = .lbfgsb_defaults(control))
}

.nloptr_defaults <- function(control)
{
    if (!is.null(control$trace)) {
        control$print_level <- control$trace
        control$trace <- NULL
    }
    if (is.null(control$algorithm)) control$algorithm <- "NLOPT_LD_SLSQP"
    if (is.null(control$xtol_rel)) control$xtol_rel <- 1e-14
    if (is.null(control$maxeval)) control$maxeval <- 1000
    if (is.null(control$xtol_abs)) control$xtol_abs <- 1e-12
    if (is.null(control$check_derivatives)) control$check_derivatives <- FALSE
    return(control)
}

.lbfgsb_defaults <- function(control)
{
    if (is.null(control$trace)) control$trace <- 0
    return(control)
}

.solnp_defaults <- function(control)
{
    if (is.null(control$trace)) control$trace <- 0
    if (is.null(control$rho)) control$rho <- 1.0
    if (is.null(control$outer.iter)) control$outer.iter <- 400
    if (is.null(control$inner.iter)) control$inner.iter <- 600
    if (is.null(control$delta)) control$delta <- 1e-7
    if (is.null(control$tol)) control$tol <- 1e-8
    return(control)
}


.solver_extract_solution <- function(sol, solver = "solnp")
{
    switch(solver,
           "solnp" = sol$values[length(sol$values)],
           "nloptr" = sol$objective,
           "optim" = sol$value)
}


.dcc_dynamic_solver <- function(solver, pars, fun, lower, upper, control, arglist){
    if (solver == "solnp") {
        sol <- try(solnp(pars = pars, fun = fun, ineqfun = arglist$inequality_cons,
                         ineqLB = -1.0, ineqUB = 0.0, LB = lower, UB = upper,
                         control = control, arglist = arglist), silent = TRUE)
    } else if (solver == "nloptr") {
        sol <- try(nloptr(x0 = pars, eval_f = fun, eval_grad_f = arglist$grad,
                          lb = lower, ub = upper, eval_g_ineq = arglist$inequality_cons,
                   eval_jac_g_ineq = arglist$inequality_jac, opts = control,
                   arglist = arglist), silent = TRUE)
    }
    return(sol)
}

.dcc_constant_solver <- function(solver, pars, fun, lower, upper, control, arglist){
    if (solver == "solnp") {
        sol <- try(solnp(pars = pars, fun = fun, ineqfun = arglist$inequality_cons,
                         LB = lower, UB = upper,
                         control = control, arglist = arglist), silent = TRUE)
    } else if (solver == "nloptr") {
        sol <- try(nloptr(x0 = pars, eval_f = fun, eval_grad_f = arglist$grad,
                          lb = lower, ub = upper, opts = control,
                          arglist = arglist), silent = TRUE)
    } else if (solver == "optim") {
        sol <- try(optim(par = pars, fn = fun, gr = arglist$grad, lower = lower, upper = upper, method = "L-BFGS-B",
                         control = control, arglist = arglist), silent = TRUE)
    }
    return(sol)
}
