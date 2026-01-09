using Adapt
using Oceananigans.Grids: prettysummary

#####
##### Abstract solver type
#####

"""
    AbstractInterfaceSolver

Abstract supertype for iterative solvers that compute interface characteristic scales
(u★, θ★, q★) for Monin-Obukhov similarity theory.
"""
abstract type AbstractInterfaceSolver end

#####
##### Solver stop criteria
#####

"""
    ConvergenceStopCriteria{FT}

Stop criteria based on convergence tolerance and maximum iterations.

Fields
======
- `tolerance`: The convergence tolerance for the solver.
- `maxiter`: The maximum number of iterations allowed.
"""
struct ConvergenceStopCriteria{FT}
    tolerance :: FT
    maxiter :: Int
end

Base.summary(c::ConvergenceStopCriteria) =
    string("ConvergenceStopCriteria(tolerance=", prettysummary(c.tolerance),
           ", maxiter=", c.maxiter, ")")

"""
    FixedIterations{I}

Stop criteria that runs a fixed number of iterations regardless of convergence.

Fields
======
- `iterations`: The number of iterations to run.
"""
struct FixedIterations{I}
    iterations :: I
end

Base.summary(f::FixedIterations) = string("FixedIterations(", f.iterations, ")")

#####
##### Fixed-point iteration solver
#####

"""
    FixedPointSolver{S} <: AbstractInterfaceSolver

Fixed-point (Picard) iteration solver for interface state computation.
This is the traditional approach that iteratively updates the characteristic
scales until convergence.

Fields
======
- `stop_criteria`: Either `ConvergenceStopCriteria` or `FixedIterations`.
"""
struct FixedPointSolver{S} <: AbstractInterfaceSolver
    stop_criteria :: S
end

"""
    FixedPointSolver(FT = Float64; tolerance = 1e-8, maxiter = 100)

Construct a `FixedPointSolver` with convergence-based stop criteria.

Keyword Arguments
=================
- `tolerance`: The convergence tolerance. Default: 1e-8.
- `maxiter`: The maximum number of iterations. Default: 100.
"""
function FixedPointSolver(FT::DataType = Float64; tolerance = 1e-8, maxiter = 100)
    stop_criteria = ConvergenceStopCriteria(convert(FT, tolerance), maxiter)
    return FixedPointSolver(stop_criteria)
end

Base.summary(::FixedPointSolver) = "FixedPointSolver"

function Base.show(io::IO, solver::FixedPointSolver)
    print(io, summary(solver), "(", summary(solver.stop_criteria), ")")
end

Adapt.adapt_structure(to, solver::FixedPointSolver) =
    FixedPointSolver(Adapt.adapt(to, solver.stop_criteria))

#####
##### Broyden's method solver
#####

"""
    BroydenSolver{N, FT, S} <: AbstractInterfaceSolver

Good Broyden's method solver for interface state computation.
Uses rank-1 Jacobian updates starting from a scaled identity matrix.

This quasi-Newton method typically converges faster than fixed-point iteration
for the nonlinear system arising in Monin-Obukhov similarity theory.

Type Parameters
===============
- `N`: Number of variables (3 for `BulkTemperature`, 4 for `SkinTemperature`).
- `FT`: Floating-point type.
- `S`: Stop criteria type.

Fields
======
- `stop_criteria`: Either `ConvergenceStopCriteria` or `FixedIterations`.
- `initial_jacobian_scale`: Scale factor for initial Jacobian J₀ = scale * I.
"""
struct BroydenSolver{N, FT, S} <: AbstractInterfaceSolver
    stop_criteria :: S
    initial_jacobian_scale :: FT
end

"""
    BroydenSolver(::Val{N}, FT = Float64; tolerance = 1e-8, maxiter = 100, initial_jacobian_scale = 1.0)

Construct a `BroydenSolver` for an N-variable system.

Arguments
=========
- `Val{N}`: The number of variables. Use `Val(3)` for `BulkTemperature`
  (solving for u★, θ★, q★) or `Val(4)` for `SkinTemperature`
  (solving for u★, θ★, q★, Tₛ).

Keyword Arguments
=================
- `tolerance`: The convergence tolerance. Default: 1e-8.
- `maxiter`: The maximum number of iterations. Default: 100.
- `initial_jacobian_scale`: Scale for initial Jacobian approximation J₀ = scale * I. Default: 1.0.
"""
function BroydenSolver(::Val{N}, FT::DataType = Float64;
                       tolerance = 1e-8,
                       maxiter = 100,
                       initial_jacobian_scale = 1.0) where N
    stop_criteria = ConvergenceStopCriteria(convert(FT, tolerance), maxiter)
    return BroydenSolver{N, FT, typeof(stop_criteria)}(stop_criteria, convert(FT, initial_jacobian_scale))
end

Base.summary(::BroydenSolver{N, FT}) where {N, FT} = "BroydenSolver{$N, $FT}"

function Base.show(io::IO, solver::BroydenSolver{N}) where N
    print(io, summary(solver), '\n',
          "├── variables: ", N, '\n',
          "├── initial_jacobian_scale: ", prettysummary(solver.initial_jacobian_scale), '\n',
          "└── stop_criteria: ", summary(solver.stop_criteria))
end

Adapt.adapt_structure(to, solver::BroydenSolver{N, FT}) where {N, FT} =
    BroydenSolver{N, FT, typeof(solver.stop_criteria)}(
        Adapt.adapt(to, solver.stop_criteria),
        solver.initial_jacobian_scale)

#####
##### Helper to determine number of variables from temperature formulation
#####

"""
    num_solver_variables(temperature_formulation)

Return `Val(N)` where N is the number of variables for the Broyden solver
based on the temperature formulation.

- `BulkTemperature`: Returns `Val(3)` for (u★, θ★, q★)
- `SkinTemperature`: Returns `Val(4)` for (u★, θ★, q★, Tₛ)
"""
num_solver_variables(::BulkTemperature) = Val(3)
num_solver_variables(::SkinTemperature) = Val(4)
