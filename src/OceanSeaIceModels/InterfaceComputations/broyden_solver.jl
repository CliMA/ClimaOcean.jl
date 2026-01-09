#####
##### GPU-compatible linear algebra helpers for Broyden's method
#####
##### All operations use tuples (NTuple) instead of arrays to avoid allocations
##### and ensure GPU compatibility. Functions are @inline for performance.
#####

#####
##### Identity matrix initialization (stored as row-major tuple)
#####

"""
    identity_matrix_inv(::Val{3}, scale)

Return a 3x3 identity matrix scaled by `scale`, stored as an NTuple{9}.
The matrix is stored in row-major order.
"""
@inline function identity_matrix_inv(::Val{3}, scale::FT) where FT
    z = zero(FT)
    return (scale, z, z,
            z, scale, z,
            z, z, scale)
end

"""
    identity_matrix_inv(::Val{4}, scale)

Return a 4x4 identity matrix scaled by `scale`, stored as an NTuple{16}.
The matrix is stored in row-major order.
"""
@inline function identity_matrix_inv(::Val{4}, scale::FT) where FT
    z = zero(FT)
    return (scale, z, z, z,
            z, scale, z, z,
            z, z, scale, z,
            z, z, z, scale)
end

#####
##### Matrix-vector multiplication
#####

"""
    mat_vec_mul(J, v, ::Val{3})

Compute J * v where J is a 3x3 matrix stored as NTuple{9} (row-major)
and v is a 3-vector stored as NTuple{3}. Returns NTuple{3}.
"""
@inline function mat_vec_mul(J, v, ::Val{3})
    v1 = J[1]*v[1] + J[2]*v[2] + J[3]*v[3]
    v2 = J[4]*v[1] + J[5]*v[2] + J[6]*v[3]
    v3 = J[7]*v[1] + J[8]*v[2] + J[9]*v[3]
    return (v1, v2, v3)
end

"""
    mat_vec_mul(J, v, ::Val{4})

Compute J * v where J is a 4x4 matrix stored as NTuple{16} (row-major)
and v is a 4-vector stored as NTuple{4}. Returns NTuple{4}.
"""
@inline function mat_vec_mul(J, v, ::Val{4})
    v1 = J[1]*v[1]  + J[2]*v[2]  + J[3]*v[3]  + J[4]*v[4]
    v2 = J[5]*v[1]  + J[6]*v[2]  + J[7]*v[3]  + J[8]*v[4]
    v3 = J[9]*v[1]  + J[10]*v[2] + J[11]*v[3] + J[12]*v[4]
    v4 = J[13]*v[1] + J[14]*v[2] + J[15]*v[3] + J[16]*v[4]
    return (v1, v2, v3, v4)
end

"""
    negative_mat_vec_mul(J, v, ::Val{N})

Compute -J * v. Used to compute the Newton step Δx = -J⁻¹ F.
"""
@inline function negative_mat_vec_mul(J, v, ::Val{3})
    result = mat_vec_mul(J, v, Val(3))
    return (-result[1], -result[2], -result[3])
end

@inline function negative_mat_vec_mul(J, v, ::Val{4})
    result = mat_vec_mul(J, v, Val(4))
    return (-result[1], -result[2], -result[3], -result[4])
end

#####
##### Vector operations
#####

"""
    dot_product(a, b, ::Val{N})

Compute the dot product a ⋅ b for N-dimensional vectors stored as NTuple.
"""
@inline dot_product(a, b, ::Val{3}) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline dot_product(a, b, ::Val{4}) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3] + a[4]*b[4]

"""
    add_tuples(a, b, ::Val{N})

Compute a + b element-wise for N-dimensional vectors stored as NTuple.
"""
@inline add_tuples(a, b, ::Val{3}) = (a[1]+b[1], a[2]+b[2], a[3]+b[3])
@inline add_tuples(a, b, ::Val{4}) = (a[1]+b[1], a[2]+b[2], a[3]+b[3], a[4]+b[4])

"""
    subtract_tuples(a, b, ::Val{N})

Compute a - b element-wise for N-dimensional vectors stored as NTuple.
"""
@inline subtract_tuples(a, b, ::Val{3}) = (a[1]-b[1], a[2]-b[2], a[3]-b[3])
@inline subtract_tuples(a, b, ::Val{4}) = (a[1]-b[1], a[2]-b[2], a[3]-b[3], a[4]-b[4])

"""
    scale_tuple(α, v, ::Val{N})

Compute α * v for scalar α and N-dimensional vector v.
"""
@inline scale_tuple(α, v, ::Val{3}) = (α*v[1], α*v[2], α*v[3])
@inline scale_tuple(α, v, ::Val{4}) = (α*v[1], α*v[2], α*v[3], α*v[4])

"""
    tuple_norm(v, ::Val{N})

Compute the Euclidean norm ||v||₂ for N-dimensional vector v.
"""
@inline tuple_norm(v, ::Val{3}) = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
@inline tuple_norm(v, ::Val{4}) = sqrt(v[1]^2 + v[2]^2 + v[3]^2 + v[4]^2)

#####
##### Transpose operation
#####

"""
    transpose_matrix(J, ::Val{3})

Return the transpose of a 3x3 matrix stored as NTuple{9} (row-major).
"""
@inline transpose_matrix(J, ::Val{3}) = (J[1], J[4], J[7],
                                         J[2], J[5], J[8],
                                         J[3], J[6], J[9])

"""
    transpose_matrix(J, ::Val{4})

Return the transpose of a 4x4 matrix stored as NTuple{16} (row-major).
"""
@inline transpose_matrix(J, ::Val{4}) = (J[1], J[5], J[9],  J[13],
                                         J[2], J[6], J[10], J[14],
                                         J[3], J[7], J[11], J[15],
                                         J[4], J[8], J[12], J[16])

#####
##### Outer product (rank-1 update)
#####

"""
    outer_product_update(J, u, v, ::Val{3})

Compute J + u ⊗ v (rank-1 update) for 3x3 matrix J and 3-vectors u, v.
Returns the updated matrix as NTuple{9}.
"""
@inline function outer_product_update(J, u, v, ::Val{3})
    return (J[1] + u[1]*v[1], J[2] + u[1]*v[2], J[3] + u[1]*v[3],
            J[4] + u[2]*v[1], J[5] + u[2]*v[2], J[6] + u[2]*v[3],
            J[7] + u[3]*v[1], J[8] + u[3]*v[2], J[9] + u[3]*v[3])
end

"""
    outer_product_update(J, u, v, ::Val{4})

Compute J + u ⊗ v (rank-1 update) for 4x4 matrix J and 4-vectors u, v.
Returns the updated matrix as NTuple{16}.
"""
@inline function outer_product_update(J, u, v, ::Val{4})
    return (J[1]  + u[1]*v[1], J[2]  + u[1]*v[2], J[3]  + u[1]*v[3], J[4]  + u[1]*v[4],
            J[5]  + u[2]*v[1], J[6]  + u[2]*v[2], J[7]  + u[2]*v[3], J[8]  + u[2]*v[4],
            J[9]  + u[3]*v[1], J[10] + u[3]*v[2], J[11] + u[3]*v[3], J[12] + u[3]*v[4],
            J[13] + u[4]*v[1], J[14] + u[4]*v[2], J[15] + u[4]*v[3], J[16] + u[4]*v[4])
end

#####
##### Good Broyden's method Jacobian inverse update
#####

"""
    broyden_update(J_inv, Δx, ΔF, ::Val{N})

Compute the Good Broyden update for the inverse Jacobian:

    J⁻¹_new = J⁻¹ + (Δx - J⁻¹ ΔF) ⊗ (Δxᵀ J⁻¹) / (Δxᵀ J⁻¹ ΔF)

where:
- J_inv is the current inverse Jacobian approximation
- Δx = x_new - x_old is the step taken
- ΔF = F_new - F_old is the change in residual

This is the "good" Broyden formula that updates J⁻¹ directly,
avoiding the need to solve a linear system at each iteration.
"""
@inline function broyden_update(J_inv, Δx, ΔF, valN::Val{N}) where N
    # Compute J⁻¹ ΔF
    J_inv_ΔF = mat_vec_mul(J_inv, ΔF, valN)

    # Compute Δxᵀ J⁻¹ (= (J⁻ᵀ Δx)ᵀ)
    J_inv_T = transpose_matrix(J_inv, valN)
    Δx_J_inv = mat_vec_mul(J_inv_T, Δx, valN)

    # Compute denominator: Δxᵀ J⁻¹ ΔF
    denom = dot_product(Δx_J_inv, ΔF, valN)

    # Guard against division by zero using ifelse (GPU-safe)
    eps_val = eps(typeof(denom))
    safe_denom = ifelse(abs(denom) < eps_val, one(denom), denom)

    # Compute numerator vector: (Δx - J⁻¹ ΔF) / denom
    diff = subtract_tuples(Δx, J_inv_ΔF, valN)
    u = scale_tuple(inv(safe_denom), diff, valN)

    # Compute rank-1 update: J⁻¹_new = J⁻¹ + u ⊗ Δx_J_inv
    return outer_product_update(J_inv, u, Δx_J_inv, valN)
end

#####
##### State <-> Tuple conversion
#####

"""
    state_to_tuple(Ψ, ::Val{3})

Extract (u★, θ★, q★) from an `InterfaceState` as a tuple.
Used for `BulkTemperature` formulation.
"""
@inline state_to_tuple(Ψ, ::Val{3}) = (Ψ.u★, Ψ.θ★, Ψ.q★)

"""
    state_to_tuple(Ψ, ::Val{4})

Extract (u★, θ★, q★, T) from an `InterfaceState` as a tuple.
Used for `SkinTemperature` formulation.
"""
@inline state_to_tuple(Ψ, ::Val{4}) = (Ψ.u★, Ψ.θ★, Ψ.q★, Ψ.T)

"""
    tuple_to_state(x, Ψ_template, ::Val{3})

Create a new `InterfaceState` with u★, θ★, q★ from tuple x,
and remaining fields from Ψ_template.
"""
@inline function tuple_to_state(x, Ψ, ::Val{3})
    FT = eltype(Ψ)
    return InterfaceState(convert(FT, x[1]),  # u★
                          convert(FT, x[2]),  # θ★
                          convert(FT, x[3]),  # q★
                          Ψ.u, Ψ.v, Ψ.T, Ψ.S, Ψ.q, Ψ.melting)
end

"""
    tuple_to_state(x, Ψ_template, ::Val{4})

Create a new `InterfaceState` with u★, θ★, q★, T from tuple x,
and remaining fields from Ψ_template.
"""
@inline function tuple_to_state(x, Ψ, ::Val{4})
    FT = eltype(Ψ)
    return InterfaceState(convert(FT, x[1]),  # u★
                          convert(FT, x[2]),  # θ★
                          convert(FT, x[3]),  # q★
                          Ψ.u, Ψ.v,
                          convert(FT, x[4]),  # T (from tuple, not template)
                          Ψ.S, Ψ.q, Ψ.melting)
end

#####
##### Convergence check
#####

"""
    broyden_converged(F, tolerance, ::Val{N})

Check if the Broyden iteration has converged based on residual norm.
Returns `true` if ||F||₂ < tolerance.
"""
@inline function broyden_converged(F, tolerance, valN::Val{N}) where N
    norm_F = tuple_norm(F, valN)
    return norm_F < tolerance
end

#####
##### Residual computation helper
#####

"""
    compute_broyden_residual(flux_formulation, Ψₛ, Ψₐ, Ψᵢ, Qᵣ, ℙₛ, ℙₐ, ℙᵢ, ::Val{N})

Compute the residual F = iterate(x) - x for Broyden's method.

The residual measures how far the current state is from being a fixed point
of the interface state iteration. Returns a tuple (F, Ψₛ_new) where F is
the residual and Ψₛ_new is the result of one iteration.
"""
@inline function compute_broyden_residual(flux_formulation,
                                          Ψₛ,
                                          atmosphere_state,
                                          interior_state,
                                          downwelling_radiation,
                                          interface_properties,
                                          atmosphere_properties,
                                          interior_properties,
                                          valN::Val{N}) where N

    # Perform one iteration
    Ψₛ_new = iterate_interface_state(flux_formulation,
                                     Ψₛ,
                                     atmosphere_state,
                                     interior_state,
                                     downwelling_radiation,
                                     interface_properties,
                                     atmosphere_properties,
                                     interior_properties)

    # Extract current and new state as tuples
    x_old = state_to_tuple(Ψₛ, valN)
    x_new = state_to_tuple(Ψₛ_new, valN)

    # Compute residual F = x_new - x_old
    F = subtract_tuples(x_new, x_old, valN)

    return F, Ψₛ_new
end
