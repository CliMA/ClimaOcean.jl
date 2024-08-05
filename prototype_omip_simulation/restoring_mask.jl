
# Build a mask that goes from 0 to 1 as a cubic function of φ between
# 50 degrees and 30 degrees and zero derivatives at 50 and 30.
function cubic_restoring_mask(x₁, x₂, y₁ = 0, y₂ = 1)
    A⁺ = [ x₁^3   x₁^2  x₁ 1
        x₂^3   x₂^2  x₂ 1
        3*x₁^2 2*x₁  1  0
        3*x₂^2 2*x₂  1  0]
            
    b⁺ = [y₁, y₂, 0, 0]
    c⁺ = A⁺ \ b⁺

    # Coefficients for the cubic mask
    c₁⁺ = c⁺[1]
    c₂⁺ = c⁺[2]
    c₃⁺ = c⁺[3]
    c₄⁺ = c⁺[4]

    c₁⁻ = - c⁺[1]
    c₂⁻ = c⁺[2]
    c₃⁻ = - c⁺[3]
    c₄⁻ = c⁺[4]

    @inline mask(λ, φ, z, t) = ifelse(φ <=  50, c₁⁺ * φ^3 + c₂⁺ * φ^2 + c₃⁺ * φ + c₄⁺, zero(eltype(φ)))

    return mask
end