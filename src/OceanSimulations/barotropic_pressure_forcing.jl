using Oceananigans.Operators: ∂xᶠᶜᶜ, ∂yᶜᶠᶜ

struct XDirection end
struct YDirection end

struct BarotropicPressureForcing{D, P}
    direction :: D
    pressure :: P
end

const XDirectionBPF = BarotropicPressureForcing{<:XDirection}
const YDirectionBPF = BarotropicPressureForcing{<:YDirection}

@inline (bpf::XDirectionBPF)(i, j, k, grid, clock, fields) = - ∂xᶠᶜᶜ(i, j, k, grid, bpf.pressure)
@inline (bpf::YDirectionBPF)(i, j, k, grid, clock, fields) = - ∂yᶜᶠᶜ(i, j, k, grid, bpf.pressure)

