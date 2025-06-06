using PythonCall

struct PyCOAREFluxFormulation{S, P, N}
    mypythonscript :: S
    np :: P
    solver_stop_criteria :: N
end

using ClimaOcean.OceanSeaIceModels.InterfaceComputations: InterfaceState, zero_interface_state
import ClimaOcean.OceanSeaIceModels.InterfaceComputations: compute_interface_state

sys = pyimport("sys")
sys.path.insert(0, "/Users/simonesilvestri/development/COARE-algorithm/Python/COARE3.6")

PyCOAREFluxFormulation() = PyCOAREFluxFormulation(pyimport("coare36vn_zrf_et"), pyimport("numpy"), nothing)

# Iterating condition for the characteristic scales solvers
@inline function compute_interface_state(flux_formulation::PyCOAREFluxFormulation,
                                         initial_interface_state,
                                         atmosphere_state,
                                         interior_state,
                                         downwelling_radiation,
                                         interface_properties,
                                         atmosphere_properties,
                                         interior_properties)

    Ψₐ = atmosphere_state
    Ψᵢ = interior_state
    Ψₛ = initial_interface_state

    np = flux_formulation.np
    uₐ = Ψₐ.u
    vₐ = Ψₐ.v
    uᵢ = Ψᵢ.u
    vᵢ = Ψᵢ.v

    𝒰 = np.array([sqrt((uₐ - uᵢ)^2 + (vₐ - vᵢ)^2)])

    rₐ = np.array([Ψₐ.r])
    zu = 10.0
    zt = 10.0
    zq = 10.0
    Sᵢ = np.array([Ψᵢ.S])
    pₐ = np.array([Ψₐ.𝒬.p / 101.325]) # Should be in mb
    Tₐ = np.array([Ψₐ.𝒬.T - 273.15])  # Should be in ᵒC
    Tᵢ = np.array([Ψᵢ.T - 273.15])    # Should be in ᵒC
    Qs = np.array([0.0])
    Qℓ = np.array([0.0])
    zi = 512

    # out = np.array([usr,tau,hsb,hlb,hbb,hsbb,hlwebb,tsr,qsr,zo,zot,zoq,Cd,Ch,Ce,L,zeta,dT_skinx,dq_skinx,dz_skin,Urf,Trf,Qrf,RHrf,UrfN,TrfN,QrfN,lw_net,sw_net,Le,rhoa,UN,U10,U10N,Cdn_10,Chn_10,Cen_10,hrain,Qs,Evap,T10,T10N,Q10,Q10N,RH10,P10,rhoa10,gust,wc_frac,Edis])

    output = flux_formulation.mypythonscript.coare36vn_zrf_et(𝒰, zu, Tₐ, zt, rₐ, zq, pₐ, Tᵢ, 0, 0, 0, 0, 1, zi, 0, Sᵢ) 
    output = PythonCall.pyconvert(Array, output[0])
    
    u★ = convert(eltype(Ψₐ.r), output[1])
    θ★ = convert(eltype(Ψₐ.r), output[8])
    q★ = convert(eltype(Ψₐ.r), output[9])
    Ψₛ = InterfaceState(u★, θ★, q★, uᵢ, vᵢ, Ψᵢ.T , Ψᵢ.S, Ψₐ.r)

    return Ψₛ
end