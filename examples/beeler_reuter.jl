# using RushLarsenSolvers
using MTKCardiacCellModels
using MTKCardiacCellModels: AlphaBetaGate, NernstPotential, OhmicCurrent
using MTKCardiacCellModels: BeelerReuterIK1, BeelerReuterIx1, BeelerReuterICa, BeelerReuterINa, BeelerReuterCalciumDynamics
using ModelingToolkit
# using Plots

import ModelingToolkit: t_nounits as t, D_nounits as D

@mtkmodel BeelerReuterModel begin
    @structural_parameters begin
        I_stim
    end
    @parameters begin
        gNaL = 0.003
        ENa = 50.0
        Cₘ = 1.0
        Caₒ   = 0.0021021513474235995, [description = "Extracellular calcium concentration"]
    end
    @variables begin
        φₘ(t)
        Caᵢ(t)
    end
    @components begin
        # Potentials
        ECa = NernstPotential(zₓ=2)
        # Channels
        cNa_leak = OhmicCurrent(g = gNaL)
        cK1 = BeelerReuterIK1()
        cx1 = BeelerReuterIx1()
        cCa = BeelerReuterICa()
        cNa = BeelerReuterINa()
        # Concentrations
        calcium = BeelerReuterCalciumDynamics()
    end
    @equations begin
        # Reversal Potential
        ECa.Xᵢ ~ Caᵢ
        ECa.Xₑ ~ Caₒ
        # Membrane equation
        D(φₘ) ~ -1 / Cₘ * (cNa_leak.Iₘ + cNa.Iₘ + cK1.Iₘ + cx1.Iₘ + cCa.Iₘ - I_stim(t))
        # Leak current
        cNa_leak.φₘ ~ φₘ
        cNa_leak.E  ~ ENa
        # x1 channel
        cx1.φₘ ~ φₘ
        # K1 channel
        cK1.φₘ ~ φₘ
        # INa channel
        cNa.φₘ ~ φₘ
        # Slow inward current
        cCa.φₘ ~ φₘ
        cCa.ECa ~ ECa.Eₓ
        # Calcium dynamics
        calcium.Caᵢ ~ Caᵢ
        calcium.iCa ~ cCa.Iₘ
    end
end

function rounded_pulse(t_c; t_start=20.0, duration=5.0, amplitude=10.0, smoothness=0.1)
    """
    Creates a smooth pulse that starts at 0, rises to amplitude, holds, then falls back to 0.

    Parameters:
    - t: time
    - t_start: when the pulse begins
    - duration: how long the pulse stays at amplitude
    - amplitude: peak value of the pulse
    - smoothness: controls the rise/fall smoothness (smaller = sharper, larger = smoother)
    """
    t_end = t_start + duration

    # Smooth step up using tanh
    rise = 0.5 * (1 + tanh((t_c - t_start) / smoothness))

    # Smooth step down using tanh
    fall = 0.5 * (1 - tanh((t_c - t_end) / smoothness))

    return amplitude * rise * fall
end

@register_symbolic rounded_pulse(t)

@mtkcompile beeler_reuter_system = BeelerReuterModel(;I_stim=rounded_pulse)

u0 = [
    1e-7,       # Caᵢ(t) - Intracellular calcium concentration [mM]
    0.975,      # sod_current₊sodium₊j₊y(t) - Sodium slow inactivation gate
    0.988,      # sod_current₊sodium₊h₊y(t) - Sodium fast inactivation gate
    0.011,      # sod_current₊sodium₊m₊y(t) - Sodium activation gate
    0.994,      # slow_current₊slow_inward₊f₊y(t) - Calcium inactivation gate
    0.003,      # slow_current₊slow_inward₊d₊y(t) - Calcium activation gate
    0.0001,     # x1_current₊x1₊y(t) - Outward current activation gate
    -84.0       # V(t) - Resting potential [mV]
]

prob = ODEProblem(beeler_reuter_system, u0, (0.0, 400.0))

using OrdinaryDiffEqTsit5, Plots
sol = solve(prob, Tsit5())
plot(sol, idxs=[8])
# rlf = RushLarsenFunction(sys)


# p = [1.0,
#     50.0,
#     4.0,
#     50.0,
#     0.003,
#     0.09,
#     0.8
# ]

# # Time span for simulation
# tspan = (0.0, 500.0)  # 2000 ms

# prob = ODEProblem(rlf, u0, tspan, p)

# sol = solve(prob, RushLarsen(), dt=0.01)