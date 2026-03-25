using RushLarsenSolvers
using MTKCardiacCellModels
using MTKCardiacCellModels: BeelerReuterModel
using ModelingToolkit
using Plots

import ModelingToolkit: t_nounits as t, D_nounits as D

function rounded_pulse(t_c; t_start=20.0, duration=5.0, amplitude=10.0, smoothness=0.1)
    t_end = t_start + duration
    rise = 0.5 * (1 + tanh((t_c - t_start) / smoothness))
    fall = 0.5 * (1 - tanh((t_c - t_end) / smoothness))
    return amplitude * rise * fall
end

@register_symbolic rounded_pulse(t)

@named beeler_reuter_system = BeelerReuterModel(I_stim=rounded_pulse)
sys = mtkcompile(beeler_reuter_system)

u0 = [
    1e-7,       # Caᵢ(t) - Intracellular calcium concentration [mM]
    0.975,      # j gate - Sodium slow inactivation
    0.988,      # h gate - Sodium fast inactivation
    0.011,      # m gate - Sodium activation
    0.994,      # f gate - Calcium inactivation
    0.003,      # d gate - Calcium activation
    0.0001,     # x1 gate - Outward current activation
    -84.0       # φₘ(t) - Resting potential [mV]
]

rlf = RushLarsenFunction(sys)

tspan = (0.0, 500.0)

p = ODEProblem(sys, u0, tspan).p
prob = ODEProblem(rlf, u0, tspan, p)

sol = solve(prob, RushLarsen(), dt=0.01)
