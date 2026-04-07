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

rlf = RushLarsenFunction(sys)

tspan = (0.0, 500.0)

base_prob = ODEProblem(sys, [], tspan)
prob = ODEProblem(rlf, base_prob.u0, tspan, base_prob.p)

sol = solve(prob, RushLarsen(), dt=0.01)
