using RushLarsenSolvers
using MTKCardiacCellModels
using MTKCardiacCellModels: AlphaBetaGate, IonChannelConductance, IonCurrent, LeakCurrent,
    SlowInwardCurrent, TimeIndependentK, TimeActivatedOutward, CalciumDynamics, StimulationSystem
using ModelingToolkit
using Plots

import ModelingToolkit: t_nounits as t, D_nounits as D

@variables V(t) Ca_i(t)
V = GlobalScope(V)
Ca_i = GlobalScope(Ca_i)

# Fast sodium current (INa) gates
@named m = AlphaBetaGate(
    α_expr=(V + 47.0) / (1.0 - exp(-0.1 * (V + 47.0))),
    β_expr=40.0 * exp(-0.056 * (V + 72.0))
)

@named h = AlphaBetaGate(
    α_expr=0.126 * exp(-0.25 * (V + 77.0)),
    β_expr=1.7 / (1.0 + exp(-0.082 * (V + 22.5)))
)

@named j = AlphaBetaGate(
    α_expr=0.055 * exp(-0.25 * (V + 78)) / (1 + exp(-0.2 * (V + 78))),
    β_expr=0.3 / (1 + exp(-0.1 * (V + 32.0)))
)

# Fast sodium channel conductance (includes leak component)
@named sodium = IonChannelConductance(g_max=4.0, gates=[m, h, j], powers=[3.0, 1.0, 1.0])
@named sod_current = IonCurrent(Vm=V, gm=sodium, Ex=50.0)

# Sodium leak current
@named na_leak = LeakCurrent(Vm=V, gL=0.003, EL=50.0)

# Slow inward current (Is) gates
@named d = AlphaBetaGate(
    α_expr=0.095 * exp(-0.01 * (V - 5.0)) / (1.0 + exp(-0.072 * (V - 5.0))),
    β_expr=0.07 * exp(-0.017 * (V + 44.0)) / (1.0 + exp(0.05 * (V + 44.0)))
)

@named f = AlphaBetaGate(
    α_expr=0.012 * exp(-0.008 * (V + 28.0)) / (1.0 + exp(0.15 * (V + 28.0))),
    β_expr=0.0065 * exp(-0.02 * (V + 30.0)) / (1.0 + exp(-0.2 * (V + 30.0)))
)

# Slow inward current (calcium)
@named slow_inward = IonChannelConductance(g_max=0.09, gates=[d, f])

@named slow_current = SlowInwardCurrent(Vm=V, Ca_i=Ca_i, gm=slow_inward)


@named k1_current = TimeIndependentK(Vm=V)

# Time-activated outward current (Ix1) gate
@named x1 = AlphaBetaGate(
    α_expr=0.0005 * exp(0.083 * (V + 50.0)) / (1.0 + exp(0.057 * (V + 50.0))),
    β_expr=0.0013 * exp(-0.06 * (V + 20.0)) / (1.0 + exp(-0.04 * (V + 20.0)))
)

@named x1_current = TimeActivatedOutward(Vm=V, gate=x1, g_x1=0.8)

@named calcium_dyn = CalciumDynamics(Ca_i=Ca_i, I_s=slow_current.Im)

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


@named stim = StimulationSystem(stim_func=rounded_pulse)


@component function BeelerReuterCompartment(; name, Vm, Cm, currents, I_stim, Ca_dyn)
    @parameters Cm = Cm
    Is = [current.Im for current in currents]

    eqs = [D(Vm) ~ -1 / Cm * (sum(Is) - I_stim)]

    System(eqs, t, [Vm, I_stim], [Cm], systems=[currents..., Ca_dyn], name=name)
end

@named com = BeelerReuterCompartment(
    Vm=V,
    Cm=1.0,
    currents=[sod_current, na_leak, slow_current, k1_current, x1_current],
    I_stim=stim.I_stim,
    Ca_dyn=calcium_dyn
)


@named beeler_reuter_system = compose(com, stim)

sys = structural_simplify(beeler_reuter_system)

rlf = RushLarsenFunction(sys)

u0 = [
    1e-7,       # Ca_i(t) - Intracellular calcium concentration [mM]
    0.0001,     # x1_current₊x1₊y(t) - Outward current activation gate
    0.994,      # slow_current₊slow_inward₊f₊y(t) - Calcium inactivation gate
    0.003,      # slow_current₊slow_inward₊d₊y(t) - Calcium activation gate
    0.975,      # sod_current₊sodium₊j₊y(t) - Sodium slow inactivation gate
    0.988,      # sod_current₊sodium₊h₊y(t) - Sodium fast inactivation gate
    0.011,      # sod_current₊sodium₊m₊y(t) - Sodium activation gate
    -84.0       # V(t) - Resting potential [mV]
]

p = [1.0,
    50.0,
    4.0,
    50.0,
    0.003,
    0.09,
    0.8
]

# Time span for simulation
tspan = (0.0, 500.0)  # 2000 ms

prob = ODEProblem(rlf, u0, tspan, p)

sol = solve(prob, RushLarsen(), dt=0.01)