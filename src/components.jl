struct GatingVariable <: ModelingToolkit.Symbolics.AbstractVariableMetadata end

ModelingToolkit.Symbolics.option_to_metadata_type(::Val{:gating}) = GatingVariable
is_gating_variable(var) = ModelingToolkit.Symbolics.hasmetadata(var, GatingVariable) ? ModelingToolkit.Symbolics.getmetadata(var, GatingVariable) : false

@component function AlphaBetaGate(; name, α_expr, β_expr)
    @variables begin
        y(t), [gating = true]
        α(t)
        β(t)
    end

    eqs = [
        α ~ α_expr
        β ~ β_expr
        D(y) ~ α * (1 - y) - β * y
    ]

    System(eqs, t, [y, α, β], [], name=name)
end

@component function SteadyStateTauGate(; name, y∞, τ)
    @variables y(t)

    @parameters begin
        y∞ = y∞
        τ = τ
    end

    eqs = [D(y) ~ (y∞ - y) / (τ)]
    System(eqs, t, [y], [y∞, τ], name=name)
end


@component function IonChannelConductance(; name, g_max, gates, powers=nothing)
    @variables g(t)
    @parameters g_max = g_max

    if isnothing(powers)
        powers = ones(Int, length(gates))
    end

    gate_product = prod(gate.y^exp for (gate, exp) in zip(gates, powers))

    eqs = [
        g ~ g_max * gate_product
    ]

    System(eqs, t, [g], [g_max], systems=[gates...], name=name)
end

"""
NernstPotential(Xₑ, Xᵢ, zₓ, T)

    Xₑ = extracellular concentration of X
    Xᵢ = extracellular concentration of X
    zₓ = charge of X
    T  = absolute temperature
"""
@component function NernstPotential(; name, Xₑ, Xᵢ, zₓ, T)
    @variables Eₓ(t)
    # Faraday and Gas constants
    @parameters F = 96.5 , R = 8.315
    eqs = [
        Eₓ ~ (zX * F) / (zₓ * R) * log(Xₑ / Xᵢ)
    ]

    System(eqs, t, [Eₓ], [F, R], systems=[gm], name=name)
end

"""
GoldmannHodgkinKatzPotential(Xₑs, Xᵢs, Pₓs, T)

    Xₑs = extracellular concentrations of X
    Xᵢs = extracellular concentrations of X
    Pₓs = permeabilities of X
    T   = absolute temperature
"""
@component function GoldmannHodgkinKatzPotential(; name, Xₑs, Xᵢs, Pₓs, T)
    @variables Eₓ(t)
    # Faraday and Gas constants
    @parameters F = 96.5 , R = 8.315
    eqs = [
        Eₓ ~ (zX * F) / R * log(sum(Xₑs .* Pₓs) / sum(Xᵢs .* Pₓs))
    ]

    System(eqs, t, [Eₓ], [F, R], systems=[gm], name=name)
end

@component function OhmicCurrent(; name, Vm, gm, Ex)
    @variables Im(t)
    @parameters Ex = Ex
    eqs = [
        Im ~ gm.g * (Vm - Ex)
    ]

    System(eqs, t, [Vm, Im], [Ex], systems=[gm], name=name)
end

"""
GoldmannHodgkinKatzCurrent(Xₑ, Xᵢ, Pₓ, zₓ, T)

    Xₑs = extracellular concentrations of X
    Xᵢs = extracellular concentrations of X
    Pₓs = permeabilities of X
    T   = absolute temperature
"""
@component function GoldmannHodgkinKatzCurrent(; name, φₘ, Xₑ, Xᵢ, Pₓ, zₓ, T)
    @variables I(t)
    # Faraday and Gas constants
    @parameters F = 96.5 , R = 8.315
    eqs = [
        I ~ Pₓ * (zₓ^2 * F^2) / (R * T) * (Xₓ - Xₑ * exp(-zₓ) * F * φₘ / (R * T)) / (1 - exp(-zₓ) * F * φₘ / (R * T))
    ]

    System(eqs, t, [I], [F, R], systems=[gm], name=name)
end

@component function LeakCurrent(; name, Vm, gL, EL)
    @variables Im(t)

    @parameters begin
        gL = gL
        EL = EL
    end

    eqs = [
        Im ~ gL * (Vm - EL)
    ]

    System(eqs, t, [Vm, Im], [EL, gL], name=name)
end

@component function CompartmentSystem(; name, Vm, Cm, currents, stim)
    @parameters Cm = Cm
    Is = (current.Im for current in currents)

    eqs = [D(Vm) ~ -1 / Cm * (sum(Is) - stim(t))]

    System(eqs, t, [Vm], [Cm], systems=[currents...], name=name)
end

# FIXME
@component function SlowInwardCurrent(; name, Vm, Ca_i, gm)
    @variables Im(t) E_s(t)

    eqs = [
        E_s ~ -82.3 - 13.0287 * log(Ca_i)
        Im ~ gm.g * (Vm - E_s)
    ]

    System(eqs, t, [Vm, Ca_i, Im, E_s], [], systems=[gm], name=name)
end

# FIXME
@component function TimeIndependentK(; name, Vm)
    @variables Im(t)

    eqs = [
        Im ~ 0.35 * (4.0 * (exp(0.04 * (Vm + 85.0)) - 1.0) /
                     (exp(0.08 * (Vm + 53.0)) + exp(0.04 * (Vm + 53.0))) +
                     0.2 * (Vm + 23.0) / (1.0 - exp(-0.04 * (Vm + 23.0))))
    ]

    System(eqs, t, [Vm, Im], [], name=name)
end

# FIXME
@component function TimeActivatedOutward(; name, Vm, gate, g_x1)
    @variables Im(t)
    @parameters g_x1 = g_x1

    eqs = [
        Im ~ g_x1 * gate.y * (exp(0.04 * (Vm + 77)) - 1) / exp(0.04 * (Vm + 35))
    ]

    System(eqs, t, [Vm, Im], [g_x1], systems=[gate], name=name)
end

# FIXME
@component function CalciumDynamics(; name, Ca_i, I_s)
    eqs = [
        D(Ca_i) ~ -10.0^-7 * I_s + 0.07 * (10.0^-7 - Ca_i)
    ]

    System(eqs, t, [Ca_i, I_s], [], name=name)
end

@component function StimulationSystem(; name, stim_func)
    @variables I_stim(t)

    eqs = [
        I_stim ~ stim_func(t)
    ]

    System(eqs, t, [I_stim], [], name=name)
end
