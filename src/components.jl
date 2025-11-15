struct GatingVariable <: ModelingToolkit.Symbolics.AbstractVariableMetadata end

ModelingToolkit.Symbolics.option_to_metadata_type(::Val{:gating}) = GatingVariable
is_gating_variable(var) = ModelingToolkit.Symbolics.hasmetadata(var, GatingVariable) ? ModelingToolkit.Symbolics.getmetadata(var, GatingVariable) : false

@component function AlphaBetaGate(; name, α_expr, β_expr)
    @variables begin
        y(t), [gating = true, bounds = (0.0, 1.0)]
        α(t)
        β(t)
    end

    eqs = [
        α ~ α_expr
        β ~ β_expr
        D(y) ~ α * (1 - y) - β * y # = α - (α + β) y
    ]

    System(eqs, t, [y, α, β], [], name=name)
end

# TODO compatibility with Rush-Larsen
#  we have α = 1/τ and β = 
@component function TauGate(; name, y∞, τ)
    @variables y(t), [bounds = (0.0, 1.0)]

    @parameters begin
        y∞ = y∞
        τ = τ
    end

    eqs = [
        D(y) ~ (y∞ - y) / (τ)
    ]

    System(eqs, t, [y], [y∞, τ], name=name)
end

# This is the lIₘit case when τ -> 0
@component function RapidTauGate(; name, y∞_expr, τ)
    @variables y(t), [bounds = (0.0, 1.0)]

    @parameters begin
        y∞ = y∞
    end

    eqs = [
        0 ~ y∞ - y,
        y∞ ~ y∞_expr(t),
    ]
    System(eqs, t, [y, y∞], [y∞, τ], name=name)
end

# @component function IonChannelConductance(; name, g_max, gates, powers=nothing)
#     @variables g(t)
#     @parameters g_max = g_max

#     if isnothing(powers)
#         powers = ones(Int, length(gates))
#     end

#     gate_product = prod(gate.y^exp for (gate, exp) in zip(gates, powers))

#     eqs = [
#         g ~ g_max * gate_product
#     ]

#     System(eqs, t, [g], [g_max], systems=[gates...], name=name)
# end

@mtkmodel NernstPotential begin
    @constants begin
        F = 96.5 # Faraday constant
        R = 8.315 # Gas constant
    end
    @parameters begin
        zₓ # Charge
        T = 310.0  # Absolute temperature
    end
    @variables begin
        Eₓ(t)
        Xₑ(t) # concentration of X on side e (,,extracellular'')
        Xᵢ(t) # concentration of X on side i (,,intracellular'')
    end
    @equations begin
        Eₓ ~ (R * T) / (zₓ * F) * log(Xₑ / Xᵢ)
    end
end

"""
GoldmannHodgkinKatzPotential(Xₑs, Xᵢs, Pₓs, T)

    Xₑs = extracellular concentrations of X
    Xᵢs = extracellular concentrations of X
    Pₓs = permeabilities of X
    T   = absolute temperature
"""
@mtkmodel GoldmannHodgkinKatzPotential begin
    @structural_parameters begin
        N # number of species
    end
    @constants begin
        F = 96.5 # Faraday constant
        R = 8.315 # Gas constant
    end
    @parameters begin
        zₓ # Charge
        T = 310.0  # Absolute temperature
    end
    @variables begin
        Eₓ(t)
        Xₑs(t)[1:N]
        Xᵢs(t)[1:N]
        Pₓs(t)[1:N]
    end
    @equations begin
        Eₓ ~ (zX * F) / R * log(sum(Xₑs .* Pₓs) / sum(Xᵢs .* Pₓs))
    end
end

@mtkmodel OhmicCurrent begin
    @parameters begin
        g
    end

    @variables begin
        Iₘ(t)
        φₘ(t)
        E(t)
    end

    @equations begin
        Iₘ ~ g * (φₘ - E)
    end
end

@mtkmodel GoldmannHodgkinKatzCurrent begin
    @constants begin
        # Faraday and Gas constants
        F = 96.5
        R = 8.315
    end
    @parameters begin
        T = 310 # absolute temperature
        zₓ # charge
        Pₓ # permability
    end
    @variables begin
        Iₘ(t)
        Xₑ(t)
        Xᵢ(t)
    end
    @equations begin
        Iₘ ~ Pₓ * (zₓ^2 * F^2) / (R * T) * (Xₓ - Xₑ * exp(-zₓ * F * φₘ) / (R * T)) / (1 - exp(-zₓ * F * φₘ) / (R * T))
    end
end


#########################


@mtkmodel BeelerReuterIK1 begin
    @variables begin
        Iₘ(t)
        ik2(t)
        φₘ(t)
    end
    @equations begin
        # This equation comes from Adrian (1969) Eqn 15.1.
        ik2 ~ 4.0 * (exp(0.04 * (φₘ + 85.0)) - 1.0) / (exp(0.08 * (φₘ + 53.0)) + exp(0.04 * (φₘ + 53.0)))
        # This is attributed to McAllister, which in turn also attributes this to Adrian (1969), but I cannot find the formula in the paper.
        Iₘ ~ 0.35 * (ik2 + 0.2 * (φₘ + 23.0) / (1.0 - exp(-0.04 * (φₘ + 23.0))))
    end
end

# TODO in their paper they do not state if they got this function from somewhere else. Double check.
function BeelerReuterRate(φₘ, C₁, C₂, C₃, C₄, C₅, C₆, C₇)
    (C₁ * exp(C₂ * (φₘ + C₃)) + C₄ * (φₘ - C₅)) / (exp(C₆ * (φₘ + C₃)) + C₇)
end

@mtkmodel BeelerReuterIx1 begin
    @variables begin
        Iₘ(t)
        φₘ(t)
    end
    @parameters begin
        g_x1 = 0.8
    end
    @components begin
        gate =  AlphaBetaGate(
            # TODO pull parameters from parameter array instead of hardcoding
            α_expr=BeelerReuterRate(φₘ, 0.0005, 0.083, 50.0, 0, 0, 0.057, 1.0),
            β_expr=BeelerReuterRate(φₘ, 0.0013, -0.06, 20.0, 0, 0, -0.04, 1.0),
        )
    end
    @equations begin
        Iₘ ~ g_x1 * gate.y * (exp(0.04 * (φₘ + 77)) - 1) / exp(0.04 * (φₘ + 35))
    end
end

@mtkmodel BeelerReuterICa begin
    @parameters begin
        g_max = 0.09
    end
    @variables begin
        φₘ(t)
        ECa(t)
        Iₘ(t)
    end
    @components begin
        d = AlphaBetaGate(
            # TODO pull parameters from parameter array instead of hardcoding
            α_expr=BeelerReuterRate(φₘ, 0.095, -0.01, -5.0, 0, 0, -0.072, 1.0),
            β_expr=BeelerReuterRate(φₘ, 0.07, -0.017, 44.0, 0, 0, 0.05, 1.0),
        )
        f = AlphaBetaGate(
            # TODO pull parameters from parameter array instead of hardcoding
            α_expr=BeelerReuterRate(φₘ, 0.012, -0.008, 28.0, 0, 0, 0.15, 1.0),
            β_expr=BeelerReuterRate(φₘ, 0.0065, -0.02, 30.0, 0, 0, -0.2, 1.0),
        )
        i = OhmicCurrent(g=g_max)
    end
    @equations begin
        i.φₘ ~ φₘ
        i.E ~ ECa
        Iₘ ~ f.y * d.y * i.Iₘ
    end
end

@mtkmodel BeelerReuterINa begin
    @parameters begin
        g_max = 4.0
        ENa = 50.0
    end
    @variables begin
        φₘ(t)
        Iₘ(t)
    end
    @components begin
        # TODO use BeelerReuterRate
        m = AlphaBetaGate(
            α_expr=(φₘ + 47.0) / (1.0 - exp(-0.1 * (φₘ + 47.0))),
            β_expr=40.0 * exp(-0.056 * (φₘ + 72.0))
        )
        h = AlphaBetaGate(
            α_expr=0.126 * exp(-0.25 * (φₘ + 77.0)),
            β_expr=1.7 / (1.0 + exp(-0.082 * (φₘ + 22.5)))
        )
        j = AlphaBetaGate(
            α_expr=0.055 * exp(-0.25 * (φₘ + 78)) / (1 + exp(-0.2 * (φₘ + 78))),
            β_expr=0.3 / (1 + exp(-0.1 * (φₘ + 32.0)))
        )
        i = OhmicCurrent(g=g_max)
    end
    @equations begin
        i.φₘ ~ φₘ
        i.E ~ ENa
        Iₘ ~ m.y^3 * h.y * j.y * i.Iₘ
    end
end

@mtkmodel BeelerReuterCalciumDynamics begin
    @parameters begin
        T     = 310.0, [description = "Absolute temperature"]
    end
    @variables begin
        Caᵢ(t)
        iCa(t)
    end
    @equations begin
        D(Caᵢ) ~ -10.0^-7 * iCa + 0.07 * (10.0^-7 - Caᵢ)
    end
end
