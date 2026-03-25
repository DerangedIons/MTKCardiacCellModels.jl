struct GatingVariable <: ModelingToolkit.Symbolics.AbstractVariableMetadata end

ModelingToolkit.Symbolics.option_to_metadata_type(::Val{:gating}) = GatingVariable
is_gating_variable(var) = ModelingToolkit.Symbolics.hasmetadata(var, GatingVariable) ? ModelingToolkit.Symbolics.getmetadata(var, GatingVariable) : false

const FARADAY = 96.5    # Faraday constant [C/mmol]
const GAS_CONST = 8.315 # Gas constant [J/(mol·K)]

# ─── Gates ────────────────────────────────────────────────────────────────────

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

# Limit case when τ -> 0
@component function RapidTauGate(; name, y∞_expr)
    @variables begin
        y(t), [bounds = (0.0, 1.0)]
        y∞(t)
    end

    eqs = [
        y∞ ~ y∞_expr,
        0 ~ y∞ - y,
    ]

    System(eqs, t, [y, y∞], [], name=name)
end

# ─── Generic components ───────────────────────────────────────────────────────

@component function NernstPotential(; name, zₓ, T=310.0)
    @parameters begin
        zₓ = zₓ
        T = T
    end
    @variables begin
        Eₓ(t)
        Xₑ(t)
        Xᵢ(t)
    end

    eqs = [
        Eₓ ~ (GAS_CONST * T) / (zₓ * FARADAY) * log(Xₑ / Xᵢ)
    ]

    System(eqs, t, [Eₓ, Xₑ, Xᵢ], [zₓ, T], name=name)
end

@doc """
    GoldmannHodgkinKatzPotential(; name, N, zₓ, T=310.0)

    Xₑs = extracellular concentrations of X
    Xᵢs = intracellular concentrations of X
    Pₓs = permeabilities of X
    T   = absolute temperature
"""
@component function GoldmannHodgkinKatzPotential(; name, N, zₓ, T=310.0)
    @parameters begin
        zₓ = zₓ
        T = T
    end
    @variables begin
        Eₓ(t)
        Xₑs(t)[1:N]
        Xᵢs(t)[1:N]
        Pₓs(t)[1:N]
    end

    eqs = [
        Eₓ ~ (zₓ * FARADAY) / GAS_CONST * log(sum(Xₑs .* Pₓs) / sum(Xᵢs .* Pₓs))
    ]

    System(eqs, t, [Eₓ, collect(Xₑs)..., collect(Xᵢs)..., collect(Pₓs)...], [zₓ, T], name=name)
end

@component function OhmicCurrent(; name, g)
    @parameters g = g
    @variables begin
        Iₘ(t)
        φₘ(t)
        E(t)
    end

    eqs = [
        Iₘ ~ g * (φₘ - E)
    ]

    System(eqs, t, [Iₘ, φₘ, E], [g], name=name)
end

@component function GoldmannHodgkinKatzCurrent(; name, zₓ, Pₓ, T=310.0)
    @parameters begin
        T = T
        zₓ = zₓ
        Pₓ = Pₓ
    end
    @variables begin
        Iₘ(t)
        φₘ(t)
        Xₑ(t)
        Xᵢ(t)
    end

    eqs = [
        Iₘ ~ Pₓ * (zₓ^2 * FARADAY^2) / (GAS_CONST * T) * (Xᵢ - Xₑ * exp(-zₓ * FARADAY * φₘ / (GAS_CONST * T))) / (1 - exp(-zₓ * FARADAY * φₘ / (GAS_CONST * T)))
    ]

    System(eqs, t, [Iₘ, φₘ, Xₑ, Xᵢ], [T, zₓ, Pₓ], name=name)
end

# ─── Beeler-Reuter model ─────────────────────────────────────────────────────

# TODO in their paper they do not state if they got this function from somewhere else. Double check.
function BeelerReuterRate(φₘ, C₁, C₂, C₃, C₄, C₅, C₆, C₇)
    (C₁ * exp(C₂ * (φₘ + C₃)) + C₄ * (φₘ - C₅)) / (exp(C₆ * (φₘ + C₃)) + C₇)
end

@component function BeelerReuterIK1(; name)
    @variables begin
        Iₘ(t)
        ik2(t)
        φₘ(t)
    end

    eqs = [
        # This equation comes from Adrian (1969) Eqn 15.1.
        ik2 ~ 4.0 * (exp(0.04 * (φₘ + 85.0)) - 1.0) / (exp(0.08 * (φₘ + 53.0)) + exp(0.04 * (φₘ + 53.0)))
        # This is attributed to McAllister, which in turn also attributes this to Adrian (1969), but I cannot find the formula in the paper.
        Iₘ ~ 0.35 * (ik2 + 0.2 * (φₘ + 23.0) / (1.0 - exp(-0.04 * (φₘ + 23.0))))
    ]

    System(eqs, t, [Iₘ, ik2, φₘ], [], name=name)
end

@component function BeelerReuterIx1(; name, g_x1=0.8)
    @variables begin
        Iₘ(t)
        φₘ(t)
    end
    @parameters g_x1 = g_x1

    # TODO pull parameters from parameter array instead of hardcoding
    φₘ_parent = ParentScope(φₘ)
    gate = AlphaBetaGate(
        name=:gate,
        α_expr=BeelerReuterRate(φₘ_parent, 0.0005, 0.083, 50.0, 0, 0, 0.057, 1.0),
        β_expr=BeelerReuterRate(φₘ_parent, 0.0013, -0.06, 20.0, 0, 0, -0.04, 1.0),
    )

    eqs = [
        Iₘ ~ g_x1 * gate.y * (exp(0.04 * (φₘ + 77)) - 1) / exp(0.04 * (φₘ + 35))
    ]

    System(eqs, t, [Iₘ, φₘ], [g_x1], systems=[gate], name=name)
end

@component function BeelerReuterICa(; name, g_max=0.09)
    @parameters g_max = g_max
    @variables begin
        φₘ(t)
        ECa(t)
        Iₘ(t)
    end

    # TODO pull parameters from parameter array instead of hardcoding
    φₘ_parent = ParentScope(φₘ)
    d = AlphaBetaGate(
        name=:d,
        α_expr=BeelerReuterRate(φₘ_parent, 0.095, -0.01, -5.0, 0, 0, -0.072, 1.0),
        β_expr=BeelerReuterRate(φₘ_parent, 0.07, -0.017, 44.0, 0, 0, 0.05, 1.0),
    )
    f = AlphaBetaGate(
        name=:f,
        α_expr=BeelerReuterRate(φₘ_parent, 0.012, -0.008, 28.0, 0, 0, 0.15, 1.0),
        β_expr=BeelerReuterRate(φₘ_parent, 0.0065, -0.02, 30.0, 0, 0, -0.2, 1.0),
    )

    eqs = [
        Iₘ ~ g_max * d.y * f.y * (φₘ - ECa)
    ]

    System(eqs, t, [φₘ, ECa, Iₘ], [g_max], systems=[d, f], name=name)
end

@component function BeelerReuterINa(; name, g_max=4.0, ENa=50.0)
    @parameters begin
        g_max = g_max
        ENa = ENa
    end
    @variables begin
        φₘ(t)
        Iₘ(t)
    end

    # TODO use BeelerReuterRate
    φₘ_parent = ParentScope(φₘ)
    m = AlphaBetaGate(
        name=:m,
        α_expr=(φₘ_parent + 47.0) / (1.0 - exp(-0.1 * (φₘ_parent + 47.0))),
        β_expr=40.0 * exp(-0.056 * (φₘ_parent + 72.0)),
    )
    h = AlphaBetaGate(
        name=:h,
        α_expr=0.126 * exp(-0.25 * (φₘ_parent + 77.0)),
        β_expr=1.7 / (1.0 + exp(-0.082 * (φₘ_parent + 22.5))),
    )
    j = AlphaBetaGate(
        name=:j,
        α_expr=0.055 * exp(-0.25 * (φₘ_parent + 78)) / (1 + exp(-0.2 * (φₘ_parent + 78))),
        β_expr=0.3 / (1 + exp(-0.1 * (φₘ_parent + 32.0))),
    )

    eqs = [
        Iₘ ~ g_max * m.y^3 * h.y * j.y * (φₘ - ENa)
    ]

    System(eqs, t, [φₘ, Iₘ], [g_max, ENa], systems=[m, h, j], name=name)
end

@component function BeelerReuterCalciumDynamics(; name)
    @variables begin
        Caᵢ(t)
        iCa(t)
    end

    eqs = [
        D(Caᵢ) ~ -1e-7 * iCa + 0.07 * (1e-7 - Caᵢ)
    ]

    System(eqs, t, [Caᵢ, iCa], [], name=name)
end

@component function BeelerReuterModel(; name, I_stim, gNaL=0.003, ENa=50.0, Cₘ=1.0, Caₒ=0.0021021513474235995)
    @parameters begin
        gNaL = gNaL
        ENa = ENa
        Cₘ = Cₘ
        Caₒ = Caₒ
    end
    @variables begin
        φₘ(t)
        Caᵢ(t)
    end

    ECa = NernstPotential(name=:ECa, zₓ=2)
    cK1 = BeelerReuterIK1(name=:cK1)
    cx1 = BeelerReuterIx1(name=:cx1)
    cCa = BeelerReuterICa(name=:cCa)
    cNa = BeelerReuterINa(name=:cNa)
    calcium = BeelerReuterCalciumDynamics(name=:calcium)

    eqs = [
        # Nernst reversal potential for calcium
        ECa.Xᵢ ~ Caᵢ
        ECa.Xₑ ~ Caₒ
        # Membrane voltage equation (leak current inlined as gNaL * (φₘ - ENa))
        D(φₘ) ~ -1 / Cₘ * (gNaL * (φₘ - ENa) + cNa.Iₘ + cK1.Iₘ + cx1.Iₘ + cCa.Iₘ - I_stim(t))
        # Channel voltage connections
        cx1.φₘ ~ φₘ
        cK1.φₘ ~ φₘ
        cNa.φₘ ~ φₘ
        cCa.φₘ ~ φₘ
        cCa.ECa ~ ECa.Eₓ
        # Calcium dynamics
        calcium.Caᵢ ~ Caᵢ
        calcium.iCa ~ cCa.Iₘ
    ]

    System(eqs, t, [φₘ, Caᵢ], [gNaL, ENa, Cₘ, Caₒ],
        systems=[ECa, cK1, cx1, cCa, cNa, calcium], name=name)
end
