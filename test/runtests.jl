using MTKCardiacCellModels
using MTKCardiacCellModels: BeelerReuterModel
using ModelingToolkit
using OrdinaryDiffEq
using Test

import ModelingToolkit: t_nounits as t, D_nounits as D

@testset "MTKCardiacCellModels.jl" begin
    @testset "Beeler-Reuter action potential" begin
        function pulse(t_c; t_start=2.0, duration=1.0, amplitude=25.0, k=10.0)
            t_end = t_start + duration
            return amplitude * (1 / (1 + exp(-k * (t_c - t_start)))) * (1 / (1 + exp(k * (t_c - t_end))))
        end
        @register_symbolic pulse(t)

        @named model = BeelerReuterModel(I_stim=pulse)
        sys = mtkcompile(model)

        u0 = unknowns(sys) .=> [
            1e-7,   # Caᵢ
            0.975,  # j gate
            0.988,  # h gate
            0.011,  # m gate
            0.994,  # f gate
            0.003,  # d gate
            0.0001, # x1 gate
            -84.0,  # φₘ
        ]
        prob = ODEProblem(sys, u0, (0.0, 500.0))
        sol = solve(prob, FBDF(); abstol=1e-8, reltol=1e-6)

        @test SciMLBase.successful_retcode(sol)
        @test all(isfinite, sol.u[end])

        # Membrane potential should depolarize (peak > 0 mV) during the action potential
        φₘ_idx = findfirst(s -> contains(string(s), "φₘ"), string.(unknowns(sys)))
        φₘ_vals = [u[φₘ_idx] for u in sol.u]
        @test maximum(φₘ_vals) > 0.0
    end
end
