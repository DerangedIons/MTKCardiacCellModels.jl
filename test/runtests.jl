using MTKCardiacCellModels
using MTKCardiacCellModels: AlphaBetaGate, TauGate, RapidTauGate,
    NernstPotential, OhmicCurrent,
    BeelerReuterIK1, BeelerReuterIx1, BeelerReuterICa, BeelerReuterINa,
    BeelerReuterCalciumDynamics, BeelerReuterModel
using ModelingToolkit
using Test

import ModelingToolkit: t_nounits as t, D_nounits as D

@testset "MTKCardiacCellModels.jl" begin
    @testset "Gates" begin
        @variables V(t)
        @test AlphaBetaGate(name=:test, α_expr=V, β_expr=1-V) isa System
        @test TauGate(name=:test, y∞=0.5, τ=1.0) isa System
        @test RapidTauGate(name=:test, y∞_expr=V) isa System
    end

    @testset "Generic components" begin
        @test NernstPotential(name=:test, zₓ=2) isa System
        @test OhmicCurrent(name=:test, g=1.0) isa System
    end

    @testset "Beeler-Reuter components" begin
        @test BeelerReuterIK1(name=:test) isa System
        @test BeelerReuterIx1(name=:test) isa System
        @test BeelerReuterICa(name=:test) isa System
        @test BeelerReuterINa(name=:test) isa System
        @test BeelerReuterCalciumDynamics(name=:test) isa System
    end

    @testset "Beeler-Reuter full model" begin
        stim(t) = 0.0
        @register_symbolic stim(t)
        @named model = BeelerReuterModel(I_stim=stim)
        sys = mtkcompile(model)
        @test sys isa System
        @test length(unknowns(sys)) == 8
    end
end
