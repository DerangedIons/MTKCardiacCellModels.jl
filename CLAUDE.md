# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MTKCardiacCellModels is a Julia package providing ModelingToolkit (MTK) v11 components for lumped cardiac cell models. It implements reusable building blocks (gates, currents, reversal potentials) that can be composed into whole-cell action potential models.

## Commands

```bash
# Run tests
julia -e "using Pkg; Pkg.test()"

# Load package in REPL (from dev directory)
julia -e "using Pkg; Pkg.activate(\".\"); using MTKCardiacCellModels"
```

## Architecture

All components live in `src/components.jl` and follow a three-layer hierarchy:

1. **Gate layer** — `AlphaBetaGate`, `TauGate`, `RapidTauGate`: Hodgkin-Huxley gating kinetics (ODEs or algebraic)
2. **Biophysical layer** — `NernstPotential`, `GoldmannHodgkinKatzPotential`, `OhmicCurrent`, `GoldmannHodgkinKatzCurrent`: generic electrochemical primitives
3. **Model layer** — `BeelerReuterIK1`, `BeelerReuterIx1`, `BeelerReuterICa`, `BeelerReuterINa`, `BeelerReuterCalciumDynamics`, `BeelerReuterModel`: Beeler-Reuter 1977 channel and whole-cell implementations

`BeelerReuterModel` composes all channel components and calcium dynamics into a single system via MTK's `systems=` parameter.

## Key Patterns

- **`@component` + `System()`**: Every component is a function decorated with `@component` that returns a `System(eqs, t, vars, params; systems=[], name)`. This is the MTK v11 pattern (not the older `@mtkmodel`/`ODESystem` pattern).
- **`ParentScope(φₘ)`**: Child components (gates) reference the parent's membrane voltage by wrapping it in `ParentScope`. This is how voltage-dependent rate expressions in gates connect to the enclosing channel's `φₘ`.
- **`GatingVariable` metadata**: Custom metadata type registered via `option_to_metadata_type(::Val{:gating})`. Gating variables are declared with `[gating = true, bounds = (0.0, 1.0)]`.
- **`t_nounits` / `D_nounits`**: The package uses unitless time and derivative operators imported from ModelingToolkit.
- **`BeelerReuterRate`**: A plain function (not a component) implementing the generic voltage-dependent rate formula used by BR77 alpha/beta expressions.
- **Parameter ↔ Variable promotion**: MTK v11's `@component` + `System()` has no built-in toggle to promote a parameter to a variable or vice versa. Instead, declare the quantity as `@variables` in both components. To fix it as a constant, connect it via an algebraic equation (e.g. `x ~ fixed_value`); to leave it free, give it a differential equation or omit the algebraic constraint. MTK's compiler eliminates constant aliases at no runtime cost (see `SciML/ModelingToolkit.jl#4402`).

## Dependencies

- **ModelingToolkit v11** (only runtime dependency)
- **OrdinaryDiffEq** and **Test** (test-only)
- Julia 1.10+
