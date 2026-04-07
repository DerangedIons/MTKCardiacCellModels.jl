# Item 1: Fix unit inconsistency in constants

## Problem
`FARADAY = 96.5` is in C/mmol but `GAS_CONST = 8.315` is in J/(mol·K) — reviewer (@termi-official) flagged the mixed units.

## Analysis
The current values are **numerically correct** — the R/F ratio produces millivolts directly because R is in J/(mol·K) and F is in C/mmol, and the factors of 1000 cancel in RT/zF. Changing one constant without the other would break the model.

## Proposed fix
Replace the comment block (components.jl lines 6-9) to explain the intentional unit pairing:
```julia
# Physical constants scaled so that RT/(zF) yields millivolts directly,
# matching the mV convention used throughout (Beeler & Reuter 1977).
# R [J/(mol·K)] paired with F [C/mmol] — the 1e3 ratio cancels, giving mV.
const FARADAY = 96.5    # Faraday constant [C/mmol]  (≈ 96 485 C/mol ÷ 1000)
const GAS_CONST = 8.315 # Gas constant [J/(mol·K)]
```

## Bonus bugs found (GHK components, currently untested)

**GHK Potential (line 102):** Formula is inverted and missing `T`. Standard GHK voltage is `E = (RT/zF) ln(...)`, but code has `(zF/R) ln(...)`.
```julia
# Current (wrong):
Eₓ ~ (zₓ * FARADAY) / GAS_CONST * log(sum(Xₑs .* Pₓs) / sum(Xᵢs .* Pₓs))
# Fixed:
Eₓ ~ (GAS_CONST * T) / (zₓ * FARADAY) * log(sum(Xₑs .* Pₓs) / sum(Xᵢs .* Pₓs))
```

**GHK Current (line 137):** Missing `φₘ` multiplier. Standard GHK current is `I = P·z²F²V/(RT) · (...)`.
```julia
# Current (wrong) — missing φₘ after / (GAS_CONST * T):
Iₘ ~ Pₓ * (zₓ^2 * FARADAY^2) / (GAS_CONST * T) * (Xᵢ - ...)
# Fixed:
Iₘ ~ Pₓ * (zₓ^2 * FARADAY^2 * φₘ) / (GAS_CONST * T) * (Xᵢ - ...)
```
