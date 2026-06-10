# ── KROME Fortran rate-expression → Julia closure ─────────────────────────────
# KROME rate strings are Fortran math in a fixed set of temperature variables
# (Tgas in K, Te in eV, lnTe, invTe, …). We translate the Fortran surface syntax
# to Julia, prepend the KROME variable preamble, and `eval` once at build time
# into a `Tgas -> k` closure (host-side codegen, like `calc_rates` building a
# table). The closures are then tabulated by `to_generic`.
#
# Supported variables match KROME's `krome_constants` / rate preamble. Reactions
# whose rate references unsupported symbols (e.g. local densities `n(...)`, dust,
# user functions) are flagged by `compile_rate` returning `nothing` so the parser
# can skip them rather than miscompile.

# KROME temperature variables available to every rate expression.
const _KROME_PREAMBLE = quote
    Te      = Tgas * 8.617343e-5          # K → eV
    lnTe    = log(Te)
    invTe   = 1.0 / Te
    T32     = Tgas / 300.0
    invT    = 1.0 / Tgas
    invsqrT = 1.0 / sqrt(Tgas)
    sqrTgas = sqrt(Tgas)
    lnTgas  = log(Tgas)
    logTgas = log10(Tgas)
    logTe   = log10(Te)
end

# symbols the expression is allowed to reference (besides Tgas + the preamble vars
# + Base math). Anything else ⇒ unsupported (density-dependent / dust / user rate).
const _ALLOWED = Set([:Tgas, :Te, :lnTe, :invTe, :T32, :invT, :invsqrT, :sqrTgas,
                      :lnTgas, :logTgas, :logTe,
                      :exp, :log, :log10, :sqrt, :abs, :max, :min, :sin, :cos,
                      :tanh, :+, :-, :*, :/, :^, :Tgas])

"""
    translate_fortran(s) -> String

Rewrite a Fortran rate expression to Julia source: `d`/`D` exponent literals
(`3.92d-13` → `3.92e-13`, `1.d0` → `1.0e0`), `**` → `^`, and the Fortran intrinsics
`dexp/dlog/dlog10/dsqrt/dabs/dmax1/dmin1` → their Julia names.
"""
function translate_fortran(s::AbstractString)
    t = String(strip(s))
    # Fortran double-precision intrinsics → Julia
    for (a, b) in ("dexp" => "exp", "dlog10" => "log10", "dlog" => "log",
                   "dsqrt" => "sqrt", "dabs" => "abs", "dmax1" => "max",
                   "dmin1" => "min", "dble" => "float")
        t = replace(t, a => b)
    end
    t = replace(t, "**" => "^")
    # d/D exponent markers in numeric literals: <digit-or-dot> d <sign?digit> → e
    t = replace(t, r"([0-9.])[dD]([+-]?[0-9])" => s"\1e\2")
    return t
end

# collect the bare symbols referenced in a parsed expression
function _symbols!(set, ex)
    if ex isa Symbol
        push!(set, ex)
    elseif ex isa Expr
        if ex.head === :call
            push!(set, ex.args[1])
            for a in ex.args[2:end]; _symbols!(set, a); end
        else
            for a in ex.args; _symbols!(set, a); end
        end
    end
    return set
end

"""
    compile_rate(s) -> (Tgas -> k)::Function or nothing

Compile a KROME rate string into a `Tgas -> rate` closure (cgs). Returns `nothing`
when the expression references symbols outside the supported temperature set
(density/dust/user-dependent rates), so callers can skip those reactions cleanly.
"""
function compile_rate(s::AbstractString)
    jl = translate_fortran(s)
    ex = try
        Meta.parse(jl)
    catch
        return nothing
    end
    syms = _symbols!(Set{Symbol}(), ex)
    for sym in syms
        sym in _ALLOWED && continue
        return nothing            # unsupported reference → skip this reaction
    end
    body = Expr(:block, _KROME_PREAMBLE, ex)
    raw = eval(Expr(:->, :Tgas, body))   # host-side codegen, once per reaction
    # Wrap so calls dispatch through invokelatest: the closure is eval'd at a
    # newer world age than its callers, so direct calls would hit a world-age
    # error otherwise. The wrapper itself is an ordinary runtime closure.
    return Tgas -> Base.invokelatest(raw, Tgas)
end
