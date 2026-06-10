# ── KROME Fortran rate-expression → Julia closure ─────────────────────────────
# KROME rate strings are Fortran math in temperature variables (Tgas, Te, lnTe,
# …) and, for the richer networks, in the per-cell total density (`ntot` /
# `Hnuclei`) and in `@var:` intermediate variables (e.g. the kl/kh/ncr/a
# density-bridging helpers for collisional dissociation). We translate the Fortran
# surface syntax to Julia, prepend a preamble defining the temperature variables,
# the `@var:` assignments (in order), and `@common:` user variables (defaulted),
# then `eval` once at build time into a `(Tgas, ntot) -> k` closure.
#
# Reactions whose rate (or whose @var dependencies) reference symbols we don't
# provide — explicit per-species densities `n(idx_X)`, dust `Tdust` tables,
# `auto` reverse rates — return `nothing` from `compile_rate` so the parser skips
# them. Density dependence is detected automatically and flagged.

# base temperature variables (KROME's rate preamble), Julia source.
const _STD_PREAMBLE = quote
    Te=Tgas*8.617343e-5; lnTe=log(Te); invTe=1.0/Te; T=Tgas; invT=1.0/Tgas
    T32=Tgas/300.0; invsqrT=1.0/sqrt(Tgas); sqrTgas=sqrt(Tgas); sqrT=sqrt(Tgas)
    lnTgas=log(Tgas); logT=log10(Tgas); logTgas=log10(Tgas); logTe=log10(Te)
    Hnuclei=ntot
end

# math intrinsics + variables the compiled rate may reference for free.
const _MATH = Set([:exp,:log,:log10,:sqrt,:abs,:max,:min,:sin,:cos,:tanh,:atan,
                   :+,:-,:*,:/,:^,:float])
const _BASEVARS = Set([:Tgas,:ntot,:Hnuclei,:Te,:lnTe,:invTe,:T,:invT,:T32,
                       :invsqrT,:sqrTgas,:sqrT,:lnTgas,:logT,:logTgas,:logTe,
                       :Av,:av,:Tdust])

"""
    translate_fortran(s) -> String

Rewrite a Fortran rate expression to Julia source: `d`/`D` exponent literals
(`3.92d-13` → `3.92e-13`), `**` → `^`, the `dexp/dlog/dsqrt/…` intrinsics, and
`get_Hnuclei(n(:))` → `ntot`.
"""
function translate_fortran(s::AbstractString)
    t = String(strip(s))
    t = replace(t, r"get_Hnuclei\s*\(\s*n\(:\)\s*\)" => "ntot")
    for (a, b) in ("dexp"=>"exp","dlog10"=>"log10","dlog"=>"log","dsqrt"=>"sqrt",
                   "dabs"=>"abs","dmax1"=>"max","dmin1"=>"min","dble"=>"float")
        t = replace(t, a => b)
    end
    t = replace(t, "**" => "^")
    # Fortran allows a bare trailing decimal point (`1.-a` = `1.0 - a`, `1.)`),
    # which Julia mis-parses; complete it to `1.0` unless a digit/exponent follows.
    t = replace(t, r"(\d)\.(?![0-9dDeE])" => s"\1.0")
    t = replace(t, r"([0-9.])[dD]([+-]?[0-9])" => s"\1e\2")
    return t
end

function _symbols!(set, ex)
    if ex isa Symbol
        push!(set, ex)
    elseif ex isa Expr
        if ex.head === :call
            push!(set, ex.args[1]); for a in ex.args[2:end]; _symbols!(set, a); end
        else
            for a in ex.args; _symbols!(set, a); end
        end
    end
    return set
end

"""
    compile_rate(s; vars=Pair[], commons=Symbol[]) -> ((Tgas,ntot)->k, density_dep) or nothing

Compile a KROME rate string into a `(Tgas, ntot) -> rate` closure (cgs) plus a
`Bool` flag of whether it depends on density. `vars` is the ordered list of
`@var:` `name => expr` definitions; `commons` the `@common:` user variables
(injected as `0.0`). Returns `nothing` if the rate references unsupported symbols.
"""
function compile_rate(s::AbstractString; vars::Vector{<:Pair} = Pair{Symbol,String}[],
                      commons::Vector{Symbol} = Symbol[])
    jl = translate_fortran(s)
    rex = try Meta.parse(jl) catch; return nothing end

    # assemble the @var assignments (translated, in order) + collect their symbols
    varnames = Set{Symbol}(); varassigns = Expr[]; allsyms = Set{Symbol}()
    for (nm, expr) in vars
        vex = try Meta.parse(translate_fortran(expr)) catch; return nothing end
        push!(varassigns, Expr(:(=), nm, vex)); push!(varnames, nm)
        _symbols!(allsyms, vex)
    end
    _symbols!(allsyms, rex)

    allowed = union(_MATH, _BASEVARS, varnames, Set(commons))
    for sym in allsyms
        sym in allowed && continue
        return nothing
    end

    commonpre = Expr(:block, (Expr(:(=), c, 0.0) for c in commons)...,
                     Expr(:(=), :Av, 0.0), Expr(:(=), :av, 0.0), Expr(:(=), :Tdust, :Tgas))
    body = Expr(:block, _STD_PREAMBLE, commonpre, varassigns..., rex)
    raw  = eval(Expr(:->, Expr(:tuple, :Tgas, :ntot), body))
    f    = (Tgas, ntot) -> Base.invokelatest(raw, Tgas, ntot)

    dd = try
        a = f(500.0, 1.0); b = f(500.0, 1.0e8)
        !(isapprox(a, b; rtol = 1e-12) || (a == 0 && b == 0))
    catch
        false
    end
    return (f, dd)
end
