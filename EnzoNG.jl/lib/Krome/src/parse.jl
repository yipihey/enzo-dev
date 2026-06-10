# ── KROME react_* file parser ─────────────────────────────────────────────────
# Reads a KROME network file into a backend-agnostic IR. The classic KROME row is
# comma-separated:
#
#     idx, R1, R2, R3, P1, P2, P3, P4, Tmin, Tmax, rate
#
# (empty reactant/product slots are blank). A `@format:` directive overrides the
# column layout for the lines that follow (modern KROME style); we honor it and
# fall back to the 11-column classic layout otherwise. Photon/cosmic-ray pseudo-
# species (`g`, `CR`, …) are dropped from the stoichiometry. Per-line temperature
# limits gate the (piecewise) rate; lines that share reactants/products with
# disjoint T-ranges sum back to the full piecewise rate in the executor.

# tokens that are NOT tracked species (photons, cosmic rays, blanks)
const _DUMMY = Set(["", "G", "GAMMA", "Γ", "CR", "CRP", "CRPHOT", "PHOTON", "DUST"])

"""
One reaction: real-species reactant/product indices + a compiled `(Tgas,ntot)→k`
rate closure and a `density_dep` flag (true if the rate varies with `ntot`).
"""
struct KromeReaction
    reactants::Vector{Int}
    products::Vector{Int}
    rate::Function
    density_dep::Bool
    raw::String
end

"""
    KromeNetwork

A parsed KROME network: the ordered `species` list (first-appearance order), the
name→index `index` map, and the `reactions`. Lower to a runnable solver with
[`to_generic`](@ref).
"""
struct KromeNetwork
    species::Vector{String}
    index::Dict{String,Int}
    reactions::Vector{KromeReaction}
    skipped::Int          # reactions dropped (unsupported density/dust/user rate)
end

const _DEFAULT_FORMAT = ["idx", "R", "R", "R", "P", "P", "P", "P", "Tmin", "Tmax", "rate"]

# parse a single Tmin/Tmax constraint field into a (lo, hi) bound contribution
function _parse_limit(field::AbstractString, role::Symbol)
    f = strip(field)
    (isempty(f) || uppercase(f) == "NONE") && return (-Inf, Inf)
    for (op, kind) in (".GE." => :lo, ".GT." => :lo, ".LE." => :hi, ".LT." => :hi,
                       ">=" => :lo, "<=" => :hi, ">" => :lo, "<" => :hi)
        if startswith(f, op)
            v = tryparse(Float64, replace(strip(f[length(op)+1:end]), "d" => "e", "D" => "e"))
            v === nothing && return (-Inf, Inf)
            return kind === :lo ? (v, Inf) : (-Inf, v)
        end
    end
    v = tryparse(Float64, replace(f, "d" => "e", "D" => "e"))
    v === nothing && return (-Inf, Inf)
    return role === :Tmin ? (v, Inf) : (-Inf, v)   # bare number: bound by column role
end

"""
    parse_krome(path) -> KromeNetwork

Parse a KROME `react_*` file at `path`. Reactions whose rate expression references
unsupported symbols (local densities, dust, user variables) are skipped and
counted in `KromeNetwork.skipped`; the rest become runnable mass-action reactions.
"""
function parse_krome(path::AbstractString)
    species = String[]; index = Dict{String,Int}()
    reactions = KromeReaction[]; skipped = 0
    fmt = copy(_DEFAULT_FORMAT)
    vars = Pair{Symbol,String}[]      # accumulated @var: name => expr (in order)
    commons = Symbol[]                # @common: user variables

    sid(tok) = begin
        s = uppercase(strip(tok))
        (s in _DUMMY) && return 0
        haskey(index, s) && return index[s]
        push!(species, s); index[s] = length(species); return index[s]
    end

    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        if startswith(line, "@format:")
            fmt = strip.(split(line[length("@format:")+1:end], ","))
            continue
        elseif startswith(line, "@var:")
            body = strip(line[length("@var:")+1:end])
            eq = findfirst('=', body)
            if eq !== nothing
                nm = Symbol(strip(body[1:eq-1])); ex = strip(body[eq+1:end])
                # later redefinitions of the same name win (drop the earlier)
                filter!(p -> p.first != nm, vars); push!(vars, nm => String(ex))
            end
            continue
        elseif startswith(line, "@common:")
            for c in strip.(split(line[length("@common:")+1:end], ","))
                isempty(c) || push!(commons, Symbol(c))
            end
            continue
        end
        startswith(line, "@") && continue          # @noTabNext, @ev, … : ignored
        cols = strip.(split(line, ","))
        length(cols) < length(fmt) && continue      # malformed / continuation
        reac = Int[]; prod = Int[]; tmin = "NONE"; tmax = "NONE"; rate = ""
        for (k, role) in enumerate(fmt)
            k > length(cols) && break
            tok = cols[k]
            if role == "R"
                s = sid(tok); s != 0 && push!(reac, s)
            elseif role == "P"
                s = sid(tok); s != 0 && push!(prod, s)
            elseif role == "Tmin"; tmin = tok
            elseif role == "Tmax"; tmax = tok
            elseif role == "rate"; rate = join(cols[k:end], ",")  # rate may contain commas
                break
            end
        end
        isempty(rate) && (skipped += 1; continue)
        compiled = compile_rate(rate; vars = vars, commons = commons)
        compiled === nothing && (skipped += 1; continue)
        base, dd = compiled
        lo1, hi1 = _parse_limit(tmin, :Tmin)
        lo2, hi2 = _parse_limit(tmax, :Tmax)
        lo = max(lo1, lo2); hi = min(hi1, hi2)
        gated = (lo == -Inf && hi == Inf) ? base :
                ((Tgas, ntot) -> (lo ≤ Tgas ≤ hi) ? base(Tgas, ntot) : 0.0)
        push!(reactions, KromeReaction(reac, prod, gated, dd, rate))
    end
    return KromeNetwork(species, index, reactions, skipped)
end
