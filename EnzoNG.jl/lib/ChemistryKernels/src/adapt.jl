# ── Adapt rules ───────────────────────────────────────────────────────────────
# KernelAbstractions converts kernel arguments to the device with Adapt.jl. Plain
# NamedTuples adapt element-wise for free, but our struct-of-arrays containers
# need explicit rules so their inner arrays become device pointers on Metal. (On
# CPU these are identity, so the rules are harmless there.)

Adapt.adapt_structure(to, s::Species) = Species(Adapt.adapt(to, getfield(s, :fields)))

Adapt.adapt_structure(to, p::PhotoRates) = PhotoRates(
    Adapt.adapt(to, p.kphHI), Adapt.adapt(to, p.kphHeI), Adapt.adapt(to, p.kphHeII),
    Adapt.adapt(to, p.kdissH2), Adapt.adapt(to, p.photogamma))

function Adapt.adapt_structure(to, rt::RateTables{T}) where {T}
    a(v) = Adapt.adapt(to, v)
    return RateTables{T,typeof(a(rt.k1))}(rt.grid,
        a(rt.k1),a(rt.k2),a(rt.k3),a(rt.k4),a(rt.k5),a(rt.k6),a(rt.k7),a(rt.k8),
        a(rt.k9),a(rt.k10),a(rt.k11),a(rt.k12),a(rt.k13),a(rt.k14),a(rt.k15),
        a(rt.k16),a(rt.k17),a(rt.k18),a(rt.k19),a(rt.k22),
        a(rt.ceHI),a(rt.ceHeI),a(rt.ceHeII),a(rt.ciHI),a(rt.ciHeI),a(rt.ciHeII),
        a(rt.ciHeIS),a(rt.reHII),a(rt.reHeII1),a(rt.reHeII2),a(rt.reHeIII),a(rt.brem),
        a(rt.vibh),a(rt.hyd01k),a(rt.h2k01),a(rt.roth),a(rt.rotl),a(rt.gphdl),a(rt.gpldl),
        rt.compa,rt.comp_xray,rt.comp_temp,rt.gammaha,rt.piHI,rt.piHeI,rt.piHeII)
end
