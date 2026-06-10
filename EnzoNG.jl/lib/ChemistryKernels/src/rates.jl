# ── Rate + cooling coefficient tables (port of calc_rates.F) ──────────────────
# calc_rates builds, ONCE per call (and only re-built when the radiation spectrum
# changes), log(T)-tabulated coefficients. It is host-side table construction --
# not a per-cell GPU kernel -- so it stays a plain precision-generic Julia
# function; the resulting tables are uploaded to the device with `to_device`.
#
# The collisional set is the canonical Abel, Anninos, Zhang & Norman (1997)
# numbering (k1..k19, k22); recombination uses Cen (1992) case-A or the Hui &
# Gnedin (1997)-style case-B fits gated by `casebrates`. Cooling coefficients are
# Black (1981)/Cen (1992). Every coefficient is divided by its code unit exactly
# as the Fortran does (`/kunit`, `/coolunit`) so the tables drop into the kernels
# unchanged.

"""
    RateTables{T}

All log(T)-tabulated coefficients consumed by [`cool1d_multi!`](@ref) and
[`solve_rate_cool!`](@ref), on a shared [`TempGrid`](@ref). Two-body rate
coefficients (`k1..k19`, `k22`) are in code units; cooling coefficients
(`ceHI`, …, `brem`) are in cooling-code units. Scalar photo/Compton terms
(`piHI`, `compa`, `gammaha`) ride along. Build with [`build_rate_tables`](@ref).

Field arrays are length-`nratec` and may be moved to a device en masse with
[`rate_tables_to_device`](@ref).
"""
struct RateTables{T,V<:AbstractVector{T}}
    grid::TempGrid{T}
    # collisional / recombination / H2 formation+destruction (Abel+97)
    k1::V; k2::V; k3::V; k4::V; k5::V; k6::V; k7::V; k8::V; k9::V; k10::V
    k11::V; k12::V; k13::V; k14::V; k15::V; k16::V; k17::V; k18::V; k19::V; k22::V
    # cooling coefficients (Black81 / Cen92)
    ceHI::V; ceHeI::V; ceHeII::V
    ciHI::V; ciHeI::V; ciHeII::V; ciHeIS::V
    reHII::V; reHeII1::V; reHeII2::V; reHeIII::V
    brem::V
    # H2 cooling (Lepp & Shull / Galli & Palla low-density limit pieces)
    vibh::V; hyd01k::V; h2k01::V; roth::V; rotl::V; gphdl::V; gpldl::V
    # scalar terms
    compa::T; comp_xray::T; comp_temp::T; gammaha::T
    piHI::T; piHeI::T; piHeII::T
end

"`tiny`-floor used by calc_rates for unpopulated table entries."
const RATE_TINY = 1e-20

"""
    build_rate_tables(T = Float64; nratec=400, temstart=1.0, temend=1e9,
                      units::ChemistryUnits, casebrates=false,
                      comp_xray=0.0, comp_temp=0.0, gammaha=0.0,
                      piHI=0.0, piHeI=0.0, piHeII=0.0) -> RateTables{T}

Construct the rate/cooling tables (port of `calc_rates.F`). `units` supplies
`kunit`/`coolunit`. `casebrates=true` selects case-B recombination for k2/k4/k6.
The optional scalars are the UV-background photo-ionization (`piHI`…) and
photoelectric-dust (`gammaha`) terms calc_rates fills from the spectrum.
"""
function build_rate_tables(::Type{T} = Float64;
                           nratec::Integer = 400, temstart = 1.0, temend = 1.0e9,
                           units::ChemistryUnits, casebrates::Bool = false,
                           comp_xray = 0.0, comp_temp = 0.0, gammaha = 0.0,
                           piHI = 0.0, piHeI = 0.0, piHeII = 0.0) where {T}
    grid = TempGrid(T; nratec = nratec, temstart = temstart, temend = temend)
    kunit    = Float64(units.kunit)
    coolunit = Float64(units.coolunit)
    tevk     = PhysConst.tevk
    kboltz   = PhysConst.kboltz
    dhuge    = 1.0e30
    tiny     = RATE_TINY

    mk() = Vector{T}(undef, nratec)
    k1=mk();k2=mk();k3=mk();k4=mk();k5=mk();k6=mk();k7=mk();k8=mk();k9=mk();k10=mk()
    k11=mk();k12=mk();k13=mk();k14=mk();k15=mk();k16=mk();k17=mk();k18=mk();k19=mk();k22=mk()
    ceHI=mk();ceHeI=mk();ceHeII=mk()
    ciHI=mk();ciHeI=mk();ciHeII=mk();ciHeIS=mk()
    reHII=mk();reHeII1=mk();reHeII2=mk();reHeIII=mk();brem=mk()
    vibh=mk();hyd01k=mk();h2k01=mk();roth=mk();rotl=mk();gphdl=mk();gpldl=mk()

    @inbounds for i in 1:nratec
        logttt = grid.logtem0 + (i - 1) * grid.dlogtem
        ttt    = exp(logttt)
        tev    = ttt / tevk
        logtev = log(tev)

        # ── collisional ionization / recombination (Abel+97 polynomial fits) ──
        if tev > 0.8
            k1[i] = exp(-32.71396786375 + 13.53655609057logtev
                - 5.739328757388logtev^2 + 1.563154982022logtev^3
                - 0.2877056004391logtev^4 + 0.03482559773736999logtev^5
                - 0.00263197617559logtev^6 + 0.0001119543953861logtev^7
                - 2.039149852002e-6*logtev^8) / kunit
            k3[i] = exp(-44.09864886561001 + 23.91596563469logtev
                - 10.75323019821logtev^2 + 3.058038757198logtev^3
                - 0.5685118909884001logtev^4 + 0.06795391233790001logtev^5
                - 0.005009056101857001logtev^6 + 0.0002067236157507logtev^7
                - 3.649161410833e-6*logtev^8) / kunit
            k4[i] = (1.54e-9*(1.0 + 0.3/exp(8.099328789667/tev))
                / (exp(40.49664394833662/tev)*tev^1.5) + 3.92e-13/tev^0.6353) / kunit
            k5[i] = exp(-68.71040990212001 + 43.93347632635logtev
                - 18.48066993568logtev^2 + 4.701626486759002logtev^3
                - 0.7692466334492logtev^4 + 0.08113042097303logtev^5
                - 0.005324020628287001logtev^6 + 0.0001975705312221logtev^7
                - 3.165581065665e-6*logtev^8) / kunit
        else
            k1[i] = tiny; k3[i] = tiny
            k4[i] = 3.92e-13/tev^0.6353 / kunit
            k5[i] = tiny
        end
        if casebrates
            k4[i] = 1.26e-14*(5.7067e5/ttt)^0.75 / kunit
        end

        # HII recombination (k2): case-A polynomial or case-B fit
        if casebrates
            k2[i] = ttt < 1.0e9 ?
                4.881357e-6*ttt^(-1.5)*(1.0 + 1.14813e2*ttt^(-0.407))^(-2.242)/kunit :
                tiny
        else
            k2[i] = ttt > 5500.0 ?
                exp(-28.61303380689232 - 0.7241125657826851logtev
                    - 0.02026044731984691logtev^2 - 0.002380861877349834logtev^3
                    - 0.0003212605213188796logtev^4 - 0.00001421502914054107logtev^5
                    + 4.989108920299513e-6*logtev^6 + 5.755614137575758e-7*logtev^7
                    - 1.856767039775261e-8*logtev^8 - 3.071135243196595e-9*logtev^9)/kunit :
                k4[i]
        end
        # HeIII recombination (k6)
        if casebrates
            k6[i] = ttt < 1.0e9 ?
                7.8155e-5*ttt^(-1.5)*(1.0 + 2.0189e2*ttt^(-0.407))^(-2.242)/kunit : tiny
        else
            k6[i] = 3.36e-10/sqrt(ttt)/(ttt/1.0e3)^0.2/(1.0 + (ttt/1.0e6)^0.7)/kunit
        end

        # ── H2 / H- formation + destruction (Abel+97) ─────────────────────────
        k7[i] = 6.77e-15*tev^0.8779 / kunit
        k8[i] = tev > 0.1 ?
            exp(-20.06913897587003 + 0.2289800603272916logtev
                + 0.03599837721023835logtev^2 - 0.004555120027032095logtev^3
                - 0.0003105115447124016logtev^4 + 0.0001073294010367247logtev^5
                - 8.36671960467864e-6*logtev^6 + 2.238306228891639e-7*logtev^7)/kunit :
            1.43e-9/kunit
        k9[i] = ttt > 6.7e3 ?
            5.81e-16*(ttt/56200.0)^(-0.6657*log10(ttt/56200.0))/kunit :
            1.85e-23*ttt^1.8/kunit
        k10[i] = 6.0e-10/kunit
        if tev > 0.3
            k13[i] = 1.0670825e-10*tev^2.012/(exp(4.463/tev)*(1.0 + 0.2472tev)^3.512)/kunit
            k11[i] = exp(-24.24914687731536 + 3.400824447095291logtev
                - 3.898003964650152logtev^2 + 2.045587822403071logtev^3
                - 0.5416182856220388logtev^4 + 0.0841077503763412logtev^5
                - 0.007879026154483455logtev^6 + 0.0004138398421504563logtev^7
                - 9.36345888928611e-6*logtev^8)/kunit
            k12[i] = 4.38e-10*exp(-102000.0/ttt)*ttt^0.35/kunit
        else
            k13[i] = tiny; k11[i] = tiny; k12[i] = tiny
        end
        k14[i] = tev > 0.04 ?
            exp(-18.01849334273 + 2.360852208681logtev - 0.2827443061704logtev^2
                + 0.01623316639567logtev^3 - 0.03365012031362999logtev^4
                + 0.01178329782711logtev^5 - 0.001656194699504logtev^6
                + 0.0001068275202678logtev^7 - 2.631285809207e-6*logtev^8)/kunit :
            tiny
        k15[i] = tev > 0.1 ?
            exp(-20.37260896533324 + 1.139449335841631logtev
                - 0.1421013521554148logtev^2 + 0.00846445538663logtev^3
                - 0.0014327641212992logtev^4 + 0.0002012250284791logtev^5
                + 0.0000866396324309logtev^6 - 0.00002585009680264logtev^7
                + 2.4555011970392e-6*logtev^8 - 8.06838246118e-8*logtev^9)/kunit :
            2.56e-9*tev^1.78186/kunit
        k16[i] = 6.5e-9/sqrt(tev)/kunit
        k17[i] = ttt > 1.0e4 ?
            4.0e-4*ttt^(-1.4)*exp(-15100.0/ttt)/kunit :
            1.0e-8*ttt^(-0.4)/kunit
        k18[i] = ttt > 617.0 ? 1.32e-6*ttt^(-0.76)/kunit : 1.0e-8/kunit
        k19[i] = 5.0e-7*sqrt(100.0/ttt)/kunit
        # 3-body H2 formation (Abel+02 low-density coefficient); k22 reused below
        k22[i] = 1.3e-32*(ttt/300.0)^(-0.38)/kunit

        # ── cooling coefficients (Black 1981 / Cen 1992) ──────────────────────
        ceHI[i]   = 7.5e-19*exp(-min(log(dhuge),118348.0/ttt))/(1.0+sqrt(ttt/1.0e5))/coolunit
        ceHeI[i]  = 9.1e-27*exp(-min(log(dhuge),13179.0/ttt))*ttt^(-0.1687)/(1.0+sqrt(ttt/1.0e5))/coolunit
        ceHeII[i] = 5.54e-17*exp(-min(log(dhuge),473638.0/ttt))*ttt^(-0.397)/(1.0+sqrt(ttt/1.0e5))/coolunit
        ciHeIS[i] = 5.01e-27*ttt^(-0.1687)/(1.0+sqrt(ttt/1.0e5))*exp(-min(log(dhuge),55338.0/ttt))/coolunit
        # collisional ionization cooling tied to the rate coefficients (Abel)
        ciHI[i]   = 2.18e-11*k1[i]*kunit/coolunit
        ciHeI[i]  = 3.94e-11*k3[i]*kunit/coolunit
        ciHeII[i] = 8.72e-11*k5[i]*kunit/coolunit
        reHII[i]   = 8.70e-27*sqrt(ttt)*(ttt/1000.0)^(-0.2)/(1.0+(ttt/1.0e6)^0.7)/coolunit
        reHeII1[i] = 1.55e-26*ttt^0.3647/coolunit
        reHeII2[i] = 1.24e-13*ttt^(-1.5)*exp(-min(log(dhuge),470000.0/ttt))*
                     (1.0+0.3*exp(-min(log(dhuge),94000.0/ttt)))/coolunit
        reHeIII[i] = 3.48e-26*sqrt(ttt)*(ttt/1000.0)^(-0.2)/(1.0+(ttt/1.0e6)^0.7)/coolunit
        brem[i]    = 1.43e-27*sqrt(ttt)*(1.1+0.34*exp(-(5.5-log10(ttt))^2/3.0))/coolunit

        # ── H2 cooling (Lepp & Shull low/high-density pieces) ─────────────────
        xx = log10(ttt/1.0e4)
        vibh[i] = 1.1e-18*exp(-min(log(dhuge),6744.0/ttt))/coolunit
        dum = ttt > 1635.0 ? 1.0e-12*sqrt(ttt)*exp(-1000.0/ttt) :
                             1.4e-13*exp((ttt/125.0)-(ttt/577.0)^2)
        hyd01k[i] = dum*exp(-min(log(dhuge), 8.152e-13/(kboltz*ttt)))/coolunit
        dum2 = 8.152e-13*(4.2/(kboltz*(ttt+1190.0)) + 1.0/(kboltz*ttt))
        h2k01[i] = 1.45e-12*sqrt(ttt)*exp(-min(log(dhuge),dum2))/coolunit
        rotl[i] = ttt > 4031.0 ? 1.38e-22*exp(-9243.0/ttt)/coolunit :
                                 10.0^(-22.9 - 0.553xx - 1.148xx^2)/coolunit
        roth[i] = ttt > 1087.0 ? 3.90e-19*exp(-6118.0/ttt)/coolunit :
                                 10.0^(-19.24 + 0.474xx - 1.247xx^2)/coolunit
        # Galli & Palla (1998) low-density H2 cooling, log10 polynomial in log10(T)
        tm = max(min(ttt, 1.0e4), 13.0); lt3 = log10(tm/1.0e3)
        gpldl[i] = 10.0^(-103.0 + 97.59lt3 - 48.05lt3^2 + 10.80lt3^3 - 0.9032lt3^4)/coolunit
        gphdl[i] = roth[i]  # high-density limit shares the rotational-H piece
    end

    return RateTables{T,Vector{T}}(grid,
        k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k22,
        ceHI,ceHeI,ceHeII,ciHI,ciHeI,ciHeII,ciHeIS,reHII,reHeII1,reHeII2,reHeIII,brem,
        vibh,hyd01k,h2k01,roth,rotl,gphdl,gpldl,
        T(2.873e-73/Float64(units.coolunit)),  # compa (Compton, Peebles 1971)·/coolunit
        T(comp_xray), T(comp_temp), T(gammaha),
        T(piHI), T(piHeI), T(piHeII))
end

"""
    rate_tables_to_device(be, rt::RateTables) -> RateTables

Move every table array of `rt` onto backend `be` (scalars unchanged). Lets the
per-cell kernels read coefficients from device memory.
"""
function rate_tables_to_device(be, rt::RateTables{T}) where {T}
    d(v) = to_device(be, v)
    return RateTables{T,typeof(d(rt.k1))}(rt.grid,
        d(rt.k1),d(rt.k2),d(rt.k3),d(rt.k4),d(rt.k5),d(rt.k6),d(rt.k7),d(rt.k8),
        d(rt.k9),d(rt.k10),d(rt.k11),d(rt.k12),d(rt.k13),d(rt.k14),d(rt.k15),
        d(rt.k16),d(rt.k17),d(rt.k18),d(rt.k19),d(rt.k22),
        d(rt.ceHI),d(rt.ceHeI),d(rt.ceHeII),d(rt.ciHI),d(rt.ciHeI),d(rt.ciHeII),
        d(rt.ciHeIS),d(rt.reHII),d(rt.reHeII1),d(rt.reHeII2),d(rt.reHeIII),d(rt.brem),
        d(rt.vibh),d(rt.hyd01k),d(rt.h2k01),d(rt.roth),d(rt.rotl),d(rt.gphdl),d(rt.gpldl),
        rt.compa,rt.comp_xray,rt.comp_temp,rt.gammaha,rt.piHI,rt.piHeI,rt.piHeII)
end
