#include "fortran.def"
#include "phys_const.def"
c=======================================================================
c  Standalone one-zone wrapper around Enzo's solve_rate_cool.F, exposing a
c  clean by-reference interface callable from Julia/C. Builds the rate
c  tables with calc_rates, then advances a single cell for `dt`. Primordial
c  only: no metals, no cloudy, no radiative transfer, no dust H2.
c=======================================================================
      subroutine enzo_chem_onezone(ispecies, nratec,
     &     d, ge, de, HI, HII, HeI, HeII, HeIII,
     &     HM, H2I, H2II, DI, DII, HDI,
     &     dt, dx, aye, redshift, temstart, temend, gamma, fh, dtoh,
     &     utem, uxyz, uaye, urho, utim, casebflag, ierr)
      implicit NONE
#include "fortran_types.def"
      INTG_PREC ispecies, nratec, casebflag, ierr
      R_PREC d, ge, de, HI, HII, HeI, HeII, HeIII
      R_PREC HM, H2I, H2II, DI, DII, HDI
      R_PREC dt, dx, aye, redshift, temstart, temend, gamma, fh, dtoh
      R_PREC utem, uxyz, uaye, urho, utim

c     local rate/cooling tables (allocated to nratec)
      INTG_PREC, parameter :: NDR = 2
      R_PREC ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec)
      R_PREC ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec)
      R_PREC ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec)
      R_PREC reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      R_PREC compa, gammahacgs, gammaha, piHI, piHeI, piHeII
      R_PREC hyd01ka(nratec), h2k01a(nratec), vibha(nratec)
      R_PREC rotha(nratec), rotla(nratec), gpldl(nratec), gphdl(nratec)
      R_PREC hdltea(nratec), hdlowa(nratec), hdcool(nratec)
      R_PREC ciecoa(nratec)
      R_PREC gaHIa(nratec), gaH2a(nratec), gaHea(nratec)
      R_PREC gaHpa(nratec), gaela(nratec), gasgra(nratec)
      R_PREC k1a(nratec), k2a(nratec), k3a(nratec), k4a(nratec)
      R_PREC k5a(nratec), k6a(nratec), k7a(nratec), k8a(nratec)
      R_PREC k9a(nratec), k10a(nratec), k11a(nratec), k12a(nratec)
      R_PREC k13a(nratec), k13dda(nratec,7), k14a(nratec), k15a(nratec)
      R_PREC k16a(nratec), k17a(nratec), k18a(nratec), k19a(nratec)
      R_PREC k20a(nratec), k21a(nratec), k22a(nratec), k23a(nratec)
      R_PREC k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec)
      R_PREC k54a(nratec), k55a(nratec), k56a(nratec)
      R_PREC k24, k25, k26, k27, k28, k29, k30, k31
      R_PREC h2dusta(nratec,NDR), ncrna(nratec)
      R_PREC ncrd1a(nratec), ncrd2a(nratec)
      R_PREC dtemstart, dtemend
      R_PREC metala(nratec,1)

c     1-cell field arrays + dummies for disabled physics
      R_PREC d1(1), e1(1), ge1(1), u1(1), v1(1), w1(1), de1(1)
      R_PREC HI1(1), HII1(1), HeI1(1), HeII1(1), HeIII1(1)
      R_PREC HM1(1), H2I1(1), H2II1(1), DI1(1), DII1(1), HDI1(1)
      R_PREC metal1(1), kphHI1(1), kphHeI1(1), kphHeII1(1)
      R_PREC kdissH2I1(1), photogamma1(1)
      R_PREC inutot(1), clPar1(1), clPar2(1), clPar3(1), clPar4(1)
      R_PREC clPar5(1), clCooling(1), clHeating(1)
      R_PREC comp_xraya, comp_temp, z_solar
      R_PREC alpha0, f3, eta1, eta2, dxa(1)
      R_PREC xe_start, xe_end, clEleFra
      R_PREC avgsighp, avgsighep, avgsighe2p, gammaha_use
      INTG_PREC iexpand, imethod, idual, imetal, imcool, idust, idim
      INTG_PREC is, js, ks, ie, je, ke, ih2co, ipiht, igammah
      INTG_PREC iradtype, nfreq, imetalregen, iradshield
      INTG_PREC iradtrans, iradcoupled, iradstep, irt_honly
      INTG_PREC ih2optical, iciecool, ithreebody, icmbTfloor, iClHeat
      INTG_PREC clGridRank, clGridDim(5), clDataSize, n_xe
      INTG_PREC threebody, ioutput, ndratec
      LOGIC_PREC casebrates

      casebrates = (casebflag .ne. 0)
      threebody = 0
      ioutput   = 0
      ndratec   = NDR
      dtemstart = 1.0d0
      dtemend   = 1.0d4
      alpha0 = 0._RKIND
      f3     = 0._RKIND
      iradtype = 0

c     --- build the rate tables (Enzo calc_rates.F) ---
      call calc_rates(nratec, aye, temstart, temend, alpha0, f3,
     &     iradtype, casebrates, threebody,
     &     utem, uxyz, uaye, urho, utim,
     &     ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa,
     &     ciHeISa, ciHeIIa, reHIIa, reHeII1a,
     &     reHeII2a, reHeIIIa, brema, compa,
     &     gammahacgs, gammaha,
     &     piHI, piHeI, piHeII,
     &     hyd01ka, h2k01a, vibha, rotha, rotla,
     &     gpldl, gphdl, hdltea, hdlowa, hdcool, ciecoa,
     &     gaHIa, gaH2a, gaHea, gaHpa, gaela, gasgra,
     &     k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &     k11a, k12a, k13a, k13dda, k14a, k15a, k16a, k17a,
     &     k18a, k19a, k20a, k21a, k22a, k23a,
     &     k24, k25, k26, k27, k28, k29, k30, k31,
     &     k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &     ndratec, dtemstart, dtemend, h2dusta,
     &     ncrna, ncrd1a, ncrd2a, ioutput)

c     --- stage the single cell (zero velocities => e == ge) ---
      d1(1)=d; ge1(1)=ge; e1(1)=ge; u1(1)=0._RKIND
      v1(1)=0._RKIND; w1(1)=0._RKIND; de1(1)=de
      HI1(1)=HI; HII1(1)=HII; HeI1(1)=HeI; HeII1(1)=HeII
      HeIII1(1)=HeIII; HM1(1)=HM; H2I1(1)=H2I; H2II1(1)=H2II
      DI1(1)=DI; DII1(1)=DII; HDI1(1)=HDI
      metal1(1)=0._RKIND; kphHI1(1)=0._RKIND; kphHeI1(1)=0._RKIND
      kphHeII1(1)=0._RKIND; kdissH2I1(1)=0._RKIND; photogamma1(1)=0._RKIND
      inutot(1)=0._RKIND; metala(:,1)=0._RKIND
      clPar1(1)=0._RKIND; clPar2(1)=0._RKIND; clPar3(1)=0._RKIND
      clPar4(1)=0._RKIND; clPar5(1)=0._RKIND
      clCooling(1)=0._RKIND; clHeating(1)=0._RKIND
      dxa(1)=dx

      iexpand=0; imethod=2; idual=1; imetal=0; imcool=0; idust=0
      idim=1; is=0; js=0; ks=0; ie=0; je=0; ke=0
      ih2co=1; ipiht=1; igammah=0
      nfreq=1; imetalregen=0; iradshield=0
      iradtrans=0; iradcoupled=0; iradstep=0; irt_honly=0
      ih2optical=0; iciecool=0; ithreebody=0; icmbTfloor=0; iClHeat=0
      clGridRank=1; clGridDim(1)=1; clGridDim(2)=1; clGridDim(3)=1
      clGridDim(4)=1; clGridDim(5)=1; clDataSize=1; n_xe=1
      comp_xraya=0._RKIND; comp_temp=0._RKIND; z_solar=0.02_RKIND
      eta1=1.0e-5_RKIND; eta2=1.0e-5_RKIND
      xe_start=0._RKIND; xe_end=0._RKIND; clEleFra=0._RKIND
      avgsighp=0._RKIND; avgsighep=0._RKIND; avgsighe2p=0._RKIND
      ierr=0

      call solve_rate_cool(d1, e1, ge1, u1, v1, w1, de1,
     &     HI1, HII1, HeI1, HeII1, HeIII1,
     &     1_IKIND, 1_IKIND, 1_IKIND, nratec, iexpand, imethod,
     &     idual, ispecies, imetal, imcool, idust, idim,
     &     is, js, ks, ie, je, ke, ih2co, ipiht, igammah,
     &     dxa, dt, aye, redshift, temstart, temend,
     &     utem, uxyz, uaye, urho, utim,
     &     eta1, eta2, gamma, fh, dtoh, z_solar,
     &     k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &     k11a, k12a, k13a, k13dda, k14a, k15a,
     &     k16a, k17a, k18a, k19a, k22a,
     &     k24, k25, k26, k27, k28, k29, k30, k31,
     &     k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &     ndratec, dtemstart, dtemend, h2dusta,
     &     ncrna, ncrd1a, ncrd2a,
     &     ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa,
     &     ciHeISa, ciHeIIa, reHIIa, reHeII1a,
     &     reHeII2a, reHeIIIa, brema, compa, gammaha,
     &     comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &     HM1, H2I1, H2II1, DI1, DII1, HDI1, metal1,
     &     hyd01ka, h2k01a, vibha, rotha, rotla,
     &     gpldl, gphdl, hdltea, hdlowa,
     &     gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &     gasgra, metala, n_xe, xe_start, xe_end,
     &     inutot, iradtype, nfreq, imetalregen,
     &     iradshield, avgsighp, avgsighep, avgsighe2p,
     &     iradtrans, iradcoupled, iradstep, ierr,
     &     irt_honly,
     &     kphHI1, kphHeI1, kphHeII1, kdissH2I1,
     &     photogamma1,
     &     ih2optical, iciecool, ithreebody, ciecoa,
     &     icmbTfloor, iClHeat,
     &     clEleFra, clGridRank, clGridDim,
     &     clPar1, clPar2, clPar3, clPar4, clPar5,
     &     clDataSize, clCooling, clHeating)

c     --- copy results back ---
      ge=ge1(1); de=de1(1)
      HI=HI1(1); HII=HII1(1); HeI=HeI1(1); HeII=HeII1(1); HeIII=HeIII1(1)
      HM=HM1(1); H2I=H2I1(1); H2II=H2II1(1)
      DI=DI1(1); DII=DII1(1); HDI=HDI1(1)
      return
      end

c     error stubs referenced by ERROR_MESSAGE in the Enzo sources
      subroutine f_error(filename, line)
      implicit none
      character*(*) filename
      integer*8 line
      write(0,*) 'enzo f_error at ', filename, line
      return
      end
      subroutine f_warning(filename, line)
      implicit none
      character*(*) filename
      integer*8 line
      return
      end
