#define WRF
#define HAIL
#define DOUBLE
#define PARTIAL

#undef WDM7

#define ICE_SPHERE
#undef NI0_HDC2004
#define NEW_NCRMIN
#define HEYMS_VTGH
#undef LH2010_PRAUT
#undef HDC2004_PRAUT
#define LB2017_PRAUT

#undef EVAP_PARTIAL
#define NEW_N0
#define WET_BULB
#undef LAMDAR_1M
#undef WDM7_LAMG
#undef WDM7_LAMR
#undef WDM7_LAMC

#define SEMI_RAIN
#define SEDI_MEAN
#undef SEDI_TEST
#define SEDI_VT_AVERAGE
#define SEDI_SHAPE
#undef SEDI_VT_BOTTOM
#define SEDI_UPDATE
#undef PVTRN_ONEHALF

#define NCACT_LEI2019
#define WYSER_REI
#undef RYAN2000_qimax
#define SIZE_LIMIT1
#undef SIZE_LIMIT2
#undef SIZE_LIMIT3
#define SIZE_LIMIT_SEDI
#define PRTOC
#define WEAK_ACCRETION
#undef PRT
#undef DEBUG
#ifndef DOUBLE
#undef PRT
#define SEMI_RAIN
#endif

#ifdef WDM7
#undef PARTIAL
#undef SEMI_RAIN
#undef SEDI_MEAN
#undef SEDI_SHAPE
#undef HEYMS_VTGH
#undef ICE_SPHERE
#define NI0_HDC2004
#undef WEAK_ACCRETION
#undef NEW_N0
#undef WET_BULB
#define HDC2004_PRAUT
#undef LB2017_PRAUT
#undef NEW_NCRMIN
#undef  NCACT_LEI2019
#undef SIZE_LIMIT1
#define WDM7_LAMG
#define WDM7_LAMR
#define WDM7_LAMC
#endif
!
#define MP_INIT_ mp_udm_init
#define MP_RUN_ mp_udm_run
#define MP_FINALIZE_ mp_udm_finalize
#define MP_ARG_INIT_ arg_table_mp_udm_init
#define MP_ARG_RUN_ arg_table_mp_udm_run
#define MODULE_NAME_ mp_udm
!>\file MODULE_NAME_.F90
!!
!! This file is the MODULE_NAME_ micorphysics shceme.
!! author Songyou Hong [hong@ucar.edu, songyouhong@gmail.com]
!!
  module MODULE_NAME_
   use module_wrf_error
   use module_mp_radar
   use ccpp_kind_types,only: kind_phys
#define IMDIM_ ncol
#define ITDIM_ ncol
#define KMDIM_ nlev
#define KTDIM_ nlev
   implicit none
   public :: MP_INIT_, MP_RUN_, MP_FINALIZE_
   public :: udm_mp_reflectivity, udm_mp_effective_radius
   logical :: if_aerosol_activated
!
!-------------------------------------------------------------------------------
! parameters for microphysical processes
!
   real(kind=kind_phys), parameter, private :: &
          dtcldcr  = 180.  ,&    !< maximum time step for sub-cycling loops  [s]
          n0r    = 8.e6    ,&    !< intercept parameter for rain           [m-4]
#ifdef NEW_N0
          n0s    = 1.e6    ,&    !< intercept parameter for snow at 0c     [m-4]
#else
          n0s    = 2.e6    ,&    !< intercept parameter for snow at 0c     [m-4]
#endif
#ifdef NEW_N0
          n0g    = 1.e6    ,&    !< intercept parameter for graupel        [m-4]
#else
          n0g    = 4.e6    ,&    !< intercept parameter for graupel        [m-4]
#endif
#ifdef HAIL
#ifdef NEW_N0
          n0h    = 1.e4    ,&    !< intercept parameter for hail           [m-4]
#else
          n0h    = 4.e4    ,&    !< intercept parameter for hail           [m-4]
#endif
#endif
#ifdef NEW_N0
          n0smax = 1.e9    ,&    !< maximum n0s                            [m-4]
          n0smin = 1.e5    ,&    !< minimum n0s                            [m-4]
          n0gmax = 1.e7    ,&    !< maximum n0g                            [m-4]
          n0gmin = 1.e4    ,&    !< minimum n0g                            [m-4]
          n0hmax = 1.e5    ,&    !< maximum n0h                            [m-4]
          n0hmin = 1.e3    ,&    !< minimum n0h                            [m-4]
#else
          n0smax = 1.e11   ,&    !< maximum n0s (t=-90c unlimited)         [m-4]
#endif
          mi0    = 1.e-12  ,&    !< initial mass of cloud ice               [kg]
#ifdef NEW_N0
          alpha_qs = .12   ,&    !< exponent factor for n0s computation
          alpha_qg = .08   ,&    !< exponent factor for n0g computation
          alpha_qh = .04   ,&    !< exponent factor for n0h computation
          beta   = .0      ,&    !< factor for mass in n0 (=0, no dependency)
          q_ref  =  1.e-4  ,&    !< a reference value for mass in n0    [kg-1kg]
#else
          alpha  = .12     ,&    !< .122 exponent factor for n0s computation
#endif
          avti   = 1.49e4  ,&    !< a constant for vt of hexagonal ice
          bvti   = 1.31    ,&    !< a constant for vt of hexagonal ice
          cxmi   = (1./11.9)**2,&!< a constant for mass of hexagonal ice
          dxmi   = 2.0     ,&    !< a constant for mass of hexagonal ice
#ifdef ICE_SPHERE
          pi_    = 3.141593,&    !< pi only for ice mass computation
          deni   = 500.    ,&    !< density of cloud ice                 [kgm-3] 
          avtis  = 2.71e3  ,&    !< a constant for vt of spherical ice
          bvtis  = 1.0     ,&    !< a constant for vt of spherical ice
          cxmis  = pi_*deni/6.,& !< a constant for mass of spherical ice
          dxmis  = 3.0     ,&    !< a constant for mass of spherical ice
#endif
          avtr   = 841.9   ,&    !< a constant for terminal velocity of rain
          bvtr   = 0.8     ,&    !< a constant for terminal velocity of rain
          avts   = 11.72   ,&    !< a constant for terminal velocity of snow
          bvts   = .41     ,&    !< a constant for terminal velocity of snow
#ifdef HEYMS_VTGH
          avtg   = 553.1   ,&    !< a constant for terminal velocity of graupel
          bvtg   = 0.97    ,&    !< a constant for terminal velocity of graupel
#else
          avtg   = 330.    ,&    !< a constant for terminal velocity of graupel
          bvtg   = 0.8     ,&    !< a constant for terminal velocity of graupel
#endif
          deng   = 500.    ,&    !< density of graupel                   [kgm-3]
#ifdef HAIL
#ifdef HEYMS_VTGH
          avth   = 168.0   ,&    !< a constant for terminal velocity of hail
          bvth   = 0.72    ,&    !< a constant for terminal velocity of hail
#else
          avth   = 285.    ,&    !< a constant for terminal velocity of hail
          bvth   = 0.8     ,&    !< a constant for terminal velocity of hail
#endif
          denh   = 912.    ,&    !< density of hail                      [kgm-3]
#endif
#ifndef DOUBLE
          lamdarmax = 8.e4 ,&    !< maximum slope parameter for rain       [m-1]
          lamdarmin = 5.e2 ,&    !< minimum slope for rain         [m-1] (5.8mm)
#endif
          lamdasmax = 1.e5 ,&    !< maximum slope for snow       (10micro) [m-1]
#ifdef WDM7_LAMG
          lamdagmax = 6.e4 ,&    !< maximum slope for graupel    (50micro) [m-1]
#else
          lamdagmax = 2.e4 ,&    !< maximum slope for graupel    (50micro) [m-1]
#endif
#ifdef HAIL
          lamdahmax = 2.e4 ,&    !< maximum slope for hail       (50micro) [m-1]
#endif
          peaut  = .55     ,&    !< collection efficiency in autoconversion
          r0     = .8e-5   ,&    !< threshold radius in autoconversion       [m]
          xncr   = 3.e8    ,&    !< number concentration of cloud droplets [m-3]
          xmyu   = 1.718e-5,&    !< the dynamic viscosity             [kgm-1s-1]
          dimax  = 500.e-6 ,&    !< maximum cloud-ice diamter                [m]
          inpmax = 500.e6  ,&    !< maximum ice nuclei particles (inp)     [m-3]
          inpmin = 0.1     ,&    !< minimum ice nuclei particles (inp)     [m-3]
          pfrz1  = 100.    ,&    !< constant in biggs freezing
          pfrz2  = 0.66    ,&    !< constant in biggs freezing
          qcmin  = 1.e-12  ,&    !< minimun values for qr, qs, and qg   [kgkg-1]
          qrmin  = 1.e-11  ,&    !< minimun values for qr, qs, and qg   [kgkg-1]
          eacrc  = 1.0     ,&    !< snow/cloud-water collection efficiency
#ifdef HAIL
          eachs  = 1.0     ,&    !< hail/snow collection efficiency
          eachg  = 0.5     ,&    !< hail/graupel collection efficiency
          cd     = 0.6     ,&    !< drag coefficient for hailstone
#endif
          dens   = 100.    ,&    !< density of snow                      [kgm-3]
          qs0    = 0.6e-3        !< threshold for snow aggretion        [kgkg-1]
!
#ifdef DOUBLE
   real(kind=kind_phys), parameter, private :: &
          drmin   = 30.e-6    ,&  !< minimum diameter for rain               [m]
          drmax   = 3000.e-6  ,&  !< maximum diameter for rain               [m]
          drcoeff = (24.)**(.3333333),& !< factor for raindrop diameter 2.89
#ifdef WDM7_LAMR
          lamdarmax = 5.0e4, &
          lamdarmin = 2.0e3, &
#else
          lamdarmax = drcoeff/(drmin)  ,& !< maximum slope for rain     [m-1] 
          lamdarmin = drcoeff/(drmax)  ,& !< minimum slope for rain    [m-1]
#endif
          nrmin     = 1.e-2   ,&  !< nr minimum                            [m-3]
          nrmax     = 300.e6  ,&  !< nr maximum                            [m-3]
          dcmin     = 1.e-6   ,&  !< minimum diameter for clouds             [m]
          dcmax     = 82.e-6  ,&  !< maximum diameter for clouds             [m]
#ifdef WDM7_LAMC
          lamdacmax = 5.0e5, &
          lamdacmin = 2.0e4, &
#else
          lamdacmax = 1./(dcmin)  ,&  !< maximum slope for cloud          [m-1]
          lamdacmin = 1./(dcmax)  ,& !< minimum slope for cloud          [m-1]
#endif
#ifdef NEW_NCRMIN
          ncmin     = 10.e6   ,&  !< nc minimum                            [m-3]
#else
          ncmin     = 1.e1   ,&  !< nc minimum                            [m-3]
#endif
          ncmax     = 30000.e6,&  !< nc maximum                            [m-3]
#ifdef HDC2004_PRAUT
          xncr0   = 5.e7     ,&  !< nc for autoconversion over ocean       [m-3]
          xncr1   = 5.e8     ,&  !< nc for autoconversion over lands       [m-3]
#endif
#ifdef  NCACT_LEI2019
          satmax  = 0.0048   ,&  !< supersatuation for ccn activation
#else
          satmax  = 1.0048   ,&  !<
#endif
          actk    = 0.6      ,&  !< a constant for ccn activation
          actr    = 1.5e-6   ,&  !< initial radius for activated cloud       [m]
          ncrk1   = 3.03e3   ,&  !< Long’s kernel coefficient           [m-3s-1]
          ncrk2   = 2.59e15  ,&  !< Long’s kernel coefficient           [m-3s-1]
          di5000  = 5000.e-6 ,&  !< diameter                                 [m]
          di2000  = 2000.e-6 ,&  !< diameter                                 [m]
          di1000  = 1000.e-6 ,&  !< diameter                                 [m]
          di600   = 600.e-6  ,&  !< diameter                                 [m]
          di100   = 100.e-6  ,&  !< diameter                                 [m]
          di82    = 82.e-6   ,&  !< diameter                                 [m]
          di15    = 15.e-6   ,&  !< diameter                                 [m]
          di5     = 5.e-6    ,&  !< diameter                                 [m]
          di2     = 2.e-6    ,&  !< diameter                                 [m]
          di1     = 1.e-6        !< diameter                                 [m]
   real(kind=kind_phys), parameter, private :: &
          ccn_0   = 50.e+6,   &  !< background ccn                         [m-3]
          ccnmax  = 20000.e+6,&  !< maximum ccn                            [m-3]
          ccnmin  = 50.e+6       !< minimum ccn                            [m-3]
#endif
!-------------------------------------------------------------------------------
! min max  effective radius (m)
!
   real(kind=kind_phys), parameter, private :: &
          recmin  = 2.51e-6  ,&    !< minimum effective radius for cloud     [m]
          recmax  = 50.e-6   ,&    !< maximum effective radius for cloud     [m]
          reimin  = 5.01e-6  ,&    !< minimum effective radius for ice       [m]
          reimax  = 125.e-6  ,&    !< maximum effective radius for ice       [m]
          resmin  = 25.e-6   ,&    !< minimum effective radius for snow      [m]
          resmax  = 999.e-6        !< maximum effective radius for snow      [m]
!-------------------------------------------------------------------------------
! semi_lagrangian option sub-stepping (greater than 100 means single-loop)
!
   real(kind=kind_phys), parameter, private :: &
          sedi_semi_cfl = 10000.0
!-------------------------------------------------------------------------------
! inline expansion table for saturation vapor pressure fsvp_water and fsvp
!
!  use constant, only : cliq_,cice_,cvap_,cpd_,hvap_,hsub_,psat_,rd_,rv_,ttp_
!
   integer, parameter                       ::  nxsvp = 7501
   real(kind=kind_phys)                     ::  c1xsvp,c2xsvp
   real(kind=kind_phys),  dimension(nxsvp)  ::  tbsvp
   integer, parameter                       ::  nxsvpw = nxsvp
   real(kind=kind_phys)                     ::  c1xsvpw,c2xsvpw
   real(kind=kind_phys),  dimension(nxsvpw) ::  tbsvpw
!
   real(kind=kind_phys), parameter, private :: &
         cliq_ = 4.1900e+3     ,&    !< specific heat of liquid water [jkg-1k-1]
         cice_ = 2.1060e+3     ,&    !< specific heat of ice water    [jkg-1k-1]
         cvap_ = 1.8700e+3     ,&    !< specific heat of vapor water  [jkg-1k-1]
         cpd_  = 1.00464e+3    ,&    !< specific heat of dry air      [jkg-1k-1]
         hvap_ = 2.5010e+6     ,&    !< latent heat for vaporization     [jkg-1]
         hsub_ = 2.8340e+6     ,&    !< latent heat for sublimation      [jkg-1]
         psat_ = 6.1078e+2     ,&    !< saturated vapor pressure at 0c      [pa]
         rd_   = 2.8704e+2     ,&    !< gas constant for dry air      [jkg-1k-1]
         rv_   = 4.6150e+2     ,&    !< gas constant for water vapor  [jkg-1k-1]
         ttp_  = 2.7316e+2           !< temperature at triple point          [k]
  real(kind=kind_phys), parameter, private :: &
         psatk = psat_,                         &
         dldt  = cvap_-cliq_,                   &
         dldti = cvap_-cice_,                   &
         xa    = -dldt/rv_,                     &
         xb    = xa+hvap_/(rv_*ttp_),           &
         xai   = -dldti/rv_,                    &
         xbi   = xai+hsub_/(rv_*ttp_)
!-------------------------------------------------------------------------------
! inline expansion table for shape parameter dependent sedimentation
!
   integer, parameter    ::  nxshape = 9991 !< 10 micro to 10 mm every 1 micro
   real(kind=kind_phys)                            ::  c1xshape,c2xshape
   real(kind=kind_phys), dimension(nxshape)        ::  tbshape
!-------------------------------------------------------------------------------
! inline expansion table for lb2017 autoconversion
!
   integer, parameter    ::  nxaut = 9991 !<  0.1 micro to 100 micro every 0.1
   integer, parameter    ::  numax = 15   !<  maximum dispersion parameter
   real(kind=kind_phys)       ::  c1xaut,c2xaut
   double precision, dimension(nxaut,numax)        ::  tbaut_qc_tot, tbaut_qc_sub
   double precision, dimension(nxaut,numax)        ::  tbaut_nc_tot, tbaut_nc_sub
#ifdef ICE_SPHERE
!-------------------------------------------------------------------------------
! apply temperature dependent ice habbit....based on aspect ratio concept 
! of Um et al. (2015, acp) and Matrosov et al. (2020, jamc)
!
   real(kind=kind_phys), parameter, private :: &
         t1_sphere = -10., &        !< threshold for spherical ice          [C]
         t2_sphere = -33.           !< threshold for hexagonal ice          [C]
#endif
!-------------------------------------------------------------------------------
! miscellaneous...
!
   real(kind=kind_phys), parameter, private :: &
          zero_0 = 0.0, &
          one_1 = 1.0
   logical, parameter, private :: flgzero = .true.
!-------------------------------------------------------------------------------
! save the microphysics constants that are initialzed in init driver
!
#ifdef DOUBLE
   real(kind=kind_phys), save :: &
             qc_ocean,qc_land,qck1,pidnc,bvtr1,bvtr2,bvtr3,bvtr4,bvtr5,        &
             bvtr6,bvtr7, bvtr2o5,bvtr3o5,                                     &
             g1pbr,g2pbr,g3pbr,g4pbr,g5pbr,g6pbr,g7pbr,                        &
             g5pbro2,g7pbro2,pi,pisq,                                          &
             pvtr,pvtrn,eacrr,pacrr,pidn0r,pidnr,                              &
             precr1,precr2,roqimax,bvts1,bvts2,                                &
             bvts3,bvts4,g1pbs,g3pbs,g4pbs,g5pbso2,                            &
             pvts,pacrs,precs1,precs2,pidn0s,xlv1,                             &
             bvtg1,bvtg2,bvtg3,bvtg4,g1pbg,g3pbg,g4pbg,                        &
             g5pbgo2,pvtg,pacrg,precg1,precg2,pidn0g,                          &
#ifdef HAIL
             g6pbgh,precg3,bvth2,bvth3,bvth4,g3pbh,g4pbh,g5pbho2,pacrh,pvth,   &
             prech1,prech2,prech3,pidn0h,                                      &
             rslopehmax,rslopehbmax,rslopeh2max,rslopeh3max,                   &
#endif
             rslopecmax,rslopec2max,rslopec3max,                               &
             rslopermax,rslopesmax,rslopegmax,                                 &
             rsloperbmax,rslopesbmax,rslopegbmax,                              &
             rsloper2max,rslopes2max,rslopeg2max,                              &
             rsloper3max,rslopes3max,rslopeg3max
#else
   real(kind=kind_phys), save :: &
          qc_ocean, qc_land, qck1, pidnc,          &
          bvtr1,bvtr2,bvtr3,bvtr4,g1pbr,           &
          g3pbr,g4pbr,g5pbro2,pvtr,eacrr,pacrr,    &
          bvtr6,g6pbr,                             &
          precr1,precr2,roqimax,bvts1,             &
          bvts2,bvts3,bvts4,g1pbs,g3pbs,g4pbs,     &
          bvtg1,bvtg2,bvtg3,bvtg4,g1pbg,           &
          g3pbg,g4pbg,g5pbgo2,pvtg,pacrg,          &
          precg1,precg2,pidn0g,                    &
          g5pbso2,pvts,pacrs,precs1,precs2,pidn0r, &
          pidn0s,xlv1,pi,pisq,               &
#ifdef HAIL
          g6pbgh,precg3,bvth2,bvth3,bvth4,g3pbh,g4pbh,g5pbho2,pacrh,pvth,   &
          prech1,prech2,prech3,pidn0h,                                      &
          rslopehmax,rslopehbmax,rslopeh2max,rslopeh3max,                   &
#endif
          rslopermax,rslopesmax,rslopegmax,        &
          rsloperbmax,rslopesbmax,rslopegbmax,     &
          rsloper2max,rslopes2max,rslopeg2max,     &
          rsloper3max,rslopes3max,rslopeg3max
#endif
!
#ifdef LB2017_PRAUT
   real(kind=kind_phys), dimension(numax), private, save :: ccg, ocg1
   real(kind=kind_phys), parameter, private :: bm_r1 = 3.0
   real(kind=kind_phys), parameter, private :: obmr1 = 1./bm_r1
   real(kind=kind_phys), save :: mrmin
#endif
!
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!! \section MP_ARG_INIT_ Argument Table
!! \htmlinclude MP_INIT_.html
!!
!-------------------------------------------------------------------------------
      subroutine MP_INIT_(den0,denr,dens,cl,cv,ccn0, &
                  nwfa2d, nwfa, nifa, if_aerosol_aware,              &
                          ids, ide, jds, jde, kds, kde,              &
                          ims, ime, jms, jme, kms, kme,              &
                          its, ite, jts, jte, kts, kte)
         implicit none
  integer,           intent(in   ) :: ids,ide, jds,jde, kds,kde ,    &
                                      ims,ime, jms,jme, kms,kme ,    &
                                      its,ite, jts,jte, kts,kte
   real, intent(in) :: den0,denr,dens,cl,cv,ccn0
   logical,          intent(in   ) :: if_aerosol_aware
!
         integer :: i, k
   real, dimension( ims:ime, kms:kme, jms:jme ), intent(inout),      &
                                                      optional ::    &
                                                            nwfa,    &
                                                            nifa
   real, dimension( ims:ime, jms:jme ), intent(in), optional ::      &
                                                          nwfa2d
!
!  local
!
   logical ::          if_aerosol_data
!
   if_aerosol_data = .true.
   if_aerosol_activated = .false.
!
! check aerosol option
!
   if ((.not.PRESENT(nwfa)) .or. (.not.PRESENT(nifa)))               &
           if_aerosol_data = .false.
!
   if_aerosol_activated = if_aerosol_aware .and. if_aerosol_data
!
         call udminit(den0,denr,dens,cl,cv,ccn0)
         call udm_funct_shape_setup
         call udm_funct_svp_setup
         call udm_funct_lb2017_setup
!
      end subroutine MP_INIT_
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
      subroutine MP_FINALIZE_()
!-------------------------------------------------------------------------------
      end subroutine MP_FINALIZE_
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
!> \section MP_ARG_RUN_ Argument Table
!! \htmlinclude MP_RUN_.html
!!
subroutine MP_RUN_(ncol, nlev, t1, q1 &
#ifdef HAIL
                   ,qc1, qr1, qi1, qs1, qg1, qh1                  &
#else 
                   ,qc1, qr1, qi1, qs1, qg1                       &
#endif
#ifdef DOUBLE
                   ,nn1,nc1,nr1                                   &
#endif
                   ,den1, p1, delz1                               &
                   ,delt, g, cpd, cpv, rd, rv, t0c                &
                   ,ep1, ep2, qmin                                &
                   ,xls, xlv0, xlf0, den0, denr                   &
                   ,cliq,cice,psat                                &
                   ,lat                                           &
                   ,area                                          &
                   ,slimsk                                        &
                   ,rain                                          &
                   ,snow                                          &
                   ,graupel                                       &
                   ,hail                                          &
#ifdef WRF
                   ,rainncv                                       &
                   ,sr                                            &
                   ,snowncv                                       &
                   ,graupelncv                                    &
                   ,hailncv                                       &
#endif
                   ,itimestep,nifa1                               &
                   ,if_aerosol_activated                          &
#ifdef wrf_chem
                   ,rainprod2d, evapprod2d                        &
#endif
                   ,errmsg,errflg                                 &
                                                                  )
!
!-------------------------------------------------------------------------------
!===============================================================================
!! The Unified Forecast System (ufs) Double Moment Microphyiscs Scheme (udm mp)
!! largely adopts ice microphys processes of wrf Single-Moment (WSM) and warm
!! microphysics from wrf Double-Moment (WDM) schemes, with bug fixes and
!! revisions based on literature and accumulated realism. udm MP utilizes
!! the in-cloud microphysics concept (Kim and Hong 2018) and physically based
!! autoconversion of Lee and Baik (2017) by solving the stochastic collection
!! equation. For grauel and hail, the corresponding termimal velocity is updated  
!! The semi-lagrangian sedimentation of Juang and Hong (2010) is reconfigured
!! for computational efficiency and numerical accuracy.
!! The scheme incorporates aerosol data interacting with cloud condensation 
!! nuclei (CCN) ice nuclei particles (INP) as options.
!!
!! All productions terms are optimized by introducing a cloud top definition.
!!
!! All units are in m.k.s. and source/sink terms in kgkg-1s-1.
!!
!!  udm mp scheme, developed by songyou hong, haiqin li, jian-wen bao (noaa)
!!         jimy dhdhia (ncar), spring, summer, fall, and winter 2021-2025
!!  wdm mp scheme, coded by kyo-sun lim and songyou hong (yonsei univ.),
!!         fall 2008
!!  wsm mp scheme, coded by songyou hong, jeong-ock jade lim (yonsei univ.),
!!         jimy dudhia, and shu-hua chen (ncar),summer 2002
!!  ncep mp scheme, coded by songyou hong, henry juang, qingyun zhao (ncep),
!!         summer 1997
!!
!!  history :
!!        03-04-2026 : CAMS aerosol option added by taeho mun (unist, Korea) 
!!
!!  references :
!!        heymsfield et al(hey, 2018) j. atmos. Sci.
!!        kim and and hong (kh, 2018) j. atmos. sci.
!!        lee and baik (lb2017, 2017) j. atmos. sci.
!!        juang and hong (jh, 2010) mon. wea. rev.
!!        lim and hong (lh, 2010) mon. wea. rev.
!!        dudhia, hong and lim (dhl, 2008) j. meteor. soc. japan
!!        hong and lim (hl, 2006) j. korean meteor. soc.
!!        hong, dudhia, and chen (hdc, 2004) mon. wea. rev.
!!        hong, juang, and zhao (hjz, 1998) mon. wea. rev.
!!
!! \section structure
!!            --- udm_mp  ---  |- slope_rain slope_snow slope_graupel slope_hail
!!                             |
!!                             |- semi_lagrangian
!!                             |
!!                             |- cldf_diag
!!                             |
!!                             |- udminit
!!                             |    |
!!                             |    ------ rgmma
!!                             |
!!                             |- adjust_number_concent
!!                             |
!!                             |- find_cloud_top
!!                             |
!!                             |- udm_mp_effective_radius
!!                             |
!!                             |- udm_mp_reflectivity
!!                             |
!!                             |- udm_funct_shape_setup
!!                             |    |
!!                             |    ------ fshapex, fshape
!!                             |
!!                             |- udm_funct_svp_setup
!!                             |    |
!!                             |    ------ fsvpx, fsvpxw, fsvp, fsvp_water
!!                             |
!!                             |- udm_funct_lb2017_setup
!!                             |    |
!!                             |    ------ funct_aut_qc_tot, funct_aut_qc_sub,
!!                             |    ------ funct_aut_nc_tot, funct_aut_nc_sub,
!!                             |   
!!                             |------ iceDeMott
!!                             |
!! \ section input and output
!!  input :
!!    delt                        - timestep
!!    g, cpd, cpv, t0c, den0,     - constant
!!    rd, rv, ep1, ep2, qmin,
!!    xls, xlv0, xlf0, denr,
!!    cliq, cice, psat
!!    ids, ide, jds, jde, kds, kde - dimension
!!    ims, ime, jms, jme, kms, kme
!!    its, ite, jts, jte, kts, kte
!!    ncloud                       - number of hydrometeor
!!    p                            - pressure, pa
!!    delz                         - depth of model layer, m
!!
!!  inout :
!!    t1                           - temperautre
!!    q1                           - specific humidity
!!    q2                           - mixing ratio of cloud, rain, ice, and snow
!!                                   qc, qr, qi, qc
!!    rain, rainncv                - precipitation
!!    sr                           - ratio of snow to rain
!!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   integer                              , intent(in   ) :: itimestep
   integer                              , intent(in   ) :: ncol, nlev
#ifndef WRF
   integer                              , intent(in   ) :: kdt
   real(kind=kind_phys)                 , intent(in   ) :: fhswr
#endif
   real(kind=kind_phys)                 , intent(in   ) :: delt
   real(kind=kind_phys)                 , intent(in   ) :: g, cpd, cpv, t0c,   &
                                                             rd, rv, ep1, ep2
   integer, dimension( ITDIM_ )         , intent(in   ) :: slimsk
   real(kind=kind_phys)                 , intent(in   ) ::                     &
                                                             qmin, xls,        &
                                                             xlv0, xlf0,       &
                                                             den0, denr,       &
                                                             cliq, cice, psat
   integer                              , intent(in   ) :: lat
   logical                              , intent(in   ) :: if_aerosol_activated
   real(kind=kind_phys), dimension( IMDIM_,KMDIM_ ), intent(in   ) :: p1
   real(kind=kind_phys), dimension( IMDIM_,KMDIM_ ), intent(in   ) :: delz1
   real(kind=kind_phys), dimension( IMDIM_,KMDIM_ ), intent(in   ) :: den1
   real(kind=kind_phys), dimension( IMDIM_,KMDIM_ ), intent(inout) :: q1
   real(kind=kind_phys), dimension( ITDIM_,KTDIM_ ), intent(inout) :: t1
   real(kind=kind_phys), dimension( IMDIM_ )       , intent(in   ) :: area
   real(kind=kind_phys), dimension( IMDIM_,KMDIM_ ), intent(inout) :: qc1, qr1
#ifdef HAIL
   real(kind=kind_phys), dimension( IMDIM_,KMDIM_ ), intent(inout) :: &
                                                           qi1, qs1, qg1, qh1
#else
   real(kind=kind_phys), dimension( IMDIM_,KMDIM_ ), intent(inout) :: &
                                                                qi1, qs1, qg1
#endif
   real(kind=kind_phys), dimension( IMDIM_,KMDIM_ ), intent(inout) :: &
                                                                nn1, nc1, nr1
   real(kind=kind_phys), dimension( IMDIM_,KMDIM_ ), intent(inout) :: &
                                                                        nifa1
   real(kind=kind_phys), dimension( IMDIM_ )          , intent(inout) :: rain
   real(kind=kind_phys), dimension( IMDIM_ ), optional, intent(inout) :: &
                                                                         snow, &
                                                                      graupel, &
                                                                         hail
#ifdef WRF
   real(kind=kind_phys), dimension( IMDIM_ )          , intent(inout) :: &
                                                                      rainncv, &
                                                                           sr
   real(kind=kind_phys), dimension( IMDIM_ ), optional, intent(inout) :: &
                                                                      snowncv, &
                                                                   graupelncv, &
                                                                      hailncv
#endif
#ifdef wrf_chem
   real(kind=kind_phys), dimension( ITDIM_ , KTDIM_ ), intent(inout)  ::       &
                                                       rainprod2d, evapprod2d
#endif
   character(len=*), intent(out)   :: errmsg
   integer,          intent(out)   :: errflg
!
! local variables
!
   real(kind=kind_phys), dimension( ITDIM_ ) ::          dxmeter
   real(kind=kind_phys), dimension( KTDIM_ , ITDIM_ ) :: &
                                                                q, &
                                                              den, &
                                                                t, &
                                                                p, &
                                                             delz
#ifdef HDC2004_PRAUT
   real(kind=kind_phys), dimension( ITDIM_ ) :: &
                                                              qc0
#endif
   real(kind=kind_phys), dimension( KTDIM_ , ITDIM_ , 2) :: &
                                                              qci
   real(kind=kind_phys), dimension( KTDIM_ , ITDIM_ , 4) :: &
                                                              qrs
   real(kind=kind_phys), dimension( KTDIM_, 4)  :: &
                                                           rslope, &
                                                          rslope2, &
                                                          rslope3, &
                                                          rslopeb
#ifdef DOUBLE
   real(kind=kind_phys), dimension( KTDIM_ , ITDIM_ , 3) :: &
                                                              ncr
   real(kind=kind_phys), dimension( KTDIM_ , ITDIM_ )    :: &
                                                             nifa
#endif
   real(kind=kind_phys), dimension( KTDIM_ ) :: &
                                                              xlv, &
                                                              xlf, &
                                                              cpm, &
                                                           rh_mul, &
                                                           rh_ice, &
                                                         qsat_mul, &
                                                         qsat_ice, &
                                                               tc, &
                                                               tw, &
                                                             dend, &
                                                           denfac, &
                                                              vtr, &
                                                              vts, &
                                                              vtg, &
#ifdef HAIL
                                                              vth, &
#endif
                                                           vtmean, &
                                                           sumice, &
                                                              vti, &
#ifdef ICE_SPHERE
                                                             vtis, &
#endif
                                                               ni, &
                                                              ni0, &
                                                               mi, &
                                                               di, &
#ifdef ICE_SPHERE
                                                              dis, &
#endif
                                                            denqr, &
                                                            denqs, &
                                                            denqg, &
#ifdef HAIL
                                                            denqh, &
#endif
                                                            denqi, &
                                                           n0sfac, &
                                                           n0gfac, &
                                                           n0hfac, &
                                                           diffus, &
                                                           viscok, &
                                                           viscod, &
                                                           ab_mul, &
                                                           ab_ice, &
                                                           venfac
#ifdef DOUBLE
   real(kind=kind_phys), dimension( KTDIM_ ) :: &
                                                          rslopec, &
                                                         rslopec2, &
                                                         rslopec3, &
                                                               dc, &
                                                               dr, &
                                                            denqn
#endif
   real(kind=kind_phys), dimension( KTDIM_ ) :: &
                                                            pigen, &
                                                            pcond, &
                                                            prevp, &
                                                            psevp, &
                                                            pgevp, &
                                                            pidep, &
                                                            psdep, &
                                                            pgdep, &
#ifdef HAIL
                                                            phdep, &
#endif
                                                            praut, &
                                                            psaut, &
                                                            pgaut, &
                                                            piacr, &
                                                            pracw, &
                                                            praci, &
                                                            pracs, &
                                                            psacw, &
                                                            psaci, &
                                                            psacr, &
#ifdef HAIL
                                                            phaut, &
                                                            phacr, &
                                                            phacs, &
                                                            phacg, &
                                                            phaci, &
                                                            phmlt, &
                                                            pheml, &
                                                            phevp, &
                                                            primh, &
                                                            pvapg, &
                                                            pvaph, &
                                                            pgwet, &
                                                            phwet, &
                                                          pgaci_w, &
                                                          phaci_w, &
                                                            nhacw, &
                                                            nhacr, &
                                                            nheml, &
#endif
                                                            pgacw, &
#ifdef HAIL
                                                            phacw, &
                                                            pracg, &
#endif
                                                            pgaci, &
                                                            pgacr, &
                                                            pgacs, &
                                                            paacw, &
                                                            psmlt, &
                                                            pgmlt, &
                                                            pseml, &
                                                            pgeml
#ifdef DOUBLE
   real(kind=kind_phys), dimension( KTDIM_ ) :: &
                                                            pcact, &
                                                            nraut, &
                                                            nracw, &
                                                            ncevp, &
                                                            nccol, &
                                                            nrcol, &
                                                            nsacw, &
                                                            ngacw, &
                                                            niacr, &
                                                            nsacr, &
                                                            ngacr, &
                                                            naacw, &
                                                            nseml, &
                                                            ngeml, &
                                                            ncact
#endif
   real(kind=kind_phys), dimension( KTDIM_ ) :: &
                                                           satrdt, &
                                                           supice, &
                                                           supsat, &
                                                             cldf
   real(kind=kind_phys) :: &
                                                           qrpath, &
                                                           qspath, &
                                                           qgpath, &
                                                           qipath, &
                                                         precip_r, &
                                                         precip_s, &
                                                         precip_g, &
                                                         precip_i
   logical, dimension( KTDIM_ ) :: &
                                                            ifsat, &
                                                            ifice
!
! for reflectivity....
!
   real(kind=kind_phys), dimension( KTDIM_ ) :: &
                                                              dbz, &
                                                            re_qc, &
                                                            re_qi, &
                                                            re_qs
#ifdef LB2017_PRAUT
!
! lb2017 autoconversion
!
   real(kind=kind_phys), dimension( KTDIM_ )   :: ncaut
   integer :: nu_c
   logical :: L_qc
   real(kind=kind_phys) :: L_c, n_c
   double precision :: n0c, n0csq, lamc
   double precision, parameter :: vt0_c = 1.09734d8
   double precision, parameter :: kc_aut = 1.35429d14
   double precision, parameter :: alpha_aut = 0.88
   double precision, dimension(0:42) :: gamma1
   double precision ::                                             &
                   qc_tot, qc_sub, qc_aut,                         &
                   nc_tot, nc_sub, nc_aut, nr_aut
#endif
!
! miscellaneous variables
!
   real(kind=kind_phys)  :: &
             rdtcld,conden, x, y, z, a, b, c, d, e,                &
             coeres, dtcld, eacrs, acrfac, egi,                    &
#ifdef HAIL
             ehi, rs0, ghw1, ghw2, ghw3, ghw4, precip_h, qhpath,   &
#endif
             qimax, roqi0, lamdar, diameter,                       &
             precip_sum, precip_ice, factor, source, htotal,       &
             hvalue, pfrzdtc, pfrzdtr,                             &
             tstepsnow, tstepgraupel, tstephail,                   &
             alpha2, delta2, delta3, dtcfl, dr_embryo,temp
#ifdef DOUBLE
   real(kind=kind_phys)  :: &
                              gfac, sfac, nfrzdtr, nfrzdtc, qnpath
   real(kind=kind_phys), dimension( KTDIM_ ) :: &
                                       qrconcr, qrcon, taucon, vtn
#endif
   integer :: &
              nstep, niter, n, lond, latd, i, j, k, loop, loops,    &
              idim, kdim, ktopini, ktopmax, ktop, ktopqc, ktopqi,   &
              ktopqr, ktopqs, ktopqg, ktopqh, ktoprh
   integer :: its, ite, kts, kte
#ifdef PRT
   integer :: iprt
   real(kind=kind_phys) :: sum1,sum2
#endif
#ifndef SEMI_RAIN
   real(kind=kind_phys), dimension( KTDIM_ ) :: & 
                            sedir, sedin, dqr, dnr, qrs_tmp, ncr_tmp
#endif
   logical :: &
                                flgcld, lqc, lqi, lqr, lqs, lqg, lqh
!----------------------------------------------------------------
! preparation ....
!
   its = 1
   ite = ncol
   kts = 1
   kte = nlev
   errmsg = ''
   errflg = 0
   idim = ite-its+1
   kdim = kte-kts+1
   ktopini = kte - 1
   dxmeter(:) = sqrt(area(:))
#ifdef DOUBLE
   dr_embryo =  1./(pi*denr/6.*di82**3)
#endif
#ifdef PRT
   lond = ncol/2 + 1
   latd = 1
#endif
!-------------------------------------------------------------------------------
! assign passing variables to local arrays
!
   t(:,:) = 0.; q(:,:) = 0.; p(:,:) =0.
   delz(:,:) = 0.; den(:,:) = 0.
   nifa(:,:) = 0. 
   ktop = kte
   do k = kts, ktop
     do i = its, ite
       t(k,i)  = t1(i,k)
       q(k,i)  = q1(i,k)
       p(k,i)  = p1(i,k)
       delz(k,i) = delz1(i,k)
       den(k,i)  = den1(i,k)
     enddo
   enddo
!-------------------------------------------------------------------------------
! padding 0 for negative values generated by dynamics - optional
!
   qci(:,:,:) = qmin
   qrs(:,:,:) = qmin
#ifdef DOUBLE
   ncr(:,:,:) = qmin
#endif
!
   do k = kts,ktop
     do i = its,ite
       if(flgzero) then
         qci(k,i,1) = max(qc1(i,k),0.0)
         qrs(k,i,1) = max(qr1(i,k),0.0)
         qci(k,i,2) = max(qi1(i,k),0.0)
         qrs(k,i,2) = max(qs1(i,k),0.0)
         qrs(k,i,3) = max(qg1(i,k),0.0)
#ifdef HAIL
         qrs(k,i,4) = max(qh1(i,k),0.0)
#endif
#ifdef DOUBLE
         ncr(k,i,1) = min(max(nn1(i,k),ccnmin),ccnmax)
         ncr(k,i,2) = max(nc1(i,k),0.0)
         ncr(k,i,3) = max(nr1(i,k),0.0)
#endif
       else
         qci(k,i,1) = qc1(i,k)
         qrs(k,i,1) = qr1(i,k)
         qci(k,i,2) = qi1(i,k)
         qrs(k,i,2) = qs1(i,k)
         qrs(k,i,3) = qg1(i,k)
#ifdef HAIL
         qrs(k,i,4) = qh1(i,k)
#endif
#ifdef DOUBLE
         ncr(k,i,1) = min(max(nn1(i,k),ccnmin),ccnmax)
         ncr(k,i,2) = nc1(i,k)
         ncr(k,i,3) = nr1(i,k)
#endif
       endif
     enddo
   enddo
!
   if (if_aerosol_activated) then 
     do i = its, ite
       do k = kts, ktop
         nifa(k,i) = min(max(nifa1(i,k),inpmin),inpmax)
       enddo
     enddo
   endif
!
!-------------------------------------------------------------------------------
! initialize the surface rain, snow, graupel and etc....
!
#ifdef WRF
   do i = its, ite
     if(present (snowncv) .and. present (snow)) snowncv(i) = 0.
     if(present (graupelncv) .and. present (graupel)) graupelncv(i) = 0.
     if(present (hailncv) .and. present (hail)) hailncv(i) = 0.
   enddo
   rainncv(:) = 0.;  sr(:) = 0.
#endif
#ifndef WRF
!
   rain(:) = 0.; snow(:) = 0.; graupel(:) = 0.
#ifdef HAIL
   hail(:) = 0.
#endif
#endif
#ifdef HDC2004_PRAUT
   qc0(:) =0.
#endif
#ifdef DOUBLE
#ifdef HDC2004_PRAUT
!
   qc0(:) = 0.0
   do i = its,ite
     if(slimsk(i)==0) then      ! water
       qc0(i) = qc_ocean
     else
       qc0(i) = qc_land
     endif
   enddo
#endif
#endif
!===============================================================================
!
! inner loop wih the time step of dtcldcr (default = 180 sec)
!
!===============================================================================
! compute the sub-cycling time steps.
!
   loops = max(ceiling(delt/dtcldcr),1)
   dtcld = delt/loops
   if(delt<=dtcldcr) dtcld = delt
   rdtcld = 1./dtcld
!
!-------------------------------------------------------------------------------
!
 inner_loop : do loop = 1, loops ! sub-stepping with cld timestep = 180s
!
!-------------------------------------------------------------------------------
!
 i_loop : do i = its, ite ! i-loop for one-dimensional code
!
!-------------------------------------------------------------------------------
! compute moisture variables
!
   qsat_mul(:) = 0.; rh_mul(:)   = 0.
   qsat_ice(:) = 0.; rh_ice(:)   = 0.
   ktopqc = 0; ktopqi = 0; ktopqr = 0; ktopqs = 0; ktopqg = 0; ktoprh = 0
#ifdef HAIL
   ktopqh = 0
#endif
!
   ktop = ktopini
   do k = kts, ktop
     qsat_mul(k) = fsvp_water(t(k,i),p(k,i))      !< mul means water in korean
     qsat_mul(k) = ep2 * qsat_mul(k) / (p(k,i) - qsat_mul(k))
     qsat_mul(k) = max(qsat_mul(k),qmin)
     rh_mul(k)   = max(q(k,i) / qsat_mul(k),qmin)
     qsat_ice(k) = fsvp(t(k,i),p(k,i))
     qsat_ice(k) = ep2 * qsat_ice(k) / (p(k,i) - qsat_ice(k))
     qsat_ice(k) = max(qsat_ice(k),qmin)
     rh_ice(k)   = max(q(k,i) / qsat_ice(k),qmin)
   enddo
!
   cpm(:) = 0.; xlv(:) = 0.; xlf(:) = 0.
   do k = kts, ktop
     cpm(k) = cpd*(1.-max(q(k,i),qmin)) + max(q(k,i),qmin)*cpv
     xlv(k) = xlv0 - xlv1*(t(k,i)-t0c)
     xlf(k) = xls - xlv(k)
   enddo
!-------------------------------------------------------------------------------
! re-define cloud top for numerical efficiencies
!
   lqc = .false.
   lqi = .false.
   lqr = .false.
   lqs = .false.
   lqg = .false.
#ifdef HAIL
   lqh = .false.
#endif
   flgcld = .false.
   call find_cloud_top(1,kdim,ktopini,qci(:,i,1),zero_0,ktopqc)
   call find_cloud_top(1,kdim,ktopini,qci(:,i,2),zero_0,ktopqi)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,1),zero_0,ktopqr)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,2),zero_0,ktopqs)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,3),zero_0,ktopqg)
#ifdef HAIL
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,4),zero_0,ktopqh)
#endif
   call find_cloud_top(1,kdim,ktopini,rh_ice(:), one_1,ktoprh)
   if(ktopqc>0.0) lqc = .true.
   if(ktopqi>0.0) lqi = .true.
   if(ktopqr>0.0) lqr = .true.
   if(ktopqs>0.0) lqs = .true.
   if(ktopqg>0.0) lqg = .true.
#ifdef HAIL
   if(ktopqh>0.0) lqh = .true.
   ktopmax = max(ktopqc,ktopqi,ktopqr,ktopqs,ktopqg,ktopqh,ktoprh)
#else
   ktopmax = max(ktopqc,ktopqi,ktopqr,ktopqs,ktopqg,ktoprh)
#endif
!-------------------------------------------------------------------------------
! early checkout
!
#ifdef HAIL
   if(lqc .or. lqi .or. lqr .or. lqs .or. lqg .or. lqh) flgcld = .true.
#else
   if(lqc .or. lqi .or. lqr .or. lqs .or. lqg) flgcld = .true.
#endif
   if((.not.flgcld) .and. ktoprh>0.0) flgcld = .true.
   if(.not.flgcld) then
     cycle i_loop
   endif
#ifdef PRT
   if(i==lond.and.lat==latd) then
     print*,' ktop r s g h c i rh ',ktopqr,ktopqs,ktopqg,ktopqh,ktopqc,ktopqi,ktoprh
     print*,' ktop lqr s g h c i flgcld ',lqr,lqs,lqg,lqh,lqc,lqi,flgcld
   endif
#endif
!-------------------------------------------------------------------------------
! initialize the local variables
!
   temp = 0.         ; hvalue = 0.       ; htotal = 0. ; source = 0.
   pcond(:) = 0.     ; pigen(:) = 0.     
   praut(:) = 0.     ; psaut(:) = 0.     ; pgaut(:) = 0.
   pidep(:) = 0.     ; psdep(:) = 0.     ; pgdep(:) = 0.
   pracw(:) = 0.     ; praci(:) = 0.     ; pracs(:) = 0.
   psacw(:) = 0.     ; psaci(:) = 0.     ; psacr(:) = 0.
   pgacw(:) = 0.     ; pgaci(:) = 0.     ; pgacr(:) = 0. ; pgacs(:) = 0.
   paacw(:) = 0.     ; piacr(:) = 0.
   psmlt(:) = 0.     ; pgmlt(:) = 0.     ; pseml(:) = 0. ; pgeml(:) = 0.
   prevp(:) = 0.     ; psevp(:) = 0.     ; pgevp(:) = 0.
   denqr(:) = 0.     ; denqs(:) = 0.     ; denqi(:) = 0.
   denqh(:) = 0.     ; denqg(:) = 0.    
   vtr(:)   = 0.     ; vts(:)   = 0.     ; vtg(:)   = 0. ; vti(:)   = 0.
   qrpath   = 0.     ; qspath   = 0.     ; qgpath   = 0. ; qipath   = 0.
   precip_r = 0.     ; precip_s = 0.     ; precip_g = 0. ; precip_i = 0.
   satrdt(:)= 0.     ; supsat(:)= 0.     ; supice(:)= 0.
   ifsat(:) = .false.; ifice(:) =.false.
   ni(:)   = 1.e3    ; mi(:)    = 0.     ; cldf(:) = 1.
   sumice(:) = 0.    ; vtmean(:) = 0. 
   n0gfac(:) = 0.    ; n0hfac(:) = 0.   ; n0sfac(:) = 0.
   tc(:)     = 0.    ; tw(:)     = 0.   ; dend(:)   = 0. ; denfac(:) = 0.
   tstepsnow = 0.    ; tstepgraupel = 0. ; tstephail = 0.
   lamdar = 0.       ; diameter = 0.     ; eacrs = 0.  ; egi =0.
   roqi0 = 0.        ; acrfac = 0.       ; ni0(:) = 0.
   qimax = 0.        ; alpha2 = 0.       ; coeres = 0.
   pfrzdtc = 0.      ; pfrzdtr = 0.
#ifdef DOUBLE
   pcact(:) = 0.     ; nsacw(:) = 0.     ;  ngacw(:) = 0.  ; naacw(:) = 0.
   niacr(:) = 0.     ; nsacr(:) = 0.     ;  ngacr(:) = 0.  ; nseml(:) = 0.
   ngeml(:) = 0.     ; nracw(:) = 0.     ;  nccol(:) = 0.  ; nrcol(:) = 0.
   ncact(:) = 0.     ; nraut(:) = 0.     ;  ncevp(:) = 0.  ; denqn(:) = 0.
   vtn(:)   = 0.     ; dc(:)  = 0.       ;  dr(:) = 0.     ; di(:)  = 0.
   qrcon(:)= 0.      ; qrconcr(:) = 0.   ;  taucon(:) = 0. 
   qnpath = 0.       ; sfac = 0.         ;  gfac = 0.
   nfrzdtc = 0.      ; nfrzdtr = 0.
#endif
#ifdef HAIL
   phdep(:)= 0.      ; phaut(:) = 0.   ; pracg(:) = 0. ; phacw(:) = 0.
   phaci(:)= 0.      ; phacr(:) = 0.   ; phacs(:) = 0. ; phacg(:) = 0.
   phmlt(:)= 0.      ; pheml(:) = 0.   ; phevp(:) = 0. ; primh(:) = 0.
   pvapg(:)= 0.      ; pvaph(:) = 0.   ; pgwet(:) = 0. ; phwet(:) = 0.
   pgaci_w(:)= 0.    ; phaci_w(:) = 0. ; qhpath = 0.   ; vth(:) = 0.
   nhacw(:) = 0.     ; nhacr(:) = 0.   ; nheml(:)= 0.  ; precip_h = 0.
#endif
#ifdef LB2017_PRAUT
   ncaut(:) = 0.     ; L_qc = .true.   ; n0c = 0.        ; n0csq = 0. 
   L_c = 0.          ; n_c = 0.        ; lamc = 0.       ; nu_c = 0
   gamma1(:) = 0.    ; qc_tot = 0.     ; qc_sub = 0.     ; qc_aut = 0.
   nc_tot = 0.       ; nc_sub = 0.     ; nc_aut = 0.     ; nr_aut = 0.
#endif
#ifdef ICE_SPHERE
   dis(:)   = 0.     ; vtis(:)   = 0.  
#endif
!-------------------------------------------------------------------------------
! compute local variables
!
   ktop = ktopmax
#ifdef DOUBLE
   do k = ktop, kts, -1
     if(qci(k,i,1) < qcmin) ncr(k,i,2) =0.
     if(qrs(k,i,1) < qrmin) ncr(k,i,3) =0.
   enddo
#endif
!
   do k = ktop, kts, -1
     tc(k) = t(k,i) - t0c
#ifdef WET_BULB
     tw(k)       = stull_wet_bulb(tc(k), rh_mul(k)*100.)  ! C
#ifdef PRT
   if(i==lond.and.lat==latd) then
     print*,' k tc tw rh',k,tc(k),tw(k),rh_mul(k)
   endif
#endif
#endif
     dend(k) = (p(k,i)/t(k,i)-den(k,i)*rv)/(rd-rv)  ! dry density
     denfac(k) = sqrt(den0/den(k,i))
   enddo
!===============================================================================
! optimizatin : A**B => exp(log(A)*(B))
!-------------------------------------------------------------------------------
! x : temperature, y: pressure, z: density
! diffus: diffusion coefficient of the water vapr ! 8.794e-5*x**1.81/y
!         diffus(x,y) = 8.794e-5 * exp(log(x)*(1.81)) / y
! viscok: kinematic viscosity(m2s-1)              ! 1.496e-6*x**1.5/(x+120.)/z
!         viscok(x,z) = 1.496e-6 * (x*sqrt(x)) /(x+120.)/z
! viscod: dynamic viscosity(kgm-1s-1)             ! 1.414e3*viscok*z
!         viscod(x,z) = 1.414e3*viscok(x,z)*z
! ab: thermodynamic terms, a and b in the denominator associated with
!             heat conduction and vapor diffusion
!         ab(a,b,c,d,e) = d*a*a/(viscod(c,d)*rv*c*c)+1./(e*diffus(c,b))
! venfac: parameter associated with the ventilation effects
!
   do k = ktop, kts, -1
     diffus(k) = (8.794e-5*exp(log(t(k,i))*(1.81)))/p(k,i)
     viscok(k) = (1.496e-6*(t(k,i)*sqrt(t(k,i))))/(t(k,i)+120.)/den(k,i)
     viscod(k) = 1.414e3*viscok(k)*den(k,i)
     ab_mul(k) = ((den(k,i)*xlv(k)*xlv(k)))/(viscod(k)*(rv*t(k,i)*t(k,i)))      &
               + 1./(qsat_mul(k)*diffus(k))
     ab_ice(k) = ((den(k,i)*xls*xls))/(viscod(k)*(rv*t(k,i)*t(k,i)))            &
               + 1./(qsat_ice(k)*diffus(k))
     venfac(k) = (exp(.3333333*log(viscok(k)/(diffus(k))))                      &
               * sqrt(sqrt(den0/(den(k,i)))))/sqrt(viscok(k))
   enddo
#ifdef SIZE_LIMIT1
!===============================================================================
!
! adjustment of nc and nr according to avaialble size of clouds and rain
!
!===============================================================================
   call adjust_number_concent(ktopqr,kdim,qrs(:,i,1),ncr(:,i,3),den(:,i),      &
        pidnr,drcoeff,qrmin,nrmin,nrmax,di1000,drmin,drmax)
   call adjust_number_concent(ktopqc,kdim,qci(:,i,1),ncr(:,i,2),den(:,i),      &
        pidnc,one_1,qcmin,ncmin,ncmax,di15,  dcmin, dcmax)
#endif
#ifdef PARTIAL
!-------------------------------------------------------------------------------
! obtain in-cloud properties
!
   ktop = max(ktopqc,ktopqi)
   call cldf_diag(kdim,t(:,i),p(:,i),q(:,i),qci(:,i,1),qci(:,i,2),             &
        dxmeter(i),cldf(:),ktop)
!
   do k = ktop, kts, -1
       ! change condensate variables to in-cloud variables
       if(cldf(k)>0.) then
         qci(k,i,1) = qci(k,i,1) / cldf(k)
         qci(k,i,2) = qci(k,i,2) / cldf(k)
         qrs(k,i,1) = qrs(k,i,1) / cldf(k)
         qrs(k,i,2) = qrs(k,i,2) / cldf(k)
         qrs(k,i,3) = qrs(k,i,3) / cldf(k)
#ifdef HAIL
         qrs(k,i,4) = qrs(k,i,4) / cldf(k)
#endif
       endif
   enddo
#endif
!-------------------------------------------------------------------------------
! update the slope parameters for microphysics computation
!
   ktop = ktopqr
   call slope_rain(kdim,ktop,qrs(:,i,1),den(:,i),denfac(:),t(:,i),             &
#ifdef DOUBLE
        ncr(:,i,3),                                                            &
#endif
        rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:),t0c,den0)
   ktop = ktopqs
   call slope_snow(kdim,ktop,qrs(:,i,2),den(:,i),denfac(:),t(:,i),             &
        rslope(:,2),rslopeb(:,2),rslope2(:,2),rslope3(:,2),vts(:),t0c,den0)
   ktop = ktopqg
   call slope_graupel(kdim,ktop,qrs(:,i,3),den(:,i),denfac(:),t(:,i),          &
        rslope(:,3),rslopeb(:,3),rslope2(:,3),rslope3(:,3),vtg(:),t0c,den0)
#ifdef HAIL
   ktop = ktopqh
   call slope_hail(kdim,ktop,qrs(:,i,4),den(:,i),denfac(:),t(:,i),             &
        rslope(:,4),rslopeb(:,4),rslope2(:,4),rslope3(:,4),vth(:),t0c,den0)
#endif
#ifdef DOUBLE
   ktop = ktopqc
   call slope_cloud(kdim,ktop,qci(:,i,1),ncr(:,i,2),den(:,i),denfac(:),        &
        t(:,i),qcmin,rslopec(:),rslopec2(:),rslopec3(:),t0c,den0)
#endif
!===============================================================================
!
! melting/freezing
!
!===============================================================================
   ktop = ktopmax
   do k = kts, ktop
     if(t(k,i)<t0c) ifice(k) = .true.
#ifdef NEW_N0
     hvalue = exp(log((max(qrs(k,i,2),q_ref)/q_ref))*(beta))
     n0sfac(k) = exp(-alpha_qs*tc(k))*hvalue
     n0sfac(k) = min(max(n0sfac(k),n0smin/n0s),n0smax/n0s)
     hvalue = exp(log((max(qrs(k,i,3),q_ref)/q_ref))*(beta))
     n0gfac(k) = exp(-alpha_qg*tc(k))*hvalue
     n0gfac(k) = min(max(n0gfac(k),n0gmin/n0g),n0gmax/n0g)
     hvalue = exp(log((max(qrs(k,i,4),q_ref)/q_ref))*(beta))
     n0hfac(k) = exp(-alpha_qh*tc(k))*hvalue
     n0hfac(k) = min(max(n0hfac(k),n0hmin/n0h),n0hmax/n0h)
#else
     n0sfac(k) = max(min(exp(-alpha*tc(k)),n0smax/n0s),1.)
#endif
   enddo
!-------------------------------------------------------------------------------
! psmlt: melting of snow [hl a33]
!       (t>t0: s->r)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
#ifdef WET_BULB
       if(tw(k)>0. .and. qrs(k,i,2)>0.) then
         coeres = rslope2(k,2)*sqrt(rslope(k,2)*rslopeb(k,2))
         psmlt(k) = viscod(k)/xlf0*tw(k)*pi/2.*n0sfac(k)                       &
#else
       if(.not.ifice(k) .and. qrs(k,i,2)>0.) then
         coeres = rslope2(k,2)*sqrt(rslope(k,2)*rslopeb(k,2))
         psmlt(k) = viscod(k)/xlf0*(t(k,i)-t0c)*pi/2.*n0sfac(k)                &
#endif
                  *(precs1*rslope2(k,2)+precs2*venfac(k)*coeres)/den(k,i)
         psmlt(k) = max(min(psmlt(k)*dtcld,qrs(k,i,2)),0.)
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! nsmlt: melting of snow [LH A27]
!       (T>T0: ->NR)
         if(qrs(k,i,2)>qrmin) then
           sfac = rslope(k,2)*n0s*n0sfac(k)/qrs(k,i,2)
           ncr(k,i,3) = min(ncr(k,i,3) + sfac*psmlt(k), nrmax)
         endif
#endif
         qrs(k,i,2) = qrs(k,i,2) - psmlt(k)
         qrs(k,i,1) = qrs(k,i,1) + psmlt(k)
         t(k,i) = t(k,i) - xlf0/cpm(k)*psmlt(k)*cldf(k)
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgmlt: melting of graupel [hl a23]
!       (t>t0: g->r)
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
#ifdef WET_BULB
       if(tw(k)>0. .and. qrs(k,i,3)>0.) then
         coeres = rslope2(k,3)*sqrt(rslope(k,3)*rslopeb(k,3))
#ifdef NEW_N0
         pgmlt(k) = viscod(k)/xlf0*tw(k)*n0gfac(k)*(precg1*rslope2(k,3)         &
#else
         pgmlt(k) = viscod(k)/xlf0*tw(k)*(precg1*rslope2(k,3)                   &
#endif
#else
       if(.not.ifice(k) .and. qrs(k,i,3)>0.) then
         coeres = rslope2(k,3)*sqrt(rslope(k,3)*rslopeb(k,3))
#ifdef NEW_N0
         pgmlt(k) = viscod(k)/xlf0*(t(k,i)-t0c)*n0gfac(k)*(precg1*rslope2(k,3) &
#else
         pgmlt(k) = viscod(k)/xlf0*(t(k,i)-t0c)*(precg1*rslope2(k,3)           &
#endif
#endif
                  + precg2*venfac(k)*coeres)/den(k,i)
         pgmlt(k) = max(min(pgmlt(k)*dtcld,qrs(k,i,3)),0.)
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! ngmlt: melting of graupel [lh a28]
!       (t>t0: ->nr)
         if(qrs(k,i,3)>qrmin) then
#ifdef NEW_N0
           gfac = rslope(k,3)*n0g*n0gfac(k)/qrs(k,i,3)
#else
           gfac = rslope(k,3)*n0g/qrs(k,i,3)
#endif
           ncr(k,i,3) = min(ncr(k,i,3) + gfac*pgmlt(k), nrmax)
         endif
#endif
         qrs(k,i,3) = qrs(k,i,3) - pgmlt(k)
         qrs(k,i,1) = qrs(k,i,1) + pgmlt(k)
         t(k,i) = t(k,i) - xlf0/cpm(k)*pgmlt(k)*cldf(k)
       endif
     enddo
   endif
#ifdef HAIL
!-------------------------------------------------------------------------------
! phmlt: melting of hail [bht a22]
!       (t>t0: qh->qr)
   if(lqh) then
     ktop = ktopqh
     do k = ktop, kts, -1
#ifdef WET_BULB
       if(tw(k)>0. .and. qrs(k,i,4)>0.) then
         coeres = rslope2(k,4)*sqrt(rslope(k,4)*rslopeb(k,4))
#ifdef NEW_N0
         phmlt(k) = viscod(k)/xlf0*tw(k)*n0hfac(k)*(prech1*rslope2(k,4)        &
#else
         phmlt(k) = viscod(k)/xlf0*tw(k)*(prech1*rslope2(k,4)                  &
#endif
#else
       if(.not.ifice(k) .and. qrs(k,i,4)>0.) then
         coeres = rslope2(k,4)*sqrt(rslope(k,4)*rslopeb(k,4))
#ifdef NEW_N0
         phmlt(k) = viscod(k)/xlf0*(t(k,i)-t0c)*n0hfac(k)*(prech1*rslope2(k,4) &
#else
         phmlt(k) = viscod(k)/xlf0*(t(k,i)-t0c)*(prech1*rslope2(k,4)           &
#endif
#endif
                  + prech2*venfac(k)*coeres)/den(k,i)
         phmlt(k) = max(min(phmlt(k)*dtcld,qrs(k,i,4)),0.)
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! nhmlt: melting of hail [lh a28]
!       (t>t0: ->nr)
         if(qrs(k,i,4)>qrmin) then
#ifdef NEW_N0
           gfac = rslope(k,4)*n0h*n0hfac(k)/qrs(k,i,4)
#else
           gfac = rslope(k,4)*n0h/qrs(k,i,4)
#endif
           ncr(k,i,3) = min(ncr(k,i,3) + gfac*phmlt(k), nrmax)
         endif
#endif
         qrs(k,i,4) = qrs(k,i,4) - phmlt(k)
         qrs(k,i,1) = qrs(k,i,1) + phmlt(k)
         t(k,i) = t(k,i) - xlf0/cpm(k)*phmlt(k)*cldf(k)
       endif
     enddo
   endif
#endif
!-------------------------------------------------------------------------------
! pimlt: instantaneous melting of cloud ice [hl a47]
!       (t>t0: i->c)
   if(lqi) then
     ktop = ktopqi
     do k = ktop, kts, -1
#ifdef WET_BULB
       if(tw(k)>0. .and. qci(k,i,2)>0.) then
#else
       if(.not.ifice(k) .and. qci(k,i,2)>0.) then
#endif
         qci(k,i,1) = qci(k,i,1) + qci(k,i,2)
         t(k,i) = t(k,i) - xlf0/cpm(k)*qci(k,i,2)*cldf(k)
         qci(k,i,2) = 0.
#ifdef DOUBLE
!---------------------------------------------------------------
! nimlt: instantaneous melting of cloud ice  [lh a18]
!        (t>t0: ->nc)
         hvalue = (den(k,i)*max(qci(k,i,2),qcmin))
         hvalue = sqrt(sqrt(hvalue*hvalue*hvalue))
         ni(k) = min(max(5.38e7*hvalue,1.e3),1.e6)
         ncr(k,i,2) = min(ncr(k,i,2) + ni(k),ncmax)
#endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pihmf: homogeneous freezing of cloud water below -40c [hl a45]
!        (t<-40c: c->i)
   if(lqc) then
     ktop = ktopqc
     do k = ktop, kts, -1
       if(tc(k)<-40. .and. qci(k,i,1)>0.) then
         qci(k,i,2) = qci(k,i,2) + qci(k,i,1)
         t(k,i) = t(k,i) + xlf(k)/cpm(k)*qci(k,i,1)*cldf(k)
         qci(k,i,1) = 0.
#ifdef DOUBLE
!---------------------------------------------------------------
! nihmf: homogeneous  freezing of cloud water below -40c [lh a17]
!        (t<-40c: nc-> )
         if(ncr(k,i,2)>0.) ncr(k,i,2) = 0.
#endif
       endif
!-------------------------------------------------------------------------------
! pihtf: heterogeneous freezing of cloud water [hl a44]
!        (t0>t>-40c: c->i)
       if(ifice(k) .and. qci(k,i,1)>qcmin) then
#ifdef DOUBLE
         hvalue = max(tc(k),-70.)
         pfrzdtc = min(pisq*pfrz1*(exp(-pfrz2*hvalue)-1.)*denr/den(k,i)        &
                   *ncr(k,i,2)*rslopec3(k)*rslopec3(k)/18.*dtcld               &
                   ,qci(k,i,1))
!---------------------------------------------------------------
! nihtf: heterogeneous  freezing of cloud water [LH A16]
!         (T0>T>-40C: NC->)
         if(ncr(k,i,2)>ncmin) then
           nfrzdtc = min(pi*pfrz1*(exp(-pfrz2*hvalue)-1.)*ncr(k,i,2)           &
                     *rslopec3(k)/6.*dtcld,ncr(k,i,2))
           ncr(k,i,2) = max(ncr(k,i,2) - nfrzdtc,0.)
         endif
#else
         hvalue = max(tc(k),-50.)
         pfrzdtc = min(pfrz1*(exp(-pfrz2*hvalue)-1.)                           &
                   *den(k,i)/denr/xncr*qci(k,i,1)*qci(k,i,1)*dtcld,qci(k,i,1))
#endif
         qci(k,i,1) = qci(k,i,1) - pfrzdtc
         qci(k,i,2) = qci(k,i,2) + pfrzdtc
         t(k,i) = t(k,i) + xlf(k)/cpm(k)*pfrzdtc*cldf(k)
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgfrz: freezing of rain water [hl a20]
!        (t<t0, r->g)
   if(lqr) then
     ktop = ktopqr
     do k = ktop, kts, -1
       if(ifice(k) .and. qrs(k,i,1)>0.) then
#ifdef DOUBLE
         hvalue = max(tc(k), -70.)
         pfrzdtr = min(140.*pisq*pfrz1*ncr(k,i,3)*denr/den(k,i)                &
                  *(exp(-pfrz2*hvalue)-1.)*rslope3(k,1)*rslope3(k,1)           &
                  *dtcld,qrs(k,i,1))
!-------------------------------------------------------------------------------
! ngfrz: freezing of rain water [lh a26]
!        (T<T0, NR-> )
         if(ncr(k,i,3)>nrmin) then
           nfrzdtr = min(4.*pi*pfrz1*ncr(k,i,3)*(exp(-pfrz2*hvalue)-1.)        &
                    *rslope3(k,1)*dtcld, ncr(k,i,3))
           ncr(k,i,3) = max(ncr(k,i,3) - nfrzdtr, 0.)
         endif
#else
         hvalue = max(tc(k), -50.)
         hvalue = rslope3(k,1)
         hvalue = hvalue*hvalue*rslope(k,1)
         pfrzdtr = min(20.*pisq*pfrz1*n0r*denr/den(k,i)                        &
                     *(exp(-pfrz2*hvalue)-1.)*hvalue*dtcld,qrs(k,i,1))
#endif
         qrs(k,i,1) = qrs(k,i,1) - pfrzdtr
         qrs(k,i,3) = qrs(k,i,3) + pfrzdtr
         t(k,i) = t(k,i) + xlf(k)/cpm(k)*pfrzdtr*cldf(k)
       endif
     enddo
   endif
#ifdef PARTIAL
!-------------------------------------------------------------------------------
! restore grid-mean variables
!
   ktop = max(ktopqc,ktopqi)
   do k = ktop, kts, -1
       ! change in-cloud condensate variables to grid-mean variables
     if(cldf(k)>0.) then
       qci(k,i,1) = qci(k,i,1) * cldf(k)
       qci(k,i,2) = qci(k,i,2) * cldf(k)
       qrs(k,i,1) = qrs(k,i,1) * cldf(k)
       qrs(k,i,2) = qrs(k,i,2) * cldf(k)
       qrs(k,i,3) = qrs(k,i,3) * cldf(k)
#ifdef HAIL
       qrs(k,i,4) = qrs(k,i,4) * cldf(k)
#endif
     endif
   enddo
#endif
#ifdef DOUBLE
#ifdef PRT
!
   ktop = ktopmax
   do k = ktop, kts, -1
       hvalue = 1.e-6
        if(i==lond.and. lat==latd.and.k==kts) print *, 'ncr 1 2 3'
        if(i==lond.and. lat==latd) print 113, k,ncr(k,i,1)*hvalue, ncr(k,i,2)*hvalue, &
              ncr(k,i,3)
113   format(1x,i5,3f8.2)
   enddo
#endif
#endif
!-------------------------------------------------------------------------------
! re-define cloud top for numerical efficiencies
!
   ktop = ktopini
   do k = kts, ktop
     qsat_mul(k) = fsvp_water(t(k,i),p(k,i))
     qsat_mul(k) = ep2 * qsat_mul(k) / (p(k,i) - qsat_mul(k))
     qsat_mul(k) = max(qsat_mul(k),qmin)
     rh_mul(k)   = max(q(k,i) / qsat_mul(k),qmin)
     qsat_ice(k) = fsvp(t(k,i),p(k,i))
     qsat_ice(k) = ep2 * qsat_ice(k) / (p(k,i) - qsat_ice(k))
     qsat_ice(k) = max(qsat_ice(k),qmin)
     rh_ice(k)   = max(q(k,i) / qsat_ice(k),qmin)
   enddo
!
   lqc = .false.
   lqi = .false.
   lqr = .false.
   lqs = .false.
   lqg = .false.
#ifdef HAIL
   lqh = .false.
#endif
   call find_cloud_top(1,kdim,ktopini,qci(:,i,1),zero_0,ktopqc)
   call find_cloud_top(1,kdim,ktopini,qci(:,i,2),zero_0,ktopqi)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,1),zero_0,ktopqr)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,2),zero_0,ktopqs)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,3),zero_0,ktopqg)
#ifdef HAIL
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,4),zero_0,ktopqh)
#endif
   call find_cloud_top(1,kdim,ktopini,rh_ice(:), one_1,ktoprh)
   if(ktopqc>0.0) lqc = .true.
   if(ktopqi>0.0) lqi = .true.
   if(ktopqr>0.0) lqr = .true.
   if(ktopqs>0.0) lqs = .true.
   if(ktopqg>0.0) lqg = .true.
#ifdef HAIL
   if(ktopqh>0.0) lqh = .true.
   ktopmax = max(ktopqc,ktopqi,ktopqr,ktopqs,ktopqg,ktopqh,ktoprh)
#else
   ktopmax = max(ktopqc,ktopqi,ktopqr,ktopqs,ktopqg,ktoprh)
#endif
#ifdef PRT
   if(i==lond.and.lat==latd) then
     print*,' ktop r s g h c i rh ',ktopqr,ktopqs,ktopqg,ktopqh,ktopqc,ktopqi,ktoprh
   endif
#endif
#ifdef DOUBLE
#ifdef SIZE_LIMIT2
!===============================================================================
!
! adjustment of nc and nr according to avaialble size of clouds and rain
!
!===============================================================================
   call adjust_number_concent(ktopqr,kdim,qrs(:,i,1),ncr(:,i,3),den(:,i),      &
        pidnr,drcoeff,qrmin,nrmin,nrmax,di1000,drmin,drmax)
   call adjust_number_concent(ktopqc,kdim,qci(:,i,1),ncr(:,i,2),den(:,i),      &
        pidnc,one_1,qcmin,ncmin,ncmax,di15,  dcmin, dcmax)
#endif
#endif
#ifdef PARTIAL
!-------------------------------------------------------------------------------
! obtain in-cloud properties
!
   call cldf_diag(kdim,t(:,i),p(:,i),q(:,i),qci(:,i,1),qci(:,i,2),             &
        dxmeter(i),cldf(:),ktop)
   do k = ktop, kts, -1
       ! change condensate variables to in-cloud variables
     if(cldf(k)>0.) then
       qci(k,i,1) = qci(k,i,1) / cldf(k)
       qci(k,i,2) = qci(k,i,2) / cldf(k)
       qrs(k,i,1) = qrs(k,i,1) / cldf(k)
       qrs(k,i,2) = qrs(k,i,2) / cldf(k)
       qrs(k,i,3) = qrs(k,i,3) / cldf(k)
#ifdef HAIL
       qrs(k,i,4) = qrs(k,i,4) / cldf(k)
#endif
     endif
   enddo
#endif
!-------------------------------------------------------------------------------
! update the slope parameters and microphysical parameters
!
   vtr(:) = 0. ; vts(:) = 0. ; vtg(:) = 0. ; vth(:) = 0. ; vti(:) = 0.
   ktop = ktopqr
   call slope_rain(kdim,ktop,qrs(:,i,1),den(:,i),denfac(:),t(:,i),             &
#ifdef DOUBLE
        ncr(:,i,3),                                                            &
#endif
        rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:),t0c,den0)
   ktop = ktopqs
   call slope_snow(kdim,ktop,qrs(:,i,2),den(:,i),denfac(:),t(:,i),             &
        rslope(:,2),rslopeb(:,2),rslope2(:,2),rslope3(:,2),vts(:),t0c,den0)
   ktop = ktopqg
   call slope_graupel(kdim,ktop,qrs(:,i,3),den(:,i),denfac(:),t(:,i),          &
        rslope(:,3),rslopeb(:,3),rslope2(:,3),rslope3(:,3),vtg(:),t0c,den0)
#ifdef HAIL
   ktop = ktopqh
   call slope_hail(kdim,ktop,qrs(:,i,4),den(:,i),denfac(:),t(:,i),             &
        rslope(:,4),rslopeb(:,4),rslope2(:,4),rslope3(:,4),vth(:),t0c,den0)
#endif
   ktop = ktopqi
   do k = ktop, kts, -1
     if(qci(k,i,2)>qcmin) then
       hvalue   = (den(k,i)*max(qci(k,i,2),qcmin))
       hvalue   = sqrt(sqrt(hvalue*hvalue*hvalue))
       ni(k)  = min(max(5.38e7*hvalue,1.e3),1.e6)
       mi(k)     = den(k,i)*qci(k,i,2)/ni(k)
       di(k)  = max(min(exp(log((mi(k)/cxmi))*(1./dxmi)),dimax), qmin)
       vti(k) = avti*exp(log(di(k))*(bvti))
#ifdef ICE_SPHERE
       hvalue = min(max((tc(k) - t2_sphere)/(t1_sphere-t2_sphere),0.),1.0)
       dis(k) = max(min(exp(log((mi(k)/cxmis))*(1./dxmis)),dimax), qmin)
       vtis(k)= avtis*exp(log(dis(k))*(bvtis))
       di(k)  = di(k)  * (1.-hvalue) + dis(k)  * hvalue
       vti(k) = vti(k) * (1.-hvalue) + vtis(k) * hvalue
#ifdef PRT
        if(i==lond.and. lat==latd) print 127, k,tc(k),di(k)*1.e6,dis(k)*1.e6,vti(k),vtis(k) 
127   format(1x,'vti',i5,5e15.5)
#endif
#endif
     endif
   enddo
!
#if defined (HAIL) && defined (SEDI_MEAN)
   ktop = max(ktopqs,ktopqg,ktopqh)
   do k = ktop, kts, -1
     sumice(k) = max( (qrs(k,i,2)+qrs(k,i,3)+qrs(k,i,4)), qmin)
     if(sumice(k)>qmin) then
       vtmean(k) = (vts(k)*qrs(k,i,2)+vtg(k)*qrs(k,i,3)                        &
                 +vth(k)*qrs(k,i,4))/sumice(k)
     else
       vtmean(k) = 0.
     endif
   enddo
#else
   ktop = max(ktopqs,ktopqg)
   do k = ktop, kts, -1
     sumice(k) = max( (qrs(k,i,2)+qrs(k,i,3)), qmin)
     if(sumice(k)>qmin) then
       vtmean(k) = (vts(k)*qrs(k,i,2)+vtg(k)*qrs(k,i,3))/sumice(k)
     else
       vtmean(k) = 0.
     endif
   enddo
#endif
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! compute the mean-volume drop diameter for raindrop distribution    [lh a10]
!
   ktop = ktopqr
   do k = ktop, kts, -1
     dr(k) = rslope(k,1)*drcoeff
   enddo
!-------------------------------------------------------------------------------
! compute the mean-volume drop diameter for cloud-droplet distribution [lh a7]
!
   ktop = ktopqc
   call slope_cloud(kdim,ktop,qci(:,i,1),ncr(:,i,2),den(:,i),denfac(:),      &
        t(:,i),qcmin,rslopec(:),rslopec2(:),rslopec3(:),t0c,den0)
   do k = ktop, kts, -1
     dc(k) = rslopec(k)
   enddo
#endif
!===============================================================================
!
! warm rain processes
!
!===============================================================================
#ifdef DOUBLE
   ktop = max(ktopqc,ktopqr)
   do k = ktop, kts, -1
     qrcon(k)  = 2.7e-2*den(k,i)*qci(k,i,1)*(1.e20/16.*rslopec2(k)             &
               *rslopec2(k)-0.4)
     qrconcr(k) = max(1.2*qrcon(k), qrmin)
   enddo
#endif
   if(lqc) then
     ktop = ktopqc
     do k = ktop, kts, -1
#ifdef DOUBLE
#ifdef LB2017_PRAUT
!-------------------------------------------------------------------------------
! praut: auto conversion rate from cloud to rain  [lb2017, eqs. 15~19]
!        (qc->qr, nc-> nr)
       if(qci(k,i,1) > qcmin) then
         L_qc = .true.
         L_c = qci(k,i,1)*den(k,i)
         n_c = max(min(ncr(k,i,2),1.e12),2.)
         nu_c = min(nint(1.e9/n_c)+2,15) 
         lamc = max(min((n_c*pidnc*8.d0*ccg(nu_c)*ocg1(nu_c)/L_c)**obmr1,      &
                lamdacmax),lamdacmin)
       else
         L_qc = .false.
         L_c = qcmin
         n_c = 2.
       endif
       if(L_qc.and.L_c>1.e-5) then
!-----------------------------------------
! for mass concentraion
!-----------------------------------------
         gamma1(0) = 1.0d0/( lamc)
         do n=1, nu_c
           gamma1(n) = gamma1(n-1)*(dble(n)/( lamc))
         enddo
         n0c = n_c/gamma1(nu_c)
         n0csq = n0c * n0c
         qc_tot = funct_aut_qc_tot(lamc,nu_c)
         qc_sub = funct_aut_qc_sub(lamc,nu_c)
!
         qc_aut =(qc_tot-alpha_aut*qc_sub)*(pidnc*8.0)*pi*vt0_c*kc_aut*n0csq
!-----------------------------------------
! for number concentraion
!-----------------------------------------
         nc_tot = funct_aut_nc_tot(lamc,nu_c)
         nc_sub = funct_aut_nc_sub(lamc,nu_c)
!
         nc_aut = (2.0*nc_tot-alpha_aut*nc_sub)*pi*vt0_c*kc_aut*n0csq
         nr_aut = (    nc_tot-alpha_aut*nc_sub)*pi*vt0_c*kc_aut*n0csq
!
         praut(k) = max(min(real(qc_aut/den(k,i)),qci(k,i,1)*rdtcld),0.)
         ncaut(k) = max(min(real(nc_aut),ncr(k,i,2)*rdtcld),0.)
         nraut(k) = min(real(nr_aut),praut(k)/mrmin)
#ifdef PRT
        if(i==lond.and. lat==latd) print 117, k,praut(k), nraut(k),qc_aut 
117   format(1x,'praut',i5,3e20.10)
#endif
       end if
#endif
#ifdef LH2010_PRAUT
!-------------------------------------------------------------------------------
! praut: auto conversion rate from cloud to rain  [lh 9]
!        (qc->qr, nc-> nr)
       if(dc(k)>di15) then
         taucon(k) = 3.7/den(k,i)/qci(k,i,1)/(0.5e6*rslopec(k)-7.5)
         taucon(k) = max(taucon(k), qmin)
         praut(k) = qrcon(k)/(taucon(k)*den(k,i))
         praut(k) = min(max(praut(k),0.),qci(k,i,1)*rdtcld)
!-------------------------------------------------------------------------------
! nraut: auto conversion rate from cloud to rain [lh a6] [cp 18 & 19]
!        (nc->nr)
         nraut(k) = dr_embryo*den(k,i)*praut(k)
         if(qrs(k,i,1)>qrconcr(k)) then
           nraut(k) = ncr(k,i,3)/qrs(k,i,1)*praut(k)
           nraut(k) = min(nraut(k),ncr(k,i,2)*rdtcld)
         endif
#ifdef PRT
         if(i==lond.and. lat==latd) print 117, k,praut(k), nraut(k),taucon(k) 
117   format(1x,'praut',i5,3e20.10)
#endif
       endif
#endif
#ifdef HDC2004_PRAUT
!-------------------------------------------------------------------------------
! praut: auto conversion rate from cloud to rain  [hdc 16]
!        (qc->qr, nc-> nr)
       if(qci(k,i,1)>qc0(i).and.ncr(k,i,2)>ncmin) then
         praut(k) = qck1*exp(log(qci(k,i,1))*(7./3.))*exp(log(ncr(k,i,2))      &
                  *(-1./3.))
         praut(k) = min(praut(k),qci(k,i,1)*rdtcld)
!-------------------------------------------------------------------------------
! nraut: auto conversion rate from cloud to rain [lh a6] [cp 18 & 19]
!        (nc->nr)
         if(qrs(k,i,1)>qrconcr(k)) then
           nraut(k) = ncr(k,i,3)/qrs(k,i,1)*praut(k)
           nraut(k) = min(nraut(k),ncr(k,i,2)*rdtcld)
         endif
       endif
#endif
#else /* single moment */
!-------------------------------------------------------------------------------
! praut: auto conversion rate from cloud to rain [hdc 16]
!        (qc->qr)
       if(qci(k,i,1)>qc_land) then
         praut(k) = qck1*exp(log(qci(k,i,1))*((7./3.)))
         praut(k) = min(praut(k),qci(k,i,1)*rdtcld)
       endif
#endif
     enddo
   endif
!
#ifdef DOUBLE
   if(lqc) then
     ktop = ktopqc
     do k = ktop, kts, -1
!-------------------------------------------------------------------------------
! nccol: self collection of cloud water             [lh a8] [cp 24 & 25]
!        (nc->)
       if(dc(k)>=di100) then
         nccol(k) = ncrk1*ncr(k,i,2)*ncr(k,i,2)*rslopec3(k)
       else
         nccol(k) = 2.*ncrk2*ncr(k,i,2)*ncr(k,i,2)*rslopec3(k)*rslopec3(k)
       endif
     enddo
   endif
#endif
!
   if(lqr) then
     ktop = ktopqr
     do k = ktop, kts, -1
       supsat(k) = max(q(k,i),qmin)-qsat_mul(k)
       satrdt(k) = supsat(k)*rdtcld
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! pracw: accretion of cloud water by rain     [lh 10] [cp 22 & 23]
!        (qc->qr)
! nracw: accretion of cloud water by rain     [lh a9]
!        (nc->)
       if(qrs(k,i,1)>=qrconcr(k)) then
         if(dr(k)>=di100) then
           nracw(k) = min(ncrk1*ncr(k,i,2)*ncr(k,i,3)*(rslopec3(k)             &
                      + 24.*rslope3(k,1)),ncr(k,i,2)*rdtcld)
           pracw(k) = min(pi/6.*(denr/den(k,i))*ncrk1*ncr(k,i,2)               &
                      *ncr(k,i,3)*rslopec3(k)*(2.*rslopec3(k)                  &
                      + 24.*rslope3(k,1)),qci(k,i,1)*rdtcld)
         else
           nracw(k) = min(ncrk2*ncr(k,i,2)*ncr(k,i,3)*(2.*rslopec3(k)          &
                      *rslopec3(k)+5040.*rslope3(k,1)                          &
                      *rslope3(k,1)),ncr(k,i,2)*rdtcld)
           pracw(k) = min(pi/6.*(denr/den(k,i))*ncrk2*ncr(k,i,2)               &
                      *ncr(k,i,3)*rslopec3(k)*(6.*rslopec3(k)                  &
                      *rslopec3(k)+5040.*rslope3(k,1)*rslope3(k,1))            &
                      ,qci(k,i,1)*rdtcld)
         endif
       endif
!-------------------------------------------------------------------------------
! nrcol: self collection of rain-drops and break-up [lh a21] [cp 24 & 25]
!        (nr->)
       if(qrs(k,i,1)>=qrconcr(k)) then
         if(dr(k)<di100) then
           nrcol(k) = 5040.*ncrk2*ncr(k,i,3)*ncr(k,i,3)*rslope3(k,1)           &
                     *rslope3(k,1)
         elseif(dr(k)>=di100 .and. dr(k)<di600) then
           nrcol(k) = 24.*ncrk1*ncr(k,i,3)*ncr(k,i,3)*rslope3(k,1)
         elseif(dr(k)>=di600 .and. dr(k)<di2000) then
           hvalue = -2.5e3*(dr(k)-di600)
           nrcol(k) = 24.*exp(hvalue)*ncrk1*ncr(k,i,3)*ncr(k,i,3)              &
                      *rslope3(k,1)
         else
           nrcol(k) = 0.
         endif
       endif
#else
!-------------------------------------------------------------------------------
! pracw: accretion of cloud water by rain [hl a40]
!        (c->r)
       if(qrs(k,i,1)>qrmin .and. qci(k,i,1)>qcmin) then
         pracw(k) = min(pacrr*rslope3(k,1)*rslopeb(k,1)                        &
                      *qci(k,i,1)*denfac(k),qci(k,i,1)*rdtcld)
       endif
#endif
!-------------------------------------------------------------------------------
! prevp: evaporation/condensation rate of rain [hdc 14]
!        (v->r or r->v)
       if(qrs(k,i,1)>0.) then
#ifdef DOUBLE
         coeres = rslope(k,1)*sqrt(rslope(k,1)*rslopeb(k,1))
         prevp(k) = (rh_mul(k)-1.)*ncr(k,i,3)*(precr1*rslope(k,1)              &
#else
         coeres = rslope2(k,1)*sqrt(rslope(k,1)*rslopeb(k,1))
         prevp(k) = (rh_mul(k)-1.)*(precr1*rslope2(k,1)                        &
#endif
                   +precr2*venfac(k)*coeres)/ab_mul(k)
         if(prevp(k)<=0.) then
           prevp(k) = max(prevp(k),-qrs(k,i,1)*rdtcld)
           prevp(k) = max(prevp(k),satrdt(k)*.5)
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! nrevp: evaporation/condensation rate of rain   [lh a14]
!        (nr->nccn)
#ifdef EVAP_PARTIAL
           hvalue = max(min(-prevp(k)*ncr(k,i,3)/qrs(k,i,1)*dtcld,ncr(k,i,3)),0.)
           ncr(k,i,1) = min(ncr(k,i,1) + hvalue,ccnmax)
           ncr(k,i,3) = max(ncr(k,i,3) - hvalue,0.)
#else
           if(prevp(k)==-qrs(k,i,1)*rdtcld) then
             ncr(k,i,1) = min(ncr(k,i,1) + ncr(k,i,3), ccnmax)
             ncr(k,i,3) = 0.
           endif
#endif
#endif
         else
           prevp(k) = min(prevp(k),satrdt(k)*.5)
         endif
       endif
     enddo
   endif
!===============================================================================
!
! cold rain processes
!
!===============================================================================
   ktop = ktopmax
   ifice(:) = .false.
   do k = ktop, kts, -1
     if(t(k,i)<t0c) ifice(k) = .true.
#ifdef NEW_N0
     hvalue = exp(log((max(qrs(k,i,2),q_ref)/q_ref))*(beta))
     n0sfac(k) = exp(-alpha_qs*tc(k))*hvalue
     n0sfac(k) = min(max(n0sfac(k),n0smin/n0s),n0smax/n0s)
     hvalue = exp(log((max(qrs(k,i,3),q_ref)/q_ref))*(beta))
     n0gfac(k) = exp(-alpha_qg*tc(k))*hvalue
     n0gfac(k) = min(max(n0gfac(k),n0gmin/n0g),n0gmax/n0g)
     hvalue = exp(log((max(qrs(k,i,4),q_ref)/q_ref))*(beta))
     n0hfac(k) = exp(-alpha_qh*tc(k))*hvalue
     n0hfac(k) = min(max(n0hfac(k),n0hmin/n0h),n0hmax/n0h)
#else
     n0sfac(k) = max(min(exp(-alpha*tc(k)),n0smax/n0s),1.)
#endif
     supsat(k) = max(q(k,i),qmin)-qsat_ice(k)
     satrdt(k) = supsat(k)*rdtcld
   enddo
!
   if(lqr .and. lqi) then
     ktop = min(ktopqr,ktopqi)
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,1)>qrmin .and. qci(k,i,2)>qcmin) then
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! praci: accretion of cloud ice by rain [hl a15] [lfo 25]
!        (t<t0: qi->qr)
           acrfac = 6.*rslope2(k,1)+4.*di(k)*rslope(k,1) + di(k)**2
           praci(k) = pi*qci(k,i,2)*ncr(k,i,3)*abs(vtr(k)-vti(k))*acrfac/4.
             ! reduce collection efficiency (suggested by B. Wilt)
           praci(k) = praci(k)*min(max(0.0,qrs(k,i,1)/qci(k,i,2)),1.)**2
           praci(k) = min(praci(k),qci(k,i,2)*rdtcld)
!-------------------------------------------------------------------------------
! piacr: accretion of rain by cloud ice [hl A19] [lfo 26]
!        (t<t0: qr->qs OR qr->qg)
           piacr(k) = pisq*avtr*ncr(k,i,3)*denr*ni(k)*denfac(k)                &
                       *g7pbr*rslope3(k,1)*rslope2(k,1)*rslopeb(k,1)           &
                       /24./den(k,i)
             ! reduce collection efficiency (suggested by B. Wilt)
           piacr(k) = piacr(k)*min(max(0.0,qci(k,i,2)/qrs(k,i,1)),1.)**2
           piacr(k) = min(piacr(k),qrs(k,i,1)*rdtcld)
!-------------------------------------------------------------------------------
! niacr: accretion of rain by cloud ice  [lh a25]
!        (t<t0: nr->)
           if(ncr(k,i,3)>nrmin) then
             niacr(k) = pi*avtr*ncr(k,i,3)*ni(k)*denfac(k)*g4pbr               &
                         *rslope2(k,1)*rslopeb(k,1)/4.
             ! reduce collection efficiency (suggested by B. Wilt)
             niacr(k) = niacr(k)*min(max(0.0,qci(k,i,2)/qrs(k,i,1)),1.)**2
             niacr(k) = min(niacr(k),ncr(k,i,3)*rdtcld)
           endif
#else
!-------------------------------------------------------------------------------
! praci: accretion of cloud ice by rain [hl a15]
!        (t<t0: i->r)
           acrfac = 2.*rslope3(k,1) + 2.*di(k)*rslope2(k,1)                    &
                   + di(k)**2*rslope(k,1)
           praci(k) = pi*qci(k,i,2)*n0r*abs(vtr(k)-vti(k))*acrfac/4.
           praci(k) = praci(k)*min(max(0.0,qrs(k,i,1)/qci(k,i,2)),1.)**2
           praci(k) = min(praci(k),qci(k,i,2)*rdtcld)
!-------------------------------------------------------------------------------
! piacr: accretion of rain by cloud ice [hl a19]
!        (t<t0: r->s or r->g)
           piacr(k) = pisq*avtr*n0r*denr*ni(k)*denfac(k)                       &
                       *g6pbr*rslope3(k,1)*rslope3(k,1)                        &
                       *rslopeb(k,1)/24./den(k,i)
           piacr(k) = piacr(k)*min(max(0.0,qci(k,i,2)/qrs(k,i,1)),1.)**2
           piacr(k) = min(piacr(k),qrs(k,i,1)*rdtcld)
#endif
         endif
       endif
     enddo
   endif
!
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
!
       if(ifice(k)) then
         if(qrs(k,i,2)>qrmin .and. qci(k,i,2)>qcmin) then
#ifdef WEAK_ACCRETION
           eacrs = exp(0.09*(tc(k)))
#else
           eacrs = exp(0.07*(tc(k)))
#endif
!-------------------------------------------------------------------------------
! psaci: accretion of cloud ice by rain [hdc 10]
!        (t<t0: i->s)
           acrfac = 2.*rslope3(k,2) + 2.*di(k)*rslope2(k,2)                    &
                   + di(k)**2*rslope(k,2)
           psaci(k) = pi*qci(k,i,2)*eacrs*n0s*n0sfac(k)                        &
                        *abs(vtmean(k)-vti(k))*acrfac/4.
           psaci(k) = psaci(k)*min(max(0.0,qrs(k,i,2)/qci(k,i,2)),1.)**2
           psaci(k) = min(psaci(k),qci(k,i,2)*rdtcld)
         endif
       endif
!-------------------------------------------------------------------------------
! psacw: accretion of cloud water by snow  [hl a7]
!        (t<t0: c->s, and t>=t0: c->r)
       if(qrs(k,i,2)>qrmin .and. qci(k,i,1)>qcmin) then
         psacw(k) = min(pacrs*n0sfac(k)*rslope3(k,2)*rslopeb(k,2)              &
                  *qci(k,i,1)*denfac(k),qci(k,i,1)*rdtcld)
             ! reduce collection efficiency (suggested by B. Wilt)
         psacw(k) = psacw(k)*min(max(0.0,qrs(k,i,2)/qci(k,i,1)),1.)**2
         psacw(k) = min(psacw(k),qci(k,i,1)*rdtcld)
       endif
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! nsacw: accretion of cloud water by snow  [lh A12]
!        (NC ->)
       if(qrs(k,i,2)>qrmin .and. qci(k,i,1)>qcmin .and. ncr(k,i,2)>ncmin) then
         nsacw(k) = min(pacrs*n0sfac(k)*rslope3(k,2)*rslopeb(k,2)              &
             ! reduce collection efficiency (suggested by B. Wilt)
                     *min(max(0.0,qrs(k,i,2)/qci(k,i,1)),1.)**2                &
                     *ncr(k,i,2)*denfac(k),ncr(k,i,2)*rdtcld)
       endif
#endif
     enddo
   endif
!
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
       if(ifice(k)) then
!-------------------------------------------------------------------------------
! pgaci: accretion of cloud ice by graupel [hl a17]
!        (t<t0: i->g)
         if(qrs(k,i,3)>qrmin .and. qci(k,i,2)>qcmin) then
#ifdef WEAK_ACCRETION
           egi = exp(0.09*(tc(k)))
#else
           egi = exp(0.07*(tc(k)))
#endif
           acrfac = 2.*rslope3(k,3) + 2.*di(k)*rslope2(k,3)                    &
                    +  di(k)**2*rslope(k,3)
#ifdef NEW_N0
           pgaci(k) = pi*egi*qci(k,i,2)*n0g*n0gfac(k)*abs(vtmean(k)-vti(k))    &
                    *acrfac/4.
#else
           pgaci(k) = pi*egi*qci(k,i,2)*n0g*abs(vtmean(k)-vti(k))*acrfac/4.
#endif
           pgaci(k) = min(pgaci(k),qci(k,i,2)*rdtcld)
         endif
       endif
!-------------------------------------------------------------------------------
! pgacw: accretion of cloud water by graupel [hl a6]
!        (t<t0: c->g, and t>=t0: c->r)
       if(qrs(k,i,3)>qrmin .and. qci(k,i,1)>qcmin) then
#ifdef NEW_N0
          pgacw(k) = min(pacrg*n0gfac(k)*rslope3(k,3)*rslopeb(k,3)             &
#else
          pgacw(k) = min(pacrg*rslope3(k,3)*rslopeb(k,3)                       &
#endif
             ! reduce collection efficiency (suggested by B. Wilt)
                    *min(max(0.0,qrs(k,i,3)/qci(k,i,1)),1.)**2                 &
                    *qci(k,i,1)*denfac(k),qci(k,i,1)*rdtcld)
       endif
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! ngacw: accretion of cloud water by graupel [lh a13]
!        (nc->
       if(qrs(k,i,3)>qrmin .and. qci(k,i,1)>qcmin .and. ncr(k,i,2)>ncmin) then
#ifdef NEW_N0
         ngacw(k) = min(pacrg*n0gfac(k)*rslope3(k,3)*rslopeb(k,3)*ncr(k,i,2)   &
#else
         ngacw(k) = min(pacrg*rslope3(k,3)*rslopeb(k,3)*ncr(k,i,2)             &
#endif
             ! reduce collection efficiency (suggested by B. Wilt)
                     *min(max(0.0,qrs(k,i,3)/qci(k,i,1)),1.)**2                &
                     *denfac(k),ncr(k,i,2)*rdtcld)
       endif
#endif
     enddo
   endif
!-------------------------------------------------------------------------------
! paacw: accretion of cloud water by averaged snow/graupel
!        (t<t0: c->g or s, and t>=t0: c->r)
   ktop = max(ktopqs,ktopqg)
   do k = ktop, kts, -1
     if(sumice(k)>qmin) then
       paacw(k) = (qrs(k,i,2)*psacw(k) + qrs(k,i,3)*pgacw(k))/sumice(k)
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! naacw: Accretion of cloud water by averaged snow/graupel
!        (nc->)
       naacw(k) = (qrs(k,i,2)*nsacw(k) + qrs(k,i,3)*ngacw(k))/sumice(k)
#endif
     endif
   enddo
!
#ifdef HAIL
   if(lqh) then
     ktop = ktopqh
     do k = ktop, kts, -1
       if(ifice(k)) then
!-------------------------------------------------------------------------------
! phaci: accretion of cloud ice by hail [bht]
!        (t<t0: i->h)
         if(qrs(k,i,4)>qrmin .and. qci(k,i,2)>qcmin) then
#ifdef WEAK_ACCRETION
           ehi = exp(0.09*(tc(k)))
#else
           ehi = exp(0.07*(tc(k)))
#endif
           acrfac = 2.*rslope3(k,4) + 2.*di(k)*rslope2(k,4)                    &
                    +  di(k)**2*rslope(k,4)
#ifdef NEW_N0
#ifdef SEDI_MEAN
           phaci(k) = pi*ehi*qci(k,i,2)*n0h*n0hfac(k)*abs(vtmean(k)-vti(k))    &
#else
           phaci(k) = pi*ehi*qci(k,i,2)*n0h*n0hfac(k)*abs(vth(k)-vti(k))    &
#endif
                    *acrfac/4.
#else
#ifdef SEDI_MEAN
           phaci(k) = pi*ehi*qci(k,i,2)*n0h*abs(vtmean(k)-vti(k))*acrfac/4.
#else
           phaci(k) = pi*ehi*qci(k,i,2)*n0h*abs(vth(k)-vti(k))*acrfac/4.
#endif
#endif
           phaci(k) = min(phaci(k),qci(k,i,2)*rdtcld)
         endif
       endif
!-------------------------------------------------------------------------------
! phacw: rccretion of cloud water by hail [bht a08]
!      (t<t0: qc->qh, AND t>=t0: qc->qr)
       if(qrs(k,i,4)>qrmin .and. qci(k,i,1)>qcmin) then
#ifdef NEW_N0
         phacw(k) = min(pacrh*n0hfac(k)*rslope3(k,4)*rslopeb(k,4)*qci(k,i,1)   &
#else
         phacw(k) = min(pacrh*rslope3(k,4)*rslopeb(k,4)*qci(k,i,1)             &
#endif
                        *min(max(0.0,qrs(k,i,4)/qci(k,i,1)),1.)**2             &
                        *denfac(k),qci(k,i,1)*rdtcld)
       endif
#ifdef DOUBLE
!
! nhacw: Accretion of cloud water by hail
!      (nc->)
       if(qrs(k,i,4)>qrmin .and. qci(k,i,1)>qcmin .and. ncr(k,i,2)>ncmin) then
#ifdef NEW_N0
         nhacw(k) = min(pacrh*n0hfac(k)*rslope3(k,4)*rslopeb(k,4)*ncr(k,i,2)   &
#else
         nhacw(k) = min(pacrh*rslope3(k,4)*rslopeb(k,4)*ncr(k,i,2)             &
#endif
                     *min(max(0.0,qrs(k,i,4)/qci(k,i,1)),1.)**2                &
                     *denfac(k),ncr(k,i,2)*rdtcld)
       endif
#endif
     enddo
   endif
#endif
!
   if(lqr) then
     ktop = ktopqr
     do k = kts, ktop
!-------------------------------------------------------------------------------
! pracs: accretion of snow by rain [hl a11]
!         (t<t0: s->g)
       if(qrs(k,i,2)>qrmin .and. qrs(k,i,1)>qrmin) then
         if(ifice(k)) then
#ifdef DOUBLE
           acrfac = 5.*rslope3(k,2)*rslope3(k,2)                               &
                  + 4.*rslope3(k,2)*rslope2(k,2)*rslope(k,1)                   &
                  + 1.5*rslope2(k,2)*rslope2(k,2)*rslope2(k,1)
           pracs(k) = pisq*ncr(k,i,3)*n0s*n0sfac(k)*abs(vtr(k)-vtmean(k))      &
                      *(dens/den(k,i))*acrfac
#else
           acrfac = 5.*rslope3(k,2)*rslope3(k,2)*rslope(k,1)                   &
                   +2.*rslope3(k,2)*rslope2(k,2)*rslope2(k,1)                  &
                   +.5*rslope2(k,2)*rslope2(k,2)*rslope3(k,1)
           pracs(k) = pisq*n0r*n0s*n0sfac(k)*abs(vtr(k)-vtmean(k))             &
                          *(dens/den(k,i))*acrfac
#endif
             ! reduce collection efficiency (suggested by B. Wilt)
           pracs(k) = pracs(k)*min(max(0.0,qrs(k,i,1)/qrs(k,i,2)),1.)**2
           pracs(k) = min(pracs(k),qrs(k,i,2)*rdtcld)
         endif
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! psacr: accretion of rain by snow [hl a10] [lfo 28]
!         (t<t0:qr->qs or qr->qg) (t>=t0: enhanced melting of snow)
         acrfac = 30.*rslope3(k,1)*rslope2(k,1)*rslope(k,2)                    &
                 +10.*rslope2(k,1)*rslope2(k,1)*rslope2(k,2)                   &
                 + 2.*rslope3(k,1)*rslope3(k,2 )
         psacr(k) = pisq*ncr(k,i,3)*n0s*n0sfac(k)*abs(vtmean(k)-vtr(k))        &
                        *(denr/den(k,i))*acrfac
#else
!-------------------------------------------------------------------------------
! psacr: accretion of rain by snow [hl a10]
!         (t<t0:r->s or r->g) (t>=t0: enhanced melting of snow)
         acrfac = 5.*rslope3(k,1)*rslope3(k,1)*rslope(k,2)                     &
                 +2.*rslope3(k,1)*rslope2(k,1)*rslope2(k,2)                    &
                 +.5*rslope2(k,1)*rslope2(k,1)*rslope3(k,2)
         psacr(k) = pisq*n0r*n0s*n0sfac(k)*abs(vtmean(k)-vtr(k))               &
                        *(denr/den(k,i))*acrfac
#endif
             ! reduce collection efficiency (suggested by B. Wilt)
         psacr(k) = psacr(k)*min(max(0.0,qrs(k,i,2)/qrs(k,i,1)),1.)**2
         psacr(k) = min(psacr(k),qrs(k,i,1)*rdtcld)
       endif
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! nsacr: accretion of rain by snow  [lh a23]
!       (t<t0:nr->)
       if(qrs(k,i,2)>qrmin .and. qrs(k,i,1)>qrmin .and. ncr(k,i,3)>nrmin) then
           acrfac = 1.5*rslope2(k,1)*rslope(k,2)                               &
                  + 1.0*rslope(k,1)*rslope2(k,2) + .5*rslope3(k,2)
           nsacr(k) = pi*ncr(k,i,3)*n0s*n0sfac(k)*abs(vtmean(k)-vtr(k))        &
                        *acrfac
             ! reduce collection efficiency (suggested by B. Wilt)
           nsacr(k) = nsacr(k)*min(max(0.0,qrs(k,i,2)/qrs(k,i,1)),1.)**2
           nsacr(k) = min(nsacr(k),ncr(k,i,3)*rdtcld)
       endif
#endif
#ifdef HAIL
!-------------------------------------------------------------------------------
! pracg: Accretion of graupel by rain [bht a17]
!         (t<t0: qg->qh)
       if(qrs(k,i,3)>qrmin .and. qrs(k,i,1)>qrmin) then
         if(ifice(k)) then
#ifdef DOUBLE
           acrfac = 5.*rslope3(k,3)*rslope3(k,3)                               &
                   +4.*rslope3(k,3)*rslope2(k,3)*rslope(k,1)                   &
                   +1.5*rslope2(k,3)*rslope2(k,3)*rslope2(k,1)
#ifdef NEW_N0
           pracg(k) = pisq*ncr(k,i,3)*n0g*n0gfac(k)*abs(vtr(k)-vtmean(k))      &
#else
           pracg(k) = pisq*ncr(k,i,3)*n0g*abs(vtr(k)-vtmean(k))                &
#endif
                       *(deng/den(k,i))*acrfac
#else
           acrfac = 5.*rslope3(k,3)*rslope3(k,3)*rslope(k,1)                   &
                   +2.*rslope3(k,3)*rslope2(k,3)*rslope2(k,1)                  &
                   +.5*rslope2(k,3)*rslope2(k,3)*rslope3(k,1)
#ifdef NEW_N0
           pracg(k) = pisq*n0r*n0g*n0gfac(k)*abs(vtr(k)-vtmean(k))             &
#else
           pracg(k) = pisq*n0r*n0g*abs(vtr(k)-vtmean(k))                       &
#endif
                       *(deng/den(k,i))*acrfac
#endif
           pracg(k) = pracg(k)*min(max(0.0,qrs(k,i,1)/qrs(k,i,3)),1.)**2
           pracg(k) = min(pracg(k),qrs(k,i,3)*rdtcld)
         endif
       endif
#endif
!-------------------------------------------------------------------------------
! pgacr: accretion of rain by graupel [hl a12]
!         (t<t0: r->g) (t>=t0: enhanced melting of graupel)
       if(qrs(k,i,3)>qrmin .and. qrs(k,i,1)>qrmin) then
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! pgacr: accretion of rain by graupel [HL A12] [LFO 42]
!         (t<t0: qr->qg) (t>=t0: enhanced melting of graupel)
         acrfac =  30.*rslope3(k,1)*rslope2(k,1)*rslope(k,3)                   &
                 + 10.*rslope2(k,1)*rslope2(k,1)*rslope2(k,3)                  &
                 +  2.*rslope3(k,1)*rslope3(k,3)
#ifdef NEW_N0
         pgacr(k) = pisq*ncr(k,i,3)*n0g*n0gfac(k)*abs(vtmean(k)-vtr(k))        &
                   *(denr/den(k,i))*acrfac
#else
         pgacr(k) = pisq*ncr(k,i,3)*n0g*abs(vtmean(k)-vtr(k))*(denr/den(k,i))  &
                       *acrfac
#endif
#else
         acrfac =  5.*rslope3(k,1)*rslope3(k,1)*rslope(k,3)                    &
                 + 2.*rslope3(k,1)*rslope2(k,1)*rslope2(k,3)                   &
                 + .5*rslope2(k,1)*rslope2(k,1)*rslope3(k,3)
         pgacr(k) = pisq*n0r*n0g*abs(vtmean(k)-vtr(k))*(denr/den(k,i))*acrfac
#endif
             ! reduce collection efficiency (suggested by B. Wilt)
         pgacr(k) = pgacr(k)*min(max(0.0,qrs(k,i,3)/qrs(k,i,1)),1.)**2
         pgacr(k) = min(pgacr(k),qrs(k,i,1)*rdtcld)
       endif
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! ngacr: accretion of rain by graupel  [lh a24]
!        (t<t0: nr->)
       if(qrs(k,i,3)>qrmin .and. qrs(k,i,1)>qrmin .and. ncr(k,i,3)>nrmin) then
         acrfac = 1.5*rslope2(k,1)*rslope(k,3)                                 &
                 +1.0*rslope(k,1)*rslope2(k,3) + .5*rslope3(k,3)
#ifdef NEW_N0
         ngacr(k) = pi*ncr(k,i,3)*n0g*n0gfac(k)*abs(vtmean(k)-vtr(k))*acrfac
#else
         ngacr(k) = pi*ncr(k,i,3)*n0g*abs(vtmean(k)-vtr(k))*acrfac
#endif
         ngacr(k) = ngacr(k)*min(max(0.0,qrs(k,i,3)/qrs(k,i,1)),1.)**2
         ngacr(k) = min(ngacr(k),ncr(k,i,3)*rdtcld)
       endif
#endif
     enddo
   endif
!
#ifdef HAIL
   if(lqh) then
     ktop = ktopqh
     do k = kts, ktop
!-------------------------------------------------------------------------------
! phacr: Accretion of rain by hail [bht a13]
!         (t<t0: qr->qh) (t>=t0: enhance melting of hail)
       if(qrs(k,i,4)>qrmin.and.qrs(k,i,1)>qrmin) then
#ifdef DOUBLE
         acrfac = 30.*rslope3(k,1)*rslope2(k,1)*rslope(k,4)                    &
                 +10.*rslope3(k,1)*rslope(k,1)*rslope2(k,4)                    &
                 + 2.*rslope3(k,1)*rslope3(k,4)
#ifdef NEW_N0
#ifdef SEDI_MEAN
         phacr(k) = pisq*ncr(k,i,3)*n0h*n0hfac(k)*abs(vtmean(k)-vtr(k))        &
#else
         phacr(k) = pisq*ncr(k,i,3)*n0h*n0hfac(k)*abs(vth(k)-vtr(k))        &
#endif
                   *(denr/den(k,i))*acrfac
#else
#ifdef SEDI_MEAN
         phacr(k) = pisq*ncr(k,i,3)*n0h*abs(vtmean(k)-vtr(k))*(denr/den(k,i))  &
#else
         phacr(k) = pisq*ncr(k,i,3)*n0h*abs(vth(k)-vtr(k))*(denr/den(k,i))  &
#endif
                     *acrfac
#endif
#else
         acrfac = 5.*rslope3(k,1)*rslope3(k,1)*rslope(k,4)                     &
                 +2.*rslope3(k,1)*rslope2(k,1)*rslope2(k,4)                    &
                 +.5*rslope2(k,1)*rslope2(k,1)*rslope3(k,4)
         phacr(k) = pisq*n0r*n0h*abs(vtmean(k)-vtr(k))*(denr/den(k,i))         &
                        *acrfac
#endif
         phacr(k) = phacr(k)*min(max(0.0,qrs(k,i,4)/qrs(k,i,1)),1.)**2
         phacr(k) = min(phacr(k),qrs(k,i,1)*rdtcld)
       endif
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! nhacr: accretion of rain by hail
!        (t<t0: nr->)
       if(qrs(k,i,4)>qrmin .and. qrs(k,i,1)>qrmin .and. ncr(k,i,3)>nrmin) then
         acrfac = 1.5*rslope2(k,1)*rslope(k,4)                                 &
                 + 1.0*rslope(k,1)*rslope2(k,4) + .5*rslope3(k,4)
#ifdef NEW_N0
#ifdef SEDI_MEAN
         nhacr(k) = pi*ncr(k,i,3)*n0h*n0hfac(k)*abs(vtmean(k)-vtr(k))*acrfac
#else
         nhacr(k) = pi*ncr(k,i,3)*n0h*n0hfac(k)*abs(vth(k)-vtr(k))*acrfac
#endif
#else
#ifdef SEDI_MEAN
         nhacr(k) = pi*ncr(k,i,3)*n0h*abs(vtmean(k)-vtr(k))*acrfac
#else
         nhacr(k) = pi*ncr(k,i,3)*n0h*abs(vth(k)-vtr(k))*acrfac
#endif
#endif
         nhacr(k) = nhacr(k)*min(max(0.0,qrs(k,i,4)/qrs(k,i,1)),1.)**2
         nhacr(k) = min(nhacr(k),ncr(k,i,3)/dtcld)
       endif
#endif
!-------------------------------------------------------------------------------
! phacs: accretion of snow by hail [bht a14]
!         (t<t0: qs->qh)
       if(qrs(k,i,4)>qrmin.and.qrs(k,i,2)>qrmin) then
         acrfac = 5.*rslope3(k,2)*rslope3(k,2)*rslope(k,4)                     &
                 +2.*rslope3(k,2)*rslope2(k,2)*rslope2(k,4)                    &
                 +.5*rslope2(k,2)*rslope2(k,2)*rslope3(k,4)
#ifdef NEW_N0
         phacs(k) = pisq*eachs*n0s*n0sfac(k)*n0h*n0hfac(k)                     &
#ifdef SEDI_MEAN
                    *abs(vtmean(k)-vtmean(k))                                  &
#else
                    *abs(vth(k)-vtmean(k))                                  &
#endif
#else
#ifdef SEDI_MEAN
         phacs(k) = pisq*eachs*n0s*n0sfac(k)*n0h*abs(vtmean(k)-vtmean(k))      &
#else
         phacs(k) = pisq*eachs*n0s*n0sfac(k)*n0h*abs(vth(k)-vtmean(k))      &
#endif
#endif
                    *(dens/den(k,i))*acrfac
         phacs(k) = min(phacs(k),qrs(k,i,2)*rdtcld)
       endif
!-------------------------------------------------------------------------------
! phacg: accretion of snow by hail [bht a15]
!         (t<t0: qg->qh)
       if(qrs(k,i,4)>qrmin.and.qrs(k,i,3)>qrmin) then
         acrfac = 5.*rslope3(k,3)*rslope3(k,3)*rslope(k,4)                     &
                 +2.*rslope3(k,3)*rslope2(k,3)*rslope2(k,4)                    &
                 +.5*rslope2(k,3)*rslope2(k,3)*rslope3(k,4)
#ifdef NEW_N0
         phacg(k) = pisq*eachg*n0g*n0gfac(k)*n0h*n0hfac(k)                     &
#ifdef SEDI_MEAN
                    *abs(vtmean(k)-vtmean(k))                                  &
#else
                    *abs(vth(k)-vtmean(k))                                  &
#endif
#else
         phacg(k) = pisq*eachg*n0g*n0h*abs(vtmean(k)-vtmean(k))                &
#endif
                     *(deng/den(k,i))*acrfac
         phacg(k) = min(phacg(k),qrs(k,i,3)*rdtcld)
       endif
!-------------------------------------------------------------------------------
! pgwet: wet growth of graupel [lfo 43]
!
       if(qrs(k,i,4)>qrmin.or.qrs(k,i,3)>qrmin) then
         rs0 = fsvp(t0c,p(k,i))
         rs0 = ep2*rs0/(p(k,i)-rs0)
         rs0 = max(rs0,qmin)
         ghw1 = den(k,i)*xlv(k)*diffus(k)*(rs0-q(k,i)) - viscod(k)*(tc(k))
         ghw2 = den(k,i)*(xlf0+cliq*(tc(k)))
         ghw3 = venfac(k)
         ghw4 = den(k,i)*(xlf0-cliq*(-tc(k))+cice*(-tc(k)))
       endif
       if(qrs(k,i,3)>qrmin) then
         if(pgaci(k)>0.0) then
#ifdef WEAK_ACCRETION
           egi = exp(0.09*(tc(k)))
#else
           egi = exp(0.07*(tc(k)))
#endif
           pgaci_w(k) = pgaci(k)/egi
         else
           pgaci_w(k) = 0.0
         endif
         pgwet(k) = ghw1/ghw2*(precg1*rslope2(k,3)                             &
                   +precg3*ghw3*exp(log(rslope(k,4))*2.75)                     &
                   +ghw4*(pgaci_w(k)+pgacs(k)))
         pgwet(k) = max(pgwet(k), 0.0)
         if(pgacw(k)+pgacr(k)<0.95*pgwet(k)) then
           pgaci(k) = 0.0
           pgacs(k) = 0.0
         endif
       endif
!-------------------------------------------------------------------------------
! phwet: wet growth of hail [lfo 43]
!
       if(qrs(k,i,4)>qrmin) then
         if(phaci(k)>0.0) then
#ifdef WEAK_ACCRETION
           ehi = exp(0.09*(tc(k)))
#else
           ehi = exp(0.07*(tc(k)))
#endif
           phaci_w(k) = phaci(k)/ehi
         else
           phaci_w(k) = 0.0
         endif
         phwet(k) = ghw1/ghw2*(prech1*rslope2(k,4)                             &
                   +prech3*ghw3*exp(log(rslope(k,4))*2.75)                     &
                   +ghw4*(phaci_w(k)+phacs(k)))
         phwet(k) = max(phwet(k), 0.0)
         if(phacw(k)+phacr(k)<0.95*phwet(k)) then
           phaci(k) = 0.0
           phacs(k) = 0.0
           phacg(k) = 0.0
         endif
       endif
     enddo
   endif
#endif
!-------------------------------------------------------------------------------
! pseml: enhanced melting of snow by accretion of water [hl a34]
!        (t>=t0: s->r)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
       if(.not.ifice(k) .and. qrs(k,i,2)>0.) then
         pseml(k) = min(max(-cliq*tc(k)*(paacw(k)+psacr(k))/xlf0,              &
                    -qrs(k,i,2)*rdtcld),0.)
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! nseml: enhanced melting of snow by accretion of water    [lh a29]
!        (t>=t0: ->nr)
         if  (qrs(k,i,2)>qrmin) then
           sfac = rslope(k,2)*n0s*n0sfac(k)/qrs(k,i,2)
           nseml(k) = -sfac*pseml(k)
         endif
#endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgeml: enhanced melting of graupel by accretion of water [hl a24]
!        (t>=t0: g->r)
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
       if(.not.ifice(k) .and. qrs(k,i,3)>0.) then
         pgeml(k) = min(max(-cliq*tc(k)*(paacw(k)+pgacr(k))/xlf0,              &
                   -qrs(k,i,3)*rdtcld),0.)
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! ngeml: enhanced melting of graupel by accretion of water [lh a30]
!         (t>=t0: -> nr)
         if  (qrs(k,i,3)>qrmin) then
#ifdef NEW_N0
           gfac = rslope(k,3)*n0g*n0gfac(k)/qrs(k,i,3)
#else
           gfac = rslope(k,3)*n0g/qrs(k,i,3)
#endif
           ngeml(k) = -gfac*pgeml(k)
         endif
#endif
       endif
     enddo
   endif
#ifdef HAIL
!-------------------------------------------------------------------------------
! pheml: enhanced melting of hail by accretion of water [bht a23]
!        (t>=t0: h->r)
   if(lqh) then
     ktop = ktopqh
     do k = ktop, kts, -1
       if(.not.ifice(k) .and. qrs(k,i,4)>0.) then
         pheml(k) = min(max(-cliq*tc(k)*(phacw(k)+phacr(k))/xlf0,              &
                   -qrs(k,i,4)*rdtcld),0.)
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! nheml: enhanced melting of graupel by accretion of water [lh a30]
!         (t>=t0: -> nr)
         if  (qrs(k,i,4)>qrmin) then
#ifdef NEW_N0
           gfac = rslope(k,4)*n0h/qrs(k,i,4)
#else
           gfac = rslope(k,4)*n0h*n0hfac(k)/qrs(k,i,4)
#endif
           nheml(k) = -gfac*pheml(k)
         endif
#endif
       endif
     enddo
   endif
#endif
! -------------------------------------------------------------------------------
! pidep: deposition/sublimation rate of ice [hdc 9]
!       (t<t0: v->i or i->v)
   if(lqi) then
     ktop = ktopqi
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qci(k,i,2)>0 .and. (.not.ifsat(k))) then
           pidep(k) = 4.*di(k)*ni(k)*(rh_ice(k)-1.)/ab_ice(k)
           supice(k) = satrdt(k) - prevp(k)
           if(pidep(k)<0.) then
             pidep(k) = max(max(pidep(k),satrdt(k)*.5),supice(k))
             pidep(k) = max(pidep(k),-qci(k,i,2)*rdtcld)
           else
             pidep(k) = min(min(pidep(k),satrdt(k)*.5),supice(k))
           endif
           if(abs(prevp(k)+pidep(k))>=abs(satrdt(k))) ifsat(k) = .true.
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! psdep: deposition/sublimation rate of snow [hdc 14]
!        (v->s or s->v)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,2)>0. .and. (.not.ifsat(k))) then
           coeres = rslope2(k,2)*sqrt(rslope(k,2)*rslopeb(k,2))
           psdep(k) = (rh_ice(k)-1.)*n0sfac(k)*(precs1*rslope2(k,2)            &
                      + precs2*venfac(k)*coeres)/ab_ice(k)
           supice(k) = satrdt(k)-prevp(k)-pidep(k)
           if(psdep(k)<0.) then
             psdep(k) = max(psdep(k),-qrs(k,i,2)*rdtcld)
             psdep(k) = max(max(psdep(k),satrdt(k)*.5),supice(k))
           else
             psdep(k) = min(min(psdep(k),satrdt(k)*.5),supice(k))
           endif
           if(abs(prevp(k)+pidep(k)+psdep(k))>=abs(satrdt(k)))                 &
             ifsat(k) = .true.
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgdep: deposition/sublimation rate of graupel [hl a21]
!        (t<t0: v->g or g->v)
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,3)>0. .and. (.not.ifsat(k))) then
           coeres = rslope2(k,3)*sqrt(rslope(k,3)*rslopeb(k,3))
#ifdef NEW_N0
           pgdep(k) = (rh_ice(k)-1.)*n0gfac(k)*(precg1*rslope2(k,3)            &
#else
           pgdep(k) = (rh_ice(k)-1.)*(precg1*rslope2(k,3)                      &
#endif
                     + precg2*venfac(k)*coeres)/ab_ice(k)
           supice(k) = satrdt(k)-prevp(k)-pidep(k)-psdep(k)
           if(pgdep(k)<0.) then
             pgdep(k) = max(pgdep(k),-qrs(k,i,3)*rdtcld)
             pgdep(k) = max(max(pgdep(k),satrdt(k)*.5),supice(k))
           else
             pgdep(k) = min(min(pgdep(k),satrdt(k)*.5),supice(k))
           endif
           if(abs(prevp(k)+pidep(k)+psdep(k)+pgdep(k))>=abs(satrdt(k)))        &
             ifsat(k) = .true.
         endif
       endif
     enddo
   endif
#ifdef HAIL
!-------------------------------------------------------------------------------
! phdep: deposition/sublimation rate of hail [bht a19]
!        (t<t0: v->h or h->v)
   if(lqh) then
     ktop = ktopqh
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,4)>0. .and. (.not.ifsat(k))) then
           coeres = rslope2(k,4)*sqrt(rslope(k,4)*rslopeb(k,4))
#ifdef NEW_N0
           phdep(k) = (rh_ice(k)-1.)*n0hfac(k)*(prech1*rslope2(k,4)            &
#else
           phdep(k) = (rh_ice(k)-1.)*(prech1*rslope2(k,4)                      &
#endif
                     + prech2*venfac(k)*coeres)/ab_ice(k)
           supice(k) = satrdt(k)-prevp(k)-pidep(k)-psdep(k)-pgdep(k)
           if(phdep(k)<0.) then
             phdep(k) = max(phdep(k),-qrs(k,i,4)*rdtcld)
             phdep(k) = max(max(phdep(k),satrdt(k)*.5),supice(k))
           else
             phdep(k) = min(min(phdep(k),satrdt(k)*.5),supice(k))
           endif
           if(abs(prevp(k)+pidep(k)+psdep(k)+pgdep(k)+phdep(k))                &
               >=abs(satrdt(k))) ifsat(k) = .true.
         endif
       endif
     enddo
   endif
#endif
!-------------------------------------------------------------------------------
! pigen: generation(nucleation) of ice from vapor [hl a50] [hdc 7-8]
!       (t<t0: v->i)
   ktop = ktoprh 
   do k = ktop, kts, -1
     if(ifice(k)) then
       if(supsat(k)>0 .and. (.not.ifsat(k))) then
         supice(k) = satrdt(k)-prevp(k)-pidep(k)-psdep(k)-pgdep(k)-phdep(k)
         if (if_aerosol_activated) then 
           ni0(k) = iceDeMott(tc(k), den(k,i), nifa(k,i))
           ni0(k) = min(ni0(k),inpmax)
           roqi0 = mi0 * ni0(k)
#ifdef DEBUG
           if(i.eq.its) print*,' it k tc nifa ni0 ni0-old',itimestep,k,tc(k), &
                     nifa(k,i)/1000.,ni0(k)/1000.,exp(-0.1*tc(k))
#endif
         else
#ifndef NI0_HDC2004
           if(slimsk(i)==0) then
             ni0(k) = 1000.*exp(-0.2*tc(k)-5.)
           else if(slimsk(i)==1) then
             ni0(k) = 1000.*exp(-0.15*tc(k)-2.5)
           else if(slimsk(i)==2) then
             ni0(k) = 1000.*exp(-0.35*tc(k)-10.)
           endif
#else
             ni0(k) = 1000.*exp(-0.1*tc(k))
#endif
           ni0(k) = min(ni0(k),inpmax)
           roqi0 = 4.92e-11*exp(log(ni0(k))*(1.33))
         endif
         pigen(k) = max(0.,(roqi0/den(k,i)-max(qci(k,i,2),0.))*rdtcld)
         pigen(k) = min(min(pigen(k),satrdt(k)),supice(k))
       endif
     endif
   enddo
!
!-------------------------------------------------------------------------------
! psaut: conversion(aggregation) of ice to snow [hdc 12]
!       (t<t0: i->s)
   if(lqi) then
     ktop = ktopqi
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qci(k,i,2)>0.) then
#ifdef RYAN2000_qimax
           if (t(k,i)>255.66) THEN
             qimax = 0.4632*10**(-6.0-0.0413*(tc(k)))*den(i,k)
           else
             qimax = 0.2316*10**(-4.0+0.0519*(tc(k)))*den(i,k)
           endif
           hvalue = 0.005*exp(0.025*(tc(k))) ![s-1]
           psaut(k) = min(max(0.,hvalue*(qci(k,i,2)-qimax)),qci(k,i,2)*rdtcld)
#else
           qimax = roqimax/den(k,i)
           psaut(k) = max(0.,(qci(k,i,2)-qimax)*rdtcld)
#endif
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgaut: conversion(aggregation) of snow to graupel [hl a4]
!        (t<t0: s->g)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,2)>0.) then
           alpha2 = 1.e-3*exp(0.09*(tc(k)))
           pgaut(k) = min(max(0.,alpha2*(qrs(k,i,2)-qs0)),qrs(k,i,2)*rdtcld)
         endif
       endif
     enddo
   endif
#ifdef HAIL
!-------------------------------------------------------------------------------
! phaut: conversion(aggregation) of graupel to hail [bht a18]
!        (t<t0: g->h)
   if(lqh) then
     ktop = ktopqh
     do k = ktop, kts, -1
       if(ifice(k)) then
         if(qrs(k,i,3)>0.) then
           alpha2 = 1.e-3*exp(0.09*(tc(k)))
           phaut(k) = min(max(0.,alpha2*(qrs(k,i,3)-qs0)),qrs(k,i,3)*rdtcld)
         endif
       endif
     enddo
   endif
#endif
!-------------------------------------------------------------------------------
! psevp: evaporation of melting snow [hl a35]
!       (t>t0: s->v)
   if(lqs) then
     ktop = ktopqs
     do k = ktop, kts, -1
       if(.not.ifice(k)) then
         if(qrs(k,i,2)>0. .and. rh_mul(k)<1.) then
           coeres = rslope2(k,2)*sqrt(rslope(k,2)*rslopeb(k,2))
           psevp(k) = (rh_mul(k)-1.)*n0sfac(k)*(precs1*rslope2(k,2)            &
                       + precs2*venfac(k)*coeres)/ab_mul(k)
           psevp(k) = min(max(psevp(k),-qrs(k,i,2)*rdtcld),0.)
         endif
       endif
     enddo
   endif
!-------------------------------------------------------------------------------
! pgevp: evaporation of melting graupel [hl a25]
!       (t>=t0: g->v)
   if(lqg) then
     ktop = ktopqg
     do k = ktop, kts, -1
       if(.not.ifice(k)) then
         if(qrs(k,i,3)>0. .and. rh_mul(k)<1.) then
           coeres = rslope2(k,3)*sqrt(rslope(k,3)*rslopeb(k,3))
#ifdef NEW_N0
           pgevp(k) = (rh_mul(k)-1.)*n0gfac(k)*(precg1*rslope2(k,3)            &
#else
           pgevp(k) = (rh_mul(k)-1.)*(precg1*rslope2(k,3)                      &
#endif
                     + precg2*venfac(k)*coeres)/ab_mul(k)
           pgevp(k) = min(max(pgevp(k),-qrs(k,i,3)*rdtcld),0.)
         endif
       endif
     enddo
   endif
#ifdef HAIL
!-------------------------------------------------------------------------------
! phevp: evaporation of melting hail [bht a20]
!       (t>=t0: h->v)
   if(lqh) then
     ktop = ktopqh
     do k = ktop, kts, -1
       if(.not.ifice(k)) then
         if(qrs(k,i,4)>0. .and. rh_mul(k)<1.) then
           coeres = rslope2(k,4)*sqrt(rslope(k,4)*rslopeb(k,4))
#ifdef NEW_N0
           phevp(k) = (rh_mul(k)-1.)*n0hfac(k)*(prech1*rslope2(k,4)            &
#else
           phevp(k) = (rh_mul(k)-1.)*(prech1*rslope2(k,4)                      &
#endif
                     + prech2*venfac(k)*coeres)/ab_mul(k)
           phevp(k) = min(max(phevp(k),-qrs(k,i,4)*rdtcld),0.)
         endif
       endif
     enddo
   endif
#endif
#ifdef PARTIAL
!===============================================================================
!
! convert the in-cloud source/sink terms to the grd-mean state
!
!===============================================================================
   ktop = max(ktopqc,ktopqi)
   do k = ktop, kts, -1
       if(cldf(k)>0.) then
           ! change in-cloud condensate variables to grid-mean variables
         qci(k,i,1) = qci(k,i,1) * cldf(k)
         qci(k,i,2) = qci(k,i,2) * cldf(k)
         qrs(k,i,1) = qrs(k,i,1) * cldf(k)
         qrs(k,i,2) = qrs(k,i,2) * cldf(k)
         qrs(k,i,3) = qrs(k,i,3) * cldf(k)
#ifdef HAIL
         qrs(k,i,4) = qrs(k,i,4) * cldf(k)
#endif
           ! get grid-mean rates
         praut(k) = praut(k) * cldf(k)
         pracw(k) = pracw(k) * cldf(k)
         prevp(k) = prevp(k) * cldf(k)
         praci(k) = praci(k) * cldf(k)
         piacr(k) = piacr(k) * cldf(k)
         psaci(k) = psaci(k) * cldf(k)
         pgaci(k) = pgaci(k) * cldf(k)
         psacw(k) = psacw(k) * cldf(k)
         pgacw(k) = pgacw(k) * cldf(k)
         paacw(k) = paacw(k) * cldf(k)
         pracs(k) = pracs(k) * cldf(k)
         pgacr(k) = pgacr(k) * cldf(k)
         pgacs(k) = pgacs(k) * cldf(k)
         pseml(k) = pseml(k) * cldf(k)
         pgeml(k) = pgeml(k) * cldf(k)
         pidep(k) = pidep(k) * cldf(k)
         psdep(k) = psdep(k) * cldf(k)
         pgdep(k) = pgdep(k) * cldf(k)
         pigen(k) = pigen(k) * cldf(k)
         psaut(k) = psaut(k) * cldf(k)
         pgaut(k) = pgaut(k) * cldf(k)
         psevp(k) = psevp(k) * cldf(k)
         pgevp(k) = pgevp(k) * cldf(k)
#ifdef HAIL
         phdep(k) = phdep(k) * cldf(k)
         phaut(k) = phaut(k) * cldf(k)
         pracg(k) = pracg(k) * cldf(k)
         phacw(k) = phacw(k) * cldf(k)
         phaci(k) = phaci(k) * cldf(k)
         phacr(k) = phacr(k) * cldf(k)
         phacs(k) = phacs(k) * cldf(k)
         pheml(k) = pheml(k) * cldf(k)
         phevp(k) = phevp(k) * cldf(k)
         primh(k) = primh(k) * cldf(k)
         pvapg(k) = pvapg(k) * cldf(k)
         pvaph(k) = pvaph(k) * cldf(k)
         pgwet(k) = pgwet(k) * cldf(k)
         phwet(k) = phwet(k) * cldf(k)
         pgaci_w(k) = pgaci_w(k) * cldf(k)
         phaci_w(k) = phaci_w(k) * cldf(k)
#endif
       endif
   enddo
#endif
!===============================================================================
!
! check mass conservation of source/sink terms and feedback to the large scale
!
!===============================================================================
   ktop = ktopmax
   do k = ktop, kts, -1
!
#ifdef PRT
       hvalue = 1.e10
       hvalue = 10.
       if(i==lond.and. lat==latd) print 111, k,praci(k)*hvalue,piacr(k)*hvalue, psaci(k)*hvalue, &
          pgaci(k)*hvalue, psacw(k)*hvalue, pgacw(k)*hvalue, paacw(k)*hvalue, pracs(k)*hvalue, pgacr(k)*hvalue, pgacs(k)*hvalue
        if(i==lond.and. lat==latd) print 111, k,pseml(k)*hvalue, pgeml(k)*hvalue, pidep(k)*hvalue, psdep(k)*hvalue, &
            pgdep(k)*hvalue, pigen(k)*hvalue, psaut(k)*hvalue,pgaut(k)*hvalue, psevp(k)*hvalue, pgevp(k)*hvalue
!        if(i==lond.and. lat==latd) print 111, k,praut(k)*hvalue, pracw(k)*hvalue, prevp(k)*hvalue, paacw(k)*hvalue, &
!       pseml(k)*hvalue, pgeml(k)*hvalue, nraut(k)*hvalue, nrcol(k)*hvalue, nseml(k)*hvalue, ngeml(k)*hvalue
111   format(1x,i5,10f8.2)
#endif
#ifdef PRT
   if(i==lond.and.lat==latd) then
     print*,' k tc nos n0g n0h  ',k,tc(k),n0sfac(k),n0gfac(k),n0hfac(k)
   endif
#endif
     delta2 = 0.
     delta3 = 0.
     if(qrs(k,i,1)<1.e-4 .and. qrs(k,i,2)<1.e-4) delta2 = 1.
     if(qrs(k,i,1)<1.e-4) delta3 = 1.
!
     if(ifice(k)) then
!
! cloud water
!
       hvalue = max(qmin,qci(k,i,1))
#ifdef HAIL
       source = (praut(k)+pracw(k)+paacw(k)+paacw(k)+phacw(k))*dtcld
#else
       source = (praut(k)+pracw(k)+paacw(k)+paacw(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         praut(k) = praut(k) * factor
         pracw(k) = pracw(k) * factor
         paacw(k) = paacw(k) * factor
#ifdef HAIL
         phacw(k) = phacw(k) * factor
#endif
       endif
!
! cloud ice
!
       hvalue = max(qmin,qci(k,i,2))
#ifdef HAIL
       source = (psaut(k)-pigen(k)-pidep(k)+praci(k)+psaci(k)+pgaci(k)         &
                + phaci(k))*dtcld
#else
       source = (psaut(k)-pigen(k)-pidep(k)+praci(k)+psaci(k)+pgaci(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         psaut(k) = psaut(k) * factor
         pigen(k) = pigen(k) * factor
         pidep(k) = pidep(k) * factor
         praci(k) = praci(k) * factor
         psaci(k) = psaci(k) * factor
         pgaci(k) = pgaci(k) * factor
#ifdef HAIL
         phaci(k) = phaci(k) * factor
#endif
       endif
!
! rain
!
       hvalue = max(qmin,qrs(k,i,1))
#ifdef HAIL
       source = (-praut(k)-prevp(k)-pracw(k)+piacr(k)+psacr(k)+pgacr(k)        &
                + phacr(k))*dtcld
#else
       source = (-praut(k)-prevp(k)-pracw(k)+piacr(k)+psacr(k)+pgacr(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         praut(k) = praut(k) * factor
         prevp(k) = prevp(k) * factor
         pracw(k) = pracw(k) * factor
         piacr(k) = piacr(k) * factor
         psacr(k) = psacr(k) * factor
         pgacr(k) = pgacr(k) * factor
#ifdef HAIL
         phacr(k) = phacr(k) * factor
#endif
       endif
!
! snow
!
       hvalue = max(qmin,qrs(k,i,2))
       source = - (psdep(k)+psaut(k)-pgaut(k)+paacw(k)+piacr(k)*delta3         &
#ifdef HAIL
                + pvapg(k)+pvaph(k)                                            &
#endif
                + praci(k)*delta3 - pracs(k)*(1.-delta2)                       &
#ifdef HAIL
                + psacr(k)*delta2 + psaci(k)-pgacs(k)-phacs(k))*dtcld
#else
                + psacr(k)*delta2 + psaci(k)-pgacs(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         psdep(k) = psdep(k) * factor
         psaut(k) = psaut(k) * factor
         pgaut(k) = pgaut(k) * factor
         paacw(k) = paacw(k) * factor
         piacr(k) = piacr(k) * factor
         praci(k) = praci(k) * factor
         psaci(k) = psaci(k) * factor
         pracs(k) = praCs(k) * factor
         psacr(k) = psacr(k) * factor
         pgacs(k) = pgacs(k) * factor
#ifdef HAIL
         pvapg(k) = pvapg(k) * factor
         pvaph(k) = pvaph(k) * factor
         phacs(k) = phacs(k) * factor
#endif
       endif
!
! graupel
!
       hvalue = max(qmin,qrs(k,i,3))
       source = - (pgdep(k)+pgaut(k)                                           &
                + piacr(k)*(1.-delta3) + praci(k)*(1.-delta3)                  &
                + psacr(k)*(1.-delta2) + pracs(k)*(1.-delta2)                  &
#ifdef HAIL
                + pgaci(k)+paacw(k)+pgacr(k)*delta2+pgacs(k)                   &
                - pracg(k)*(1.-delta2)-phacg(k)-phaut(k)                       &
                - pvapg(k)+primh(k))*dtcld
#else
                + pgaci(k)+paacw(k)+pgacr(k)+pgacs(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         pgdep(k) = pgdep(k) * factor
         pgaut(k) = pgaut(k) * factor
         piacr(k) = piacr(k) * factor
         praci(k) = praci(k) * factor
         psacr(k) = psacr(k) * factor
         pracs(k) = pracs(k) * factor
         paacw(k) = paacw(k) * factor
         pgaci(k) = pgaci(k) * factor
         pgacr(k) = pgacr(k) * factor
         pgacs(k) = pgacs(k) * factor
#ifdef HAIL
         phaut(k) = phaut(k) * factor
         pracg(k) = pracg(k) * factor
         phacg(k) = phacg(k) * factor
         pvapg(k) = pvapg(k) * factor
         primh(k) = primh(k) * factor
#endif
       endif
#ifdef HAIL
!
!     hail
!
       hvalue = max(qmin,qrs(k,i,4))
       source = -(phdep(k)+phaut(k)                                            &
                +pgacr(k)*(1.-delta2)+pracg(k)*(1.-delta2)                     &
                +phacw(k)+phacr(k)+phaci(k)+phacs(k)                           &
                +phacg(k)-pvaph(k)-primh(k))*dtcld
       if (source>hvalue) then
         factor = hvalue/source
         phdep(k) = phdep(k) * factor
         phaut(k) = phaut(k) * factor
         pracg(k) = pracg(k) * factor
         pgacr(k) = pgacr(k) * factor
         phacw(k) = phacw(k) * factor
         phaci(k) = phaci(k) * factor
         phacr(k) = phacr(k) * factor
         phacs(k) = phacs(k) * factor
         phacg(k) = phacg(k) * factor
         pvaph(k) = pvaph(k) * factor
         primh(k) = primh(k) * factor
       endif
#endif
#ifdef DOUBLE
!
!     cloud
!
       hvalue = max(qmin,ncr(k,i,2))
#ifdef HAIL
#ifdef PRAUT_LB2017
       source = (ncaut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k)                  &
#else
       source = (nraut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k)                  &
#endif
               +nhacw(k))*dtcld
#else
#ifdef PRAUT_LB2017
       source = (ncaut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k))*dtcld
#else
       source = (nraut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k))*dtcld
#endif
#endif
       if (source>hvalue) then
         factor = hvalue/source
#ifdef PRAUT_LB2017
         ncaut(k) = ncaut(k) * factor
#else
         nraut(k) = nraut(k) * factor
#endif
         nccol(k) = nccol(k) * factor
         nracw(k) = nracw(k) * factor
         naacw(k) = naacw(k) * factor
#ifdef HAIL
         nhacw(k) = nhacw(k) * factor
#endif
       endif
!
!     rain
!
       hvalue = max(qmin,ncr(k,i,3))
#ifdef HAIL
       source = (-nraut(k)+nrcol(k)+niacr(k)+nsacr(k)+ngacr(k)                 &
                 +nhacr(k))*dtcld
#else
       source = (-nraut(k)+nrcol(k)+niacr(k)+nsacr(k)+ngacr(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         nraut(k) = nraut(k) * factor
         nrcol(k) = nrcol(k) * factor
         niacr(k) = niacr(k) * factor
         nsacr(k) = nsacr(k) * factor
         ngacr(k) = ngacr(k) * factor
#ifdef HAIL
         nhacr(k) = nhacr(k) * factor
#endif
       endif
#endif
!
! update
!
#ifdef HAIL
       htotal = -(prevp(k)+psdep(k)+pgdep(k)+pigen(k)+pidep(k)+phdep(k))
#else
       htotal = -(prevp(k)+psdep(k)+pgdep(k)+pigen(k)+pidep(k))
#endif
       q(k,i) = q(k,i) + htotal*dtcld
#ifdef HAIL
       htotal = (praut(k)+pracw(k)+paacw(k)+paacw(k)+phacw(k))*dtcld
#else
       htotal = (praut(k)+pracw(k)+paacw(k)+paacw(k))*dtcld
#endif
       qci(k,i,1) = max(qci(k,i,1)-htotal, 0.) 
       htotal =    (praut(k)+pracw(k)                                          &
#ifdef HAIL
                   +prevp(k)-piacr(k)-pgacr(k)-psacr(k)-phacr(k))*dtcld
#else
                   +prevp(k)-piacr(k)-pgacr(k)-psacr(k))*dtcld
#endif
       qrs(k,i,1) = max(qrs(k,i,1)+htotal, 0.)
       htotal =    (psaut(k)+praci(k)                                          &
#ifdef HAIL
                   +phaci(k)                                                   &
#endif
                   +psaci(k)+pgaci(k)-pigen(k)-pidep(k))*dtcld
       qci(k,i,2) = max(qci(k,i,2)-htotal, 0.)

       htotal =    (psdep(k)+psaut(k)+paacw(k)                                 &
                   -pgaut(k)+piacr(k)*delta3 + praci(k)*delta3                 &
#ifdef HAIL
                   +pvapg(k)+pvaph(k)-phacs(k)                                 &
#endif
                   +psaci(k)-pgacs(k)-pracs(k)*(1.-delta2) + psacr(k)*delta2)  &
                   *dtcld
       qrs(k,i,2) = max(qrs(k,i,2)+htotal, 0.)
       htotal =    (pgdep(k)+pgaut(k)+piacr(k)*(1.-delta3)                     &
                   +praci(k)*(1.-delta3) + psacr(k)*(1.-delta2)                &
                   +pracs(k)*(1.-delta2) + pgaci(k)+paacw(k)                   &
#ifndef HAIL
                   +pgacr(k)+pgacs(k))*dtcld
#else
                   +pgacr(k)*delta2+pgacs(k)+primh(k)                          &
                   -pracg(k)*(1.-delta2)-phacg(k)-phaut(k)                     &
                   -pvapg(k))*dtcld
       qrs(k,i,3) = max(qrs(k,i,3)+htotal, 0.)
       htotal =    (phdep(k)+phaut(k)                                          &
                   +pgacr(k)*(1.-delta2)+pracg(k)*(1.-delta2)                  &
                   +phacw(k)+phacr(k)+phaci(k)+phacs(k)                        &
                   +phacg(k)-pvaph(k)-primh(k))                                &
                   *dtcld
       qrs(k,i,4) = max(qrs(k,i,4)+htotal, 0.)
#endif
#ifdef DOUBLE
#ifdef PRAUT_LB2017
       htotal = (-ncaut(k)-nccol(k)-nracw(k)                                   &
#else
       htotal = (-nraut(k)-nccol(k)-nracw(k)                                   &
#endif
#ifdef HAIL
                   -naacw(k)-naacw(k)-nhacw(k))*dtcld
#else
                   -naacw(k)-naacw(k))*dtcld
#endif
       ncr(k,i,2) = max(ncr(k,i,2)+htotal, 0.)
       htotal =    (nraut(k)-nrcol(k)-niacr(k)                                 &
#ifdef HAIL
                   -nsacr(k)-ngacr(k)-nhacr(k))*dtcld
#else
                   -nsacr(k)-ngacr(k))*dtcld
#endif
       ncr(k,i,3) = max(ncr(k,i,3)+htotal, 0.)
#endif
#ifndef HAIL
       hvalue = - xls*(psdep(k)+pgdep(k)+pidep(k)+pigen(k))                    &
               - xlv(k)*prevp(k) - xlf(k)*(piacr(k)+paacw(k)                   &
               + paacw(k)+pgacr(k)+psacr(k))
#else
       hvalue = - xls*(psdep(k)+pgdep(k)+phdep(k)+pidep(k)+pigen(k))           &
               - xlv(k)*prevp(k) - xlf(k)*(piacr(k)+paacw(k)+phacw(k)          &
               + paacw(k)+pgacr(k)+psacr(k)+phacr(k))
#endif
       t(k,i) = t(k,i) - hvalue/cpm(k)*dtcld
!
     else    /* if t > 0c */
!
! cloud water
!
       hvalue = max(qmin,qci(k,i,1))
#ifndef HAIL
       source = (praut(k)+pracw(k)+paacw(k)+paacw(k))*dtcld
#else
       source = (praut(k)+pracw(k)+paacw(k)+paacw(k)-phacw(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         praut(k) = praut(k) * factor
         pracw(k) = pracw(k) * factor
         paacw(k) = paacw(k) * factor
#ifdef HAIL
         phacw(k) = phacw(k) * factor
#endif
       endif
!
! rain
!
       hvalue = max(qmin,qrs(k,i,1))
#ifndef HAIL
       source = (-paacw(k)-praut(k)+pseml(k)+pgeml(k)-pracw(k)                 &
                 -paacw(k)-prevp(k))*dtcld
#else
       source = (-prevp(k)-praut(k)+pseml(k)+pgeml(k)+pheml(k)                 &
                 -pracw(k)-paacw(k)-paacw(k)-phacw(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         praut(k) = praut(k) * factor
         prevp(k) = prevp(k) * factor
         pracw(k) = pracw(k) * factor
         paacw(k) = paacw(k) * factor
         pseml(k) = pseml(k) * factor
         pgeml(k) = pgeml(k) * factor
#ifdef HAIL
         phacw(k) = phacw(k) * factor
         pheml(k) = pheml(k) * factor
#endif
       endif
!
! snow
!
       hvalue = max(qmin,qrs(k,i,2))
#ifndef HAIL
       source=(pgacs(k)-pseml(k)-psevp(k))*dtcld
#else
       source=(pgacs(k)+phacs(k)-pseml(k)-psevp(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         pgacs(k) = pgacs(k) * factor
         pseml(k) = pseml(k) * factor
         psevp(k) = psevp(k) * factor
#ifdef HAIL
         phacs(k) = phacs(k) * factor
#endif
       endif
!
! graupel
!
       hvalue = max(qmin,qrs(k,i,3))
#ifndef HAIL
       source=-(pgacs(k)+pgevp(k)+pgeml(k))*dtcld
#else
       source=-(pgacs(k)+pgevp(k)+pgeml(k)-phacg(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         pgacs(k) = pgacs(k) * factor
         pgevp(k) = pgevp(k) * factor
         pgeml(k) = pgeml(k) * factor
#ifdef HAIL
         phacg(k) = phacg(k) * factor
#endif
       endif
#ifdef HAIL
!
!     hail
!
       hvalue = max(qrmin,qrs(k,i,4))
       source=-(phacs(k)+phacg(k)+phevp(k)+pheml(k))*dtcld
       if (source>hvalue) then
         factor = hvalue/source
         phacs(k) = phacs(k)*factor
         phacg(k) = phacg(k)*factor
         phevp(k) = phevp(k)*factor
         pheml(k) = pheml(k)*factor
       endif
#endif
#ifdef DOUBLE
!
!     cloud
!
       hvalue = max(qmin,ncr(k,i,2))
#ifndef HAIL
#ifdef PRAUT_LB2017
       source = (+ncaut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k))*dtcld
#else
       source = (+nraut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k))*dtcld
#endif
#else
#ifdef PRAUT_LB2017
       source = (+ncaut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k)+nhacw(k))*dtcld
#else
       source = (+nraut(k)+nccol(k)+nracw(k)+naacw(k)+naacw(k)+nhacw(k))*dtcld
#endif
#endif
       if (source>hvalue) then
         factor = hvalue/source
#ifdef PRAUT_LB2017
         ncaut(k) = ncaut(k) * factor
#else
         nraut(k) = nraut(k) * factor
#endif
         nccol(k) = nccol(k) * factor
         nracw(k) = nracw(k) * factor
         naacw(k) = naacw(k) * factor
#ifdef HAIL
         nhacw(k) = nhacw(k) * factor
#endif
       endif
!
!     rain
!
       hvalue = max(qmin,ncr(k,i,3))
#ifndef HAIL
       source = (-nraut(k)+nrcol(k)-nseml(k)-ngeml(k))*dtcld
#else
       source = (-nraut(k)+nrcol(k)-nseml(k)-ngeml(k)-nheml(k))*dtcld
#endif
       if (source>hvalue) then
         factor = hvalue/source
         nraut(k) = nraut(k) * factor
         nrcol(k) = nrcol(k) * factor
         nseml(k) = nseml(k) * factor
         ngeml(k) = ngeml(k) * factor
#ifdef HAIL
         nheml(k) = nheml(k) * factor
#endif
       endif
#endif
!
! update
!
#ifdef HAIL
       htotal = -(prevp(k)+psevp(k)+pgevp(k)+phevp(k))
#else
       htotal = -(prevp(k)+psevp(k)+pgevp(k))
#endif
       q(k,i) = q(k,i) + htotal*dtcld
#ifdef HAIL
       htotal =   (praut(k)+pracw(k)+paacw(k)+paacw(k)+phacw(k))*dtcld
#else
       htotal =   (praut(k)+pracw(k)+paacw(k)+paacw(k))*dtcld
#endif
       qci(k,i,1) = max(qci(k,i,1)-htotal,0.)

       htotal =    (praut(k)+pracw(k)                                          &
#ifdef HAIL
                   +prevp(k)+paacw(k)+paacw(k)+phacw(k)                        &
                   -pseml(k)-pgeml(k)-pheml(k))*dtcld
#else
                   +prevp(k)+paacw(k)+paacw(k)-pseml(k)-pgeml(k))
#endif
       qrs(k,i,1) = max(qrs(k,i,1)+htotal, 0.)
#ifdef HAIL
       htotal =    (psevp(k)-pgacs(k)-phacs(k)+pseml(k))*dtcld
#else
       htotal =    (psevp(k)-pgacs(k)+pseml(k))*dtcld
#endif
       qrs(k,i,2) = max(qrs(k,i,2)+htotal, 0.)
#ifdef HAIL
       htotal =    (pgacs(k)+pgevp(k)+pgeml(k)-phacg(k))*dtcld
#else
       htotal =    (pgacs(k)+pgevp(k)+pgeml(k))*dtcld
#endif
       qrs(k,i,3) = max(qrs(k,i,3)+htotal, 0.)
#ifdef HAIL
       htotal =    (phacs(k)+phacg(k)+phevp(k)+pheml(k))*dtcld
       qrs(k,i,4) = max(qrs(k,i,4)+htotal, 0.)
#endif
#ifdef DOUBLE
#ifdef PRAUT_LB2017
       htotal =    (-ncaut(k)-nccol(k)-nracw(k)                                &
#else
       htotal =    (-nraut(k)-nccol(k)-nracw(k)                                &
#endif
#ifdef HAIL
                   -naacw(k)-naacw(k)-nhacw(k))*dtcld
#else
                   -naacw(k)-naacw(k))*dtcld
#endif
       ncr(k,i,2) = min(max(ncr(k,i,2)+htotal,0.),ncmax)
#ifdef HAIL
       htotal =    (nraut(k)-nrcol(k)+nseml(k)+ngeml(k)+nheml(k))*dtcld
#else
       htotal =    (nraut(k)-nrcol(k)+nseml(k)+ngeml(k))*dtcld
#endif
       ncr(k,i,3) = min(max(ncr(k,i,3) + htotal,0.),nrmax)
#endif /* DOUBLE */
#ifdef HAIL
       hvalue = -xlv(k)*(prevp(k)+psevp(k)+pgevp(k)+phevp(k))                  &
                -xlf(k)*(pseml(k)+pgeml(k)+pheml(k))
#else
       hvalue = -xlv(k)*(prevp(k)+psevp(k)+pgevp(k))                           &
                -xlf(k)*(pseml(k)+pgeml(k))
#endif
       t(k,i) = t(k,i) - hvalue/cpm(k)*dtcld
     endif
   enddo
#ifdef DOUBLE
#ifdef SIZE_LIMIT3
!===============================================================================
!
! adjustment of nc and nr according to avaialble size of clouds and rain
!
!===============================================================================
   call adjust_number_concent(ktopqr,kdim,qrs(:,i,1),ncr(:,i,3),den(:,i),      &
        pidnr,drcoeff,qrmin,nrmin,nrmax,di1000,drmin,drmax)
   call adjust_number_concent(ktopqc,kdim,qci(:,i,1),ncr(:,i,2),den(:,i),      &
        pidnc,one_1,qcmin,ncmin,ncmax,di15,  dcmin, dcmax)
#endif
!===============================================================================
!
! ccn activaiton and condensation/evaporation of clouds
!
!===============================================================================
#ifdef PRTOC
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,1),zero_0,ktopqr)
   if(ktopqr>0.0) lqr = .true.
   if(lqr) then
     ktop = ktopqr
     call slope_rain(kdim,ktop,qrs(:,i,1),den(:,i),denfac(:),t(:,i),           &
        ncr(:,i,3),                                                            &
        rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:),t0c,den0)
     do k = ktop, kts, -1
       if(qrs(k,i,1)>0.) then
         dr(k) = rslope(k,1)*drcoeff
!-------------------------------------------------------------------------------
! nrtoc: conversion from rain to cloud  [lh a14]
!        (nr->nc)
         if(dr(k)<=di82) then
           ncr(k,i,2) = min(ncr(k,i,2) + ncr(k,i,3), ncmax)
           ncr(k,i,3) = 0.
!-------------------------------------------------------------------------------
! prtoc: conversion from rain to cloud  [lh a15]
!        (qr->qc)
           qci(k,i,1) = qci(k,i,1) + qrs(k,i,1)
           qrs(k,i,1) = 0.
         endif
       endif
     enddo
   endif
#endif
!-------------------------------------------------------------------------------
! pcact: qv -> qc                                 [lh 8]  [kk 14]
! ncact: nccn -> nc                               [lh a2] [kk 12]
!
   ktop = ktopini
   do k = ktop, kts, -1
     qsat_mul(k) = fsvp_water(t(k,i),p(k,i))
     qsat_mul(k) = ep2 * qsat_mul(k) / (p(k,i) - qsat_mul(k))
     qsat_mul(k) = max(qsat_mul(k),qmin)
     rh_mul(k)   = max(q(k,i) / qsat_mul(k),qmin)
   enddo
!
   call find_cloud_top(1,kdim,ktopini,rh_mul(:), one_1,ktoprh)
!
   ktop = ktoprh
   do k = ktop, kts, -1
     if(rh_mul(k)>1.) then
       hvalue = ((rh_mul(k)-1.)/satmax)
       hvalue = min(1.,exp(log(hvalue)*actk))
#ifdef  NCACT_LEI2019
       htotal = ((ncr(k,i,1)+ncr(k,i,2))*hvalue - ncr(k,i,2))*rdtcld
#else
       htotal = ((ncr(k,i,1) + ncr(k,i,2))*min(1.,(rh_mul(k)/satmax)**actk)    &
               - ncr(k,i,2))*rdtcld
#endif
       ncact(k) = max(0.,htotal)
       ncact(k) = min(ncact(k),max(ncr(k,i,1),0.)*rdtcld)
       pcact(k) = min(4.*pi*denr*exp(log(actr)*3)*ncact(k)/                    &
                     (3.*den(k,i)),max(q(k,i),0.)*rdtcld)
       q(k,i) = max(q(k,i) - pcact(k)*dtcld,0.)
       qci(k,i,1) = max(qci(k,i,1) + pcact(k)*dtcld,0.)
       ncr(k,i,1) = max(ncr(k,i,1) - ncact(k)*dtcld, ccnmin)
       ncr(k,i,2) = max(ncr(k,i,2) + ncact(k)*dtcld, ncmin)
       t(k,i) = t(k,i) + pcact(k)*xlv(k)/cpm(k)*dtcld
     endif
   enddo
#endif
!-------------------------------------------------------------------------------
! pcond: condensational/evaporational rate of cloud water [hl a46]
! (qv -> qc, qc->qv)
   ktop = ktopini
   do k = ktop, kts, -1
     qsat_mul(k) = fsvp_water(t(k,i),p(k,i))
     qsat_mul(k) = ep2 * qsat_mul(k) / (p(k,i) - qsat_mul(k))
     qsat_mul(k) = max(qsat_mul(k),qmin)
   enddo
   call find_cloud_top(1,kdim,ktopini,qci(:,i,1),zero_0,ktopqc)
   call find_cloud_top(1,kdim,ktopini,rh_mul(:), one_1,ktoprh)
   ktop = max(ktopqc,ktoprh)
   do k = ktop, kts, -1
     hvalue = ((max(q(k,i),qmin)-(qsat_mul(k))) /(1.+(xlv(k))*(xlv(k))         &
              /(rv*(cpm(k)))*(qsat_mul(k)) /((t(k,i))*(t(k,i)))))
     pcond(k) = min(max(hvalue*rdtcld,0.),max(q(k,i),0.)*rdtcld)
     if(qci(k,i,1)>0. .and. hvalue<0.) then
       pcond(k) = max(hvalue,-qci(k,i,1))*rdtcld
#ifdef DOUBLE
!-------------------------------------------------------------------------------
! ncevp: evpration of cloud number concentration  [lh A3]
! (nc->nccn)
#ifdef EVAP_PARTIAL
       hvalue = max(min(-pcond(k)*ncr(k,i,2)/qci(k,i,1)*dtcld,ncr(k,i,2)),0.)
       ncr(k,i,1) = min(ncr(k,i,1) + hvalue,ccnmax)
       ncr(k,i,2) = max(ncr(k,i,2) - hvalue,0.)
#else
       if(pcond(k)==-qci(k,i,1)*rdtcld) then
         ncr(k,i,1) = ncr(k,i,1) + ncr(k,i,2)
         ncr(k,i,2) = 0.
       endif
#endif
#endif
     endif
     q(k,i) = q(k,i) - pcond(k)*dtcld
     qci(k,i,1) = max(qci(k,i,1)+pcond(k)*dtcld,0.)
     t(k,i) = t(k,i) + pcond(k)*xlv(k)/cpm(k)*dtcld
   enddo
!===============================================================================
!
! sedimentation of hydrometeors
!
!===============================================================================
! re-define cloud top for numerical efficiencies
!
   lqc = .false.
   lqi = .false.
   lqr = .false.
   lqs = .false.
   lqg = .false.
#ifdef HAIL
   lqh = .false.
#endif
   flgcld = .false.
   call find_cloud_top(1,kdim,ktopini,qci(:,i,1),zero_0,ktopqc)
   call find_cloud_top(1,kdim,ktopini,qci(:,i,2),zero_0,ktopqi)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,1),zero_0,ktopqr)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,2),zero_0,ktopqs)
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,3),zero_0,ktopqg)
#ifdef HAIL
   call find_cloud_top(1,kdim,ktopini,qrs(:,i,4),zero_0,ktopqh)
#endif
   call find_cloud_top(1,kdim,ktopini,rh_ice(:), one_1,ktoprh)
   if(ktopqc>0.0) lqc = .true.
   if(ktopqi>0.0) lqi = .true.
   if(ktopqr>0.0) lqr = .true.
   if(ktopqs>0.0) lqs = .true.
   if(ktopqg>0.0) lqg = .true.
#ifdef HAIL
   if(ktopqh>0.0) lqh = .true.
   ktopmax = max(ktopqc,ktopqi,ktopqr,ktopqs,ktopqg,ktopqh,ktoprh)
#else
   ktopmax = max(ktopqc,ktopqi,ktopqr,ktopqs,ktopqg,ktoprh)
#endif
!
! fall out for rain
!
   if(lqr) then
     vtr(:) = 0. ; vtn(:) = 0.
     ktop = ktopqr
#if (defined DOUBLE) && (defined SIZE_LIMIT_SEDI)
     call adjust_number_concent(ktopqr,kdim,qrs(:,i,1),ncr(:,i,3),den(:,i),    &
          pidnr,drcoeff,qrmin,nrmin,nrmax,di1000,drmin,drmax)
#endif
     call slope_rain(kdim,ktop,qrs(:,i,1),den(:,i),denfac(:),t(:,i),           &
#ifdef DOUBLE
        ncr(:,i,3),                                                            &
#endif
        rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:),t0c,den0)
#ifdef PRT
       sum1 = 0.0
       do k = ktop, kts, -1
         if(i==lond.and.lat==latd) then
           print 102, n,k,vtr(k),vtn(k),qrs(k,i,1)*1000,ncr(k,i,3)/1000
           sum1 = sum1 + delz(k,i)*denqr(k)
         endif
       enddo
#endif
#ifdef DOUBLE
     nstep = 1
     do k = ktop, kts, -1
       vtn(k) = pvtrn*rslopeb(k,1)*denfac(k)
       hvalue = max(vtr(k),vtn(k))
       nstep = max(nstep,ceiling(dtcld/delz(k,i)*hvalue))
#ifdef SEDI_VT_BOTTOM
       if(qrs(k,i,1)<=qrmin) vtr(k) = vtr(k+1)
       if(ncr(k,i,3)<=nrmin) vtn(k) = vtn(k+1)
#endif
     enddo
#ifdef  SEDI_VT_AVERAGE
!
     if(ktop>2) then
       do k = ktop-1, kts+1, -1
         hvalue = (2.*vtr(k)+vtr(k+1)+vtr(k-1))*0.25
         vtr(k) = hvalue
       enddo
     endif
#endif /* end of SEDI_AVERAGE */
#ifdef SEDI_SHAPE
!
     do k = ktop, kts, -1
       if(qrs(k,i,1)>qrmin) then
         diameter = drcoeff*rslope(k,1)
         hvalue = fshape(diameter)
         vtn(k) = vtr(k)/hvalue
       else
         vtn(k) = 0.
       endif
     enddo
#else
     if(ktop>2) then
       do k = ktop-1, kts+1, -1
         hvalue = (2.*vtn(k)+vtn(k+1)+vtn(k-1))*0.25
         vtn(k) = hvalue
       enddo
     endif
#endif /* end of SEDI_SHAPE */
#else /* undef DOUBLE */
     nstep = max(1,ceiling(maxval(dtcld/delz(:,i)*vtr(:))))
#ifdef SEDI_VT_BOTTOM
     do k = ktop, kts, -1
       if(qrs(k,i,1)<=qrmin) vtr(k) = vtr(k+1)
     enddo
#endif
#ifdef  SEDI_VT_AVERAGE
     if(ktop>2) then
       do k = ktop-1, kts+1, -1
         vtr(k) = (2*vtr(k)+vtr(k+1)+vtr(k-1))/4.
       enddo
     endif
#endif
#endif
!
#ifdef SEMI_RAIN
     niter = ceiling(nstep/sedi_semi_cfl)
     dtcfl = dtcld/niter
!
     do n = 1, niter
       do k = ktop, kts, -1
#ifdef DOUBLE
         if(qrs(k,i,1)>0.) then
           denqr(k) = dend(k)*qrs(k,i,1)
           denqn(k) = dend(k)*ncr(k,i,3)
         else
           denqr(k) = 0.0
           denqn(k) = 0.0
           vtr(k) = 0.0
           vtn(k) = 0.0
         endif
#else
         denqr(k) = dend(k)*qrs(k,i,1)
#endif
       enddo
#ifdef PRT
       sum1 = 0.0
           print*,' before n k vtr vtn qr nr '
       do k = ktop, kts, -1
         if(i==lond.and.lat==latd) then
           print 102, n,k,vtr(k),vtn(k),qrs(k,i,1)*1000,ncr(k,i,3)/1000
           sum1 = sum1 + delz(k,i)*denqr(k)
         endif
       enddo
#endif
!
!-------------------------------------------------------------------------------
       call semi_lagrangian(1,kdim,max(ktop,2),denfac(:),t(:,i),               &
            delz(:,i),vtr(:),denqr(:),qrpath,dtcfl,lat,i)
#ifdef DOUBLE
       call semi_lagrangian(1,kdim,max(ktop,2),denfac(:),t(:,i),               &
            delz(:,i),vtn(:),denqn(:),qnpath,dtcfl,lat,i)
       do k = ktop, kts, -1
           qrs(k,i,1) = max(denqr(k)/dend(k),0.)
           ncr(k,i,3) = max(denqn(k)/dend(k),0.)
       enddo
#else
       do k = ktop, kts, -1
         qrs(k,i,1) = max(denqr(k)/dend(k),0.)
       enddo
#endif
#ifdef PRT
        if(i==lond.and.lat==latd) print*,' qrpath qnpath ',qrpath, qnpath
       sum2 = 0.0
           print*,' after n k vtr vtn qr nr '
       do k = ktop, kts, -1
         if(i==lond.and.lat==latd) then
           print 102, n,k,vtr(k),vtn(k),qrs(k,i,1)*1000,ncr(k,i,3)/1000
           sum2 = sum2 + delz(k,i)*denqr(k)
         endif
       enddo
       if(i==lond.and. lat==latd) write(85,101) sum1,sum2,sum1-sum2,qrpath
101   format(1x,4e17.10)
#endif
!
       precip_r = qrpath/dtcld + precip_r   ! [kgm-2s-1]
!
     enddo
#else /* no semi_rain */
     niter = nstep
     dtcfl = dtcld/niter
!
     do n = 1, niter
       do k = ktop, kts, -1
         denqr(k) = dend(k)*qrs(k,i,1)
         denqn(k) = ncr(k,i,3)
       enddo
#ifdef PRT
           print*,' before n k vtr vtn qr nr '
       do k = ktop, kts, -1
         if(i==lond.and.lat==latd) then
           print 102, n,k,vtr(k),vtn(k),qrs(k,i,1)*1000,ncr(k,i,3)/1000
         endif
       enddo
#endif
!
       k = ktop
         sedir(k) = denqr(k)*vtr(k)/delz(k,i)/niter
         sedin(k) = denqn(k)*vtn(k)/delz(k,i)/niter
         qrs(k,i,1) = max(qrs(k,i,1)-sedir(k)*dtcld/dend(k),0.)
         ncr(k,i,3) = max(ncr(k,i,3)-sedin(k)*dtcld,0.)
       do k = ktop-1,kts,-1
         sedir(k) = denqr(k)*vtr(k)/delz(k,i)/niter
         sedin(k) = denqn(k)*vtn(k)/delz(k,i)/niter
         dqr(k) = min(sedir(k)*dtcld/dend(k),qrs(k,i,1))
         dqr(k+1) = min(sedir(k+1)*delz(k+1,i)/delz(k,i)*dtcld/dend(k)         &
                    ,qrs(k+1,i,1))
         dnr(k) = min(sedin(k)*dtcld,ncr(k,i,3))
         dnr(k+1) = min(sedin(k+1)*delz(k+1,i)/delz(k,i)*dtcld,ncr(k+1,i,3))
         qrs(k,i,1) = max(qrs(k,i,1)-dqr(k)+dqr(k+1),0.)
         ncr(k,i,3) = max(ncr(k,i,3)-dnr(k)+dnr(k+1),0.)
       enddo
#ifdef PRT
       if(i==lond.and.lat==latd) print*,' qrpath, qnpath ',qrpath, qnpath
       sum2 = 0.0
#endif
#ifdef PRT
           print*,' after n k vtr vtn qr nr 4mmnr'
       do k = ktop, kts, -1
         if(i==lond.and.lat==latd) then
           print 103, n,k,vtr(k),vtn(k),qrs(k,i,1)*1000,ncr(k,i,3)/1000,3.e4*den(k,i)*qrs(k,i,1)
         endif
       enddo
#endif
!
       precip_r = sedir(1)/dend(1)*delz(1,i) + precip_r   ! [kgm-2s-1]
#ifdef PRT
       if(i==lond.and.lat==latd) print*,' precip_r ',precip_r
101   format(1x,4e17.10)
#endif
!
#ifdef SEDI_UPDATE
       call slope_rain(kdim,ktop,qrs(:,i,1),den(:,i),denfac(:),t(:,i),         &
            ncr(:,i,3),                                                        &
            rslope(:,1),rslopeb(:,1),rslope2(:,1),rslope3(:,1),vtr(:),t0c,den0)
       do k = ktop, kts, -1
         vtn(k) = pvtrn*rslopeb(k,1)*denfac(k)
         if(ncr(k,i,3)<=0.0) vtn(k) = 0.0
       enddo
#endif
     enddo
#endif
   endif
#ifdef PRT
 102 format(1x,2i5,4f10.3)
 103 format(1x,2i5,5f10.3)
#endif
#ifdef SEDI_TEST
       do k = kts, kte
       qrs1(i,k,1) = qrs(k,i,1)
#ifdef DOUBLE
       ncr1(i,k,3) = ncr(k,i,3)
#endif
       enddo
   return
#endif
!-------------------------------------------------------------------------------
! fall out for precipitating ice : mass-weighted combined sedimentation (dhl )
!
#ifdef HAIL
   if(lqs .or. lqg .or. lqh) then
#else
   if(lqs .or. lqg) then
#endif
     vts(:) = 0. ; vtg(:) = 0. ; vth(:) = 0. ; vtmean(:) = 0.
     ktop =ktopqs
     call slope_snow(kdim,ktop,qrs(:,i,2),den(:,i),denfac(:),t(:,i),           &
          rslope(:,2),rslopeb(:,2),rslope2(:,2),rslope3(:,2),vts(:),t0c,den0)
     ktop =ktopqg
     call slope_graupel(kdim,ktop,qrs(:,i,3),den(:,i),denfac(:),t(:,i),        &
          rslope(:,3),rslopeb(:,3),rslope2(:,3),rslope3(:,3),vtg(:),t0c,den0)
#ifdef HAIL
     ktop =ktopqh
     call slope_hail(kdim,ktop,qrs(:,i,4),den(:,i),denfac(:),t(:,i),           &
          rslope(:,4),rslopeb(:,4),rslope2(:,4),rslope3(:,4),vth(:),t0c,den0)
#endif
!
#if defined (HAIL) && defined (SEDI_MEAN)
     ktop = max(ktopqs,ktopqg,ktopqh)
     do k = ktop, kts, -1
       sumice(k) = max( (qrs(k,i,2) + qrs(k,i,3) + qrs(k,i,4)), qmin)
       if(sumice(k)>qmin) then
         vtmean(k) = (vts(k)*qrs(k,i,2) + vtg(k)*qrs(k,i,3)                    &
                   + vth(k)*qrs(k,i,4))/sumice(k)
       else
         vtmean(k) = 0.
       endif
     enddo
#else
     ktop = max(ktopqs,ktopqg)
     do k = ktop, kts, -1
       sumice(k) = max( (qrs(k,i,2) + qrs(k,i,3)), qmin)
       if(sumice(k)>qmin) then
         vtmean(k) = (vts(k)*qrs(k,i,2) + vtg(k)*qrs(k,i,3))/sumice(k)
       else
         vtmean(k) = 0.
       endif
     enddo
#endif
#ifdef  SEDI_VT_AVERAGE
     if(ktop>2) then
       do k = ktop-1, kts+1, -1
         vtmean(k) = (2*vtmean(k)+vtmean(k+1)+vtmean(k-1))/4.
       enddo
     endif
#endif
!
     nstep = max(1,ceiling(maxval(dtcld/delz(:,i)*vtmean(:))))
     niter = ceiling(nstep/sedi_semi_cfl)
     dtcfl = dtcld/niter
!
     do n = 1, niter
!
       ktop =ktopqs
       do k = ktop, kts, -1
         denqs(k) = dend(k)*qrs(k,i,2)
       enddo
!-------------------------------------------------------------------------------
       call semi_lagrangian(1,kdim,max(ktop,2),denfac(:),t(:,i),               &
            delz(:,i),vtmean(:),denqs(:),qspath,dtcfl,lat,i)
       do k = kts, ktop
         qrs(k,i,2) = max(denqs(k)/dend(k),0.)
       enddo
       ktop =ktopqg
       do k = ktop, kts, -1
         denqg(k) = dend(k)*qrs(k,i,3)
       enddo
       call semi_lagrangian(1,kdim,max(ktop,2),denfac(:),t(:,i),               &
            delz(:,i),vtmean(:),denqg(:),qgpath,dtcfl,lat,i)
       do k = kts, ktop
         qrs(k,i,3) = max(denqg(k)/dend(k),0.)
       enddo
#ifdef HAIL
       ktop =ktopqh
       do k = ktop, kts, -1
         denqh(k) = dend(k)*qrs(k,i,4)
       enddo
       call semi_lagrangian(1,kdim,max(ktop,2),denfac(:),t(:,i),               &
#ifdef SEDI_MEAN
            delz(:,i),vtmean(:),denqh(:),qhpath,dtcfl,lat,i)
#else
            delz(:,i),vth(:),denqh(:),qhpath,dtcfl,lat,i)
#endif
       do k = kts, ktop
         qrs(k,i,4) = max(denqh(k)/dend(k),0.)
       enddo
#endif
!
       precip_s = qspath/dtcld + precip_s   ! [kgm-2s-1]
       precip_g = qgpath/dtcld + precip_g   ! [kgm-2s-1]
#ifdef HAIL
       precip_h = qhpath/dtcld + precip_h   ! [kgm-2s-1]
#endif
     enddo
   endif
!-------------------------------------------------------------------------------
! fall out for cloud ice
!
   if(lqi) then
     ktop = ktopqi
     do k = ktop, kts, -1
       hvalue = (den(k,i)*max(qci(k,i,2),qcmin))
       hvalue = exp(log(hvalue)*0.75)
       ni(k) = min(max(5.38e7*hvalue,1.e3),1.e6)
       if(qci(k,i,2)<=0.0) then
         vti(k) = 0.
       else
         mi(k)  = den(k,i)*qci(k,i,2)/ni(k)
         di(k)  = max(min(exp(log((mi(k)/cxmi))*(1./dxmi)),dimax), qmin)
         vti(k) = avti*exp(log(di(k))*(bvti))
#ifdef ICE_SPHERE
         hvalue = min(max((tc(k) - t2_sphere)/(t1_sphere-t2_sphere),0.),1.0)
         dis(k) = max(min(exp(log((mi(k)/cxmis))*(1./dxmis)),dimax), qmin)
         vtis(k)= avtis*exp(log(dis(k))*(bvtis))
         di(k)  = di(k)  * (1.-hvalue) + dis(k)  * hvalue
         vti(k) = vti(k) * (1.-hvalue) + vtis(k) * hvalue
#endif
       endif
     enddo
#ifdef  SEDI_VT_AVERAGE
     if(ktop>2) then
       do k = ktop-1, kts+1, -1
         vti(k) = (2*vti(k)+vti(k+1)+vti(k-1))/4.
       enddo
     endif
#endif
!
     nstep = max(1,ceiling(maxval(dtcld/delz(:,i)*vti(:))))
     niter = ceiling(nstep/sedi_semi_cfl)
     dtcfl = dtcld/niter
!
     do n = 1, niter
       do k = ktop, kts, -1
         denqi(k) = dend(k)*qci(k,i,2)
       enddo
!
       call semi_lagrangian(1,kdim,max(ktop,2),denfac(:),t(:,i),               &
            delz(:,i),vti(:),denqi(:),qipath,dtcfl,lat,i)
       do k = kts, ktop
         qci(k,i,2) = max(denqi(k)/dend(k),0.)
       enddo
!
       precip_i = qipath/dtcld + precip_i   ! [kgm-2s-1]
     enddo
   endif
!===============================================================================
!
! precip (den*qrsi*dz/dt) : [kgm-2s-1] ==> rain ( precip/denr*dt*1000 ) : [mm ]
! for wrf unit is mm, whereas it is m for ufs
!
!===============================================================================
#ifdef WRF
   hvalue = 1000.
#else
   hvalue = 1.
#endif
#ifdef HAIL
   precip_sum = precip_r + precip_s + precip_i + precip_g + precip_h
   precip_ice = precip_s + precip_i + precip_h
#else
   precip_sum = precip_r + precip_s + precip_i + precip_g
   precip_ice = precip_s + precip_i
#endif
!
   if(precip_sum>0.) then
     hvalue =  precip_sum/denr*dtcld*hvalue
#ifdef WRF
     rainncv(i) = hvalue + rainncv(i)
#endif
     rain(i)    = hvalue + rain(i)
   endif
!
   if(precip_ice>0.) then
     hvalue = precip_ice/denr*dtcld*hvalue
     tstepsnow   = hvalue + tstepsnow
#ifdef WRF
     if ( present (snowncv) .and. present (snow)) then
       snowncv(i) = hvalue + snowncv(i)
#else
     if ( present (snow)) then
#endif
       snow(i)    = hvalue + snow(i)
      endif
   endif
!
   if(precip_g>0.) then
     hvalue = precip_g/denr*dtcld*hvalue
     tstepgraupel   = hvalue + tstepgraupel
#ifdef WRF
     if ( present (graupelncv) .and. present (graupel)) then
       graupelncv(i) = hvalue + graupelncv(i)
#else
     if ( present (graupel)) then
#endif
       graupel(i)    = hvalue + graupel(i)
      endif
   endif
#ifdef HAIL
!
   if(precip_h>0.) then
     hvalue = precip_h/denr*dtcld*hvalue
     tstephail   = hvalue + tstephail
#ifdef WRF
     if ( present (hailncv) .and. present (hail)) then
       hailncv(i) = hvalue + hailncv(i)
#else
     if ( present (hail)) then
#endif
       hail(i)    = hvalue + hail(i)
      endif
   endif
#endif
!
#ifdef WRF
   if ( present (snow) ) then
     if(precip_sum>0.) sr(i) = snowncv(i)/(rainncv(i) + qmin)
   else
     if(precip_sum>0.) sr(i) = tstepsnow/(rainncv(i) + qmin)
   endif
#endif
!-------------------------------------------------------------------------------
!
   enddo i_loop                  ! i-loops
!
!===============================================================================
!
   enddo inner_loop                    ! dtcldcr- loops
!
!===============================================================================
!
#ifndef WRF
!-------------------------------------------------------------------------------
!  compute reflectivity and effective radius 
!
   if (mod(kdt,(int(fhswr)/int(delt)))==1) then
   do i = its, ite
     dbz(:) = 0.; re_qc(:) = 0.; re_qi(:) = 0.; re_qs(:) = 0.
     call udm_mp_reflectivity (q(:,i), qrs(:,i,1), qrs(:,i,2), qrs(:,i,3),     &
#ifdef DOUBLE
                         ncr(:,i,3),                                           &
#else
                         t(:,i),                                               &
#endif
#ifdef HAIL
                         qrs(:,i,4),                                           &
#endif
                         t(:,i), p(:,i), dbz, kts, kte, i, i)
     do k = kts, kte
       refl_10cm(i,k) = max(-35., dbz(k))
     enddo
!
     call udm_mp_effective_radius (t(:,i), qci(:,i,1), qci(:,i,2), qrs(:,i,2), &
                         den(:,i), qmin, t0c,                                  &
#ifdef DOUBLE
                         ncr(:,i,2),                                           &
#else
                         t(:,i),                                               &
#endif
                         re_qc, re_qi, re_qs, kts, kte, i, i)
     do k = kts, kte
       re_cloud(i,k) = max(recmin, min(re_qc(k), recmax))*1.e6
       re_ice(i,k)   = max(reimin, min(re_qi(k), reimax))*1.e6
       re_snow(i,k)  = max(resmin, min(re_qs(k), resmax))*1.e6
     enddo
   enddo
   endif
#endif
!-------------------------------------------------------------------------------
! assign local to passing variables
!
   ktop = kte
   do k = kts, ktop
     do i = its, ite
       t1(i,k) = t(k,i)
       q1(i,k) = q(k,i)
#ifdef DOUBLE
       nn1(i,k) = ncr(k,i,1)
       nc1(i,k) = ncr(k,i,2)
       nr1(i,k) = ncr(k,i,3)
#endif
     enddo
   enddo
!
   do k = kts, ktop
     do i = its, ite
       qc1(i,k) = qci(k,i,1)
       qr1(i,k) = qrs(k,i,1)
       qi1(i,k) = qci(k,i,2)
       qs1(i,k) = qrs(k,i,2)
       qg1(i,k) = qrs(k,i,3)
#ifdef HAIL
       qh1(i,k) = qrs(k,i,4)
#endif
     enddo
   enddo
!
#ifdef wrf_chem
      rainprod2d = praut+pracw+praci+psaci+pgaci+psacw+pgacw+paacw+psaut
      evapprod2d = -(prevp+psevp+pgevp+psdep+pgdep)
#endif
!
   end subroutine MP_RUN_
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_rain(kdim, ktop, qrs, den, denfac, t,                      &
#ifdef DOUBLE
                         ncr,                                                  &
#endif
                         rslope, rslopeb, rslope2, rslope3, vt, t0c, den0)
!-------------------------------------------------------------------------------
   implicit none
!
   integer,              intent(in   ) :: kdim, ktop
   real(kind=kind_phys), intent(in   ) :: t0c, den0 
   real(kind=kind_phys), dimension(kdim), &
        intent(in   ) ::                                                       &
                                                                          qrs, &
#ifdef DOUBLE
                                                                          ncr, &
#endif
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), dimension(kdim), &
        intent(inout   ) ::                                                    &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt
  real(kind=kind_phys)       ::  lamda, lamdar1m, lamdar2m, x, y, z, hvalue
  integer :: k
!-------------------------------------------------------------------------------
!  size distributions: (x=mixing ratio, y=air density):
#ifdef DOUBLE
   lamdar2m(x,y,z)= exp(log(((pidnr*z)/(x*y)))*((.33333333)))
#endif
   lamdar1m(x,y)=   exp(log(pidn0r/(x*y))*(0.25))      ! (pidn0r/(x*y))**.25
!
   do k = 1,ktop
#ifdef DOUBLE
     if(qrs(k)<=qrmin .or. ncr(k)<=nrmin ) then
#else
     if(qrs(k)<=qrmin) then
#endif
       rslope(k) = rslopermax
       rslopeb(k) = rsloperbmax
       rslope2(k) = rsloper2max
       rslope3(k) = rsloper3max
     else
#ifdef DOUBLE
       lamda = min(max(lamdar2m(qrs(k),den(k),ncr(k)),lamdarmin),lamdarmax)
#ifdef LAMDAR_1M
        if(qrs(k)<=0.1e-3)then
          hvalue = lamdar1m(qrs(k),den(k))
          lamda = max(lamda,hvalue)
        endif
#endif
#else
       lamda = min(max(lamdar1m(qrs(k),den(k)),lamdarmin),lamdarmax)
#endif
       rslope(k) = 1./lamda
       rslopeb(k) = exp(log(rslope(k))*(bvtr))
       rslope2(k) = rslope(k)*rslope(k)
       rslope3(k) = rslope2(k)*rslope(k)
     endif
     vt(k) = pvtr*rslopeb(k)*denfac(k)
   enddo
end subroutine slope_rain
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_graupel(kdim, ktop, qrs, den, denfac, t, rslope,           &
                         rslopeb, rslope2, rslope3, vt, t0c, den0)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer,              intent(in   ) :: kdim, ktop
   real(kind=kind_phys), intent(in   ) :: t0c, den0
   real(kind=kind_phys), dimension( kdim),                                     &
        intent(in   ) ::                                                       &
                                                                          qrs, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), dimension(kdim),                                       &
        intent(inout   ) ::                                                    &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt
#ifdef NEW_N0
   real(kind=kind_phys), dimension( kdim ) ::                                  &
                                                                       n0gfac
#endif
  real(kind=kind_phys)       ::  lamda, lamdag, x, y, z, tc, hvalue
  integer :: k
!----------------------------------------------------------------
#ifdef NEW_N0
!-------------------------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
   lamdag(x,y,z)= exp(log((pidn0g*z)/(x*y))*(0.25))  ! (pidn0g*z/(x*y))**.25
!
   n0gfac(:) = 0.0
   do k = 1, ktop
     tc = t(k) - t0c
     hvalue = exp(log((max(qrs(k),q_ref)/q_ref))*(beta))
     n0gfac(k) = exp(-alpha_qg*tc)*hvalue
     n0gfac(k) = min(max(n0gfac(k),n0gmin/n0g),n0gmax/n0g)
     if(qrs(k)<=qrmin)then
       rslope(k) = rslopegmax
       rslopeb(k) = rslopegbmax
       rslope2(k) = rslopeg2max
       rslope3(k) = rslopeg3max
     else
       lamda = min(lamdag(qrs(k),den(k),n0gfac(k)),lamdagmax)
       rslope(k) = 1./lamda
       rslopeb(k) = exp(log(rslope(k))*(bvtg))
       rslope2(k) = rslope(k)*rslope(k)
       rslope3(k) = rslope2(k)*rslope(k)
     endif
     vt(k) = pvtg*rslopeb(k)*denfac(k)
   enddo
#else
!  size distributions: (x=mixing ratio, y=air density):
   lamdag(x,y)=   exp(log(pidn0g/(x*y))*(0.25))      ! (pidn0g/(x*y))**.25
!
   do k = 1,ktop
     if(qrs(k)<=qrmin)then
       rslope(k) = rslopegmax
       rslopeb(k) = rslopegbmax
       rslope2(k) = rslopeg2max
       rslope3(k) = rslopeg3max
     else
       lamda = min(lamdag(qrs(k),den(k)),lamdagmax)
       rslope(k) = 1./lamda
       rslopeb(k) = exp(log(rslope(k))*(bvtg))
       rslope2(k) = rslope(k)*rslope(k)
       rslope3(k) = rslope2(k)*rslope(k)
     endif
     vt(k) = pvtg*rslopeb(k)*denfac(k)
   enddo
#endif
end subroutine slope_graupel
#ifdef HAIL
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_hail(kdim, ktop, qrs, den, denfac, t, rslope,              &
                         rslopeb, rslope2, rslope3, vt, t0c, den0)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
   integer,              intent(in   ) :: kdim, ktop
   real(kind=kind_phys), intent(in   ) :: t0c, den0
   real(kind=kind_phys), dimension(kdim),                                      &
        intent(in   ) ::                                                       &
                                                                          qrs, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), dimension(kdim),                                       &
        intent(inout   ) ::                                                    &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt
#ifdef NEW_N0
   real(kind=kind_phys), dimension(kdim ) ::                                   &
                                                                       n0hfac
#endif
  real(kind=kind_phys)       ::  lamda, lamdah, x, y, z, tc, hvalue
  integer ::  k
#ifdef NEW_N0
!-------------------------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
   lamdah(x,y,z)= exp(log((pidn0h*z)/(x*y))*(0.25))  ! (pidn0h*z/(x*y))**.25
!
   n0hfac(:) = 0.0
     do k = 1, ktop
       tc = t(k) - t0c
       hvalue = exp(log((max(qrs(k),q_ref)/q_ref))*(beta))
       n0hfac(k) = exp(-alpha_qh*tc)*hvalue
       n0hfac(k) = min(max(n0hfac(k),n0hmin/n0h),n0hmax/n0h)
       if(qrs(k)<=qrmin)then
         rslope(k) = rslopehmax
         rslopeb(k) = rslopehbmax
         rslope2(k) = rslopeh2max
         rslope3(k) = rslopeh3max
       else
         lamda = min(lamdah(qrs(k),den(k),n0hfac(k)),lamdahmax)
         rslope(k) = 1./lamda
         rslopeb(k) = exp(log(rslope(k))*(bvth))
         rslope2(k) = rslope(k)*rslope(k)
         rslope3(k) = rslope2(k)*rslope(k)
       endif
       vt(k) = pvth*rslopeb(k)*denfac(k)
     enddo
#else
!-------------------------------------------------------------------------------
!  size distributions: (x=mixing ratio, y=air density):
   lamdah(x,y)=   exp(log(pidn0h/(x*y))*(0.25))      ! (pidn0h/(x*y))**.25
!
     do k = 1,ktop
       if(qrs(k)<=qrmin)then
         rslope(k) = rslopehmax
         rslopeb(k) = rslopehbmax
         rslope2(k) = rslopeh2max
         rslope3(k) = rslopeh3max
       else
         lamda = min(lamdah(qrs(k),den(k)),lamdahmax)
         rslope(k) = 1./lamda
         rslopeb(k) = exp(log(rslope(k))*(bvth))
         rslope2(k) = rslope(k)*rslope(k)
         rslope3(k) = rslope2(k)*rslope(k)
       endif
       vt(k) = pvth*rslopeb(k)*denfac(k)
     enddo
#endif
end subroutine slope_hail
#endif
#ifdef DOUBLE
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine slope_cloud(kdim, ktop, qrs, ncr, den, denfac, t, qmin,          &
                          rslope, rslope2, rslope3, t0c, den0)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer,              intent(in   ) :: kdim, ktop
   real(kind=kind_phys), intent(in   ) :: t0c, den0
   real(kind=kind_phys), dimension(kdim),                                      &
        intent(in   ) ::                                                       &
                                                                          qrs, &
                                                                          ncr, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
  real(kind=kind_phys), dimension(kdim),                                       &
        intent(inout   ) ::                                                    &
                                                                       rslope, &
                                                                      rslope2, &
                                                                      rslope3
  real(kind=kind_phys)       ::  lamda, lamdac, x, y, z, tc, qmin
  integer ::  k
!-------------------------------------------------------------------------------
!  size distributions: (x=mixing ratio, y=air density):
   lamdac(x,y,z)= exp(log(((pidnc*z)/(x*y)))*((.33333333)))
!
   do k = 1,ktop
     if(qrs(k)<=qmin .or. ncr(k)<=ncmin )then
       rslope(k) = rslopecmax
       rslope2(k) = rslopec2max
       rslope3(k) = rslopec3max
     else
       lamda = min(max(lamdac(qrs(k),den(k),ncr(k)),lamdacmin),lamdacmax)
       rslope(k) = 1./lamda
       rslope2(k) = rslope(k)*rslope(k)
       rslope3(k) = rslope2(k)*rslope(k)
     endif
   enddo
end subroutine slope_cloud
#endif
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine slope_snow(kdim, ktop, qrs, den, denfac, t, rslope, rslopeb,        &
                            rslope2, rslope3, vt, t0c, den0)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer,              intent(in   ) :: kdim, ktop
   real(kind=kind_phys), intent(in   ) :: t0c, den0
   real(kind=kind_phys), dimension(kdim),                                      &
         intent(in   ) ::                                                      &
                                                                          qrs, &
                                                                          den, &
                                                                       denfac, &
                                                                            t
   real(kind=kind_phys), dimension(kdim),                              &
         intent(inout   ) ::                                                   &
                                                                       rslope, &
                                                                      rslopeb, &
                                                                      rslope2, &
                                                                      rslope3, &
                                                                           vt
   real(kind=kind_phys), dimension( kdim ) ::                                   &
                                                                       n0sfac
   real(kind=kind_phys)       ::  lamda, lamdas, x, y, z, tc, hvalue
   integer ::  k
!-------------------------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
   lamdas(x,y,z)= exp(log((pidn0s*z)/(x*y))*(0.25))  ! (pidn0s*z/(x*y))**.25
!
   n0sfac(:) = 0.0
   do k = 1, ktop
     tc = t(k) - t0c
!
! n0s: intercept parameter for snow [m-4] [hdc 6]
!
#ifdef NEW_N0
     hvalue = exp(log((max(qrs(k),q_ref)/q_ref))*(beta))
     n0sfac(k) = exp(-alpha_qs*tc)*hvalue
     n0sfac(k) = min(max(n0sfac(k),n0smin/n0s),n0smax/n0s)
#else
     n0sfac(k) = max(min(exp(-alpha*tc),n0smax/n0s),1.)
#endif
     if(qrs(k)<=qrmin)then
       rslope(k) = rslopesmax
       rslopeb(k) = rslopesbmax
       rslope2(k) = rslopes2max
       rslope3(k) = rslopes3max
     else
       lamda = min(lamdas(qrs(k),den(k),n0sfac(k)),lamdasmax)
       rslope(k) = 1./lamda
       rslopeb(k) = exp(log(rslope(k))*(bvts))
       rslope2(k) = rslope(k)*rslope(k)
       rslope3(k) = rslope2(k)*rslope(k)
     endif
     vt(k) = pvts*rslopeb(k)*denfac(k)
   enddo
end subroutine slope_snow
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine semi_lagrangian(im, km, ktop, denfacl, tkl, dzl,                 &
                                  wwl, rql, precip, dt, lat, lon)
!-------------------------------------------------------------------------------
!
! this routine is a semi-lagrangain forward advection for hydrometeors
! with mass conservation and positive definite advection
!
!  input :
!    im, km   - dimension
!    denfacl  - sqrt(den0/den)
!    tkl      - temperature in k
!    dzl      - depth of model layer in meter
!    wwl      - terminal velocity at model layer in m/s
!    dt       - time step
!    ktop     - top layer for computing
!
!  inout :
!    precip   - precipitation
!    rql      - cloud density*mixing ratio
!
! author: hann-ming henry juang <henry.juang@noaa.gov>
!         song-you hong <songyou.hong@noaa.gov>
! reference: Juang and  Hong, 2010: Forward semi-Lagrangian advection
!         with mass conservation and positive definiteness for falling
!         hydrometeors. Mon.  Wea. Rev., 138, 1778-1791
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer               , intent(in   ) :: im, km, ktop
   integer               , intent(in   ) :: lat, lon
   real(kind=kind_phys)                  , intent(in   ) :: dt
   real(kind=kind_phys), dimension(im,km), intent(in   ) :: tkl
   real(kind=kind_phys), dimension(im,km), intent(in   ) :: dzl
   real(kind=kind_phys), dimension(km), intent(in   )    :: denfacl
   real(kind=kind_phys), dimension(km), intent(in   )    :: wwl
   real(kind=kind_phys), dimension(km), intent(inout)    :: rql
   real(kind=kind_phys), intent(inout)                   :: precip
!
!  local variables
!
   real(kind=kind_phys), dimension(km)   ::                 &
                                                        dz, &
                                                        ww, &
                                                        qq, &
                                                        wd, &
                                                        wa, &
                                                    denfac, &
                                                        tk, &
                                                        qn
   real(kind=kind_phys), dimension(km+1) ::             wi, &
                                                        zi, &
                                                       qmi, &
                                                       qpi, &
                                                       dza, &
                                                        qa
   real(kind=kind_phys), dimension(km+2) ::             za
!
   integer               :: i, k, n, m, kk, kb, kt, lond, latd
   real(kind=kind_phys)  ::                                           &
                            tl, tl2, qql, dql, qqd                   ,&
                            th, th2, qqh, dqh                        ,&
                            zsum, qsum, dim, dip                     ,&
                            zsumt, qsumt, zsumb, qsumb               ,&
                            allold, allnew, zz, dzamin, cflmax, decfl
   real(kind=kind_phys)  ::      tmp
!
   real(kind=kind_phys), parameter       ::                           &
                            cfac = 0.05, fa1 = 9./16., fa2 = 1./16.
!
#if defined (DEBUG) || defined (PRT)
   integer  :: iprt
!
!===============================================================================
!     set iprt = 0 for no unit fort.84 output
!
!     lond = im/2 + 1
   lond = 101
   latd = 1
#else
   lond = 101
   latd = 1
#endif
!
   zi(:) = 0.    ; wi(:) = 0.   ; qa(:) = 0.
   wa(:) = 0.    ; wd(:) = 0.   ; dza(:) = 0.
   precip = 0.0  ; tmp = 0.0
!
!-------------------------------------------------------------------------------
!
   semi_loop : do i = 1,im
#if defined (DEBUG) || defined (PRT)
     iprt = 0
     if(lon==lond.and.lat==latd) iprt = 1
#endif
!
! assign local variables
     dz(:) = dzl(i,:)
     tk(:) = tkl(i,:)
     denfac(:) = denfacl(:)
     qq(:) = rql(:)
     ww(:) = wwl(:)
!
! skip for no precipitation for all layers
     allold = 0.0
     do k = 1,ktop
       allold = allold + qq(k)
     enddo
     if(allold<=0.0) then
       cycle semi_loop
     endif
!
! compute interface value
     zi(1)=0.0
     do k = 1,ktop
       zi(k+1) = zi(k)+dz(k)
     enddo
!
! save departure wind
     wd(:) = ww(:)
     n=1
     wi(1) = ww(1)
     wi(ktop+1) = ww(ktop)
     do k=2,ktop
       wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
     enddo
!
! 3rd order interpolation to get wi
     wi(1) = ww(1)
     wi(2) = 0.5*(ww(2)+ww(1))
     do k=3,ktop-1
       wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
     enddo
     wi(ktop) = 0.5*(ww(ktop)+ww(ktop-1))
     wi(ktop+1) = ww(ktop)
!
! terminate of top of raingroup
     do k=2,ktop
       if( ww(k)==0.0 ) wi(k)=ww(k-1)
     enddo
!
! diffusivity of wi
     do k=ktop,1,-1
       decfl = (wi(k+1)-wi(k))*dt/dz(k)
       if( decfl > cfac ) then
         wi(k) = wi(k+1) - cfac*dz(k)/dt
       endif
     enddo
!
! compute arrival point
     do k=1,ktop+1
       za(k) = zi(k) - wi(k)*dt
     enddo
     za(ktop+2) = zi(ktop+1)   !hmhj
!
     do k=1,ktop+1  !hmhj
       dza(k) = za(k+1)-za(k)
     enddo
!
! compute deformation at arrival point
     do k=1,ktop
       tmp = qq(k)*dz(k)/dza(k)
       qa(k) = tmp
     enddo
     qa(ktop+1) = 0.0
!
! estimate value at arrival cell interface with monotone
     do k=2,ktop
       dip=(qa(k+1)-qa(k))/(dza(k+1)+dza(k))
       dim=(qa(k)-qa(k-1))/(dza(k-1)+dza(k))
       if( dip*dim<=0.0 ) then
         qmi(k)=qa(k)
         qpi(k)=qa(k)
       else
         qpi(k)=qa(k)+0.5*(dip+dim)*dza(k)
         qmi(k)=2.0*qa(k)-qpi(k)
         if( qpi(k)<0.0 .or. qmi(k)<0.0 ) then
           qpi(k) = qa(k)
           qmi(k) = qa(k)
         endif
       endif
     enddo
     qpi(1)=qa(1)
     qmi(1)=qa(1)
     qmi(ktop+1)=qa(ktop+1)
     qpi(ktop+1)=qa(ktop+1)
#ifdef DEBUG
      if(iprt==1) then
        write(84,103)
      endif
103   format('  k','        zi','        dz','        ww','        wi','        za',&
             '       dza','   qq*1000','   qa*1000')
      if(iprt==1) then
      do k=1,ktop
        write(84,104) k,zi(k),dz(k),ww(k),wi(k),za(k),dza(k),qq(k)*1.e3,qa(k)*1.e3
      enddo
      endif
104   format(i3,8f10.3)
      if(iprt==1) then
        write(84,102)
      endif
102   format('   k','  kt','  kb',' t-b','    qsum','    zsum', &
             ' qpi(kt)',' qpi(kb)',' qmi(kt)',' qmi(kb)',' qq*1000',' qn*1000')
#endif
!
! interpolation to regular point
     qn = 0.0
     kb=1
     kt=1
     intp : do k=1,ktop
            kb=max(kb-1,1)
            kt=max(kt-1,1)
! find kb and kt
            if( zi(k)>=za(ktop+1) ) then
              exit intp
            else
              find_kb : do kk=kb,ktop
                        if( zi(k)<=za(kk+1) ) then
                          kb = kk
                          exit find_kb
                        else
                          cycle find_kb
                        endif
              enddo find_kb
              find_kt : do kk=kt,ktop+2    !hmhj
                        if( zi(k+1)<=za(kk) ) then
                          kt = kk
                          exit find_kt
                        else
                          cycle find_kt
                        endif
              enddo find_kt
              kt = kt - 1
!
! compute q with piecewise parabolic method
              if( kt==kb ) then
                tl=(zi(k)-za(kb))/dza(kb)
                th=(zi(k+1)-za(kb))/dza(kb)
                tl2=tl*tl
                th2=th*th
                qqd=0.5*(qpi(kb)-qmi(kb))
                qqh=qqd*th2+qmi(kb)*th
                qql=qqd*tl2+qmi(kb)*tl
               qn(k) = (qqh-qql)/(th-tl)
              else if( kt>kb ) then
                tl=(zi(k)-za(kb))/dza(kb)
                tl2=tl*tl
                qqd=0.5*(qpi(kb)-qmi(kb))
                qql=qqd*tl2+qmi(kb)*tl
                dql = qa(kb)-qql
                zsum  = (1.-tl)*dza(kb)
                qsum  = dql*dza(kb)
                if( kt-kb>1 ) then
                do m=kb+1,kt-1
                  zsum = zsum + dza(m)
                  qsum = qsum + qa(m) * dza(m)
                enddo
                endif
                th=(zi(k+1)-za(kt))/dza(kt)
                th2=th*th
                qqd=0.5*(qpi(kt)-qmi(kt))
                dqh=qqd*th2+qmi(kt)*th
                zsum  = zsum + th*dza(kt)
                qsum  = qsum + dqh*dza(kt)
                qn(k) = qsum/zsum
              endif
#ifdef DEBUG
      if(iprt==1) then
        write(84,101) k,kt,kb,kt-kb,qsum*1000.,zsum,qpi(kt)*1000.,qpi(kb)*1000., &
                      qmi(kt)*1000.,qpi(kb)*1000.,qq(k)*1000.,qn(k)*1000.
      endif
101   format(4i4,8f8.3)
#endif
               cycle intp
             endif
!
     enddo intp
!
! rain out
     sum_precip: do k=1,ktop
                   if( za(k)<0.0 .and. za(k+1)<=0.0 ) then          !hmhj
                      precip = precip + qa(k)*dza(k)
                      cycle sum_precip
                   else if ( za(k)<0.0 .and. za(k+1)>0.0 ) then    !hmhj
                      th = (0.0-za(k))/dza(k)               !hmhj
                      th2 = th*th                           !hmhj
                      qqd = 0.5*(qpi(k)-qmi(k))             !hmhj
                      qqh = qqd*th2+qmi(k)*th               !hmhj
                      precip = precip + qqh*dza(k)    !hmhj
                      exit sum_precip
                   endif
                   exit sum_precip
     enddo sum_precip
!
! replace the new values
     rql(:) = qn(:)
!
! ----------------------------------
   enddo semi_loop
!
end subroutine semi_lagrangian
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine cldf_diag(kdim, t, p, q, qc, qi, dx,cldf, ktop)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
!  abstract: computes the cloudiness from cloud ice and cloud water. 
!
!  references : Park et al. (2016, m.w.r.) Park and Hong (2024, meteor.)
!
   integer,                     intent(in   ) :: kdim, ktop
   real(kind=kind_phys),        intent(in   ) :: dx
   real(kind=kind_phys), dimension(kdim), intent(in   ) :: t, p, q
   real(kind=kind_phys), dimension(kdim), intent(in   ) :: qc
   real(kind=kind_phys), dimension(kdim), intent(in   ) :: qi
   real(kind=kind_phys), dimension(kdim), intent(out  ) :: cldf
   ! local
   integer :: i, k, kk
   real(kind=kind_phys), parameter :: cldmin = 50., cldmax = 100.
   real(kind=kind_phys), parameter :: clddiff = cldmax - cldmin
   real(kind=kind_phys), parameter :: cldf_min = 0.5
   real(kind=kind_phys) :: cv_w_min,cv_i_min,cvf_min
   real(kind=kind_phys) :: cv_w_max,cv_i_max,cvf_max
   real(kind=kind_phys) :: dxkm
!
   do k = 1,ktop
     cldf(k) = 0.
   enddo
! diagnostic method (kiaps)
   do k = 1,ktop
     cv_w_min = 4.82*(max(0.,qc(k))*1000.)**0.94/1.04
     cv_w_max = 5.77*(max(0.,qc(k))*1000.)**1.07/1.04
     cv_i_min = 4.82*(max(0.,qi(k))*1000.)**0.94/0.96
     cv_i_max = 5.77*(max(0.,qi(k))*1000.)**1.07/0.96
     cvf_min = cv_i_min+cv_w_min
     cvf_max = cv_i_max+cv_w_max
     dxkm    = dx/1000.
     cldf(k)= ((dxkm-cldmin)*(cvf_max-cvf_min)+clddiff*cvf_min)/clddiff
   enddo
   do k = 1,ktop
     cldf(k)=min(1.,max(0.,cldf(k)))
     if(qc(k)+qi(k)<1.e-6) then
       cldf(k) = 0.
     endif
     if(cldf(k)<0.01) cldf(k) = 0.
     if(cldf(k)>0.99) cldf(k) = 1.
     if(cldf(k)>=0.01 .and. cldf(k)<cldf_min) cldf(k) = cldf_min ! min
   enddo
!
end subroutine cldf_diag
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine udminit(den0,denr,dens,cl,cv,ccn0)
!-------------------------------------------------------------------------------
  implicit none
!-------------------------------------------------------------------------------
   real(kind=kind_phys), intent(in) :: den0,denr,dens,cl,cv,ccn0
   integer :: n
!
   pi = 4.*atan(1.)
   pisq = pi * pi
   xlv1 = cl-cv
!
#ifdef DOUBLE
#ifdef HDC2004_PRAUT
   qc_land   = 4./3.*pi*denr*r0**3.*xncr1/den0  ! 500 cm-3, --- 0.837 gkg-1
   qc_ocean  = 4./3.*pi*denr*r0**3.*xncr0/den0  !  50 cm-3, --- 0.0837 gkg-1
   qck1 = .104*9.8*peaut/(denr)**(1./3.)/xmyu*den0**(4./3.) !  4706.08203
#endif
#else
   qc_land   = 4./3.*pi*denr*r0**3.*xncr/den0  ! xncr 300 cm-3 --- 0.502  gkg-1
   qc_ocean  = qc_land
   qck1 = .104*9.8*peaut/(xncr*denr)**(1./3.)/xmyu*den0**(4./3.) ! 7.03
#endif
   pidnc = pi*denr/6.
!
   bvtr1 = 1.+bvtr
#ifdef DOUBLE
   bvtr2 = 2.+bvtr
#else
   bvtr2 = 2.5+.5*bvtr
#endif
   bvtr3 = 3.+bvtr
   bvtr4 = 4.+bvtr
   bvtr6 = 6.+bvtr
   g1pbr = rgmma(bvtr1)
   g3pbr = rgmma(bvtr3)
   g4pbr = rgmma(bvtr4)            ! 17.837825
   g6pbr = rgmma(bvtr6)
#ifdef DOUBLE
   bvtr5 = 5.+bvtr
   bvtr7 = 7.+bvtr
   bvtr2o5 = 2.5+.5*bvtr
   bvtr3o5 = 3.5+.5*bvtr
   g2pbr = rgmma(bvtr2)
   g5pbr = rgmma(bvtr5)
   g7pbr = rgmma(bvtr7)
   g5pbro2 = rgmma(bvtr2o5)
   g7pbro2 = rgmma(bvtr3o5)
   pvtr = avtr*g5pbr/24.
#ifdef PVTRN_ONEHALF
   pvtrn = avtr*g2pbr*1.5
#else
   pvtrn = avtr*g2pbr
#endif
#else
   g5pbro2 = rgmma(bvtr2)          ! 1.8273
   pvtr = avtr*g4pbr/6.
#endif
   eacrr = 1.0
   pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
   pidn0r =  pi*denr*n0r
#ifdef DOUBLE
   pidnr  = 4.*pi*denr
   precr1 = 2.*pi*1.56
   precr2 = 2.*pi*.31*avtr**.5*g7pbro2
#else
   precr1 = 2.*pi*n0r*.78
   precr2 = 2.*pi*n0r*.31*avtr**.5*g5pbro2
#endif
   roqimax = 2.08e22*dimax**8
!
   bvts1 = 1.+bvts
   bvts2 = 2.5+.5*bvts
   bvts3 = 3.+bvts
   bvts4 = 4.+bvts
   g1pbs = rgmma(bvts1)    !.8875
   g3pbs = rgmma(bvts3)
   g4pbs = rgmma(bvts4)    ! 12.0786
   g5pbso2 = rgmma(bvts2)
   pvts = avts*g4pbs/6.
   precs1 = 4.*n0s*.65
   precs2 = 4.*n0s*.44*avts**.5*g5pbso2
   pidn0r =  pi*denr*n0r
   pidn0s =  pi*dens*n0s
   pacrs = pi*n0s*avts*g3pbs*.25*eacrc
!
   bvtg1 = 1.+bvtg
   bvtg2 = 2.5+.5*bvtg
   bvtg3 = 3.+bvtg
   bvtg4 = 4.+bvtg
   g1pbg = rgmma(bvtg1)
   g3pbg = rgmma(bvtg3)
   g4pbg = rgmma(bvtg4)
   pacrg = pi*n0g*avtg*g3pbg*.25
   g5pbgo2 = rgmma(bvtg2)
   pvtg = avtg*g4pbg/6.
   precg1 = 2.*pi*n0g*.78
   precg2 = 2.*pi*n0g*.31*avtg**.5*g5pbgo2
   pidn0g =  pi*deng*n0g
!
   rslopermax = 1./lamdarmax
   rslopesmax = 1./lamdasmax
   rslopegmax = 1./lamdagmax
   rsloperbmax = rslopermax ** bvtr
   rslopesbmax = rslopesmax ** bvts
   rslopegbmax = rslopegmax ** bvtg
   rsloper2max = rslopermax * rslopermax
   rslopes2max = rslopesmax * rslopesmax
   rslopeg2max = rslopegmax * rslopegmax
   rsloper3max = rsloper2max * rslopermax
   rslopes3max = rslopes2max * rslopesmax
   rslopeg3max = rslopeg2max * rslopegmax
#ifdef DOUBLE
!
   rslopecmax  = 1./lamdacmax
   rslopec2max = rslopecmax * rslopecmax
   rslopec3max = rslopec2max * rslopecmax
#endif
#ifdef HAIL
!
   g6pbgh = rgmma(2.75)
   precg3 = 2.*pi*n0g*.31*g6pbgh*sqrt(sqrt(4.*deng/3./cd))
   bvth2 = 2.5+.5*bvth
   bvth3 = 3.+bvth
   bvth4 = 4.+bvth
   g3pbh = rgmma(bvth3)
   g4pbh = rgmma(bvth4)
   g5pbho2 = rgmma(bvth2)
   pacrh = pi*n0h*avth*g3pbh*.25
   pvth = avth*g4pbh/6.
   prech1 = 2.*pi*n0h*.78
   prech2 = 2.*pi*n0h*.31*avth**.5*g5pbho2
   prech3 = 2.*pi*n0h*.31*g6pbgh*sqrt(sqrt(4.*denh/3./cd))
   pidn0h = pi*denh*n0h
!
   rslopehmax = 1./lamdahmax
   rslopehbmax = rslopehmax ** bvth
   rslopeh2max = rslopehmax * rslopehmax
   rslopeh3max = rslopeh2max * rslopehmax
#endif
#ifdef LB2017_PRAUT
!..LB17
   mrmin = pi/6.*denr*di82**3.
   do n = 1, numax
     ccg(n) = rgmma(bm_r1+n+1.)
     ocg1(n) = 1./rgmma(1.+n)
   enddo
#endif
#ifdef WRF
!
!..set these variables needed for computing radar reflectivity.  these
!.. get used within radar_init to create other variables used in the
!.. radar module.
!
   xam_r = pi*denr/6.
   xbm_r = 3.
   xmu_r = 2.
   xam_s = pi*dens/6.
   xbm_s = 3.
   xmu_s = 0.
   xam_g = pi*deng/6.
   xbm_g = 3.
   xmu_g = 0.
#ifdef HAIL
   xam_h = pi*denh/6.
   xbm_h = 3.
   xmu_h = 0.
#endif
!-------------------------------------------------------------------------------
   call radar_init
!-------------------------------------------------------------------------------
#endif
!
end subroutine udminit
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
real function rgmma(x)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!     rgmma function:  use infinite product form
   real(kind=kind_phys) :: euler
   parameter (euler=0.577215664901532)
   real(kind=kind_phys) :: x, y
   integer :: i
   if(x==1.)then
     rgmma=0.
       else
     rgmma=x*exp(euler*x)
     do i=1,10000
       y=float(i)
       rgmma=rgmma*(1.000+x/y)*exp(-x/y)
     enddo
     rgmma=1./rgmma
   endif
end function rgmma
#ifdef DOUBLE
!-------------------------------------------------------------------------------
   subroutine adjust_number_concent(ktop,kdim,qq,nn,den,piconst,drconst,       &
              qqmin,nnmin,nnmax,di0,dimin,dimax)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer , intent(in ) :: ktop, kdim
   real(kind=kind_phys), dimension(kdim) , intent(inout) :: qq,nn,den
   real(kind=kind_phys), intent(in ) :: piconst,drconst,qqmin,nnmin,nnmax
   real(kind=kind_phys), intent(in ) :: di0,dimin,dimax
! local variables
   real(kind=kind_phys) :: lamdar, diameter
   integer :: k
!-------------------------------------------------------------------------------
   do k = ktop, 1, -1
     if(qq(k)>qqmin) then
       if(nn(k) <= nnmin) then
         lamdar = drconst/di0
         nn(k) = den(k)*qq(k)*exp(log(lamdar)*3.)/piconst
       endif
       lamdar = exp(log(((piconst*nn(k))/(den(k)*qq(k))))*((.33333333)))
       diameter = drconst/lamdar
       if (diameter > dimax) then
         lamdar = drconst/dimax
         nn(k) = den(k)*qq(k)*exp(log(lamdar)*3.)/piconst
       elseif (diameter < dimin) then
         lamdar = drconst/dimin
         nn(k) = den(k)*qq(k)*exp(log(lamdar)*3.)/piconst
       endif
       nn(k) = min(nn(k),nnmax)
     else
       qq(k) = 0.0
       nn(k) = 0.0
     endif
   enddo
!
  return
end subroutine adjust_number_concent
#endif
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine find_cloud_top(im,km,ktop,qq,hvalue,ktopout)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer               , intent(in   ) :: im, km, ktop
   real(kind=kind_phys), dimension(km)   , intent(in   ) :: qq
   real(kind=kind_phys)                  , intent(in   ) :: hvalue
   integer               , intent(inout) :: ktopout
   integer                               :: i,k
!
! do i = 1,im
    ktopout = 0
    find_qrtop : do k = ktop,1, -1
      if(qq(k)>hvalue) then
        ktopout = k
        exit find_qrtop
      else
        cycle find_qrtop
      endif
    enddo find_qrtop
! enddo
  return
end subroutine find_cloud_top
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine udm_mp_effective_radius (t, qc, qi, qs, rho, qmin, t0c,             &
#ifdef DOUBLE
                                nc,                                            &
#else
                             dummy,                                            &
#endif
                                re_qc, re_qi, re_qs, kts, kte, ii, jj)
!-------------------------------------------------------------------------------
!  compute radiation effective radii of cloud water, ice, and snow for
!  single-moment microphysics.
!  these are entirely consistent with microphysics assumptions, not
!  constant or otherwise ad hoc as is internal to most radiation
!  schemes.
!  coded and implemented by soo ya bae, kiaps, january 2015.
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!..sub arguments
   integer, intent(in) :: kts, kte, ii, jj
   real(kind=kind_phys), intent(in) :: qmin
   real(kind=kind_phys), intent(in) :: t0c
   real(kind=kind_phys), dimension( kts:kte ), intent(in)::  t
   real(kind=kind_phys), dimension( kts:kte ), intent(in)::  qc
#ifdef DOUBLE
   real(kind=kind_phys), dimension( kts:kte ), intent(in)::  nc
#else
   real(kind=kind_phys), dimension( kts:kte ), intent(in)::  dummy
#endif
   real(kind=kind_phys), dimension( kts:kte ), intent(in)::  qi
   real(kind=kind_phys), dimension( kts:kte ), intent(in)::  qs
   real(kind=kind_phys), dimension( kts:kte ), intent(in)::  rho
   real(kind=kind_phys), dimension( kts:kte ), intent(inout):: re_qc
   real(kind=kind_phys), dimension( kts:kte ), intent(inout):: re_qi
   real(kind=kind_phys), dimension( kts:kte ), intent(inout):: re_qs
!..local variables
   integer:: i,k
   integer :: nu_c
   real(kind=kind_phys), dimension( kts:kte ):: ni
   real(kind=kind_phys), dimension( kts:kte ):: rqc
   real(kind=kind_phys), dimension( kts:kte ):: rnc
   real(kind=kind_phys), dimension( kts:kte ):: rqi
   real(kind=kind_phys), dimension( kts:kte ):: rni
   real(kind=kind_phys), dimension( kts:kte ):: rqs
   real(kind=kind_phys) :: hvalue
   real(kind=kind_phys) :: lamdac
   real(kind=kind_phys) :: tc, n0sfac, lamdas
   real(kind=kind_phys) :: di      ! diameter of ice in m
#ifdef WYSER_REI
   real(kind=kind_phys) :: bfactor, bfactor2, bfactor3
#endif
   double precision :: lamc
   logical :: has_qc, has_qi, has_qs
!..minimum microphys values
   real(kind=kind_phys), parameter :: r1 = 1.e-12
   real(kind=kind_phys), parameter :: r2 = 1.e-6
!..mass power law relations:  mass = am*d**bm
   real(kind=kind_phys), parameter :: bm_r = 3.0
   real(kind=kind_phys), parameter :: obmr = 1.0/bm_r
   real(kind=kind_phys), parameter :: nc0  = 3.e8
#ifdef WYSER_REI
   real(kind=kind_phys), parameter :: rqi0  = 50.e-3   ! 50 g m-3
#endif
!-------------------------------------------------------------------------------
   has_qc = .false.
   has_qi = .false.
   has_qs = .false.
   do k = kts, kte
     ! for cloud
     rqc(k) = max(r1, qc(k)*rho(k))
#ifdef DOUBLE
     rnc(k) = max(R2, nc(k)*rho(k))
#endif
#ifdef DOUBLE
     if (rqc(k)>R1 .and. rnc(k)>R2) has_qc = .true.
#else
     if (rqc(k)>r1) has_qc = .true.
#endif
     if (rqc(k)>r1) has_qc = .true.
     ! for ice
     rqi(k) = max(r1, qi(k)*rho(k))
     hvalue = (rho(k)*max(qi(k),qmin))
     hvalue = sqrt(sqrt(hvalue*hvalue*hvalue))
     ni(k) = min(max(5.38e7*hvalue,1.e3),1.e6)
     rni(k)= max(r2, ni(k)*rho(k))
     if (rqi(k)>r1 .and. rni(k)>r2) has_qi = .true.
     ! for snow
     rqs(k) = max(r1, qs(k)*rho(k))
     if (rqs(k)>r1) has_qs = .true.
   enddo
   if (has_qc) then
     do k=kts,kte
#ifdef DOUBLE
       if (rqc(k)<=R1 .or. rnc(k)<=R2) CYCLE
       lamc = (pidnc*nc(k)/rqc(k))**obmr
       re_qc(k) = max(recmin,min(0.5*(1./lamc),recmax))
#else
       if (rqc(k)<=r1) cycle
       lamdac   = (pidnc*nc0/rqc(k))**obmr
       re_qc(k) =  max(recmin,min(1.5*(1.0/lamdac),recmax))
#endif
     enddo
   endif
  if (has_qi) then
     do k=kts,kte
       if (rqi(k)<=r1 .or. rni(k)<=r2) cycle
#ifdef WYSER_REI
       hvalue = t0c - t(k)
       bfactor = -2.0 + 1.0e-3*hvalue*sqrt(hvalue)*log10(rqi(k)/rqi0)
       bfactor2 = bfactor*bfactor
       bfactor3 = bfactor2*bfactor
       hvalue = 377.4 + 203.3*bfactor+ 37.91*bfactor2 + 2.3696*bfactor3
       re_qi(k) = max(reimin,min(hvalue*1.e-6,reimax))
#else
       di = 11.9*sqrt(rqi(k)/ni(k))
       re_qi(k) = max(reimin,min(0.75*0.163*di,reimax))
#endif
     enddo
   endif
   if (has_qs) then
     do k=kts,kte
       if (rqs(k)<=r1) cycle
       tc = t(k) - t0c
#ifdef NEW_N0
       n0sfac = max(min(exp(-alpha_qs*tc),n0smax/n0s),n0smin/n0s)
#else
       n0sfac = max(min(exp(-alpha*tc),n0smax/n0s),1.)
#endif
       lamdas = sqrt(sqrt(pidn0s*n0sfac/rqs(k)))
       re_qs(k) = max(resmin,min(0.5*(1./lamdas), resmax))
     enddo
   endif
end subroutine udm_mp_effective_radius
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine udm_mp_reflectivity (qv1d, qr1d, qs1d, qg1d,                        &
#ifdef DOUBLE
                          nr1d,                                                &
#else
                         dummy,                                                &
#endif
#ifdef HAIL
                          qh1d,                                                &
#endif
            t1d, p1d, dbz, kts, kte, ii, jj)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!..sub arguments
   integer, intent(in):: kts, kte, ii, jj
   real(kind=kind_phys), dimension(kts:kte), intent(in):: t1d, p1d
   real(kind=kind_phys), dimension(kts:kte), intent(in):: qv1d, qr1d, qs1d, qg1d
#ifdef HAIL
   real(kind=kind_phys), dimension(kts:kte), intent(in):: qh1d
#endif
#ifdef DOUBLE
   real(kind=kind_phys), dimension(kts:kte), intent(in):: nr1d
#else
   real(kind=kind_phys), dimension(kts:kte), intent(in):: dummy
#endif
   real(kind=kind_phys), dimension(kts:kte), intent(inout):: dbz
!..local variables
   real(kind=kind_phys), dimension(kts:kte):: temp, pres, qv, rho
   real(kind=kind_phys), dimension(kts:kte):: rr, nr, rs, rg, rh
   real(kind=kind_phys):: temp_c
!
   double precision, dimension(kts:kte):: ilamr, ilams, ilamg, ilamh
   double precision, dimension(kts:kte):: n0_r, n0_s, n0_g, n0_h
   double precision:: lamr, lams, lamg, lamh
   logical, dimension(kts:kte):: l_qr, l_qs, l_qg, l_qh
!
   real(kind=kind_phys), dimension(kts:kte):: ze_rain, ze_snow, ze_graupel
   double precision:: fmelt_s, fmelt_g
#ifdef HAIL
   real(kind=kind_phys), dimension(kts:kte):: ze_hail
   double precision:: fmelt_h
#endif
!
   integer:: i, k, k_0, kbot, n
   logical:: melti
!
   double precision:: cback, x, eta, f_d
   real(kind=kind_phys), parameter:: r=287.
!
!+---initialize the local variables 
!+
      temp(:) = 0.  ; pres(:) = 0.  ; qv(:) = 0.   ; rho(:) = 0.
      rr(:) = 0. ;  nr(:) = 0. ;  rs(:) = 0. ;  rg(:) = 0. ;  rh(:) = 0.
      temp_c = 0.
!
      ilamr(:) = 0. ; ilams(:) = 0. ; ilamg(:) = 0. ; ilamh(:) = 0.
      n0_r(:) = 0.;  n0_s(:) = 0. ; n0_g(:) = 0. ; n0_h(:) = 0.
      lamr = 0. ;  lams = 0. ;  lamg = 0. ;  lamh = 0.
      l_qr = .false. ;  l_qs = .false. ;  l_qg = .false. ;  l_qh = .false.
!
      ze_rain(:) = 0. ;  ze_snow(:) = 0. ;  ze_graupel(:) = 0.
      fmelt_s = 0. ; fmelt_g = 0. 
#ifdef HAIL
      ze_hail(:) = 0. ; fmelt_h = 0.
#endif
!
      melti = .false.
!
      cback = 0. ; x = 0. ; x = 0. ; eta = 0. ; f_d = 0.
!
   do k = kts, kte
     dbz(k) = -35.0
   enddo
!   return
!+---+-----------------------------------------------------------------+
!..put column of data into local arrays.
!+---+-----------------------------------------------------------------+
   do k = kts, kte
     temp(k) = t1d(k)
     temp_c = min(-0.001, temp(k)-273.15)
     qv(k) = max(1.e-10, qv1d(k))
     pres(k) = p1d(k)
     rho(k) = 0.622*pres(k)/(r*temp(k)*(qv(k)+0.622))
!
     if (qr1d(k) > 1.e-9) then
       rr(k) = qr1d(k)*rho(k)
#ifdef DOUBLE
       nr(k) = nr1d(k)*rho(k)
       lamr = (xam_r*xcrg(3)*xorg2*nr(k)/rr(k))**xobmr
       N0_r(k) = nr(k)*xorg2*lamr**xcre(2)
#else
       n0_r(k) = n0r
       lamr = (xam_r*xcrg(3)*n0_r(k)/rr(k))**(1./xcre(1))
#endif
       ilamr(k) = 1./lamr
       l_qr(k) = .true.
     else
       rr(k) = 1.e-12
#ifdef DOUBLE
       nr(k) = 1.e-12
#endif
       l_qr(k) = .false.
     endif
!
     if (qs1d(k) > 1.e-9) then
       rs(k) = qs1d(k)*rho(k)
#ifdef NEW_N0
       n0_s(k) = min(n0smax, n0s*exp(-alpha_qs*temp_c))
#else
       n0_s(k) = min(n0smax, n0s*exp(-alpha*temp_c))
#endif
       lams = (xam_s*xcsg(3)*n0_s(k)/rs(k))**(1./xcse(1))
       ilams(k) = 1./lams
       l_qs(k) = .true.
     else
       rs(k) = 1.e-12
       l_qs(k) = .false.
     endif
!
     if (qg1d(k) > 1.e-9) then
       rg(k) = qg1d(k)*rho(k)
#ifdef NEW_N0
       n0_g(k) = min(n0gmax, n0g*exp(-alpha_qg*temp_c))
#else
       n0_g(k) = n0g
#endif
       lamg = (xam_g*xcgg(3)*n0_g(k)/rg(k))**(1./xcge(1))
       ilamg(k) = 1./lamg
       l_qg(k) = .true.
     else
       rg(k) = 1.e-12
       l_qg(k) = .false.
     endif
#ifdef HAIL
     if (qh1d(k) > 1.e-9) then
       rh(k) = qh1d(k)*rho(k)
#ifdef NEW_N0
       n0_h(k) = min(n0hmax, n0h*exp(-alpha_qh*temp_c))
#else
       n0_h(k) = n0h
#endif
       lamh = (xam_h*xchg(3)*n0_h(k)/rh(k))**(1./xche(1))
       ilamh(k) = 1./lamh
       l_qh(k) = .true.
     else
       rh(k) = 1.e-12
       l_qh(k) = .false.
     endif
#endif
   enddo
!
!+---+-----------------------------------------------------------------+
!..locate k-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
   melti = .false.
   k_0 = kts
   do k = kte-1, kts, -1
#ifdef HAIL
      if ( (temp(k)>273.15) .and. L_qr(k)                                      &
                     .and. (L_qs(k+1).or.L_qg(k+1).or.l_qh(k+1)) ) then
#else
      if ( (temp(k)>273.15) .and. L_qr(k)                                      &
                     .and. (L_qs(k+1).or.L_qg(k+1)) ) then
#endif
         k_0 = max(k+1, k_0)
         melti=.true.
         goto 195
      endif
   enddo
 195  continue
!+---+-----------------------------------------------------------------+
!..assume rayleigh approximation at 10 cm wavelength. rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. integrations of m(d)*m(d)*n(d)*dd.
!+---+-----------------------------------------------------------------+
   do k = kts, kte
      ze_rain(k) = 1.e-22
      ze_snow(k) = 1.e-22
      if (l_qr(k)) ze_rain(k) = n0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
      if (l_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/pi)*(6.0/pi)               &
                              * (xam_s/900.0)*(xam_s/900.0)                    &
                              * n0_s(k)*xcsg(4)*ilams(k)**xcse(4)
      ze_graupel(k) = 1.e-22
      if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)            &
     &                           * (xam_g/900.0)*(xam_g/900.0)                 &
     &                           * N0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
#ifdef HAIL
      ze_hail(k) = 1.e-22
      if (L_qh(k)) ze_hail(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)               &
     &                           * (xam_h/900.0)*(xam_h/900.0)                 &
     &                           * N0_h(k)*xchg(4)*ilamh(k)**xche(4)
#endif
  enddo
!+---+-----------------------------------------------------------------+
!..special case of melting ice (snow/graupel) particles.  assume the
!.. ice is surrounded by the liquid water.  fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. uses code from uli blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+
   if (melti .and. k_0>=kts+1) then
     do k = k_0-1, kts, -1
!..reflectivity contributed by melting snow
       if (l_qs(k) .and. l_qs(k_0) ) then
        fmelt_s = max(0.005d0, min(1.0d0-rs(k)/rs(k_0), 0.99d0))
        eta = 0.d0
        lams = 1./ilams(k)
        do n = 1, nrbins
           x = xam_s * xxds(n)**xbm_s
           call rayleigh_soak_wetgraupel (x,dble(xocms),dble(xobms),           &
                 fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar,           &
                 cback, mixingrulestring_s, matrixstring_s,                    &
                 inclusionstring_s, hoststring_s,                              &
                 hostmatrixstring_s, hostinclusionstring_s)
           f_d = n0_s(k)*xxds(n)**xmu_s * dexp(-lams*xxds(n))
           eta = eta + f_d * cback * simpson(n) * xdts(n)
        enddo
        ze_snow(k) = sngl(lamda4 / (pi5 * k_w) * eta)
       endif
!..Reflectivity contributed by melting graupel
          if (L_qg(k) .and. L_qg(k_0) ) then
           fmelt_g = MAX(0.05d0, MIN(1.0d0-rg(k)/rg(k_0), 0.99d0))
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = xam_g * xxDg(n)**xbm_g
              call rayleigh_soak_wetgraupel (x, DBLE(xocmg), DBLE(xobmg),      &
     &              fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar,        &
     &              CBACK, mixingrulestring_g, matrixstring_g,                 &
     &              inclusionstring_g, hoststring_g,                           &
     &              hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)*xxDg(n)**xmu_g * DEXP(-lamg*xxDg(n))
              eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif
#ifdef HAIL
!..Reflectivity contributed by melting hail
          if (L_qh(k) .and. L_qh(k_0) ) then
           fmelt_h = MAX(0.05d0, MIN(1.0d0-rh(k)/rh(k_0), 0.99d0))
           eta = 0.d0
           lamh = 1./ilamh(k)
           do n = 1, nrbins
              x = xam_h * xxDh(n)**xbm_h
              call rayleigh_soak_wetgraupel (x, DBLE(xocmh), DBLE(xobmh),      &
     &              fmelt_h, melt_outside_h, m_w_0, m_i_0, lamda_radar,        &
     &              CBACK, mixingrulestring_h, matrixstring_h,                 &
     &              inclusionstring_h, hoststring_h,                           &
     &              hostmatrixstring_h, hostinclusionstring_h)
              f_d = N0_h(k)*xxDh(n)**xmu_h * DEXP(-lamh*xxDh(n))
              eta = eta + f_d * CBACK * simpson(n) * xdth(n)
           enddo
           ze_hail(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif
#endif
     enddo
   endif
   do k = kte, kts, -1
#ifdef HAIL
     dbz(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k)+ze_hail(k))*1.d18)
#else
     dbz(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
#endif
   enddo
end subroutine udm_mp_reflectivity
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine udm_funct_shape_setup
!-------------------------------------------------------------------------------
!
!  abstract: computes the ratio of fall velocity for mass over number
!    concentations, vm/vn, as a function of variable shape parameter.
!    the shape paramter is a diagnosed with respect to diamater of
!    drops. exact shape dependent ratio is calculated in subprogram
!    fshapex. the current implementation computes a table with a length
!    of 9991 for diameter ranging from 10 micro to 10 mm.
!
!  usage : call udm_funct_shape_setup
!
!  history : first development                        songyou hong    2023-11-27
!
!  subprograms called :
!    (fshapex)  - inlinable function to compute velocity ratio
!
!-------------------------------------------------------------------------------
   implicit none
!
   real   , parameter  ::  xmin = 10.e-6, xmax = 10.e-3
   integer             ::  jx
   real                ::  xinc,x,d
!
   xinc = (xmax-xmin)/(nxshape-1)
   c1xshape = 1.-xmin/xinc
   c2xshape = 1./xinc
!
   do jx = 1,nxshape
     x = xmin+(jx-1)*xinc
     d = x
     tbshape(jx) = fshapex(d)
   enddo
!
   return
   end subroutine udm_funct_shape_setup
!-------------------------------------------------------------------------------
   real function fshapex(d)
!-------------------------------------------------------------------------------
!
!  abstract : exactly compute velocity ratio from diameter.
!    vratio = gamma(1+myu)*gamma(4+bvtr+myu)/gamma(1+bvtr+myu)/gamma(4+myu)
!    where bvtr is the power parameter in terminal velocity, and myu is
!    the dignosed shape parameter,which is a function of drop size.
!    myu = 11.8*[1000D**(1/3)-0.7)**2 + 2, from Milbrandt and
!    Mctaggart-Cowan (JAS, 2010). D is the unit in  meters.
!    this function should be expanded inline in the calling routine.
!
!  usage :   shape=fshapex(d)
!
!  input :
!    d        - real diameter in meters
!
!  output :
!    fshapex    - real velocity ratio, Vm/Vn
!
!-------------------------------------------------------------------------------
!
   implicit none
!
! passing variable
!
   real   , intent(in   )  ::  d
!
! local variable
!
   real                    ::  shape_rain
!
     shape_rain =  min(11.8*(1000.*d-0.7)**2 + 2.,10.)
     fshapex = max(rgmma(4.+shape_rain+bvtr)*rgmma(1.+shape_rain)              &
             /rgmma(1.+shape_rain+bvtr)/rgmma(4.+shape_rain),1.)
!
   return
   end function fshapex
!-------------------------------------------------------------------------------
   real function fshape(d)
!-------------------------------------------------------------------------------
!
!  abstract : compute velocity ratio from the diameter of drops.
!    a linear interpolation is done between values in a lookup table
!    computed in funct_shape_setup.
!    input values outside table range are reset to table extrema.
!
!  usage :   shape=fshape(d)
!
!  input :
!    d        - real diameter in meters
!
!  output :
!    fshape     - real velocity ratio
!
!-------------------------------------------------------------------------------
   implicit none
!
   real   , intent(in   )            ::  d
!
   integer  ::  jx
   real     ::  xj
!
   xj = min(max(c1xshape+c2xshape*d,1.),real(nxshape))
   jx = min(xj,nxshape-1.)
   fshape = tbshape(jx)+(xj-jx)*(tbshape(jx+1)-tbshape(jx))
!
   return
   end function fshape
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine udm_funct_lb2017_setup
!-------------------------------------------------------------------------------
!
!  abstract: computes the production of qc to qr and nc to nr in the
!    lb207 autoconversion physics in terms of lamda_r, slope parameter of rain.
!    the current implementation computes a table with a length of 100 for
!    diameter ranging from 0.1 micro [10**7 m-1] to 100 micro [10**4 m-1].
!
!  usage : udm_call funct_lb2017_setup
!
!  history : first development                        songyou hong    2024-11-27
!
!  subprograms called :
!    (funct_qc_totx)  - inlinable function to compute qc_total
!    (funct_qc_subx)  - inlinable function to compute qc_sub
!    (funct_nc_totx)  - inlinable function to compute nc_total
!    (funct_nc_subx)  - inlinable function to compute nc_sub
!
!-------------------------------------------------------------------------------
   implicit none
!
   real   , parameter  ::  xmin = 1.e+4, xmax = 1.e+7
   integer             ::  jx, nu
   real                ::  xinc,x,d
!
   double precision    ::  funct_aut_qc_totx,funct_aut_qc_subx
   double precision    ::  funct_aut_nc_totx,funct_aut_nc_subx
!
   funct_aut_qc_totx = 0. ; funct_aut_qc_subx = 0.
   funct_aut_nc_totx = 0. ; funct_aut_nc_subx = 0.
   tbaut_qc_tot(:,:) = 0. ; tbaut_qc_sub(:,:) = 0.
   tbaut_nc_tot(:,:) = 0. ; tbaut_nc_sub(:,:) = 0.
!
   xinc = (xmax-xmin)/(nxaut-1)
   c1xaut = 1.-xmin/xinc
   c2xaut = 1./xinc
!
   do nu = 1,numax
   do jx = 1,nxaut
     x = xmin+(jx-1)*xinc
     d = x
     call comp_aut_table(d,nu,funct_aut_qc_totx,funct_aut_qc_subx,             &
                           funct_aut_nc_totx,funct_aut_nc_subx)
     tbaut_qc_tot(jx,nu) = funct_aut_qc_totx
     tbaut_qc_sub(jx,nu) = funct_aut_qc_subx
     tbaut_nc_tot(jx,nu) = funct_aut_nc_totx
     tbaut_nc_sub(jx,nu) = funct_aut_nc_subx
   enddo
   enddo
!
   return
   end subroutine udm_funct_lb2017_setup
!-------------------------------------------------------------------------------
   subroutine comp_aut_table(lamc,nu_c,qc_tot,qc_sub,nc_tot,nc_sub)
!-------------------------------------------------------------------------------
!
!  abstract : exactly compute conversoin rate for qc to qr, and nc and nr.
!    this function should be expanded inline in the calling routine.
!    see equations 15~19 of Lee and Baik (2017)
!
!  usage :   call comp_aut_table
!
!  input :
!    lamc        - slope parameter for rain
!
!  output :
!    qc_tot, qc_sub, nc_tot, nc_sub
!
!-------------------------------------------------------------------------------
!
   implicit none
!
! passing variable
!
   real   , intent(in   )    ::  lamc
   double precision, intent(inout  )   ::  qc_tot,qc_sub,nc_tot,nc_sub
   double precision, parameter         :: &
                  r_star    = 50.0d-6,    &  !< from 40 micro-m in lb2017
                  a_aut     = 0.21421d0,  &
                  b_aut     = -1.11347d4
   integer,parameter           :: n_gamma = 42
   double precision, dimension(0:n_gamma) :: gamma1, gamma2,       &
                          exp1, exp2, exp3, err1, err2
   double precision :: term, temp11, temp12, temp21, temp22
   double precision, dimension(10), parameter :: ak1 =             &
  (/a_aut, 1+a_aut, 1-2*a_aut, -2+a_aut, -1+2*a_aut,               &
         2-a_aut, -1-2*a_aut, -2+a_aut, 1+a_aut, 1.0d0       /)
   double precision, dimension( 7), parameter :: ak2 =             &
  (/a_aut, 1+a_aut, 1-2*a_aut, -2-2*a_aut, -2+a_aut, 1+a_aut,1.0d0/)
   integer :: n, nn, nu_c
!
   gamma1(:) = 0. ; gamma2(:) = 0. 
   exp1(:) = 0.   ; exp2(:) = 0. ; exp3(:) = 0. ; err1(:) = 0.; err2(:) = 0.
   term = 0.      ; temp11 = 0.0 ; temp12 = 0.  ; temp21 = 0. ; temp22 = 0.
!
         gamma1(0) = 1.0d0/(lamc)
         gamma2(0) = 1.0d0/(2.0d0*lamc)
         exp1(0) = 1.0d0
         err1(0) = 1.0d0
         exp2(0) = 1.0d0
         err2(0) = 1.0d0
!
         do n = 1, n_gamma
           gamma1(n) = gamma1(n-1)*(dble(n)/(lamc))
           gamma2(n) = gamma2(n-1)*(dble(n)/(2.0d0*lamc))
           exp1(n) = exp1(n-1) *(( lamc*r_star)/dble(n))
           exp2(n) = exp2(n-1) *((2.0d0*lamc*r_star)/dble(n))
         end do
!
         exp3 = 1.0d0/(lamc*gamma1)
!
         do n = 1, n_gamma
           err1(n) = 0.0d0
           err2(n) = 0.0d0
           do nn = 0, n
             err1(n) = err1(n) + exp1(nn)
             err2(n) = err2(n) + exp2(nn)
           end do
           err1(n) = 1.0d0 - exp(- lamc*r_star)*err1(n)
           err2(n) = 1.0d0 - exp(-2.0d0*lamc*r_star)*err2(n)
         end do
!-----------------------------------------
! for mass concentraion
!-----------------------------------------
         qc_tot = 0.0d0
         qc_sub = 0.0d0
!
         do n = 1, 10
           temp11 = 0.0d0
           temp12 = 0.0d0
           do nn = 0, nu_c+n
             term = exp3(nn)*gamma2(nu_c+10-n+nn)
             temp11 = temp11 + term
             temp12 = temp12 + term*err2(nu_c+10-n+nn)
           end do
           temp21 = 0.0d0
           temp22 = 0.0d0
           do nn = 0, nu_c+n
             term = exp3(nn)*gamma2(nu_c+11-n+nn)
             temp21 = temp21 + term
             temp22 = temp22 + term*err2(nu_c+11-n+nn)
           end do
!
           qc_tot = qc_tot + ak1(n)*gamma1(nu_c+n)*(gamma1(nu_c+10-n)          &
                           - temp11 + b_aut*(gamma1(nu_c+11-n) -temp21))
           qc_sub = qc_sub + ak1(n)*gamma1(nu_c+n)*(gamma1(nu_c+10-n)          &
                           *err1(nu_c+10-n)-temp12 +b_aut*(gamma1(nu_c+11-n)   &
                           *err1(nu_c+11-n)-temp22))
         end do
!-----------------------------------------
! for number concentraion
!-----------------------------------------
         nc_tot = 0.0d0
         nc_sub = 0.0d0
!
         do n = 1, 7
           temp11 = 0.0d0
           temp12 = 0.0d0
           do nn=0, nu_c+n
             term = exp3(nn)*gamma2(nu_c+7-n+nn)
             temp11 = temp11 + term
             temp12 = temp12 + term*err2(nu_c+7-n+nn)
           end do
           temp21 = 0.0d0
           temp22 = 0.0d0
           do nn = 0, nu_c+n
             term = exp3(nn)*gamma2(nu_c+8-n+nn)
             temp21 = temp21 + term
             temp22 = temp22 + term*err2(nu_c+8-n+nn)
           end do
!
           nc_tot = nc_tot + ak2(n)*gamma1(nu_c+n)*(gamma1(nu_c+7-n)           &
                           - temp11 + b_aut*(gamma1(nu_c+8-n) - temp21))
           nc_sub = nc_sub + ak2(n)*gamma1(nu_c+n)*(gamma1(nu_c+7-n)           &
                           *err1(nu_c+7-n)-temp12 +b_aut*(gamma1(nu_c+8-n)     &
                           *err1(nu_c+8-n)-temp22))
         end do
!
   return
   end subroutine comp_aut_table
!-------------------------------------------------------------------------------
   double precision function funct_aut_qc_tot(d,nu)
!-------------------------------------------------------------------------------
!
!  abstract :
!    a linear interpolation is done between values in a lookup table
!    computed in funct_lb2017_setup.
!    input values outside table range are reset to table extrema.
!
!  usage :   funct_aut_qc_tot(d)
!
!  input :
!    d        - slope parameter
!
!  output :
!    qc_tot   - real total conversion rate of qc
!
!-------------------------------------------------------------------------------
   implicit none
!
   double precision      , intent(in   )            ::  d
   integer   , intent(in   )            :: nu
!
   integer  ::  jx
   real     ::  xj
!
   xj = min(max(c1xaut+c2xaut*d,1.d0),real(nxaut))
   jx = min(xj,nxaut-1.d0)
   funct_aut_qc_tot = tbaut_qc_tot(jx,nu)+(xj-jx)*(tbaut_qc_tot(jx+1,nu)       &
                     -tbaut_qc_tot(jx,nu))
!
   return
   end function funct_aut_qc_tot
!-------------------------------------------------------------------------------
   double precision function funct_aut_qc_sub(d,nu)
!-------------------------------------------------------------------------------
!
!  abstract :
!    a linear interpolation is done between values in a lookup table
!    computed in funct_lb2017_setup.
!    input values outside table range are reset to table extrema.
!
!  usage :   funct_aut_qc_sub(d)
!
!  input :
!    d        - slope parameter
!
!  output :
!    qc_tot   - real sub conversion rate of qc
!
!-------------------------------------------------------------------------------
   implicit none
!
   double precision      , intent(in   )            ::  d
   integer   , intent(in   )            :: nu
!
   integer  ::  jx
   real     ::  xj
!
   xj = min(max(c1xaut+c2xaut*d,1.d0),real(nxaut))
   jx = min(xj,nxaut-1.d0)
   funct_aut_qc_sub = tbaut_qc_sub(jx,nu)+(xj-jx)*(tbaut_qc_sub(jx+1,nu)       &
                     -tbaut_qc_sub(jx,nu))
!
   return
   end function funct_aut_qc_sub
!-------------------------------------------------------------------------------
   double precision function funct_aut_nc_tot(d,nu)
!-------------------------------------------------------------------------------
!
!  abstract :
!    a linear interpolation is done between values in a lookup table
!    computed in funct_lb2017_setup.
!    input values outside table range are reset to table extrema.
!
!  usage :   funct_aut_nc_tot(d)
!
!  input :
!    d        - slope parameter
!
!  output :
!    qc_tot   - real sub conversion rate of nc
!
!-------------------------------------------------------------------------------
   implicit none
!
   double precision      , intent(in   )            ::  d
   integer   , intent(in   )            :: nu
!
   integer  ::  jx
   real     ::  xj
!
   xj = min(max(c1xaut+c2xaut*d,1.d0),real(nxaut))
   jx = min(xj,nxaut-1.d0)
   funct_aut_nc_tot = tbaut_nc_tot(jx,nu)+(xj-jx)*(tbaut_nc_tot(jx+1,nu)       &
                     -tbaut_nc_tot(jx,nu))
!
   return
   end function funct_aut_nc_tot
!-------------------------------------------------------------------------------
   double precision function funct_aut_nc_sub(d,nu)
!-------------------------------------------------------------------------------
!
!  abstract :
!    a linear interpolation is done between values in a lookup table
!    computed in funct_lb2017_setup.
!    input values outside table range are reset to table extrema.
!
!  usage :   funct_aut_nc_sub(d)
!
!  input :
!    d        - slope parameter
!
!  output :
!    nc_sub   - real sub conversion rate of nc
!
!-------------------------------------------------------------------------------
   implicit none
!
   double precision      , intent(in   )            ::  d
   integer   , intent(in   )            :: nu
!
   integer  ::  jx
   real     ::  xj
!
   xj = min(max(c1xaut+c2xaut*d,1.d0),real(nxaut))
   jx = min(xj,nxaut-1.d0)
   funct_aut_nc_sub = tbaut_nc_sub(jx,nu)+(xj-jx)*(tbaut_nc_sub(jx+1,nu)       &
                     -tbaut_nc_sub(jx,nu))
!
   return
   end function funct_aut_nc_sub
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine udm_funct_svp_setup
!-------------------------------------------------------------------------------
!
!  abstract: compute saturation vapor pressure table for the table lookup
!    function fsvp.exact saturation vapor pressures are calculated in
!    subprogram fsvpx.
!    the current implementation computes a table with a length
!    of 7501 for temperatures ranging from 180. to 330. kelvin.
!
!  history log:
!    1991-05-07 iredell         made into inlinable function
!    1996-02-19 song-you hong   ice effect & increased range and accuracy
!    2009-10-01 jung-eun kim    f90 format with standard physics modules
!    2010-07-01 myung-seo koo   dimension allocatable with namelist input
!    2023-11-27 song-you hong   clean-up and implemented on ufs model
!
!  usage : call funct_svp_setup
!
!  subprograms called :
!    (fsvpx)  - inlinable function to compute saturation vapor pressure
!
!-------------------------------------------------------------------------------
   implicit none
!
   real   , parameter  ::  xmin = 180.0,                                       &
                           xmax = 330.0
   integer             ::  jx
   real                ::  xinc,x,t
!
   xinc = (xmax-xmin)/(nxsvp-1)
   c1xsvp = 1.-xmin/xinc
   c2xsvp = 1./xinc
   c1xsvpw = c1xsvp
   c2xsvpw = c2xsvp
!
   do jx = 1,nxsvp
     x = xmin+(jx-1)*xinc
     t = x
     tbsvp(jx) = fsvpx(t)
     tbsvpw(jx) = fsvpxw(t)
   enddo
!
   return
   end subroutine udm_funct_svp_setup
!-------------------------------------------------------------------------------
   real function fsvpx(t)
!-------------------------------------------------------------------------------
!
!  abstract : exactly compute saturation vapor pressure from temperature.
!    the clausius-clapeyron equation is integrated from the triple point
!    to get the formula
!       svp=psatk*(tr**xa)*exp(xb*(1.-tr))
!    where tr is ttp/t and other values are physical constants
!    this function should be expanded inline in the calling routine.
!
!  usage :   svp=fsvpx(t)
!
!  input :
!    t        - real temperature in kelvin
!
!  output :
!    fsvpx    - real saturation vapor pressure in pascale (pa)
!
!-------------------------------------------------------------------------------
   implicit none
!
! passing variable
!
   real   , intent(in   )  ::  t
!
! local variable
!
   real                    ::  tr
!
   tr = ttp_/t
   if (t >= ttp_) then
     fsvpx = psatk*(tr**xa)*exp(xb*(1.-tr))
   else
     fsvpx = psatk*(tr**xai)*exp(xbi*(1.-tr))
   endif
!
   return
   end function fsvpx
!-------------------------------------------------------------------------------
   real function fsvpxw(t)
!-------------------------------------------------------------------------------
!
!  abstract : same to fxvps but water only
!
!-------------------------------------------------------------------------------
   implicit none
!
   real   , intent(in   )  ::  t
   real                    ::  tr
!
   tr = ttp_/t
   fsvpxw = psatk*(tr**xa)*exp(xb*(1.-tr))
!
   return
   end function fsvpxw
!-------------------------------------------------------------------------------
   real function fsvp(t,p)
!-------------------------------------------------------------------------------
!
!  abstract : compute saturation vapor pressure from the temperature.
!    a linear interpolation is done between values in a lookup table
!    computed in funct_svp_setup.
!
!  usage :   svp=fsvp(t)
!
!  input :
!    t        - real temperature in kelvin
!    p        - real pressure in pa (optional)
!
!  output :
!    fsvp     - real saturation vapor pressure in pascals (pa)
!
!-------------------------------------------------------------------------------
   implicit none
!
   real   , intent(in   )            ::  t
   real   , intent(in   ), optional  ::  p
!
   integer  ::  jx
   real     ::  xj
!
   xj = min(max(c1xsvp+c2xsvp*t,1.),real(nxsvp))
   jx = min(xj,nxsvp-1.)
   fsvp = tbsvp(jx)+(xj-jx)*(tbsvp(jx+1)-tbsvp(jx))
!
   if (present(p)) fsvp = min(fsvp,0.99*p)
!
   return
   end function fsvp
!-------------------------------------------------------------------------------
   real function fsvp_water(t,p)
!-------------------------------------------------------------------------------
!
!  abstract : same to fsvp but water only
!
!-------------------------------------------------------------------------------
   implicit none
!
   real   , intent(in   )            ::  t
   real   , intent(in   ), optional  ::  p
!
   integer  ::  jx1
   real     ::  xj1
!
   xj1 = min(max(c1xsvpw+c2xsvpw*t,1.),real(nxsvpw))
   jx1 = min(xj1,nxsvpw-1.)
   fsvp_water = tbsvpw(jx1)+(xj1-jx1)*(tbsvpw(jx1+1)-tbsvpw(jx1))
!
   if (present(p)) fsvp_water = min(fsvp_water,0.99*p)
!
   return
   end function fsvp_water
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
      real function iceDeMott(tc, den, nifa) 
      implicit none
      REAL, INTENT(IN):: tc, nifa, den
      REAL:: xni, nifa_cc, hvalue
!..Local vars
      REAL, PARAMETER:: rho_not0 = 101325./(287.05*273.15)
!     tc(k) = t(k,i) - t0c
!     if(t(k,i)<t0c) ifice(k) = .true.
      xni = 0.0
      nifa_cc = 0.0
      hvalue = 0.0
      hvalue = min(tc, 0.)
      nifa_cc = MAX(0.5, nifa*rho_not0/den*1.E-6) ! n(aer,0.5) (unit:#/m3) to (unit:#/scm3)
      xni = (5.94e-5*(-hvalue)**3.33)                                 & ! [DeMott, 2010]
                 * (nifa_cc**((-0.0264*(hvalue))+0.0033)) ! to n(IN,Tk) (unit:std/L at Tk)
      xni = xni*den/rho_not0*1000. ! (unit:std/L at Tk) to (unit:#/m3)
      iceDeMott = MAX(inpmin, xni)
      end FUNCTION iceDeMott 
#ifdef WET_BULB
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
      real function stull_wet_bulb(t, rh)
        real, intent(in) :: t, rh
        real :: tw_val
        tw_val = t * atan(0.151977 * sqrt(rh + 8.313659)) +                    &
                 atan(t + rh) - atan(rh - 1.676331) +                          &
                 0.00391838 * (rh**1.5) * atan(0.023101 * rh) - 4.686035
        stull_wet_bulb = tw_val
      end function stull_wet_bulb
#endif
!-------------------------------------------------------------------------------
!
end module MODULE_NAME_
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
