module constants
    use kinds


    implicit none
    save

!-----------------------------------------------------------------------
!
!     physical constants (all in cgs units except for those in MKS)
!
!-----------------------------------------------------------------------

    real (kind=dbl_kind), parameter ::     &
   &  grav      = 980.6_dbl_kind,          & ! gravit. accel. (cm/s**2)
   &  omega     = 7.292123625e-5_dbl_kind, & ! angular vel. of Earth 1/s
   &  radius    = 6370.0e5_dbl_kind,       & ! radius of Earth (cm)
   &  rho_sw    = 1.026_dbl_kind,          & ! density of sea water (g/cm^3)
   &  rho_fw    = 1.0_dbl_kind,            & ! density of pure water(g/cm^3)
   &  cp_water  = 4.180e7_dbl_kind,        & ! specific heat of water
   &  cp_sw     = 3.996e7_dbl_kind,        & ! specific heat salt water
   &  sound     = 1.5e5_dbl_kind,          & ! speed of sound (cm/s)
   &  rho0      = 1.03_dbl_kind,           & ! average water density
   &  vonkar    = 0.4_dbl_kind,            & ! von Karman constant
   &  cp_air    = 1005.0_dbl_kind,         & ! heat capacity of air (J/kg/K)
   &  rho_air   = 1.2_dbl_kind,            & ! ambient air density (kg/m^3)
   &  emissivity         = 1.0_dbl_kind,   &
   &  stefan_boltzmann   = 567.0e-10_dbl_kind, & !  W/m^2/K^4
   &  latent_heat_vapor  = 2.5e6_dbl_kind, & ! latent heat of vapor. (erg/g)
   &  latent_heat_fusion = 3.34e9_dbl_kind,& ! latent heat of fusion (erg/g)
   &  ocn_ref_salinity   = 34.7_dbl_kind     ! (psu)

!-----------------------------------------------------------------------
!
!     numbers
!
!-----------------------------------------------------------------------

    character (char_len) :: char_blank          ! empty character string


    real (kind=dbl_kind), parameter ::  &
   &  c0   = 0.0_dbl_kind, &
   &  c1   = 1.0_dbl_kind, &
   &  c1p5 = 1.5_dbl_kind, &
   &  c2   = 2.0_dbl_kind, &
   &  c3   = 3.0_dbl_kind, &
   &  c4   = 4.0_dbl_kind, &
   &  c5   = 5.0_dbl_kind, &
   &  c8   = 8.0_dbl_kind, &
   &  c10  = 10.0_dbl_kind, &
   &  c16  = 16.0_dbl_kind, &
   &  c1000= 1000.0_dbl_kind, &
   &  p33  = c1/c3,           &
   &  p5   = 0.5_dbl_kind, &
   &  p25  = 0.25_dbl_kind, &
   &  p125 = 0.125_dbl_kind, &
   &  p001 = 0.001_dbl_kind, &
   &  eps  = 1.0e-10_dbl_kind, &
   &  eps2  = 1.0e-20_dbl_kind

    real (kind=dbl_kind) :: pi, pih, pi2         ! pi, pi/2 and 2pi

!-----------------------------------------------------------------------
!
!     conversion factors
!
!-----------------------------------------------------------------------

    real (kind=dbl_kind), parameter ::     &
   &  T0_Kelvin     = 273.16_dbl_kind      & ! zero point for Celcius
   &, mpercm        = .01_dbl_kind         & ! meters per cm
   &, cmperm        = 100._dbl_kind        & ! cm per meter
   &, salt_to_ppt   = 1000._dbl_kind       & ! salt (g/g) to ppt
   &, ppt_to_salt   = 1.e-3_dbl_kind       & ! salt ppt to g/g
   &, mass_to_Sv    = 1.0e-12_dbl_kind     & ! mass flux to Sverdrups
   &, heat_to_PW    = 4.186e-15_dbl_kind   & ! heat flux to Petawatts
   &, salt_to_Svppt = 1.0e-9_dbl_kind      & ! salt flux to Sv*ppt
   &, salt_to_mmday = 3.1536e+5_dbl_kind     ! salt to water (mm/day)

    real (kind=dbl_kind) :: radian           ! degree-radian conversion

!-----------------------------------------------------------------------
!
!   Flux conversions:
!
!    velocity       flux = windstress / rho
!    temperature    flux = surface heat flux / rho / cp
!    salinity       flux = surface freshwater flux
!                          * ocean reference salinity
!                          / rho(freshwater)
!
!   Units:
!
!    windstress           (N/m^2)
!    velocity       flux  (cm^2/s^2)
!
!    heat           flux  (W/m^2)
!    temperature    flux  (C*cm/s)
!
!    freshwater     flux  (Kg/m^2/s)
!    salinity       flux  (msu*cm/s)
!
!
!      convert windstress to velocity flux:
!      ---------------------------------------
!      windstress in (N/m^2) = (kg/s^2/m) = 10(g/s^2/cm) = 10(dyn/cm^2)
!      assume here that density of seawater rho = 1 (g/cm^3)
!
!      vel_flux   = windstress / rho
!
!      vel_flux (cm^2/s^2)
!                 = windstress (N/m^2)
!                 * 10 (g/s^2/cm)/(N/m^2)
!                 / [1 (g/cm^3)]
!
!                 = windstress (N/m^2)
!                 * momentum_factor (cm^2/s^2)/N/m^2)
!
!      ==>  momentum_factor = 10
!
!
!      convert heat and solar flux to temperature flux:
!      -----------------------------------------------
!      heat_flux in (W/m^2) = (J/s/m^2) = 1000(g/s^3)
!      density of seawater rho_sw in (g/cm^3)
!      specific heat of seawater cp_sw in (erg/g/C) = (cm^2/s^2/C)
!
!      temp_flux = heat_flux / (rho_sw*cp_sw)
!
!      temp_flux (C*cm/s)
!                 = heat_flux (W/m^2)
!                 * 1000 (g/s^3)/(W/m^2)
!                 / [(rho_sw*cp_sw) (g/cm/s^2/C)]
!
!                 = heat_flux (W/m^2)
!                 * hflux_factor (C*cm/s)/(W/m^2)
!
!      ==>  hflux_factor = 1000/(rho_sw*cp_sw)
!
!
!      convert net fresh water flux to salt flux (in model units):
!      ----------------------------------------------------------
!      ocean reference salinity in (o/oo=psu)
!      density of freshwater rho_fw = 1.0 (g/cm^3)
!      h2o_flux in (kg/m^2/s) = 0.1 (g/cm^2/s)
!
!      salt_flux  = - h2o_flux * ocn_ref_salinity / rho_fw
!
!      salt_flux (msu*cm/s)
!                 = - h2o_flux (kg/m^2/s)
!                 * ocn_ref_salinity (psu)
!                 * 1.e-3 (msu/psu)
!                 * 0.1 (g/cm^2/s)/(kg/m^2/s)
!                 / 1.0 (g/cm^3)
!
!                 = - h2o_flux (kg/m^2/s)
!                 * ocn_ref_salinity (psu)
!                 * fwflux_factor (cm/s)(msu/psu)/(kg/m^2/s)
!
!      ==>  fwflux_factor = 1.e-4
!
!      salt_flux(msu*cm/s) = h2oflux(kg/m^2/s) * salinity_factor
!
!      ==> salinity_factor = - ocn_ref_salinity(psu) * fwflux_factor
!
!      convert salt flux to salt flux (in model units):
!      ----------------------------------------------------------
!      density of freshwater rho_fw = 1.0 (g/cm^3)
!      salt_flux_kg in (kg/m^2/s) = 0.1 (g/cm^2/s)
!
!      salt_flux  = - h2o_flux * ocn_ref_salinity / rho_fw
!
!      salt_flux (msu*cm/s)
!                 = salt_flux_kg (kg/m^2/s)
!                 * 0.1 (g/cm^2/s)/(kg/m^2/s)
!                 / 1.0 (g/cm^3)
!
!                 = salt_flux_kg (kg/m^2/s)
!                 * sflux_factor (msu*cm/s)/(kg/m^2/s)
!
!      ==>  slux_factor = 0.1
!
!-----------------------------------------------------------------------

    real (kind=dbl_kind), parameter ::                           &
   &  momentum_factor  = 10.0_dbl_kind                           &
   &, hflux_factor = 1000.0_dbl_kind / (rho_sw * cp_sw)          &
   &, fwflux_factor =  1.e-4_dbl_kind                            &
   &, salinity_factor =  - ocn_ref_salinity * fwflux_factor      &
   &, sflux_factor = 0.1_dbl_kind

!***********************************************************************

    contains

!***********************************************************************

    subroutine init_constants

    implicit none

!-----------------------------------------------------------------------
!
!     This subroutine initializes constants that are best defined
!     at run time (e.g. pi).
!
!-----------------------------------------------------------------------

    integer (kind=int_kind) :: n, ierr

!-----------------------------------------------------------------------

    if ( ierr .NE. 0 ) then
       print *,'ERROR in int_constants: could not allocate unity arrays'
       stop 999
    endif
    
    pi  = c4*atan(c1)
    pi2 = c2*pi
    pih = p5*pi

    radian = 180.0_dbl_kind/pi

    do n=1,char_len
      char_blank(n:n) = ' '
    end do

!-----------------------------------------------------------------------

    end subroutine init_constants

!***********************************************************************

    end module constants

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
