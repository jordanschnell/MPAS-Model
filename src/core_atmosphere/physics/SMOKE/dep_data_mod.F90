!>\file  dep_data_mod.F90
!! This file contains data for the dry deposition modules.
module dep_data_mod

  use mpas_kind_types
  use mpas_smoke_init

  implicit none

  public :: aero_dry_dep_init

  real(RKIND), parameter :: max_dep_vel = 0.005                   ! m/s (may need to set per species)
  real(RKIND), parameter :: dep_ref_hgt = 2.0                     ! Meters 
  real(RKIND), parameter :: pi = 3.1415926536
!  3*PI
  REAL(RKIND), PARAMETER :: threepi=3.0*pi
  real(RKIND), parameter :: gravity =  9.81
! mean gravitational acceleration [ m/sec**2 ]
  REAL(RKIND), PARAMETER :: grav=9.80622
  real(RKIND), parameter :: boltzmann = 1.3807e-16
! universal gas constant [ J/mol-K ]
  REAL(RKIND), PARAMETER :: rgasuniv=8.314510
! Avogadro's Constant [ 1/mol ]
  REAL, PARAMETER :: avo=6.0221367E23
  ! Boltzmann's Constant [ J / K ]i\
  REAL(RKIND), PARAMETER :: boltz=rgasuniv/avo
  real(RKIND), parameter :: Cb = 2., Cim = 0.4, alpha = 0.8, Cin = 2.5, vv = 0.8
  real(RKIND), parameter :: A_for = 0.1 ! forest
  real(RKIND), parameter :: A_grs = 0.2 ! grass
  real(RKIND), parameter :: A_wat = 100. ! water
  real(RKIND), parameter :: eps0_for = 0.8*0.01 ! forest
  real(RKIND), parameter :: eps0_grs = 0.4*0.01 ! grass
  real(RKIND), parameter :: eps0_wat = 0.6*0.01 ! water
  REAL(RKIND), PARAMETER :: one3=1.0/3.0
  REAL(RKIND), PARAMETER :: two3=2.0/3.0
!  SQRT( 2 )
  REAL(RKIND), PARAMETER :: sqrt2=1.4142135623731
!  SQRT( PI )
  REAL(RKIND), PARAMETER :: sqrtpi=1.7724539
  REAL(RKIND) :: karman = 0.4                             ! von Karman constant
  REAL(RKIND), PARAMETER :: conmin= 1.E-16
  REAL(RKIND), PARAMETER :: pirs=3.14159265358979324
  REAL(RKIND), PARAMETER :: f6dpi=6.0/pirs
  REAL(RKIND), PARAMETER :: f6dpim9=1.0E-9*f6dpi
!  starting standard surface temperature [ K ]
  REAL(RKIND), PARAMETER :: sigma1 = 1.8
  REAL(RKIND), PARAMETER :: mean_diameter1 = 4.e-8
  REAL(RKIND), PARAMETER :: fact_wfa = 1.e-9*6.0/pirs*exp(4.5*log(sigma1)**2)/mean_diameter1**3
! lowest particle diameter ( m )   
  REAL(RKIND), PARAMETER :: dgmin=1.0E-09
! lowest particle density ( Kg/m**3 )
  REAL(RKIND), PARAMETER :: densmin=1.0E03

! Arrays to hold density and diameters of aerosols
  real(RKIND), dimension(1000), save :: aero_dens, aero_diam 

  contains

   subroutine aero_dry_dep_init

      implicit none

      aero_dens(:) = 1.E3_RKIND
      aero_diam(:) = 1.E-6_RKIND

    ! Aerosol densities (kg/m3)
      if (p_smoke_fine>0)     aero_dens(p_smoke_fine)       = 1.4E3_RKIND
      if (p_smoke_coarse>0)   aero_dens(p_smoke_coarse)     = 1.4E3_RKIND
      if (p_dust_fine>0)      aero_dens(p_dust_fine)        = 2.6E3_RKIND 
      if (p_unspc_fine>0)     aero_dens(p_unspc_fine)       = 2.6E3_RKIND
      if (p_dust_coarse>0)    aero_dens(p_dust_coarse)      = 2.6E3_RKIND
      if (p_unspc_coarse>0)   aero_dens(p_unspc_coarse)     = 2.6E3_RKIND
      if (p_ssalt_fine>0)     aero_dens(p_ssalt_fine)       = 2.2E3_RKIND
      if (p_ssalt_coarse>0)   aero_dens(p_ssalt_coarse)     = 2.2E3_RKIND
!
      if (p_polp_all>0)       aero_dens(p_polp_all)         = 1.2E3_RKIND
      if (p_polp_tree>0)      aero_dens(p_polp_tree)        = 1.2E3_RKIND
      if (p_polp_grass>0)     aero_dens(p_polp_grass)       = 1.0E3_RKIND
      if (p_polp_weed>0)      aero_dens(p_polp_weed)        = 1.3E3_RKIND
!
      if (p_pols_all>0)       aero_dens(p_pols_all)         = 1.425E3_RKIND
      if (p_pols_tree>0)      aero_dens(p_pols_tree)        = 1.425E3_RKIND
      if (p_pols_grass>0)     aero_dens(p_pols_grass)       = 1.425E3_RKIND
      if (p_pols_weed>0)      aero_dens(p_pols_weed)        = 1.425E3_RKIND
    ! Aerosol diameters (m)
      if (p_smoke_fine>0)    aero_diam(p_smoke_fine)        = 4E-8_RKIND
      if (p_smoke_coarse>0)  aero_diam(p_smoke_coarse)      = 10E-6_RKIND
      if (p_dust_fine>0)     aero_diam(p_dust_fine)         = 1E-6_RKIND
      if (p_unspc_fine>0)    aero_diam(p_unspc_fine)        = 1E-7_RKIND
      if (p_dust_coarse>0)   aero_diam(p_dust_coarse)       = 4.5E-6_RKIND
      if (p_unspc_coarse>0)  aero_diam(p_unspc_coarse)      = 4.5E-6_RKIND
      if (p_ssalt_fine>0)    aero_diam(p_ssalt_fine)        = 6.32E-7_RKIND
      if (p_ssalt_coarse>0)  aero_diam(p_ssalt_coarse)      = 5.632E-6_RKIND
!
      if (p_polp_all>0)      aero_diam(p_polp_all)          = 30E-6_RKIND
      if (p_polp_tree>0)     aero_diam(p_polp_tree)         = 35E-6_RKIND
      if (p_polp_grass>0)    aero_diam(p_polp_grass)        = 35E-6_RKIND
      if (p_polp_weed>0)     aero_diam(p_polp_weed)         = 20E-6_RKIND
!
      if (p_pols_all>0)      aero_diam(p_pols_all)          = 1.5E-7_RKIND
      if (p_pols_tree>0)     aero_diam(p_pols_tree)         = 1.5E-7_RKIND
      if (p_pols_grass>0)    aero_diam(p_pols_grass)        = 1.5E-7_RKIND
      if (p_pols_weed>0)     aero_diam(p_pols_weed)         = 1.5E-7_RKIND

   end subroutine aero_dry_dep_init

end module dep_data_mod
