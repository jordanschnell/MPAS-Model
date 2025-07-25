!>\file  rad_data_mod.F90
!! This file contains data for the radiation modules.
module rad_data_mod

  use mpas_kind_types
  use mpas_smoke_init

  implicit none

  public :: aero_rad_init

! Arrays to hold scattering and absorbing efficiencies of aerosols
  real(RKIND), dimension(1000), save :: sc_eff, ab_eff

  contains

   subroutine aero_rad_init
      implicit none

      sc_eff(:) = 0._RKIND
      ab_eff(:) = 0._RKIND
! Mass scattering efficencies from Hand and Malm, (2007)
! Maybe use IMPROVE instead?
      if (p_smoke_fine>0)     sc_eff(p_smoke_fine)       = 3.9_RKIND
      if (p_smoke_coarse>0)   sc_eff(p_smoke_coarse)     = 2.6_RKIND
      if (p_dust_fine>0)      sc_eff(p_dust_fine)        = 3.3_RKIND
      if (p_unspc_fine>0)     sc_eff(p_unspc_fine)       = 3.3_RKIND
      if (p_dust_coarse>0)    sc_eff(p_dust_coarse)      = 0.7_RKIND
      if (p_unspc_coarse>0)   sc_eff(p_unspc_coarse)     = 0.7_RKIND
      if (p_ssalt_fine>0)     sc_eff(p_ssalt_fine)       = 4.5_RKIND
      if (p_ssalt_coarse>0)   sc_eff(p_ssalt_coarse)     = 1.0_RKIND
!
      if (p_polp_all>0)       sc_eff(p_polp_all)         = 2.6_RKIND
      if (p_polp_tree>0)      sc_eff(p_polp_tree)        = 2.6_RKIND
      if (p_polp_grass>0)     sc_eff(p_polp_grass)       = 2.6_RKIND
      if (p_polp_weed>0)      sc_eff(p_polp_weed)        = 2.6_RKIND
!
      if (p_pols_all>0)       sc_eff(p_pols_all)         = 3.9_RKIND 
      if (p_pols_tree>0)      sc_eff(p_pols_tree)        = 3.9_RKIND
      if (p_pols_grass>0)     sc_eff(p_pols_grass)       = 3.9_RKIND
      if (p_pols_weed>0)      sc_eff(p_pols_weed)        = 3.9_RKIND
! Mass Absorbing efficiencies from 
      if (p_smoke_fine>0)    ab_eff(p_smoke_fine)        = 0.5_RKIND 
      if (p_smoke_coarse>0)  ab_eff(p_smoke_coarse)      = 0.0_RKIND
      if (p_dust_fine>0)     ab_eff(p_dust_fine)         = 0.0_RKIND
      if (p_unspc_fine>0)    ab_eff(p_unspc_fine)        = 0.0_RKIND
      if (p_dust_coarse>0)   ab_eff(p_dust_coarse)       = 0.0_RKIND
      if (p_unspc_coarse>0)  ab_eff(p_unspc_coarse)      = 0.0_RKIND
      if (p_ssalt_fine>0)    ab_eff(p_ssalt_fine)        = 0.0_RKIND
      if (p_ssalt_coarse>0)  ab_eff(p_ssalt_coarse)      = 0.0_RKIND
!
      if (p_polp_all>0)      ab_eff(p_polp_all)          = 0.0_RKIND
      if (p_polp_tree>0)     ab_eff(p_polp_tree)         = 0.0_RKIND
      if (p_polp_grass>0)    ab_eff(p_polp_grass)        = 0.0_RKIND
      if (p_polp_weed>0)     ab_eff(p_polp_weed)         = 0.0_RKIND
!
      if (p_pols_all>0)      ab_eff(p_pols_all)          = 0.0_RKIND
      if (p_pols_tree>0)     ab_eff(p_pols_tree)         = 0.0_RKIND
      if (p_pols_grass>0)    ab_eff(p_pols_grass)        = 0.0_RKIND
      if (p_pols_weed>0)     ab_eff(p_pols_weed)         = 0.0_RKIND

      !write(*,*) 'JLS, radiation parameters initialized',sc_eff,ab_eff

   end subroutine aero_rad_init

end module rad_data_mod
