!>\file  module_tactic_sna.F90
!! This file contains the MPAS-Aerosols/RRFS SNA module

module module_tactic_sna
!
!  This module developed by Jordan Schnell (CIRES/NOAA GSL) following 
!  Druge et al., (2019) - https://doi.org/10.5194/acp-19-3707-2019
!  For serious questions contact jordan.schnell@noaa.gov
!
  use mpas_kind_types
  use mpas_smoke_init
  use mpas_smoke_config, only : pi

  implicit none

  private

  public :: mpas_smoke_tactic_sna_driver

contains

  subroutine mpas_smoke_tactic_sna_driver (                                          &
                           dt, chem, num_chem,                                       &
                           relhum, t_phy, dz8w, rho_phy,                             &
                           nifa, nwfa, hno3_bkgd,                                    &
                           ids,ide, jds,jde, kds,kde,                                &
                           ims,ime, jms,jme, kms,kme,                                &
                           its,ite, jts,jte, kts,kte                                 )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: num_chem,                          &
                                  ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte
   REAL(RKIND),INTENT(IN) :: dt
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN) :: dz8w, rho_phy, t_phy, relhum, nifa, nwfa
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(INOUT) :: hno3_bkgd
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT) :: chem
                                                                               
  ! local
   INTEGER :: i,j,k,n
   REAL(RKIND) :: TA, TS, TN, DRH, Kp, RH1, p1, p2, p3, pc, TA_star, GAMMA, RSA
   REAL(RKIND) :: NH4SO4, TA_rem, NH3SO4, NH4NO3, NH3_u, NH3_eq, HNO3_eq, tau
   REAL(RKIND) :: nh3_m, nh4_a_fine_m, so4_a_fine_m, hno3_m, no3_a_fine_m
   REAL(RKIND) :: Dg, lambda, Kn, khno3
  ! parameters
   REAL(RKIND), PARAMETER :: Ac = 6.022E23_RKIND
   REAL(RKIND), PARAMETER :: Ra = 8.314_RKIND
   REAL(RKIND), PARAMETER :: da = 4.5e-10_RKIND
   REAL(RKIND), PARAMETER :: ma = 0.029_RKIND
   REAL(RKIND), PARAMETER :: mhno3 = 0.063_RKIND
   REAL(RKIND), PARAMETER :: Dw = 0.5E-6_RKIND
   REAL(RKIND), PARAMETER :: mw_nh3 = 17.031_RKIND   ! g/mol
   REAL(RKIND), PARAMETER :: mw_nh4 = 18.04_RKIND   ! g/mol
   REAL(RKIND), PARAMETER :: mw_so4 = 96.06_RKIND   ! g/mol
   REAL(RKIND), PARAMETER :: mw_no3 = 62.0049_RKIND ! g/mol
   REAL(RKIND), PARAMETER :: mw_hno3 = 63.01_RKIND  ! g/mol

   do j = jts,jte
   do k = kts,kte
   do i = its,ite
    
    ! Convert to molar concentrations    
      nh3_m        = chem(i,k,j,p_nh3)        * rho_phy(i,k,j)/mw_nh3
      nh4_a_fine_m = chem(i,k,j,p_nh4_a_fine) * rho_phy(i,k,j)/mw_nh4
      so4_a_fine_m = chem(i,k,j,p_so4_a_fine) * rho_phy(i,k,j)/mw_so4
      hno3_m       = hno3_bkgd(i,k,j)         * rho_phy(i,k,j)/mw_hno3
      no3_a_fine_m = chem(i,k,j,p_no3_a_fine) * rho_phy(i,k,j)/mw_no3

    ! Determine total molar ammonia, nitrate, sulfate
      TA = nh3_m + nh4_a_fine_m
      TS = so4_a_fine_m
      TN = hno3_m + no3_a_fine_m

    ! Deterimine the sulfate regime
      RSA = TS / TA
      if ( RSA .gt. 1._RKIND) then
         GAMMA = 1.0_RKIND
      elseif ( RSA .gt. 0.5_RKIND .and. RSA .le. 1. ) then
         GAMMA = 1.5_RKIND
      else
         GAMMA = 2.0_RKIND
      endif

    ! Calculate the deliquescence relative humidity
      DRH = exp(723.7_RKIND / t_phy(i,k,j) + 1.6954_RKIND)

    ! Calculate some constants
      RH1 = 1._RKIND - relhum(i,k,j)
      p1  = exp(-135.94_RKIND + 8763._RKIND /t_phy(i,k,j) + 19.12_RKIND*log(t_phy(i,k,j)))
      p2  = exp(-122.65_RKIND + 9969._RKIND /t_phy(i,k,j) + 16.22_RKIND*log(t_phy(i,k,j)))
      p3  = exp(-182.61_RKIND + 13875._RKIND/t_phy(i,k,j) + 24.46_RKIND*log(t_phy(i,k,j)))
      pc  = (p1 - p2*RH1 + p3*RH1**2._RKIND)*RH1**1.75_RKIND

    ! Deterimine the humidity dependent equilibrium constant 
      Kp = exp(118.87_RKIND - 24084._RKIND/t_phy(i,k,j) - 6.025_RKIND*log(t_phy(i,k,j)))
      if ( relhum(i,k,j) .ge. DRH ) then
         Kp = KP * pc
      endif
 
    ! Moles of NH3 after SO4 neutralization
      TA_star = MAX(0._RKIND,TA - GAMMA*TS)

    ! Set the initial conditions
      TA_rem = TA
 
    ! Total NH4 associated with exisiting SO4
      NH4SO4 = min(TA_rem,GAMMA*TS)    

    ! Where the amount reomoved is = 
      TA_rem = max(0._RKIND,TA_rem-NH4SO4)

    ! Calculate
      NH3SO4 = MAX(0._RKIND,NH4SO4 - nh4_a_fine_m)
    
    ! Equilibrium Concentration 
      if ( TN*TA_star .le. Kp ) then 
         NH4NO3 = 0._RKIND
      else
         NH4NO3 = 0.5_RKIND *  (TA_star + TN - sqrt((TA_star+TN)**2._RKIND - 4._RKIND*(TN*TA_star - Kp)))
      endif

    ! Updated NH3 from any SO4 neuralization
      NH3_u = max(0._RKIND,nh3_m - NH3SO4)

    ! Equillibrium concentrations
      NH3_eq = max(0._RKIND,TA_star - NH4NO3)
      HNO3_eq = max(0._RKIND,TN - NH4NO3)

    ! Determine the time step
      Dg     = (3._RKIND / (8._RKIND * Ac * rho_phy(i,k,j) * da**2._RKIND)) * &
               (  ((ma*Ra*t_phy(i,k,j))/(2._RKIND*pi)) * ((ma+mhno3)/mhno3) )**0.5_RKIND
      lambda = (3._RKIND*Dg)/sqrt((8._RKIND*Ra*t_phy(i,k,j))/(pi*mhno3))
      Kn     = 2._RKIND * lambda / Dw
      khno3  = (2._RKIND * pi * Dw * Dg) / &
               (1._RKIND + ((4._RKIND*Kn)/(3._RKIND*lambda)) * &
               (1._RKIND - (0.47_RKIND*lambda)/(1._RKIND + Kn)))
!      if ( nifa(i,k,j) * khno3 .eq. 0. ) then
!         tau = 1.e10_RKIND * khno3
!         if ( tau .eq. 0._RKIND) tau = 1.e10
!      else
!         tau    = 1._RKIND / ( nifa(i,k,j) * khno3)
!      endif
      tau = 10._RKIND * 60._RKIND

    ! Update the chemistry array, converting back to kg/kg
      hno3_bkgd(i,k,j)         = mw_hno3/rho_phy(i,k,j) * &
                                              (hno3_m - (1._RKIND - exp(dt/tau)) * &
                                              (hno3_m - HNO3_eq))
      chem(i,k,j,p_nh3)        = mw_nh3/rho_phy(i,k,j) * &
                                              (NH3_u - (1._RKIND - exp(dt/tau)) * (NH3_u - NH3_eq))
      chem(i,k,j,p_nh4_a_fine) = mw_nh4/rho_phy(i,k,j) * &
                                             (TA_star - chem(i,k,j,p_nh3)*rho_phy(i,k,j)/mw_nh3)
      chem(i,k,j,p_no3_a_fine) = mw_no3/rho_phy(i,k,j) * (TN - hno3_m)
      chem(i,k,j,p_so4_a_fine) = (mw_nh4+mw_so4)/rho_phy(i,k,j)*NH4SO4 - chem(i,k,j,p_nh4_a_fine)
      
!      hno3_bkgd(i,k,j)         = max(0._RKIND,mw_hno3/rho_phy(i,k,j) * &
!                                              (hno3_m - (1._RKIND - exp(dt/tau)) * &
!                                              (hno3_m - HNO3_eq)))
!      chem(i,k,j,p_nh3)        = max(0._RKIND,mw_nh3/rho_phy(i,k,j) * &
!                                              (NH3_u - (1._RKIND - exp(dt/tau)) * (NH3_u - NH3_eq)))
!      chem(i,k,j,p_nh4_a_fine) = max(0._RKIND,mw_nh4/rho_phy(i,k,j) * &
!                                             (TA_star - chem(i,k,j,p_nh3)*rho_phy(i,k,j)/mw_nh3))
!      chem(i,k,j,p_no3_a_fine) = max(0._RKIND,mw_no3/rho_phy(i,k,j) * (TN - hno3_m))
!      chem(i,k,j,p_so4_a_fine) = (max(0._RKIND,mw_nh4+mw_so4)/rho_phy(i,k,j)*NH4SO4 - chem(i,k,j,p_nh4_a_fine))

   enddo ! i
   enddo ! k
   enddo ! j

  end subroutine mpas_smoke_tactic_sna_driver

end module module_tactic_sna
