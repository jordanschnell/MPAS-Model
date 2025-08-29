!>\file  module_wildfire_smoke_emissions.F90
!! This file contains the MPAS-Aerosols/RRFS wildfire emission module

module module_rwc_emissions
!
!  This module developed by Johana Romero-Alvarez and Jordan Schnell (NOAA GSL)
!  For serious questions contact johana.romero-alvarez@noaa.gov
!
  use mpas_kind_types
  use mpas_smoke_init

  implicit none

  private

  public :: mpas_smoke_rwc_emis_driver

contains


  subroutine mpas_smoke_rwc_emis_driver(dt,gmt,julday,kemit,                         &
                           xlat,xlong, chem,num_chem,dz8w,t_phy,rho_phy,             &
                           rwc_emis_scale_factor,                                    &
                           online_rwc_emis,RWC_denominator,RWC_annual_sum,           &
                           RWC_annual_sum_smoke_fine, RWC_annual_sum_smoke_coarse,   &
                           RWC_annual_sum_unspc_fine, RWC_annual_sum_unspc_coarse,   &
                           e_ant_out,  num_e_ant_out,         &
                           index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
                           index_e_ant_in_smoke_fine, index_e_ant_in_smoke_coarse,   &
                           index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
                           index_e_ant_out_smoke_fine, index_e_ant_out_smoke_coarse, &
                           ids,ide, jds,jde, kds,kde,                                &
                           ims,ime, jms,jme, kms,kme,                                &
                           its,ite, jts,jte, kts,kte                                 )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: julday, num_chem, kemit,           &
                                  ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte,         &
                                  num_e_ant_out,       &
           index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
           index_e_ant_in_smoke_fine, index_e_ant_in_smoke_coarse,   &
           index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
           index_e_ant_out_smoke_fine, index_e_ant_out_smoke_coarse, &
           online_rwc_emis

   REAL(RKIND), INTENT(IN    ) :: dt,gmt,rwc_emis_scale_factor

   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: xlat,xlong
   REAL(RKIND),DIMENSION(ims:ime,1:kemit,jms:jme),INTENT(IN) :: RWC_annual_sum_smoke_fine,RWC_annual_sum_smoke_coarse, &
                                                        RWC_annual_sum_unspc_fine,RWC_annual_sum_unspc_coarse, &
                                                        RWC_annual_sum
   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: RWC_denominator
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN) :: dz8w,rho_phy,t_phy
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_e_ant_out),INTENT(INOUT) :: e_ant_out
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT)     :: chem
                                                                               
  ! local
   INTEGER :: i,j,k,n
   REAL(RKIND) :: conv_aer, conv_gas, emis, t_phy_f, frac

   REAL(RKIND), PARAMETER :: rwc_t_thresh = 283.15_RKIND ! [ 50 F]

! For now, just the surface
   k = kts

   if (online_rwc_emis .eq. 1 ) then  

   do j = jts, jte
   do i = its, ite
     ! Is it cold enough to emit wood burning emissions?
      if ( t_phy(i,k,j) .lt. rwc_t_thresh ) then
        ! Conversion factor for aerosol emissions (ug/m2/s) --> ug/kg
         conv_aer = dt / (rho_phy(i,k,j) *  dz8w(i,k,j))
        ! Conversion factor for gas phase emissions (mol/m2/s) --> ppm/ppm
         conv_gas = 60._RKIND * 1.E6_RKIND * 4.828E-4_RKIND * dt / ( rho_phy(i,k,j) * dz8w(i,k,j) )
        ! Convert temperature to Fahrenheit
         t_phy_f = 9._RKIND/5._RKIND * (t_phy(i,k,j)-273.15_RKIND) + 32._RKIND
        ! Calculate the fraction of total emisisons based on the linear equation: TODO, make coefficients namelist?
         frac = (42.12_RKIND - 0.79_RKIND*t_phy_f) / RWC_denominator(i,j)
         if (p_smoke_fine   .gt. 0 .and. index_e_ant_out_smoke_fine .gt. 0 ) then
            emis = rwc_emis_scale_factor * conv_aer * frac * RWC_annual_sum_smoke_fine(i,1,j)
            e_ant_out(i,k,j,index_e_ant_out_smoke_fine) = e_ant_out(i,k,j,index_e_ant_out_smoke_fine) + emis
         endif
         if (p_smoke_coarse   .gt. 0 .and. index_e_ant_out_smoke_coarse .gt. 0 ) then
            emis = rwc_emis_scale_factor * conv_aer * frac * RWC_annual_sum_smoke_coarse(i,1,j)
            e_ant_out(i,k,j,index_e_ant_out_smoke_coarse) = e_ant_out(i,k,j,index_e_ant_out_smoke_coarse) + emis
         endif
      endif
   enddo
   enddo        
   
   endif

  end subroutine mpas_smoke_rwc_emis_driver

end module module_rwc_emissions
