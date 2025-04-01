!>\file  module_wildfire_smoke_emissions.F90
!! This file contains the MPAS-Aerosols/RRFS wildfire emission module

module module_anthro_emissions
!
!  This module developed by Johana Romero-Alvarez and Jordan Schnell (NOAA GSL)
!  For serious questions contact johana.romero-alvarez@noaa.gov
!
  use mpas_kind_types
  use mpas_smoke_init

  implicit none

  private

  public :: mpas_smoke_anthro_emis_driver

contains


  subroutine mpas_smoke_anthro_emis_driver(dt,gmt,julday,kemit,                      &
                           xlat,xlong, chem,num_chem,dz8w,t_phy,rho_phy,             &    
                           e_ant_in, e_ant_out, num_e_ant_in, num_e_ant_out,         &
                           index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
                           index_e_ant_in_smoke_fine, index_e_ant_in_smoke_coarse,   &
                           index_e_ant_in_dust_fine,  index_e_ant_in_dust_coarse,    &
                           index_e_ant_in_no3_a_fine, index_e_ant_in_so4_a_fine,     &
                           index_e_ant_in_nh4_a_fine,                                &
                           index_e_ant_in_so2, index_e_ant_in_nh3,                   & 
                           index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
                           index_e_ant_out_smoke_fine, index_e_ant_out_smoke_coarse, &
                           index_e_ant_out_dust_fine, index_e_ant_out_dust_coarse,   &
                           index_e_ant_out_no3_a_fine, index_e_ant_out_so4_a_fine,   &
                           index_e_ant_out_nh4_a_fine,                               &
                           index_e_ant_out_so2, index_e_ant_out_nh3,                 &
                           ids,ide, jds,jde, kds,kde,                                &
                           ims,ime, jms,jme, kms,kme,                                &
                           its,ite, jts,jte, kts,kte                                 )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: julday, num_chem, kemit,           &
                                  ids,ide, jds,jde, kds,kde,         &
                                  ims,ime, jms,jme, kms,kme,         &
                                  its,ite, jts,jte, kts,kte,         &
                                  num_e_ant_in, num_e_ant_out,       &
           index_e_ant_in_unspc_fine, index_e_ant_in_unspc_coarse,   &
           index_e_ant_in_smoke_fine, index_e_ant_in_smoke_coarse,   &
           index_e_ant_in_dust_fine,  index_e_ant_in_dust_coarse,    &
           index_e_ant_in_no3_a_fine, index_e_ant_in_so4_a_fine,     &
           index_e_ant_in_nh4_a_fine,                                &
           index_e_ant_in_so2, index_e_ant_in_nh3,                   &
           index_e_ant_out_unspc_fine, index_e_ant_out_unspc_coarse, &
           index_e_ant_out_smoke_fine, index_e_ant_out_smoke_coarse, &
           index_e_ant_out_dust_fine, index_e_ant_out_dust_coarse,   &
           index_e_ant_out_no3_a_fine, index_e_ant_out_so4_a_fine,   &
                           index_e_ant_out_nh4_a_fine,                               &
                           index_e_ant_out_so2, index_e_ant_out_nh3

   REAL(RKIND), INTENT(IN    ) :: dt,gmt

   REAL(RKIND),DIMENSION(ims:ime,jms:jme),INTENT(IN) :: xlat,xlong
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme),INTENT(IN) :: dz8w,rho_phy,t_phy
   REAL(RKIND),DIMENSION(ims:ime,1:kemit,jms:jme,1:num_e_ant_in), INTENT(IN)    :: e_ant_in
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_e_ant_out),INTENT(INOUT) :: e_ant_out
   REAL(RKIND),DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem), INTENT(INOUT)     :: chem
                                                                               
  ! local
   INTEGER :: i,j,k,n
   REAL(RKIND) :: conv_aer, conv_gas

   REAL(RKIND), PARAMETER :: rwc_t_thresh = 283.15 ! [ 50 F]


   do j = jts,jte
   do k = kts, kemit
   do i = its,ite
!  
!     NEMO emissions are in g/s (per grid cell) == 1 grid cell / 1e6 m2
!     Want to convert to ug/m2/s
!     scaling factor of 1e-6 for /m2
!     scaling factor of 1e6  for g -> ug
      conv_aer = dt / (rho_phy(i,k,j) *  dz8w(i,k,j))
      conv_gas = 4.828e-4_RKIND/rho_phy(i,k,j)*dt/(dz8w(i,k,j) * 60._RKIND)
!
!      if (p_smoke_fine   .gt. 0) chem(i,k,j,p_smoke_fine)   = chem(i,k,j,p_smoke_fine)   + conv_aer*e_ant_in(i,k,j,index_e_ant_in_smoke_fine)
!      if (p_smoke_coarse .gt. 0) chem(i,k,j,p_smoke_coarse) = chem(i,k,j,p_smoke_coarse) + conv_aer*e_ant_in(i,k,j,index_e_ant_in_smoke_coarse)
      if (p_unspc_fine   .gt. 0) chem(i,k,j,p_unspc_fine)   = chem(i,k,j,p_unspc_fine)   + conv_aer*e_ant_in(i,k,j,index_e_ant_in_unspc_fine)
      if (p_unspc_coarse .gt. 0) chem(i,k,j,p_unspc_coarse) = chem(i,k,j,p_unspc_coarse) + conv_aer*e_ant_in(i,k,j,index_e_ant_in_unspc_coarse)
      if (p_dust_fine    .gt. 0) chem(i,k,j,p_dust_fine)    = chem(i,k,j,p_dust_fine)    + conv_aer*e_ant_in(i,k,j,index_e_ant_in_dust_fine)
      if (p_dust_coarse  .gt. 0) chem(i,k,j,p_dust_coarse)  = chem(i,k,j,p_dust_coarse)  + conv_aer*e_ant_in(i,k,j,index_e_ant_in_dust_coarse)
      if (p_no3_a_fine   .gt. 0) chem(i,k,j,p_no3_a_fine)   = chem(i,k,j,p_no3_a_fine)   + conv_aer*e_ant_in(i,k,j,index_e_ant_in_no3_a_fine)
      if (p_so4_a_fine   .gt. 0) chem(i,k,j,p_so4_a_fine)   = chem(i,k,j,p_so4_a_fine)   + conv_aer*e_ant_in(i,k,j,index_e_ant_in_so4_a_fine)
      if (p_nh4_a_fine   .gt. 0) chem(i,k,j,p_nh4_a_fine)   = chem(i,k,j,p_nh4_a_fine)   + conv_aer*e_ant_in(i,k,j,index_e_ant_in_nh4_a_fine)
!
      if (p_nh3          .gt. 0) chem(i,k,j,p_nh3)          = chem(i,k,j,p_nh3)          + conv_gas*e_ant_in(i,k,j,index_e_ant_in_nh3)
      if (p_so2          .gt. 0) chem(i,k,j,p_so2)          = chem(i,k,j,p_so2)          + conv_gas*e_ant_in(i,k,j,index_e_ant_in_so2)
!
   enddo ! i
   enddo ! k
   enddo ! j


!   do j = jts, jte
!   do i = its, ite
!      if ( t_phy(i,kts,j) .lt. rwc_t_thresh ) then
!         emis = (42.12_RKIND - 0.79_RKIND*t_min(i,kts,j)) / total_rwc_emis(i,j)
!      if (p_smoke_fine   .gt. 0) chem(i,k,j,p_smoke_fine)   = chem(i,k,j,p_smoke_fine)   + conv_aer*e_ant_in(i,k,j,index_e_ant_in_smoke_fine)
!      if (p_smoke_coarse .gt. 0) chem(i,k,j,p_smoke_coarse) = chem(i,k,j,p_smoke_coarse) + conv_aer*e_ant_in(i,k,j,index_e_ant_in_smoke_coarse)
!      endif
!   enddo
!   enddo        

  end subroutine mpas_smoke_anthro_emis_driver

end module module_anthro_emissions
