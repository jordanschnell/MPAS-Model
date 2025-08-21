!>\file  module_add_emiss_burn.F90
!! This file adds the biomass burning emissions to the smoke field.

module module_add_emiss_burn
!RAR: significantly modified for the new BB emissions
  use mpas_kind_types
  use mpas_smoke_config
  use mpas_smoke_init
CONTAINS
  subroutine add_emis_burn(dtstep,dz8w,rho_phy,pi,ebb_min,          &
                           chem,num_chem,julday,gmt,xlat,xlong,     &
                           fire_end_hr, peak_hr,time_int,           &
                           coef_bb_dc, fire_hist, hwp, hwp_prevd,   &
                           swdown,ebb_dcycle, ebu,ebu_coarse,ebu_ch4, &
                           fire_type,&
                           q_vap, add_fire_moist_flux,              &
                           sc_factor, aod3d,                        &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte                )

   IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: julday, num_chem,                 &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte

   real(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, 1:num_chem ),                 &
         INTENT(INOUT ) ::                                   chem

   real(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         INTENT(INOUT ) ::                                   ebu, ebu_coarse, ebu_ch4, q_vap ! SRB: added q_vap

   real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: xlat,xlong, swdown
   real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp, peak_hr, fire_end_hr !RAR: Shall we make fire_end integer?
   real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(INOUT)  :: coef_bb_dc    ! RAR:
   real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp_prevd
   real(RKIND), INTENT(IN) :: sc_factor
   real(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: dz8w,rho_phy  !,rel_hum
   real(RKIND), INTENT(IN) ::  dtstep, gmt
   real(RKIND), INTENT(IN) ::  time_int, pi, ebb_min       ! RAR: time in seconds since start of simulation
   INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: fire_type
   integer, INTENT(IN) ::  ebb_dcycle     ! RAR: this is going to be namelist dependent, ebb_dcycle=means 
   real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: fire_hist
   real, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT) ::  aod3d
!>--local 
   logical, intent(in)  :: add_fire_moist_flux
   integer :: i,j,k,n,m
   integer :: icall=0
   real(RKIND) :: conv_gas, conv_rho, conv, dm_smoke, dm_smoke_coarse, dm_ch4, dc_hwp, dc_gp, dc_fn, ext2 !daero_num_wfa, daero_num_ifa !, lu_sum1_5, lu_sum12_14
   INTEGER, PARAMETER :: kfire_max=51    ! max vertical level for BB plume rise
   real(RKIND), PARAMETER :: ef_h2o=324.22  ! Emission factor for water vapor ! TODO, REFERENCE
   real(RKIND), PARAMETER :: sc_me= 4.0, ab_me=0.5     ! m2/g, scattering and absorption efficiency for smoke

! For Gaussian diurnal cycle

     if (mod(int(time_int),1800) .eq. 0) then
        icall = 0
     endif

     ext2= sc_me + ab_me
     do j=jts,jte
      do i=its,ite
       do k=kts,kfire_max
          if (ebu(i,k,j)<ebb_min) cycle

           if (ebb_dcycle==1) then
            conv= dtstep/(rho_phy(i,k,j)* dz8w(i,k,j))
            conv_gas = dtstep * 0.02897 / (rho_phy(i,k,j)* dz8w(i,k,j)) 
           elseif (ebb_dcycle==2) then
            conv= coef_bb_dc(i,j)*dtstep/(rho_phy(i,k,j)* dz8w(i,k,j))
           endif
            
            
           if ( p_smoke_fine > 0 ) then
              dm_smoke= conv*ebu(i,k,j)
              chem(i,k,j,p_smoke_fine) = chem(i,k,j,p_smoke_fine) + dm_smoke
              chem(i,k,j,p_smoke_fine) = MIN(MAX(chem(i,k,j,p_smoke_fine),epsilc),5.e+3_RKIND)        
              aod3d(i,k,j)= 1.e-6* ext2* chem(i,k,j,p_smoke_fine)*rho_phy(i,k,j)*dz8w(i,k,j)
           endif

           if ( p_smoke_coarse > 0 ) then 
              dm_smoke_coarse= conv*ebu_coarse(i,k,j)
              chem(i,k,j,p_smoke_coarse) = chem(i,k,j,p_smoke_coarse) + dm_smoke_coarse
              chem(i,k,j,p_smoke_coarse) = MIN(MAX(chem(i,k,j,p_smoke_coarse),epsilc),5.e+3_RKIND)        
           endif
           if ( p_ch4 > 0 ) then 
              dm_ch4= conv_gas*ebu_ch4(i,k,j)
              chem(i,k,j,p_ch4) = chem(i,k,j,p_ch4) + dm_ch4
              chem(i,k,j,p_ch4) = MIN(MAX(chem(i,k,j,p_ch4),epsilc),5.e+3_RKIND)        
           endif

           ! SRB: Modifying Water Vapor content based on Emissions
           if (add_fire_moist_flux) then
             q_vap(i,k,j) = q_vap(i,k,j) + (dm_smoke * ef_h2o * 1.e-9)  ! kg/kg:used 1.e-9 as dm_smoke is in ug/kg
             q_vap(i,k,j) = MIN(MAX(q_vap(i,k,j),0._RKIND),1.e+3_RKIND)
           endif

           if ( dbg_opt .and. (k==kts .OR. k==kfire_max) .and. (icall .le. n_dbg_lines) ) then
          endif
       enddo
       icall = icall + 1
      enddo
     enddo

    END subroutine add_emis_burn

END module module_add_emiss_burn

