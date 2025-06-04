!> Debugged and corrected module_wildfire_smoke_emissions

module module_wildfire_smoke_emissions
    use mpas_kind_types
    use mpas_smoke_init, only : p_smoke_fine, p_smoke_coarse
    use mpas_smoke_config
  
    implicit none
    private
    public :: compute_emission_factors, calculate_smoke_emissions, diurnal_cycle
  
  
  contains
  
  !-------------------------------------------------------------------------------
  ! Subroutine: compute_emission_factors
  ! Description:
  !   Computes a 2D emission factor map (EFs_map) based on vegetation fractions
  !   and ecoregions. Emission factors are calculated using a combination of
  !   flaming/smoldering coefficients and beta scaling. If bb_beta = 1, beta_kaiser
  !   is used only for Canadian forest/shrubland types; otherwise, beta_wooster is applied.
  !   EF_east factors are applied to Grasslands and Croplands based on region masks.
  !-------------------------------------------------------------------------------
    subroutine compute_emission_factors(EFs_map, vegfrac, bb_beta,                       &
                                        eco_id, efs_smold, efs_flam, efs_rsmold,         &
                                        fmc_avg, wetness,                                &
                                        nlcat, dt, gmt,                                  &
                                        ids, ide, jds, jde, kds, kde,                    &
                                        ims, ime, jms, jme, kms, kme,                    &
                                        its, ite, jts, jte, kts, kte)
  
      implicit none
      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN) :: its, ite, jts, jte, kts, kte
      REAL(RKIND), INTENT(IN) :: dt, gmt
      INTEGER, INTENT(IN) :: bb_beta, nlcat
      INTEGER, INTENT(IN), DIMENSION(ims:ime,jms:jme) :: eco_id
      REAL(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN) :: efs_smold, efs_flam, efs_rsmold, wetness, fmc_avg
      REAL(RKIND), INTENT(IN), DIMENSION(ims:ime, nlcat,jms:jme) :: vegfrac !TODO, check if correct way of using it
  
      !local variables
      INTEGER :: i, j, n, eco
      REAL(RKIND) :: ef, frac, beta_val, sum_frac
  
      !>-- Emissions factors parameters:frst,hwd,mxd,shrb,shrb_grs,grs,org_soil_mix,org_soil_woody
      REAL(RKIND), DIMENSION(1:8), PARAMETER :: EF_FLM = (/19.0, 9.4,  14.6, 9.3,  10.7, 13.3, 10.8, 11.8/) !!PM2.5 EFs for flaming US FS
      REAL(RKIND), DIMENSION(1:8), PARAMETER :: EF_SML = (/28.0, 37.7, 17.6, 36.6, 36.7, 38.4, 20.7, 23.3/) !PM2.5 EFs for smoldering US FS
      !>-- Emissions factors parameters GEUFE:IGBP 10, 12 14
      REAL(RKIND), DIMENSION(1:3), PARAMETER :: EF_east = (/9.29, 11.91, 8.20/) !!EFs eastern US GEUFE

      !CHARACTER(len=25), DIMENSION(20), PARAMETER :: land_cats = (/ &
      !                           "Evergreen_Nleaf_Frst", "Evergreen_Bleaf_Frst", "Deciduous_NleafFrst", "Deciduous_Bleaf_Frst", &
      !                           "Mixed_Frst", "Closed_Shrublands", "Open_Shrublands", "Woody_Savannas", "Savannas", "Grasslands", &
      !                           "Permanent_Wetlands", "Croplands", "Urban_Lands", "Cropland_NVeg", "snow_ice", "barren", &
      !                           "water", "Wooded_Tundra", "Mixed_Tundra", "Barren_Tundra" /)

      !>-- Biomass Combustion Coefficient (Beta Factor, β)
      !>-- Biomass Combustion Coefficient: SA,SAOS,AG,AGOS,EF,EFOS,PEAT
      REAL(RKIND), DIMENSION(1:7), PARAMETER :: beta_kaiser = (/0.78,0.26,0.29,0.13,0.49,1.55,5.87/) !Kaiser et al., 2012
      REAL(RKIND), PARAMETER :: beta_wooster = 0.38 !! Wooster et al. 2005
      REAL(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(OUT) :: EFs_map

      EFs_map = 0.0_RKIND
  
      DO j = jts, jte
        DO i = its, ite
          eco = eco_id(i,j)
          sum_frac = 0.0
          DO n = 1, 20
            sum_frac = sum_frac + vegfrac(i,j,n)
          END DO
          IF (sum_frac > 1.01_RKIND) THEN
            PRINT *, 'WARNING: vegfrac sum > 1 at (', i, ',', j, ') = ', sum_frac
          END IF
  
          DO n = 1, 20
            frac = vegfrac(i,j,n)
  
            SELECT CASE (n)
            CASE (1)  ! Evergreen_Nleaf_Frst
              ef = (0.75 * EF_FLM(1) + 0.25 * EF_SML(1))
              IF (bb_beta == 1 .AND. (eco == 2 .OR. eco == 3 .OR. eco == 4)) THEN
                beta_val = beta_kaiser(6)
              ELSE
                beta_val = beta_wooster
              END IF
  
            CASE (2) ! Evergreen_Bleaf_Frst
              ef = (0.75 * EF_FLM(1) + 0.25 * EF_SML(1))
              IF (bb_beta == 1 .AND. eco == 2) THEN
                beta_val = beta_kaiser(6)
              ELSE IF (bb_beta == 1 .AND. eco == 5) THEN
                beta_val = beta_kaiser(5)
              ELSE
                beta_val = beta_wooster
              END IF
  
            CASE (3) ! Deciduous_NleafFrst
              ef = (0.80 * EF_FLM(2) + 0.20 * EF_SML(2))
              IF (bb_beta == 1 .AND. eco == 2) THEN
                beta_val = beta_kaiser(6)
              ELSE IF (bb_beta == 1 .AND. eco == 5) THEN
                beta_val = beta_kaiser(5)
              ELSE
                beta_val = beta_wooster
              END IF
  
            CASE (4) ! Deciduous_Bleaf_Frst
              ef = (0.80 * EF_FLM(2) + 0.20 * EF_SML(2))
              IF (bb_beta == 1 .AND. eco == 2) THEN
                beta_val = beta_kaiser(6)
              ELSE IF (bb_beta == 1 .AND. eco == 5) THEN
                beta_val = beta_kaiser(5)
              ELSE
                beta_val = beta_wooster
              END IF
  
            CASE (5) ! Mixed_Frst
              ef = (0.85 * EF_FLM(3) + 0.15 * EF_SML(3))
              IF (bb_beta == 1 .AND. eco == 2) THEN
                beta_val = beta_kaiser(6)
              ELSE IF (bb_beta == 1 .AND. eco == 5) THEN
                beta_val = beta_kaiser(5)
              ELSE
                beta_val = beta_wooster
              END IF
  
            CASE (6,7) ! Shrublands
              ef = (0.95 * EF_FLM(4) + 0.05 * EF_SML(4))
              IF (bb_beta == 1 .AND. (eco == 2 .OR. eco == 3 .OR. eco == 4)) THEN
                beta_val = beta_kaiser(2)
              ELSE
                beta_val = beta_wooster
              END IF
  
            CASE (8) ! Woody Savannas
              ef = (0.95 * EF_FLM(5) + 0.05 * EF_SML(5))
              beta_val = beta_wooster
  
            CASE (9) ! Savannas
              ef = (0.95 * EF_FLM(6) + 0.05 * EF_SML(6))
              beta_val = beta_wooster
  
            CASE (10) ! Grasslands
              IF (eco == 10) THEN
                ef = EF_east(1)
              ELSE
                ef = EF_FLM(6)
              END IF
              beta_val = beta_wooster

            CASE (11) !Permanent_Wetlands
              ef = 18.9 
              beta_val = beta_wooster       
  
            CASE (12) ! Croplands
              IF (eco == 8 .OR. eco == 9) THEN
                ef = EF_east(2)
              ELSE
                ef = 8.2
              END IF
              beta_val = beta_wooster
  
            CASE (14) ! Cropland_NVeg
              IF (eco == 8 .OR. eco == 9) THEN
                ef = EF_east(3)
              ELSE
                ef = 8.2
              END IF
              beta_val = beta_wooster
  
            CASE (18) ! Wooded_Tundra
              ef = (0.7 * EF_FLM(8) + 0.3 * EF_SML(8))
              beta_val = beta_wooster
  
            CASE (19) ! Mixed_Tundra
              ef = (0.7 * EF_FLM(7) + 0.3 * EF_SML(7))
              beta_val = beta_wooster
  
            CASE DEFAULT
              ef = 0.0_RKIND
              beta_val = 0.0_RKIND
            END SELECT
  
            EFs_map(i,j) = EFs_map(i,j) + (frac * ef * beta_val)
          END DO
        END DO
      END DO
  
    end subroutine compute_emission_factors
  
  !-------------------------------------------------------------------------------
  ! Subroutine: calculate_smoke_emissions
  ! Description:
  !   Calculates smoke emissions (ebu) based on Fire Radiative Energy (fre_avg),
  !   emission factors (EFs_map), and target grid area. Emissions are estimated 
  !   for each temporal block, in units of micrograms per second per square meter.
  !   Only grid points with positive FRP are considered; emissions for zero/non-burning
  !   grid points are set to zero. 
  !   fre_avg declared always as a 3-D dummy
  !-------------------------------------------------------------------------------
    subroutine calculate_smoke_emissions(dt, julday, nlcat, EFs_map, fre_avg,       &
                                         ebb_dcycle, area, nblocks, ktau,           &
                                         bb_input_prevh, ebu,                       &
                                         ids, ide, jds, jde, kds, kde,              &
                                         ims, ime, jms, jme, kms, kme,              &
                                         its, ite, jts, jte, kts, kte)
  
      implicit none
      INTEGER, INTENT(IN) :: julday, nlcat
      INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde
      INTEGER, INTENT(IN) :: ims, ime, jms, jme, kms, kme
      INTEGER, INTENT(IN) :: its, ite, jts, jte, kts, kte
       ! Timestep, day, constants
      REAL(RKIND), INTENT(IN) :: dt
      integer,intent(in):: ktau
      INTEGER, INTENT(IN) :: ebb_dcycle, nblocks, bb_input_prevh

      REAL(RKIND), INTENT(IN),   DIMENSION(ims:ime, jms:jme, nblocks)   :: fre_avg
      REAL(RKIND), INTENT(IN),   DIMENSION(ims:ime, jms:jme)            :: EFs_map
      REAL(RKIND), INTENT(IN),   DIMENSION(ims:ime, jms:jme)            :: area !we need it for first v level 
      REAL(RKIND), INTENT(INOUT),DIMENSION(ims:ime, kms:kme, jms:jme)   :: ebu

      !Local variables
      INTEGER :: i, j, blk, hour_int
      REAL(RKIND) ::  hour_tmp
      REAL(RKIND), PARAMETER :: fg_to_ug = 1.0e6
      REAL(RKIND), PARAMETER :: to_s = 3600.0
       
     ! for EBB = 2, blk > 1
      IF ( ebb_dcycle .eq. 2 ) then
         ! Current integration hour
         hour_int = FLOOR(ktau*dt/3600.)
         ! Reset if > 24 and determine how many bb_input_prevh(s) have passed
         hour_tmp = MOD(hour_int, 24) / bb_input_prevh
         ! Make it an integer + 1
         blk = FLOOR(hour_tmp) + 1
         IF (blk < 1) blk = 1
         IF (blk > 24 / bb_input_prevh) blk = 24 / bb_input_prevh
      ELSE
         ! EBB = 1 has only one time
         blk = 1
      ENDIF
    
      !if ebb_dcycle == 2 then ebb_dc1 do not use blocks? should be 24
      DO j = jts, jte
        DO i = its, ite
           IF (fre_avg(i,j,blk) > 0.0_RKIND) THEN
              ebu(i,kts,j) = (fre_avg(i,j,blk) * EFs_map(i,j) * fg_to_ug) / ( area(i,j) * to_s) !we need the area of the grids
           ELSE
              ebu(i,kts,j) = 0.0_RKIND
           END IF
        END DO
      END DO
  
    end subroutine calculate_smoke_emissions
  
  !peak_hr only needed at first time step -> from wrapper?
  !-------------------------------------------------------------------------------
  ! Subroutine: diurnal_cycle
  ! Apply fire diurnal cycle emission scaling depending on fire type and time.
  ! Handles multiple fire age "blocks" (nblocks) using coef_bb_dc and fire_hist.
  !-------------------------------------------------------------------------------
  ! Apply diurnal cycle depending on fire type
subroutine diurnal_cycle(  dtstep,dz8w,rho_phy,pi,ebb_min,            &
                           chem,num_chem,julday,gmt,xlat,xlong,     &
                           fire_end_hr,peak_hr,time_int,coef_bb_dc, &
                           fire_hist,hwp,hwp_avg,hwp_prev_day,      &  !I think fire_hist replaced sc_factor
                           vegfrac, eco_id, nblocks,                 &
                           lu_nofire, lu_qfire, lu_sfire,           &
                           swdown,ebb_dcycle,ebu,fire_type,  &
                           q_vap, add_fire_moist_flux,              &
                           plume_beta_qv, hwp_alpha,                &
                           ids,ide, jds,jde, kds,kde,               &
                           ims,ime, jms,jme, kms,kme,               &
                           its,ite, jts,jte, kts,kte                )

  implicit none

  INTEGER,      INTENT(IN   ) ::  julday,num_chem,                  &
                                  ids,ide, jds,jde, kds,kde,        &
                                  ims,ime, jms,jme, kms,kme,        &
                                  its,ite, jts,jte, kts,kte

  REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem   ! shall we set num_chem=1 here?

  REAL(RKIND), DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         INTENT(INOUT ) ::                                   ebu, q_vap ! SRB: added q_vap
  INTEGER,      INTENT(IN), DIMENSION(ims:ime,jms:jme) :: eco_id
  REAL(RKIND),  INTENT(IN), DIMENSION(ims:ime,jms:jme,20) :: vegfrac      
 
  REAL(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp_prev_day ! JLS --- TODO
  REAL(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(INOUT)    :: fire_hist
 
  real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: xlat,xlong, swdown
  real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp, peak_hr, fire_end_hr !RAR: Shall we make fire_end integer?
  real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(OUT)    :: lu_nofire, lu_qfire, lu_sfire, coef_bb_dc
  real(RKIND), DIMENSION(ims:ime,jms:jme), INTENT(IN)     :: hwp_avg
  real(RKIND), DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: dz8w,rho_phy  !,rel_hum
  real(RKIND), INTENT(IN) ::  dtstep, gmt
  real(RKIND), INTENT(IN) ::  time_int, pi, ebb_min       ! RAR: time in seconds since start of simulation
  INTEGER, DIMENSION(ims:ime,jms:jme), INTENT(OUT) :: fire_type
  integer, INTENT(IN) ::  ebb_dcycle     ! RAR: this is going to be namelist dependent, ebb_dcycle=means 
  real(RKIND), INTENT(IN) :: plume_beta_qv ! SRB: namelist scaling factor for fire qv fluxes

  INTEGER, INTENT(IN) ::  nblocks ! e.g. 1, 4, or 24
  
  !>--local 
  logical, intent(in)  :: add_fire_moist_flux
  integer :: i,j,k,n,m, blk
  integer :: icall=0
  real(RKIND) :: conv_rho, conv, dm_smoke, dc_hwp, dc_gp, dc_fn !daero_num_wfa, daero_num_ifa !, lu_sum1_5, lu_sum12_14

  INTEGER, PARAMETER :: kfire_max=51    ! max vertical level for BB plume rise
  real(RKIND), PARAMETER :: ef_h2o=324.22  ! Emission factor for water vapor
  ! Constants for the fire diurnal cycle calculation ! JLS - needs to be
  ! defined below due to intent(in) of pi
  real(RKIND) :: coef_con
  real(RKIND) :: timeq, fire_age, age_hr, dt1,dt2,dtm         ! For BB emis. diurnal cycle calculation
  
  ! For Gaussian diurnal cycle
  real(RKIND), INTENT(IN) :: hwp_alpha  ! to scale up the wildfire emissions, Jordan please make this a namelist option
  real(RKIND), PARAMETER :: rinti=2.1813936e-8, ax2=3400., const2=130., &
                           coef2=10.6712963e-4, cx2=7200., timeq_max=3600.*24.
  !>-- Fire parameters: Fores west, Forest east, Shrubland, Savannas, Grassland, Cropland
  real(RKIND), dimension(1:5), parameter :: avg_fire_dur   = (/8.9, 4.2, 3.3, 3.0, 1.4/)
  real(RKIND), dimension(1:5), parameter :: sigma_fire_dur = (/8.7, 6.0, 5.5, 5.2, 2.4/)
  ! For fire diurnal cycle calculation
  !real(kind_phys), parameter :: avgx1=-2.0, sigmx1=0.7, C1=0.083  ! Ag fires
  !real(kind_phys), parameter :: avgx2=-0.1, sigmx2=0.8, C2=0.55   ! Grass fires, slash burns
  real(RKIND), parameter :: avgx1=0.,  sigmx1=2.2, C1=0.2  ! Ag fires
  real(RKIND), parameter :: avgx2=0.5, sigmx2=0.8, C2=1.1   ! Grass fires, slash burns

  timeq= gmt*3600._RKIND + real(time_int,4)
  timeq= mod(timeq,timeq_max)
  coef_con=1._RKIND/((2._RKIND*pi)**0.5)

  !peak_hr     = 0.
  fire_type   = 0
  lu_qfire    = 0.
  lu_sfire    = 0.
  lu_nofire   = 0.

  !RAR: change this to the fractional LU type; fire_type: 0- no fires, 1- Ag
  ! or urban fires, 2- prescribed fires in wooded area, 3- wildfires
  if (ebb_dcycle==2) then
    do j=jts,jte
      do i=its,ite
        if (ebu(i,kts,j)<ebb_min) then
           fire_type(i,j) = 0
           lu_nofire(i,j) = 1.0
        else
          ! cropland, urban, cropland/natural mosaic, barren and sparsely
          ! vegetated and non-vegetation areas:
          lu_qfire(i,j) = lu_nofire(i,j) + vegfrac(i,12,j) + vegfrac(i,13,j) + vegfrac(i,14,j) + vegfrac(i,16,j)
          ! Savannas and grassland fires, these fires last longer than the Ag fires:
          lu_sfire(i,j) = lu_qfire(i,j) + vegfrac(i,8,j) + vegfrac(i,9,j) + vegfrac(i,10,j)
          if (lu_nofire(i,j)>0.95) then ! no fires
            fire_type(i,j) = 0
          else if (lu_qfire(i,j)>0.9) then   ! Ag. and urban fires
            fire_type(i,j) = 1
          else if (eco_id(i,j)==8) then
            fire_type(i,j) = 2    ! slash burn and wildfires in the east, eastern temperate forest ecosystem
            ! SRB: Eastern wildland fires if FRE>1E6MJ and eco region is 8 and the fire is not older than 7 hrs
            if ( ebu(i,kts,j)*3600*(1/0.416) .ge. 1.E6 .and. fire_end_hr(i,j) .le. 8 ) then
              fire_type(i,j) = 4
            endif 
          else if (lu_sfire(i,j)>0.8) then
            fire_type(i,j) = 3    ! savanna and grassland fires
          else
            fire_type(i,j) = 4    ! potential wildfires     
          end if
        end if
      end do
    end do
  endif


  if (ebb_dcycle==2) then 
    !do blk = 1, nblocks
     do j=jts,jte
       do i=its,ite
        fire_age= MAX(0.01_RKIND,time_int/3600. + (fire_end_hr(i,j)-1.0))  !One hour delay is due to the latency of the RAVE files, hours

          SELECT CASE ( fire_type(i,j) )   !Ag, urban fires, bare land etc.
          CASE (1)
             ! these fires will have exponentially decreasing diurnal cycle,
             !coef_bb_dc(i,j) = coef_con*1._kind_phys/(sigma_fire_dur(5) *fire_age) *                          &
             !                exp(- ( log(fire_age) - avg_fire_dur(5))**2 /(2._kind_phys*sigma_fire_dur(5)**2 ))
            coef_bb_dc(i,j)= C1/(sigmx1* fire_age)* exp(- (log(fire_age) - avgx1)**2 /(2.*sigmx1**2 ) )

             IF ( dbg_opt .AND. time_int<5000.) then
               WRITE(6,*) 'i,j,peak_hr(i,j),fire_type(i,j) ',i,j,peak_hr(i,j),fire_type(i,j)
               WRITE(6,*) 'coef_bb_dc(i,j) ',coef_bb_dc(i,j)
             END IF

          CASE (2)    ! Savanna and grassland fires, or fires in the eastern US 
            !  coef_bb_dc(i,j) = coef_con*1._kind_phys/(sigma_fire_dur(4) *fire_age) *                          &
            !                  exp(- ( log(fire_age) - avg_fire_dur(4))**2 /(2._kind_phys*sigma_fire_dur(4)**2 ))
            coef_bb_dc(i,j)= C2/(sigmx2* fire_age)* exp(- (log(fire_age) - avgx2)**2 /(2.*sigmx2**2 ) )

              IF ( dbg_opt .AND. time_int<5000.) then
                WRITE(6,*) 'i,j,peak_hr(i,j),fire_type(i,j) ',i,j,peak_hr(i,j),fire_type(i,j)
                WRITE(6,*) 'coef_bb_dc(i,j) ',coef_bb_dc(i,j)
              END IF

          CASE (3,4)    ! wildfires 
             IF (swdown(i,j)<.1 .AND. fire_age> 6. .AND. fire_hist(i,j)>0.75) THEN
                 fire_hist(i,j)= 0.75_RKIND
             ENDIF
             IF (swdown(i,j)<.1 .AND. fire_age> 18. .AND. fire_hist(i,j)>0.5) THEN
                 fire_hist(i,j)= 0.5_RKIND
             ENDIF
             IF (swdown(i,j)<.1 .AND. fire_age> 36. .AND. fire_hist(i,j)>0.25) THEN
                 fire_hist(i,j)= 0.25_RKIND
             ENDIF
             ! this is based on hwp, hourly or instantenous TBD
             dc_hwp= hwp_alpha * hwp(i,j)/ MAX(10._RKIND,hwp_prev_day(i,j)) + (1-hwp_alpha) * (hwp(i,j)/MAX(10._RKIND,hwp_avg(i,j)))
             dc_hwp= MAX(0._RKIND,dc_hwp)
             dc_hwp= MIN(20._RKIND,dc_hwp)

       
             coef_bb_dc(i,j) = fire_hist(i,j)* dc_hwp     ! RAR: scaling factor is applied to the forest fires only, except the eastern US

             IF ( dbg_opt .AND. time_int<5000.) then
               WRITE(6,*) 'i,j,fire_hist(i,j),peak_hr(i,j) ', i,j,fire_hist(i,j),peak_hr(i,j)
               WRITE(6,*) 'dc_gp,dc_hwp,dc_fn,coef_bb_dc(i,j) ',dc_gp,dc_hwp,dc_fn,coef_bb_dc(i,j)
             END IF

          CASE DEFAULT
          END SELECT
       !enddo
      enddo
    enddo
  endif
end subroutine diurnal_cycle

end module module_wildfire_smoke_emissions
