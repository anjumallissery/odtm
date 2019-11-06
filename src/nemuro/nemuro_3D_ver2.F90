
module nemuro_mod

use mpp_mod, only : mpp_npes, mpp_pe, mpp_error, stdout, FATAL, WARNING, NOTE, mpp_init
use mpp_mod, only : mpp_exit, mpp_max, mpp_sum, mpp_sync, mpp_root_pe
use mpp_mod, only: mpp_clock_id, mpp_clock_begin, mpp_clock_end
use fms_mod,  only : field_exist, field_size, read_data, fms_init, fms_end
use fms_io_mod, only : register_restart_field, restart_file_type, save_restart, restore_state
use fms_io_mod, only : open_namelist_file, open_file, close_file, file_exist

use mpp_domains_mod, only : domain2d, domain1d, mpp_define_layout, mpp_define_domains
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_domain_components, mpp_update_domains
use mpp_domains_mod, only : mpp_get_data_domain, CGRID_SW, BITWISE_EXACT_SUM, mpp_global_sum
use diag_manager_mod, only : diag_manager_init, register_diag_field, register_static_field
use diag_manager_mod, only : diag_axis_init, send_data, diag_manager_end
use diag_data_mod, only : FILL_VALUE
use data_override_mod, only : data_override_init, data_override
use time_manager_mod, only : set_calendar_type, NO_CALENDAR, JULIAN, NOLEAP, date_to_string
use time_manager_mod, only : time_type, set_time, set_date, operator(+), assignment(=), operator(/)
use time_manager_mod, only : print_date, print_time, set_ticks_per_second, increment_date, operator(>=)
use time_manager_mod, only : operator(-), get_date

use size_mod, only : lm, rkmh

implicit none

private

integer, parameter :: ntracer_nm=15
integer, public :: ntracers=0
character(len=16), parameter :: trnms(ntracer_nm) =["PS","PL","ZS","ZL","ZP","NO3","NH4", &
                                "POM","DOM","SIOH4","Opal","Ca","CaCO3","TCO2","TALK"]
real, parameter :: tr_ini(ntracer_nm)= [0.1e-6,0.1e-6,0.1e-6,0.1e-6,0.1e-6,5.0e-6,0.5e-6, &
                    0.1e-6,0.1e-6,2.0e-6,0.1e-6,4.16e-6,1.0e-6,2.0e-3,2.0e-3]
real, public, allocatable :: tr(:,:,:,:,:) ! tracers 
real, public, allocatable :: tr_read(:,:,:,:,:) ! tracers 

integer, allocatable :: id_tr(:)

logical :: initialized=.false., use_this_module=.true.

integer :: isc,jsc,iec,jec
integer :: isd,jsd,ied,jed
integer :: levs

type(restart_file_type) :: restart_nemuro
character (len=32) :: tr_clim_file='INPUT/tr_clim.nc'
character (len=32) :: restart_file='nemuro_restart.nc'

public :: get_num_tracers, initialize_nemuro, diag_out_nemuro, save_restart_nemuro
public :: nemuro_driver

contains

subroutine initialize_nemuro(Time,domain,levs_in,id_lon,id_lat,id_depth)
    type(time_type), intent(in) :: Time
    type(domain2d), intent(in) :: domain
    integer, intent(in) :: id_lon, id_lat, id_depth
    integer, intent(in) :: levs_in 
    
    integer :: unit, ntr, id_restart, dimz(4), i, j, k, nt

    namelist/nemuro_nml/use_this_module

    call mpp_get_compute_domain(domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(domain, isd, ied, jsd, jed)
    levs=levs_in

    unit = open_namelist_file()
    read(unit,nml=nemuro_nml)
    close(unit)

    if (.not.use_this_module) then
        ntracers = 0
        initialized=.true.
        call mpp_error(NOTE, "-------NOT USING NEMURO!-------")
        return
    endif

    ntracers = ntracer_nm
    call mpp_error(NOTE, "-------USING NEMURO!-------")

    allocate(tr(isd:ied,jsd:jed,levs,2,ntracers))
    allocate(tr_read(isc:iec,jsc:jec,levs,12,ntracers))
    allocate(id_tr(ntracers))

    do ntr = 1, ntracers
        id_restart = register_restart_field(restart_nemuro, restart_file, trim(trnms(ntr)), &
                      tr(:,:,:,1,ntr), tr(:,:,:,2,ntr),domain=domain)
    end do
    tr(:,:,:,:,:) = 0.0
 
    do ntr = 1, ntracers
        if (.not. field_exist(tr_clim_file, trim(trnms(ntr)))) then
          call mpp_error(WARNING, 'field '//trim(trnms(ntr))//' not found in '//trim(tr_clim_file))
          tr_read(:,:,:,:,ntr) = tr_ini(ntr)
          cycle
        end if
          
        call field_size(tr_clim_file, trim(trnms(ntr)), dimz)

        do nt = 1, lm-1
            call read_data(tr_clim_file, trim(trnms(ntr)), tr_read(:,:,:,nt,ntr), timelevel=nt, domain=domain)
        enddo
    end do

    if (file_exist('INPUT/'//trim(restart_file))) then
        call restore_state(restart_nemuro)
    else
        call mpp_error(NOTE,'Model starting from initial state')
        do ntr = 1, ntracers
            do i=isc,iec
                do j=jsc,jec
                    if (rkmh(i,j)/=1.) cycle 
                    do k=1,levs
                        tr(i,j,k,1,ntr) = tr_read(i,j,k,1,ntr)
                        tr(i,j,k,2,ntr) = tr_read(i,j,k,1,ntr)
                    enddo
                enddo
            enddo
        end do
    endif 

    do ntr = 1, ntracers
        id_tr(ntr) = register_diag_field('nemuro', trim(trnms(ntr)), (/id_lon,id_lat,id_depth/), &
                init_time=Time, long_name=trim(trnms(ntr)), units='mol/L', &
                missing_value=FILL_VALUE)
    end do

    call save_restart(restart_nemuro,'initial') 

    initialized = .true.
end subroutine initialize_nemuro

subroutine save_restart_nemuro(timestamp)
    character(len=*), optional :: timestamp
    if(.not.initialized) call mpp_error(FATAL,"nemuro not initialized!!!")
    if (.not.use_this_module) return
    call save_restart(restart_nemuro,timestamp)
end subroutine save_restart_nemuro

subroutine diag_out_nemuro(Time)
    type(time_type), intent(in) :: Time
    logical :: used
    integer :: ntr

    if(.not.initialized) call mpp_error(FATAL,"nemuro not initialized!!!")

    if(.not.use_this_module) return

    do ntr = 1, ntracers
        used = send_data(id_tr(ntr),tr(isc:iec,jsc:jec,:,2,ntr),Time)
    end do

    return 
end subroutine diag_out_nemuro


integer function get_num_tracers()
    if(.not.initialized) call mpp_error(FATAL,"nemuro not initialized!!!")
    get_num_tracers = ntracers
end function get_num_tracers


subroutine nemuro_driver(Time, T0, S0, sfx, dt, dz)
  type(time_type), intent(in) :: Time
  real, intent(in), dimension(isd:, jsd:, :) :: T0, S0
  real, intent(in), dimension(isc:, jsc:) :: sfx
  real, intent(in) :: dt, dz(:)
  real, dimension(isc:iec, jsc:jec, 1:levs) :: POM_F_depth, Opal_F_depth, CaCO3_F_depth, &
                                          RnewS, RnewL
   
  real :: rlight(levs)
  integer :: i, j, k

  if(.not.initialized) call mpp_error(FATAL,"nemuro not initialized!!!")

  if(.not.use_this_module) return

  do i = isc, iec
    do j = jsc, jec
      ! calculate rlight
      if (rkmh(i,j)/=1.) cycle 
      call calc_rlight(sfx(i,j), tr(i,j,:,2,1), tr(i,j,:,2,2), &
           dz(:), rlight(:))
      do k = 1, levs
        call nemuro(tr(i,j,k,2,1), tr(i,j,k,2,2), tr(i,j,k,2,3), &
                    tr(i,j,k,2,4), tr(i,j,k,2,5), tr(i,j,k,2,6), &
                    tr(i,j,k,2,7), tr(i,j,k,2,8), tr(i,j,k,2,9), &
                    tr(i,j,k,2,10), tr(i,j,k,2,11), tr(i,j,k,2,12), &
                    tr(i,j,k,2,13), tr(i,j,k,2,14), tr(i,j,k,2,15), &
                    T0(i,j,k), S0(i,j,k), &
                    dt, rlight(k), k, levs, dz(k), POM_F_depth(i,j,k), &
                    Opal_F_depth(i,j,k), CaCO3_F_depth(i,j,k), RnewS(i,j,k), &
                    RnewL(i,j,k))
      end do
    end do
  end do

end subroutine nemuro_driver


subroutine calc_rlight(sfx, PS, PL, dz, rlight)
  real, intent(in) :: sfx, PS(:), PL(:), dz(:)
  real, intent(out) :: rlight(:)
  real :: PAR, rsum, rkappa
 
  integer :: k, kk 
  real, parameter :: Alpha1    =       0.04  ! Light Dissip Coeff. of sea water                       [/m]
  real, parameter :: Alpha2    =       0.04/1e-6 ! Self Shading Coeff.                                    [1/mol m]
       
  rsum = 0.0
        
  PAR  = sfx * 0.45
  do k = 1,levs
    do kk =1,k
      rkappa = Alpha1 + Alpha2 * (PS(kk) + PL(kk))
      rsum = rsum + rkappa * dz(kk)
    end do
    rlight(k) = PAR * exp(-1.0d0 * rsum)
  end do

end subroutine calc_rlight 


subroutine nemuro ( PS, PL, ZS, ZL, ZP, NO3, NH4, POM, DOM, SIOH4, Opal, &
                       Ca, CaCO3, TCO2, TALK, T0, S0, &
                       dt, rlight, k, km, dz, POM_F_depth, &
                       Opal_F_depth, CaCO3_F_depth, RnewS, RnewL)  


  implicit none

  real, intent(inout) :: PS    ! Phytoplankton small 
  real, intent(inout) :: PL    ! Phytoplankton large
  real, intent(inout) :: ZS    ! Zooplankton Small
  real, intent(inout) :: ZL    ! Zooplankton Large
  real, intent(inout) :: ZP    ! Zooplanktton Predatory
  real, intent(inout) :: NO3   ! Nitrate
  real, intent(inout) :: NH4   ! Ammonium
  real, intent(inout) :: POM   ! Particular Organic Matter
  real, intent(inout) :: DOM   ! Dissolved Organic Matter
  real, intent(inout) :: SIOH4 ! Silicic acid
  real, intent(inout) :: Opal  ! Opal
  real, intent(inout) :: Ca    ! Calcium
  real, intent(inout) :: CaCO3 ! Calcium Carbonate
  real, intent(inout) :: TCO2  ! T. Carbondioxide
  real, intent(inout) :: TALK  ! Total Alkalinity

  real, intent(in) :: T0 ! Temperature        
  real, intent(in) :: S0 ! Salinity
  real, intent(in) :: dt ! time-step
  real, intent(in) :: rlight 
  integer, intent(in) :: k, km 
  real, intent(in) :: dz
  real, intent(out) :: POM_F_depth, Opal_F_depth, CaCO3_F_depth, RnewS, &
                       RnewL

!local variables           
  real :: pomin, opalin, caco3in, no3in,sioh4in


!Small Phytoplankton - Tra01
  real :: PS_Photosynthesis
  real :: PS_Respiration
  real :: PS_Extracellular_Excretion
  real :: PS_Mortality
  real :: PS_Grazing_by_ZS
  real :: PS_Grazing_by_ZL

!Large Phytoplankton -Tra02
  real :: PL_Photosynthesis
  real :: PL_Respiration
  real :: PL_Extracellular_Excretion
  real :: PL_Mortality
  real :: PL_Grazing_by_ZL
  real :: PL_Grazing_by_ZP

!Small Zooplankton - Tra03
  real :: ZS_Excretion
  real :: ZS_Egestion
  real :: ZS_Mortality
  real :: ZS_Predating_by_ZL
  real :: ZS_Predating_by_ZP

!Large Zooplankton - Tra04
  real :: ZL_Excretion
  real :: ZL_Egestion
  real :: ZL_Mortality
  real :: ZL_Predating_by_ZP

!Predatory Zooplankton - Tra05       
  real :: ZP_Excretion
  real :: ZP_Egestion
  real :: ZP_Mortality

!Nitrate (NO3) - Tra06               
  real :: Nitrification

!Ammonium - Tra07
  real :: DOM_Remineralization
  real :: POM_Remineralization

!POM - Tra08
  real :: POM_Decomposition_to_DOM
  real :: Sinking_POM 

!Silicate - Tra10
 real :: Opal_Decomposition
 real :: Si_PL_Shell_Formation

!Opal - Tra11
 real :: Si_ZL_Egestion
 real :: Si_ZP_Egestion
 real :: Si_PL_Mortality
 real :: Sinking_Opal

! Ca - Tra12
 real :: CaCO3_Decomposition
 real :: Ca_PS_Shell_Formation
 real :: Ca_ZS_Shell_Formation

! CaCO3 - Tra13
 real :: Ca_ZS_Egestion 
 real :: Ca_ZL_Egestion
 real :: Ca_ZP_Egestion
 real :: Ca_PS_Mortality
 real :: Ca_ZS_Mortality
 real :: Sinking_CaCO3

!  CO2 - Tra14     
 real :: CO2_Air_Sea_Gas_Exchange

!*****************************************************************

  real :: PS_f, PL_f, ZS_f, ZL_f, ZP_f, NO3_f, NH4_f, POM_f, DOM_f, &
                  SIOH4_f, Opal_f, Ca_f, CaCO3_f, TCO2_f, TALK_f 
  real :: ExpPON, ExpOpal, ExcNO3, ExcSiOH4
  real :: ExpCaCO3,ExcT0,ExcS0

  integer          :: kk   
! real(8) Monod, expPsi
  real             :: rtemp1, rtemp2
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       Biological Parameters.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 real, parameter :: D2S       =       86400.0           ! days in seconds                   
 real, parameter :: PSmax     =       200.0             !  Max C:Chla ratio for PS,                              [g:g]             
 real, parameter :: PLmax     =       120.0             ! Max C:Chla ratio for PL,                               [g:g]            
 real, parameter :: KPPS      =       95.0*D2S          ! Photoassimilation coefficient for PS                   [(g:g)day] 
 real, parameter :: KPPL      =       70.0*D2S          ! photoassimilation coefficient for PL                   [(g:g)day] 
 real, parameter :: PSmin     =       30.0              ! min C:Chla ratio for PS,                               [g:g]
 real, parameter :: PLmin     =       15.0              ! min C:Chla ratio for PL,                               [g:g]

 real, parameter :: VmaxS     =       0.4/D2S           ! PS Maximum Photosynthetic rate @0degC                  [/s] 
 real, parameter :: KNO3S     =       1.0e-6            ! PS Half satuation constant for Nitrate                 [molN/l] 
 real, parameter :: KNH4S     =       0.1e-6            ! PS Half satuation constant for Ammonium                [molN/l] 
 real, parameter :: PsiS      =       1.5/1e-6          ! PS Ammonium Inhibition Coefficient                     [l/molN]         
 real, parameter :: KS        =       0.0693            ! PS Temp. Coeff. for Photosynthetic Rate                [/degC]  
 real, parameter :: IoptS     =       104.7             ! PS Optimum Light Intensity                             [W/m2]    
 real, parameter :: MPS0      =       0.0585/1e-6/D2S   ! PS Mortality Rate @0degC                               [l/mol/s]
 real, parameter :: KMPS      =       0.0693            ! PS Temp. Coeff. for Mortality                          [/degC]
 real, parameter :: RPS0      =       0.03/D2S          ! PS Respiration Rate at @0degC                          [/s]
 real, parameter :: KRS       =       0.0519            ! PS Temp. Coeff. for Respiration                        [/degC]
 real, parameter :: GammaS    =       0.135             ! PS Ratio of Extracell. Excret. to Photo.               [(nodim)]
 real, parameter :: VmaxL     =       0.8/D2S           ! PL Maximum Photosynthetic rate @0degC                  [/s] 
 real, parameter :: KNO3L     =       3.0e-6            ! PL Half satuation constant for Nitrate                 [molN/l]
 real, parameter :: KNH4L     =       0.3e-6            ! PL Half satuation constant for Ammonium                [molN/l]
 real, parameter :: KSIL      =       6.0e-6            ! PL Half satuation constant for Silicate                [molSi/l]
 real, parameter :: PsiL      =       1.5/1e-6          ! PL Ammonium Inhibition Coefficient                     [l/molN]
 real, parameter :: KL        =       0.0693            ! PL Temp. Coeff. for Photosynthetic Rate                [/degC]
 real, parameter :: IoptL     =       104.7             ! PL Optimum Light Intensity                             [W/m2]
 real, parameter :: MPL0      =       0.029/1e-6/D2S    ! PL Mortality Rate @0degC                               [l/mol/s]
 real, parameter :: KMPL      =       0.0693            ! PL Temp. Coeff. for Mortality                          [/degC]
 real, parameter :: RPL0      =       0.0300/D2S        ! PL Respiration Rate at @0degC                          [/s]
 real, parameter :: KRL       =       0.0519            ! PL Temp. Coeff. for Respiration                        [/degC]
 real, parameter :: GammaL    =       0.135             ! PL Ratio of Extracell. Excret. to Photo.               [(nodim)]
        
 real, parameter :: GRmaxS    =       0.4/D2S           ! ZS Maximum Rate of Grazing PS @0degC                   [/s]
 real, parameter :: KGS       =       0.0693            ! ZS Temp. Coeff. for Grazing                            [/degC]
 real, parameter :: LamdaS    =       1.4/1e-6          ! ZS Ivlev constant                                      [l/molN]
 real, parameter :: PSZS      =       0.043e-6          ! ZS Threshold Value for Grazing PS                      [molN/l]
 real, parameter :: GammaZS   =       0.7               ! ZS Assimilation Efficiency                             [(nodim)]
 real, parameter :: BettaZS   =       0.3               ! ZS Growth Efficiency                                   [(nodim)]
 real, parameter :: MZS0      =       0.0585/1e-6/D2S   ! ZS Mortality Rate @0degC                               [l/mol/s]
 real, parameter :: KZS       =       0.0693            ! ZS Temp. Coeff. for Mortality                          [/degC]

 real, parameter :: GRmaxLPS  =       0.1/D2S           ! ZL Maximum Rate of Grazing PS @0degC                   [/s]
 real, parameter :: GRmaxLPL  =       0.4/D2S           ! ZL Maximum Rate of Grazing PL @0degC                   [/s]
 real, parameter :: GRmaxLZS  =       0.4/D2S           ! ZL Maximum Rate of Grazing ZS @0degC                   [/s]
 real, parameter :: KGL       =       0.0693            ! ZL Temp. Coeff. for Grazing                            [/degC]
 real, parameter :: LamdaL    =       1.4/1e-6          ! ZL Ivlev constant                                      [l/molN]
 real, parameter :: PSZL      =       0.04e-6           ! ZL Threshold Value for Grazing PS                      [molN/l]
 real, parameter :: PLZL      =       0.04e-6           ! ZL Threshold Value for Grazing PL                      [molN/l]
 real, parameter :: ZSZL      =       0.04e-6           ! ZL Threshold Value for Grazing ZS                      [molN/l]
 real, parameter :: GammaZL   =       0.7               ! ZL Assimilation Efficiency                             [(nodim)]
 real, parameter :: BettaZL   =       0.3               ! ZL Growth Efficiency                                   [(nodim)]
 real, parameter :: MZL0      =       0.0585/1e-6/D2S   ! ZL Mortality Rate @0degC                               [l/mol/s]
 real, parameter :: KZL       =       0.0693            ! ZL Temp. Coeff. for Mortality                          [/degC]
       
 real, parameter :: GRmaxPPL  =       0.4/D2S           ! was 0.1! ZP Maximum Rate of Grazing PS @0degC,         [/s]
                                                        ! "Grazing PL"(in kishi study.) 
                                                        ! it is no t frm the paper. 
                                                        ! from A7station or PAPA station
 real, parameter :: GRmaxPZS  =       0.4/D2S           ! was 0.2 ZP Maximum Rate of Grazing PL @0degC 
                                                        ! "Predating ZS" (in kishi)                              [/s] 
 real, parameter :: GRmaxPZL  =       0.4/D2S           ! was 0.2 ZP Maximum Rate of Grazing ZS @0degC 
                                                        ! "predatitng ZL"(in kishi)                              [/s] 
 real, parameter :: KGP       =       0.0693            ! ZP Temp. Coeff. for Grazing !/predation                [/degC]
 real, parameter :: LamdaP    =       1.5/1e-6          ! ZP Ivlev constant                                      [l/molN]
 real, parameter :: PLZP      =       0.04e-6           ! ZP Threshold Value for Grazing PL                      [molN/l]
                                                        ! "", but it was PSZP, we dont have that why???
 real, parameter :: ZSZP      =       0.04e-6           ! ZP Threshold Value for Grazing ZS                      [molN/l]
 real, parameter :: ZLZP      =       0.04e-6           ! ZP Threshold Value for Grazing ZL                      [molN/l]
 real, parameter :: GammaZP   =       0.7               ! ZP Assimilation Efficiency                             [(nodim)]
 real, parameter :: BettaZP   =       0.3               ! ZP Growth Efficiency                                   [(nodim)]
 real, parameter :: MZP0      =       0.0500/1e-6/D2S   ! ZP Mortality Rate @0degC                               [/s] !l/mol/s unit was error
 real, parameter :: KZP       =       0.0693            ! ZP Temp. Coeff. for Mortality                          [/degC]
 real, parameter :: PusaiPL   =       4.605e6           ! ZP Preference Coeff. for PL                            [l/molN]
 real, parameter :: PusaiZS   =       3.010e6           ! ZP Preference Coeff. for ZS                            [l/molN]

 real, parameter :: NNit0     =       0.03/D2S          ! NH4 Nitrification Rate @0degC                          [/s]
 real, parameter :: KNit      =       0.0693            ! NH4 Temp. coefficient for Nitrification                [/degC]

 real, parameter :: SPOM      =       40.0/D2S          ! Part. Org. Matt. Sinking velocity                      [m/sec]
 real, parameter :: VPA0      =       0.10/D2S          ! PON Decomp. Rate to Ammonium @0degC                    [/s]         
 real, parameter :: KPA       =       0.0693            ! PON Temp. Coeff. for Decomp. to Ammon.                 [/degC]
 real, parameter :: VPD0      =       0.10/D2S          ! PON Decomp. Rate to DON @0degC                         [/s]
 real, parameter :: KPD       =       0.0693            ! PON Temp. Coeff. for Decomp. to DON                    [/degC]
 real, parameter :: VDA0      =       0.2/D2S           ! 0.02/D2S  ! DON Decomp. Rate to Ammonium @0degC        [/s] 
 real, parameter :: KDA       =       0.0693            ! DON Temp. Coeff. for Decomp. to Ammon.                 [/degC]
 real, parameter :: SOpal     =       40.0/D2S          ! (100m/day in org. nemuro)! Sinking Velocity of Opal    [m/sec]
 real, parameter :: VOpal     =       0.10/D2S          ! Opal Decomp. Rate to Silicate @0degC                   [/s]
 real, parameter :: KOpal     =       0.0693            ! Opal Temp. Coeff. for Decomp.to Silicate               [/degC]
 real, parameter :: ExcvNO3   =       1e-6              ! exchange velocity of no3                               [m/s] 
 real, parameter :: ExcvSIOH4 =       1e-6              ! exchange velocity of sioh4                             [m/s] 

 real, parameter :: RSiN      =       2.0               ! Si/N ratio                                             [Nodim]

 real, parameter :: ExcTime   =     3.3d0/(100.0d0*D2S) ! Exch. Coeff. between Sur-Deep                          [/s]
 real, parameter :: TNO3d     =       10.0d-6           ! Nitrate Concentraion in the Deep Layer                 [molN/l]
 real, parameter :: TSiOH4d   =       10.0d-6           ! Silicate Concentraion in the Deep Layer                [molSi/l]
      
 real, parameter :: SCaCO3    =       40.0/D2S          ! Calcium Carbonate Sinking Velocity                     [m/s]   
 real, parameter :: VCaCO3    =       0.050/D2S         ! Decomposition rate of CaCO3 at 0 degC                  [/s]    
 real, parameter :: kCaCO3    =       0.0693            ! Tem. Coefficient for CaCO3 decomposition               [/degC] 

 real, parameter :: RCN       =       6.625             ! Stoichiomery of carbon to Nitrogen                     [Nodim]  
       
 real, parameter :: Rcoco     =       0.1               ! Ratio of Cocolithophoroids to small phytoplankton      [Nodim]  
 real, parameter :: Rfora     =       0.1               ! Ratio of foraminifera to small phytoplankton           [Nodim]  
 real, parameter :: RCco      =       0.5               ! Ratio of Inorganic C to total C in Cocolithophoroids   [Nodim]  
 real, parameter :: RCfo      =       0.5               ! Ratio of Inorganic C to total C in Foraminifera        [Nodim]  


  real :: z_c, grad_POM, grad_Opal, grad_CaCO3
  integer :: kzc, kzcp

       
 
  pomin = pom
  opalin = opal 
  caco3in = caco3 
  no3in = no3
  sioh4in = sioh4
        

! NEMURO EQUATIONS FOR STATE VARIABLES

                                                   
        PS_Photosynthesis = VmaxS *( Monod (NO3, KNO3S)*expPsi(PsiS, NH4) + & 
                  Monod (NH4, KNH4S))                 * & 
                 exp(KS*T0)*(rlight/IoptS)*exp(1.0d0 - rlight/IoptS) * PS
        

        RnewS = Monod(Monod(NO3, KNO3S)* & 
          expPsi(PsiS, NH4) ,Monod(NH4, KNH4S))
       
        rtemp1 = Monod (NO3, KNO3L)*expPsi(PsiL, NH4) + Monod (NH4, KNH4L)
        rtemp2 = Monod(SIOH4,KSIL)

        PL_Photosynthesis = VmaxL * min (rtemp1, rtemp2)* & 
         exp(KL*T0)*(rlight/IoptL)*exp(1.0 - rlight/IoptL) * PL

        
        RnewL = Monod(Monod(NO3, KNO3L)* expPsi(PsiL, NH4) ,Monod(NH4, KNH4L))
         
        
        PS_Respiration = RPS0 * exp(KRS*T0) * PS

        PL_Respiration = RPL0 * exp(KRL*T0) * PL

        PS_Extracellular_Excretion = GammaS * PS_Photosynthesis

        PL_Extracellular_Excretion = GammaL * PL_Photosynthesis
        
        PS_Mortality = MPS0 * exp(KMPS*T0) * PS**2

        PL_Mortality = MPL0 * exp(KMPL*T0) * PL**2

        ZS_Mortality = MZS0 * exp(KZS*T0) * ZS**2

        ZL_Mortality = MZL0 * exp(KZL*T0) * ZL**2

        ZP_Mortality = MZP0 * exp(KZP*T0) * ZP**2


        rtemp1 = 1.0 - exp(LamdaS*(PSZS - PS))
        PS_Grazing_by_ZS = GRmaxS * max (0.0d0, rtemp1 ) * exp (KGS*T0) * ZS

        rtemp1 = 1.0 - exp(LamdaL*(PSZL - PS))
        PS_Grazing_by_ZL = GRmaxLPS * max (0.0d0, rtemp1 ) * exp (KGL*T0) * ZL

        rtemp1 = 1.0 - exp(LamdaL*(PLZL - PL))
        PL_Grazing_by_ZL = GRmaxLPL * max (0.0d0, rtemp1 ) * exp (KGL*T0) * ZL

        rtemp1 = 1.0 - exp(LamdaL*(ZSZL - ZS))
        ZS_Predating_by_ZL = GRmaxLZS * &
        max (0.0d0, rtemp1 ) * &
        exp (KGL*T0) * ZL      

        rtemp1 = 1.0 - exp(LamdaP*(PLZP - PL))
        PL_Grazing_by_ZP = GRmaxPPL * &
        max (0.0d0, rtemp1 ) * &
        expPsi( PusaiPL, ZS+ZL)*exp(KGP*T0)*ZP

        rtemp1 = 1.0 - exp(LamdaP*(ZSZP - ZS))
        ZS_Predating_by_ZP = GRmaxPZS *  &
        max (0.0d0, rtemp1 ) * &
        expPsi( PusaiZS, ZL)*exp(KGP*T0)*ZP

        rtemp1 = 1.0 - exp(LamdaP*(ZLZP - ZL))
        ZL_Predating_by_ZP = GRmaxPZL * &
        max (0.0d0, rtemp1 ) * &
        exp (KGP*T0) * ZP

        ZS_Excretion = (GammaZS - BettaZS)* PS_Grazing_by_ZS
        
        ZL_Excretion = (GammaZL - BettaZL)* (PS_Grazing_by_ZL +  &
                        PL_Grazing_by_ZL + ZS_Predating_by_ZL )

        
        ZP_Excretion = (GammaZP - BettaZP)* (PL_Grazing_by_ZP + &
                       ZS_Predating_by_ZP + ZL_Predating_by_ZP )

        ZS_Egestion =  (1.0d0 - GammaZS) * PS_Grazing_by_ZS

        
        ZL_Egestion = (1.0d0 - GammaZL) * (PS_Grazing_by_ZL + &
                       PL_Grazing_by_ZL + ZS_Predating_by_ZL)

        ZP_Egestion = (1.0d0 - GammaZP) * (PL_Grazing_by_ZP + &
                       ZS_Predating_by_ZP + ZL_Predating_by_ZP)

        
        POM_Remineralization =  VPA0*exp(KPA*T0)*POM

        POM_Decomposition_to_DOM = VPD0*exp(KPD*T0)*POM

        DOM_Remineralization = VDA0*exp(KDA*T0)*DOM

        Opal_Decomposition = VOpal*exp(KOpal*T0)*Opal

        CaCO3_Decomposition =VCaCO3*exp(kCaCO3*T0)*CaCO3
        
        

        Nitrification = NNit0 * exp(KNit*T0)*NH4
        
        
        Si_PL_Shell_Formation = (PL_Photosynthesis - PL_Respiration &
                                - PL_Extracellular_Excretion)*RSiN


        Si_PL_Mortality = PL_Mortality * RSiN
        
        Si_ZL_Egestion = PL_Grazing_by_ZL * RSiN

        Si_ZP_Egestion = PL_Grazing_by_ZP * RSiN


        Ca_PS_Shell_Formation =( PS_Photosynthesis - PS_Respiration &
                                - PS_Extracellular_Excretion) * RCN * Rcoco * RCco  
        
        
        Ca_ZS_Shell_Formation = (PS_Grazing_by_ZS) * BettaZS * RCN &
                                * Rfora * RCfo
         
        

        Ca_PS_Mortality = (PS_Mortality) * RCN * Rcoco * RCco 

        Ca_ZS_Mortality = (ZS_Mortality) * RCN * Rfora * RCfo 
        
        Ca_ZS_Egestion = (PS_Grazing_by_ZS) * RCN *  Rcoco * RCco 

        Ca_ZL_Egestion = (PS_Grazing_by_ZL) * RCN * Rcoco * RCco  &
                         + (ZS_Predating_by_ZL) * RCN * Rfora * RCfo 

        Ca_ZP_Egestion = (ZS_Predating_by_ZP) * RCN * Rfora * RCfo                

        PS_f = PS_Photosynthesis - PS_Respiration - &
               PS_Extracellular_Excretion - PS_Mortality - &
               PS_Grazing_by_ZS  - PS_Grazing_by_ZL


        PL_f = PL_Photosynthesis - PL_Respiration - &
               PL_Extracellular_Excretion - PL_Mortality -  &
               PL_Grazing_by_ZL - PL_Grazing_by_ZP

        ZS_f = PS_Grazing_by_ZS - ZS_Excretion - &
               ZS_Egestion - ZS_Mortality - &
               ZS_Predating_by_ZL - ZS_Predating_by_ZP


        ZL_f = PS_Grazing_by_ZL + PL_Grazing_by_ZL + &
               ZS_Predating_by_ZL - ZL_Excretion - &
               ZL_Egestion - ZL_Mortality -  &
               ZL_Predating_by_ZP
        
        ZP_f = PL_Grazing_by_ZP + ZS_Predating_by_ZP + &
               ZL_Predating_by_ZP - ZP_Excretion - &
               ZP_Egestion - ZP_Mortality

        NO3_f = Nitrification  - &
                (PS_Photosynthesis - PS_Respiration)*RnewS - &
                (PL_Photosynthesis - PL_Respiration)*RnewL

        NH4_f =  (ZS_Excretion + ZL_Excretion + ZP_Excretion) *1.0    + &
                 (DOM_Remineralization + POM_Remineralization)*1.0 -  &
                 Nitrification -   &
                 (PS_Photosynthesis - PS_Respiration)*(1.0d0 - RnewS) - &
                 (PL_Photosynthesis - PL_Respiration)*(1.0d0 - RnewL)
                
                
        POM_f = PS_Mortality + PL_Mortality + ZS_Mortality + &
                ZL_Mortality + ZP_Mortality + ZS_Egestion + &
                ZL_Egestion + ZP_Egestion - POM_Remineralization - & 
                POM_Decomposition_to_DOM

        DOM_f = PS_Extracellular_Excretion + PL_Extracellular_Excretion + &
                POM_Decomposition_to_DOM - DOM_Remineralization 
        
        SIOH4_f = Opal_Decomposition - Si_PL_Shell_Formation


        Opal_f = Si_ZL_Egestion + Si_ZP_Egestion + Si_PL_Mortality - &
                 Opal_Decomposition

        
        Ca_f = CaCO3_Decomposition - Ca_PS_Shell_Formation - &
               Ca_ZS_Shell_Formation        


        CaCO3_f = Ca_ZS_Egestion + Ca_ZL_Egestion + Ca_ZP_Egestion + &
                  Ca_PS_Mortality + Ca_ZS_Mortality - CaCO3_Decomposition 
        

        CO2_Air_Sea_Gas_Exchange = 0.0

        TCO2_f = ((NO3_f + NH4_f) * RCN) +  Ca_f &
                 + CO2_Air_Sea_Gas_Exchange 
       
       
        TALK_f = (2.0 * Ca_f) - NO3_f + NH4_f    
        

        
! DOUBT IN EXPORT and
!        ExpPON   = SPOM  * pom_grad(k)                         
!        ExpOpal  = SOpal * opal_grad(k)                         
!        ExpCaCO3 = SCaCO3  * caco3_grad(k)                     
!*****************************************************************************
        
!        ExcNO3   = no3_grad(k)                                
                                                                !no3_grad(k) is not
                                                                !simply a gradient term it is gradient * wvelocity term.same in
                                                                !sioh4_grad(k)
!        ExcSiOH4 = sioh4_grad(k)                                !(ExcvSIOH4) * sioh4_grad(k)
!        ExcT0    = t0_grad(k)
!        ExcS0    = s0_grad(k)

!***********************old nemuro exchange and
!        ExpPON   =   1.0*SPOM  /100.0 * POM
!         ExpOpal  =   1.0*SOpal / 100.0 * Opal
!         ExpCaCO3 =   1.0*SCaCO3 /100.0 * CaCO3 
!**********************************************************************************************************
!         ExcNO3   =  ExcTime * ( TNO3d   - NO3   )
!         ExcSiOH4 = ExcTime * ( TSiOH4d - SIOH4 )
!!**********************************************************************************************************
!       go to 101 !PROFILE METHOD
         Z_c = 55.0  !meter
         kzc = 51    !11         
         kzcp = kzc + 1
        
         if(k .le. kzc) then
            grad_POM = pomin/dz
            grad_Opal = opalin/dz
            grad_CaCO3 = caco3in/dz
            
            POM_f = POM_f - (grad_POM * SPOM)
            Opal_f = Opal_f - (grad_Opal * SOpal )

!            POM_F_depth = POM_F_depth  + max(0.0, POM) * dz 
!            Opal_F_depth = Opal_F_depth + max(0.0,Opal) * dz

!            POM_F_Z(k)   = 0.0
!            Opal_F_Z(k)  = 0.0
!           CaCO3_F_Z(k) = 0.0
         endif

!         POM_F_depth = max(0.0, POM_F_depth)
!         print *, "Opal_F_depth=", Opal_F_depth 
!         Opal_F_depth = max(0.0, Opal_F_depth)
        
!         if((k.ge.kzc).and.(k.le.km)) then
!         POM_F_Z(k) = POM_F_depth *  ((k*dz)/Z_c)**(-0.9) 
!
!         Opal_F_Z(k) = Opal_F_depth *  ((k*dz)/Z_c)**(-0.9) !a=0.9 ref the papr for
!!        Opal_F_Z(k) =  r_Si * Opal_F_depth * exp(-((k*dz)-(Z_c))/d) 
!!        CaCO3_F_Z(k) = r_Cp * CaCO3_F_depth * exp(-((k*dz)-(Z_c))/d) 
!!         POM_F_Z,Opal_F_Z,CaCO3_F_Z = flux of pom opal and caco3 at a
!         
!        end if
         
!         if(k.ge.kzcp) then
!         POM_f = POM_f + ((POM_F_Z(k-1)-POM_F_Z(k))/dz)
!         Opal_f = Opal_f + ((Opal_F_Z(k-1)-Opal_F_Z(k))/dz)
!!        Opal_f = Opal_f - ((Opal_F_Z(k-1)-Opal_F_Z(k))/dz)/(86400*10)
!!        CaCO3_f = CaCO3_f - ((CaCO3_F_Z(k)-CaCO3_F_Z(k-1))/dz)
!         end if

!101     continue


!**************************************************************
!         NO3_f   = NO3_f - ExcNO3  
!         SIOH4_f = SIOH4_f - ExcSiOH4 
!         POM_f   = POM_f   - ExpPON *0.0 
!         Opal_f  = Opal_f  - ExpOpal *0.0
!         CaCO3_f = CaCO3_f - ExpCaCO3 *0.0
        
        PS = max(0.0, PS + PS_f*dt)
        PL = max(0.0, PL + PL_f*dt)
        ZS = max(0.0, ZS + ZS_f*dt)
        ZL = max(0.0, ZL + ZL_f*dt)
        ZP = max(0.0, ZP + ZP_f*dt)
        NO3 = max (0.0 , NO3 + NO3_f*dt)
        NH4 = max (0.0, NH4 + NH4_f*dt)
        POM = max (0.0, POM + POM_f*dt)
        DOM = max(0.0, DOM + DOM_f*dt)
        SIOH4 = max(0.0, SIOH4 + SIOH4_f*dt)
        Opal = max(0.0, Opal + Opal_f*dt)
        Ca = Ca + Ca_f*dt  
        CaCO3 = CaCO3 + CaCO3_f * dt
        TCO2 = TCO2 + TCO2_f * dt               
        TALK =TALK + TALK_f * dt                
!        T0 = T0 - (ExcT0 * dt)                 
!        S0 = S0 - (ExcS0 * dt)                 

        return
  end subroutine nemuro


  function Monod (a, b)
    real :: a, b, Monod
    Monod = a/(a+b)
  end function
  
  function expPsi (a, b)
    real :: expPsi, a, b
    expPsi = exp(-1.0*a*b)
  end function

end module nemuro_mod 
