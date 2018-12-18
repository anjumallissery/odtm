
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

integer, public :: ntracers=15
character(len=16), allocatable :: trnms(:)
real, public, allocatable :: tr(:,:,:,:,:) ! tracers 
real, public, allocatable :: tr_read(:,:,:,:,:) ! tracers 

integer, allocatable :: id_tr(:)

logical :: initialized=.false., use_this_module=.false.

integer :: isc,jsc,iec,jec
integer :: isd,jsd,ied,jed
integer :: levs

type(restart_file_type) :: restart_nemuro
character (len=32) :: tr_clim_file='INPUT/tr_clim.nc'
character (len=32) :: restart_file='nemuro_restart.nc'

public :: get_num_tracers, initialize_nemuro, diag_out_nemuro, save_restart_nemuro

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

    call mpp_error(NOTE, "-------USING NEMURO!-------")

    allocate(trnms(ntracers))
    trnms(:)=(/"PS","PL","ZS","ZL","ZP","NO3","NH4","POM","DOM","SIOH4","Opal","Ca","CaCO3","TCO2","TALK"/)
    allocate(tr(isd:ied,jsd:jed,levs,2,ntracers))
    allocate(tr_read(isc:iec,jsc:jec,levs,12,ntracers))
    allocate(id_tr(ntracers))

    do ntr = 1, ntracers
        id_restart = register_restart_field(restart_nemuro, restart_file, trim(trnms(ntr)), &
                      tr(:,:,:,1,ntr), tr(:,:,:,2,ntr),domain=domain)
    end do
    tr(:,:,:,:,:) = 0.0
 
    do ntr = 1, ntracers
        if (.not. field_exist(tr_clim_file, trim(trnms(ntr)))) &
            call mpp_error(FATAL, 'field '//trim(trnms(ntr))//' not found in '//trim(tr_clim_file))
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

    do ntr = 1, ntracers
        used = send_data(id_tr(ntr),tr(isc:iec,jsc:jec,:,2,ntr),Time)
    end do

    return 
end subroutine diag_out_nemuro


integer function get_num_tracers()
    if(.not.initialized) call mpp_error(FATAL,"nemuro not initialized!!!")
    get_num_tracers = ntracers
end function get_num_tracers


!        subroutine nemuro_driver()
!
!                call nemuro_3D
!
!        end subroutine nemuro_driver
!

!ccccc CHECK include old3_nemuro.h_ver3 cccccccccccccccccccccccccc
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c      This module calculates the Nemuro biogeochemistry model
!c      as in Yamanaka et al., (2004), J. Oceanography, Vol. 60, 227-241
!c
!c
!c
!c
!c
!c
!c      Model coding start on:  04-Dec-2014
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!subroutine nemuro_3D ( PS,PL,ZS,ZL,ZP,NO3,NH4,POM,DOM,SIOH4,Opal,Ca,CaCO3,TCO2,TALK,&
!       T0,S0,UVEL0,VVEL0, loop, dt, rlight,k,km,pomin,opalin,caco3in,&
!       no3in,sioh4in,&
!       wvelocity,tempin,saltin,uvelin,vvelin,t0_grad, s0_grad,&
!       uvel0_grad, vvel0_grad,no3_grad,POM_F_depth,Opal_F_depth,&
!       CaCO3_F_depth,no3_obss,sioh4_obss,&
!       RnewS,RnewL,NO3_lim,NH4_lim,N_lim,Si_lim,diff_N_and_Si_lim,&
!       no3_relax_obs,sioh4_relax_obs,sfx)  !PSCChla,PLCChla,mlight)
!
!
!  implicit none
!           
!!ccccccccccccccccccccccccccccccccccccccccccccccccc
!!c      Define local variables for each Tracer  c
!!ccccccccccccccccccccccccccccccccccccccccccccccccc
!  real(8) ZS, ZL, ZP
!  real(8) NO3, NH4, POM, DOM, SIOH4, Opal, CaCO3
!  real(8) Ca, TCO2, TALK    
!  real(8) T0, S0, UVEL0, VVEL0
!
!        real(8) ExpPON, ExpOpal, ExcNO3, ExcSiOH4
!        real(8) ExpCaCO3,ExcT0,ExcS0,ExcUVEL0,ExcVVEL0
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Small Phyto Plankton, Tracer = 3        c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) PS_Photosynthesis
!        real(8) PS_Respiration
!        real(8) PS_Extracellular_Excretion
!        real(8) PS_Mortality
!        real(8) PS_Grazing_by_ZS
!        real(8) PS_Grazing_by_ZL
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Large Phyto Plankton, Tracer = 4        c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) PL_Photosynthesis
!        real(8) PL_Respiration
!        real(8) PL_Extracellular_Excretion
!        real(8) PL_Mortality
!        real(8) PL_Grazing_by_ZL
!        real(8) PL_Grazing_by_ZP
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Small Zoo Plankton, Tracer = 5          c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) ZS_Excretion
!        real(8) ZS_Egestion
!        real(8) ZS_Mortality
!        real(8) ZS_Predating_by_ZL
!        real(8) ZS_Predating_by_ZP
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Large Zoo Plankton, Tracer = 6          c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) ZL_Excretion
!        real(8) ZL_Egestion
!        real(8) ZL_Mortality
!        real(8) ZL_Predating_by_ZP
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Predatory Zoo Pankton, Tracer = 7       c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) ZP_Excretion
!        real(8) ZP_Egestion
!        real(8) ZP_Mortality
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Nitrate (NO3), Tracer = 8               c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) Nitrification
!        real(8) RnewS
!        real(8) RnewL
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Amonia (NH4), Tracer = 9                c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) DOM_Remineralization
!        real(8) POM_Remineralization
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Purt. Org. Matt (POM), Tracer = 10      c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) POM_Decomposition_to_DOM
!        real(8) Sinking_POM 
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Diss. Org. Matt. (DOM), Tracer = 11     c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Silicate, (Si(OH)4), Tracer = 12        c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) Opal_Decomposition
!        real(8)Si_PL_Shell_Formation
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!!       Opal, Tracer = 13                       c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) Si_ZL_Egestion
!        real(8) Si_ZP_Egestion
!        real(8) Si_PL_Mortality
!        real(8) Sinking_Opal
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc  
!!        Ca, Tracer=14                          c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) CaCO3_Decomposition
!        real(8) Ca_PS_Shell_Formation
!        real(8) Ca_ZS_Shell_Formation
!
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc 
!!       CaCO3, Tracer=15                        c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) Ca_ZS_Egestion 
!        real(8) Ca_ZL_Egestion
!        real(8) Ca_ZP_Egestion
!        real(8) Ca_PS_Mortality
!        real(8) Ca_ZS_Mortality
!        real(8) Sinking_CaCO3
!
!
!!cccccccccccccccccccccccccccccccccccccccccccccccc !added by ANJU/22-06-2016 -it is not complete.. 
!!      CO2, Tracer=16                           c
!!cccccccccccccccccccccccccccccccccccccccccccccccc
!        real(8) CO2_Air_Sea_Gas_Exchange
!
!
!!ccccccccccccccccccccccccccccccccccccccccccccccccc !added by ANJU/22-06-2016 -it is not complete.. check in the paper and complete it.
!!      Alkalinity, Tracer=17                     c
!!ccccccccccccccccccccccccccccccccccccccccccccccccc
!        
!
!
!
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!       Biological Parameters.
!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        
!        real(8), parameter :: D2S       =       86400.0         ! days in seconds                                !refer murthus paper!!!
!        real(8), parameter :: PSmax     =       200.0           !  Max C:Chla ratio for PS, [g:g]             
!        real(8), parameter :: PLmax     =       120.0           ! Max C:Chla ratio for PL, [g:g]            
!        real(8), parameter :: KPPS      =       95.0*D2S                ! Photoassimilation coefficient for PS [(g:g)day] 
!        real(8), parameter :: KPPL      =       70.0*D2S                ! photoassimilation coefficient for PL [(g:g)day] !day or secnd refer above paper
!        real(8), parameter :: PSmin     =       30.0            ! min C:Chla ratio for PS, [g:g]
!        real(8), parameter :: PLmin     =       15.0            ! min C:Chla ratio for PL, [g:g]
!
!        real(8), parameter :: VmaxS     =       0.3/D2S   !0.4/D2S         ! PS Maximum Photosynthetic rate @0degC   [/s] !delete
!        real(8), parameter :: KNO3S     =       1.0e-6    !1.0e-6               ! PS Half satuation constant for Nitrate  [molN/l] 
!        real(8), parameter :: KNH4S     =       0.1e-6    !0.1e-6               ! PS Half satuation constant for Ammonium [molN/l] !
!        real(8), parameter :: PsiS      =       1.5/1e-6        ! PS Ammonium Inhibition Coefficient      [l/molN]         !
!        real(8), parameter :: KS        =       0.0693          ! PS Temp. Coeff. for Photosynthetic Rate [/degC]  !
!        real(8), parameter :: IoptS     =       104.7           ! PS Optimum Light Intensity              [W/m2]    !
!        real(8), parameter :: MPS0      =       0.0585/1e-6/D2S ! PS Mortality Rate @0degC                [l/mol/s]
!        real(8), parameter :: KMPS      =       0.0693          ! PS Temp. Coeff. for Mortality           [/degC]
!        real(8), parameter :: RPS0      =       0.03/D2S        ! PS Respiration Rate at @0degC           [/s]
!        real(8), parameter :: KRS       =       0.0519          ! PS Temp. Coeff. for Respiration         [/degC]
!        real(8), parameter :: GammaS    =       0.135           ! PS Ratio of Extracell. Excret. to Photo.[(nodim)]
!
!        real(8), parameter :: VmaxL     =       2.0/D2S   !0.8/D2S              ! PL Maximum Photosynthetic rate @0degC   [/s] !delet
!        real(8), parameter :: KNO3L     =       3.0e-6    !3.0e-6               ! PL Half satuation constant for Nitrate  [molN/l]
!        real(8), parameter :: KNH4L     =       0.3e-6    !0.3e-6               ! PL Half satuation constant for Ammonium [molN/l]
!        real(8), parameter :: KSIL      =       6.0e-6    !6.0e-6               ! PL Half satuation constant for Silicate [molSi/l]
!        real(8), parameter :: PsiL      =       1.5/1e-6        ! PL Ammonium Inhibition Coefficient      [l/molN]
!        real(8), parameter :: KL        =       0.0693          ! PL Temp. Coeff. for Photosynthetic Rate [/degC]
!        real(8), parameter :: IoptL     =       104.7           ! PL Optimum Light Intensity              [W/m2]
!        real(8), parameter :: MPL0      =       0.029/1e-6/D2S  ! PL Mortality Rate @0degC                [l/mol/s]
!        real(8), parameter :: KMPL      =       0.0693          ! PL Temp. Coeff. for Mortality           [/degC]
!        real(8), parameter :: RPL0      =       0.0300/D2S      ! PL Respiration Rate at @0degC           [/s]
!        real(8), parameter :: KRL       =       0.0519          ! PL Temp. Coeff. for Respiration         [/degC]
!        real(8), parameter :: GammaL    =       0.135           ! PL Ratio of Extracell. Excret. to Photo.[(nodim)]
!        
!        real(8), parameter :: GRmaxS    =       0.4/D2S         ! ZS Maximum Rate of Grazing PS @0degC    [/s]
!        real(8), parameter :: KGS       =       0.0693          ! ZS Temp. Coeff. for Grazing             [/degC]
!        real(8), parameter :: LamdaS    =       1.4/1e-6        ! ZS Ivlev constant                       [l/molN]
!        real(8), parameter :: PSZS      =       0.043e-6        ! ZS Threshold Value for Grazing PS       [molN/l]
!        real(8), parameter :: GammaZS   =       0.7             ! ZS Assimilation Efficiency              [(nodim)]
!        real(8), parameter :: BettaZS   =       0.3             ! ZS Growth Efficiency                    [(nodim)]
!        real(8), parameter :: MZS0      =       0.0585/1e-6/D2S ! ZS Mortality Rate @0degC                [l/mol/s]
!        real(8), parameter :: KZS       =       0.0693          ! ZS Temp. Coeff. for Mortality           [/degC]
!
!        real(8), parameter :: GRmaxLPS  =       0.1/D2S         ! ZL Maximum Rate of Grazing PS @0degC    [/s]
!        real(8), parameter :: GRmaxLPL  =       0.4/D2S         ! ZL Maximum Rate of Grazing PL @0degC    [/s]
!        real(8), parameter :: GRmaxLZS  =       0.4/D2S         ! ZL Maximum Rate of Grazing ZS @0degC    [/s]
!        real(8), parameter :: KGL       =       0.0693          ! ZL Temp. Coeff. for Grazing             [/degC]
!        real(8), parameter :: LamdaL    =       1.4/1e-6        ! ZL Ivlev constant                       [l/molN]
!        real(8), parameter :: PSZL      =       0.04e-6         ! ZL Threshold Value for Grazing PS       [molN/l]
!        real(8), parameter :: PLZL      =       0.04e-6         ! ZL Threshold Value for Grazing PL       [molN/l]
!        real(8), parameter :: ZSZL      =       0.04e-6         ! ZL Threshold Value for Grazing ZS       [molN/l]
!        real(8), parameter :: GammaZL   =       0.7             ! ZL Assimilation Efficiency              [(nodim)]
!        real(8), parameter :: BettaZL   =       0.3             ! ZL Growth Efficiency                    [(nodim)]
!        real(8), parameter :: MZL0      =       0.0585/1e-6/D2S ! ZL Mortality Rate @0degC                [l/mol/s]
!        real(8), parameter :: KZL       =       0.0693          ! ZL Temp. Coeff. for Mortality           [/degC]
!        
!        real(8), parameter :: GRmaxPPL  =       0.4/D2S !was 0.1        ! ZP Maximum Rate of Grazing PS @0degC, !"Grazing PL"(in kishi study.)[/s] !it is not frm the paper. from A7station or PAPA station
!        real(8), parameter :: GRmaxPZS  =       0.4/D2S !was 0.2 ZP Maximum Rate of Grazing PL @0degC !"Predating ZS" (in kishi) [/s] !""
!        real(8), parameter :: GRmaxPZL  =       0.4/D2S !was 0.2 ZP Maximum Rate of Grazing ZS @0degC !"predatitng ZL"(in kishi) [/s] !""
!        real(8), parameter :: KGP       =       0.0693  ! ZP Temp. Coeff. for Grazing !/predation   [/degC]!""
!        real(8), parameter :: LamdaP    =       1.5/1e-6! ZP Ivlev constant                       [l/molN]!""
!        real(8), parameter :: PLZP      =       0.04e-6 ! ZP Threshold Value for Grazing PL       [molN/l]!"", but it was PSZP, we dont have that why???
!        real(8), parameter :: ZSZP      =       0.04e-6 ! ZP Threshold Value for Grazing ZS       [molN/l]
!        real(8), parameter :: ZLZP      =       0.04e-6 ! ZP Threshold Value for Grazing ZL       [molN/l]
!        real(8), parameter :: GammaZP   =       0.7     ! ZP Assimilation Efficiency              [(nodim)]
!        real(8), parameter :: BettaZP   =       0.3     ! ZP Growth Efficiency                    [(nodim)]
!        real(8), parameter :: MZP0      =       0.0500/1e-6/D2S ! ZP Mortality Rate @0degC                [/s] !l/mol/s unit was error
!        real(8), parameter :: KZP       =       0.0693          ! ZP Temp. Coeff. for Mortality           [/degC]
!        real(8), parameter :: PusaiPL   =       4.605e6         ! ZP Preference Coeff. for PL             [l/molN]
!        real(8), parameter :: PusaiZS   =       3.010e6         ! ZP Preference Coeff. for ZS             [l/molN]
!
!        real(8), parameter :: Alpha1    =       0.04            ! Light Dissip Coeff. of sea water        [/m]
!        real(8), parameter :: Alpha2    =       0.04/1e-6       ! Self Shading Coeff.                     [1/mol m]
!        
!        real(8), parameter :: NNit0     =       0.03/D2S        ! NH4 Nitrification Rate @0degC           [/s]
!        real(8), parameter :: KNit      =       0.0693          ! NH4 Temp. coefficient for Nitrification [/degC]
!
!        real(8), parameter :: SPOM      =       40.0/D2S        ! Part. Org. Matt. Sinking velocity       [m/sec]
!        real(8), parameter :: VPA0      =       0.10/D2S        ! PON Decomp. Rate to Ammonium @0degC     [/s]         !
!        real(8), parameter :: KPA       =       0.0693          ! PON Temp. Coeff. for Decomp. to Ammon.  [/degC]
!        real(8), parameter :: VPD0      =       0.10/D2S        ! PON Decomp. Rate to DON @0degC          [/s]
!        real(8), parameter :: KPD       =       0.0693          ! PON Temp. Coeff. for Decomp. to DON     [/degC]
!        real(8), parameter :: VDA0      =       0.2/D2S      !0.02/D2S  ! DON Decomp. Rate to Ammonium @0degC     [/s] !
!        real(8), parameter :: KDA       =       0.0693          ! DON Temp. Coeff. for Decomp. to Ammon.  [/degC]
!        real(8), parameter :: SOpal     =       40.0/D2S!(100m/day in org. nemuro)      ! Sinking Velocity of Opal                [m/sec]
!        real(8), parameter :: VOpal     =       0.10/D2S        ! Opal Decomp. Rate to Silicate @0degC    [/s]
!        real(8), parameter :: KOpal     =       0.0693          ! Opal Temp. Coeff. for Decomp.to Silicate[/degC]
!        real(8), parameter :: ExcvNO3   =       1e-6            ! exchange velocity of no3                [m/s] !by anju doubt in it
!        real(8), parameter :: ExcvSIOH4 =       1e-6            ! exchange velocity of sioh4              [m/s] !by anju doubt in it
!
!        real(8), parameter :: RSiN      =       2.0             ! Si/N ratio                              [Nodim]
!
!        real(8), parameter :: ExcTime   =       3.3d0 / (100.0d0*D2S) ! Exch. Coeff. between Sur-Deep     [/s]
!        real(8), parameter :: TNO3d     =       10.0d-6         ! Nitrate Concentraion in the Deep Layer  [molN/l]
!        real(8), parameter :: TSiOH4d   =       10.0d-6         ! Silicate Concentraion in the Deep Layer [molSi/l]
!        
!        real(8), parameter :: SCaCO3    =       40.0/D2S           ! Calcium Carbonate Sinking Velocity   [m/s]   added by ANJU/22-06-2016
!        real(8), parameter :: VCaCO3    =       0.050/D2S           ! Decomposition rate of CaCO3 at 0 degC [/s]     ! added by ANJU/22-06-2016
!        real(8), parameter :: kCaCO3    =       0.0693          !Tem. Coefficient for CaCO3 decomposition [/degC]  ! added by ANJU/22-06-2016
!
!        real(8), parameter :: RCN       =       6.625           !Stoichiomery of carbon to Nitrogen       [Nodim]  ! added by ANJU/22-06-2016
!        
!        real(8), parameter :: Rcoco     =       0.1  !Ratio of Cocolithophoroids to small phytoplankton   [Nodim]  ! added by ANJU/22-06-2016
!        real(8), parameter :: Rfora     =       0.1  !Ratio of foraminifera to small phytoplankton        [Nodim]  ! added by ANJU/22-06-2016
!        real(8), parameter :: RCco      =       0.5  !Ratio of Inorganic C to total C in Cocolithophoroids[Nodim]  ! added by ANJU/22-06-2016
!        real(8), parameter :: RCfo      =       0.5  !Ratio of Inorganic C to total C in Foraminifera     [Nodim]  ! added by ANJU/22-06-2016
!       
!        real, intent(inout) :: PS   ! Phytoplankton small 
!        real, intent(inout) :: PL   ! Phytoplankton large
! 
!        integer :: kk   
!        real dz !real 8 or 4 
!!       real(8) Monod, expPsi
!        real(8) rkappa, rlight, rsum, sfx, day_night
!        real(8) rtemp1, rtemp2
!        real w1,w2,no31,no32,sioh41,sioh42,t01,t02,s01,s02,uvel01,uvel02
!        real vvel01,vvel02 
!        real pom_grad(51),opal_grad(51),caco3_grad(51),no3_grad(51)
!        real sioh4_grad(51)
!        real t0_grad(51),s0_grad(51),uvel0_grad(51),vvel0_grad(51) 
!        real pomin(51),opalin(51),caco3in(51),no3in(51),sioh4in(51) 
!        real pomin_sm(51),opalin_sm(51),caco3in_sm(51),no3in_sm(51)
!        real sioh4in_sm(51) 
!        real wvelocity(51)
!        real tempin_sm(51),saltin_sm(51),uvelin_sm(51),vvelin_sm(51)
!        real tempin(51),saltin(51),uvelin(51),vvelin(51) 
!        real POM_F_depth,Opal_F_depth,CaCO3_F_depth,a
!        real POM_F_Z(51),Opal_F_Z(51),CaCO3_F_Z(51),R
!        real r_CP,d
!        real w_diff,wb
!        real delta,grad_POM, grad_Opal, grad_CaCO3 
!        real POM_ave
!        real no3_obss,sioh4_obss
!        real relaxation_NO3,relaxation_SIOH4      
!        real NO3_lim,NH4_lim,N_lim,Si_lim,diff_N_and_Si_lim
!        real no3_relax_obs,sioh4_relax_obs
!!        real(8) NLLGR,ter1,PSCCh,PSC,PSCChla(51)
!!        real(8) PLCCh,PLC,PLCChla(51),mlight
!        integer loop, dt
!        integer k,km,kmm,Z_c, kzc, kzcp                                 !depth upto which flux calculating
!        
!        real(8) PS_f, PL_f, ZS_f, ZL_f, ZP_f, NO3_f, NH4_f, POM_f, DOM_f,
!     &             SIOH4_f, Opal_f, Ca_f, CaCO3_f, TCO2_f, TALK_f 
!        
!        kmm = km -3 
!        dz  = 5.0
!        a   = 0.9   
!       
!             
!                                     !constant used in POM flux calculation refr more....
!        R= 0.08                           !0.07 also possible ref paper
!        r_CP = 117.0                      ! 106 is also possible. check ref. paper..
!!        r_Si = 2.0 !From NEmuro          
!        d = 3500.0   
!                                                   ! meter
!        PS_Photosynthesis = VmaxS *( Monod (NO3, KNO3S)*expPsi(PsiS, NH4) + & 
!                  Monod (NH4, KNH4S))                 * & 
!                 exp(KS*T0)*(rlight/IoptS)*exp(1.0d0 - rlight/IoptS) * PS
!        
!
!        RnewS = Monod(Monod(NO3, KNO3S)* & 
!          expPsi(PsiS, NH4) ,Monod(NH4, KNH4S))
!       
!        rtemp1 = Monod (NO3, KNO3L)*expPsi(PsiL, NH4) + Monod (NH4, KNH4L)
!        rtemp2 = Monod(SIOH4,KSIL)
!
!        PL_Photosynthesis = VmaxL * min (rtemp1, rtemp2)* & 
!         exp(KL*T0)*(rlight/IoptL)*exp(1.0 - rlight/IoptL) * PL
!
!
!!!!DELETE THIS****************NO SILICA LIMITATION CASE
!!         PL_Photosynthesis = VmaxL *
!!     &  rtemp1 * 
!!     &  exp(KL*T0)*(rlight/IoptL)*exp(1.0 - rlight/IoptL) * PL
!!*****************************
!
!        
!!************************************************************************************************
!                ! depth profile of C:Chla ratio for PS
!!        NLLGR = VmaxS* exp(KS*T0) * ( Monod (NO3, KNO3S)*expPsi(PsiS,
!!     &  NH4) + Monod (NH4, KNH4S))
!               
!              
!!        PSC = PSmax - (KPPS * NLLGR)
!        
!!        ter1 = (log(sfx)-log(mlight))/4.605
!!        PSCCh = (PSC - (PSC - PSmin))* ter1
!        !PSCCh = PSC - ((PSC - PSmin)* ter1)
!
!
!!        PSCChla(k) = PSCCh
!        
!        
!               !depth profile of C:Chla ratio for PL
!        
!!        NLLGR = VmaxL * exp(KL*T0) * min (rtemp1,rtemp2)
!!        PLC  =PLmax - (KPPL * NLLGR)
!!        PLCCh = (PLC -(PLC - PLmin))*ter1
!!        !PLCCh = PLC - ((PLC - PLmin)* ter1)
!!        PLCChla(k) = PLCCh
!        
!        
!!********************************************************************************************************                        
!
!        NO3_lim =  Monod (NO3, KNO3L)*expPsi(PsiL, NH4)         !nitrate lim.
!        NH4_lim= Monod (NH4, KNH4L)                             !ammonium limitation
!        N_lim= rtemp1                                           !nitrogen limitation
!                
!
!        NO3_lim =  Monod (NO3, KNO3L)*expPsi(PsiL, NH4)         !nitrate lim.
!        NH4_lim= Monod (NH4, KNH4L)                             !ammonium limitation
!        N_lim= rtemp1                                           !nitrogen limitation
!        Si_lim = rtemp2                                         !silicate limitation
!        
!
!        NO3_lim =  Monod (NO3, KNO3L)*expPsi(PsiL, NH4)!nitrate lim.
!        NH4_lim= Monod (NH4, KNH4L)!ammonium limitation
!        N_lim= rtemp1 ! nitrogen limitation
!        Si_lim = rtemp2 !silicate limitation
!        diff_N_and_Si_lim = rtemp1-rtemp2 ! difference b/w N and Si 
!        !limitation
!        
!        
!!        if((loop.le.200).and.(k.eq.1))then
!!        print*, loop,k,NH4_lim
!!        end if
!        
!        RnewL = Monod(Monod(NO3, KNO3L)* expPsi(PsiL, NH4) ,Monod(NH4, KNH4L))
!         
!        
!        PS_Respiration = RPS0 * exp(KRS*T0) * PS
!
!        PL_Respiration = RPL0 * exp(KRL*T0) * PL
!
!        PS_Extracellular_Excretion = GammaS * PS_Photosynthesis
!
!        PL_Extracellular_Excretion = GammaL * PL_Photosynthesis
!        
!        PS_Mortality = MPS0 * exp(KMPS*T0) * PS**2
!
!        PL_Mortality = MPL0 * exp(KMPL*T0) * PL**2
!
!        ZS_Mortality = MZS0 * exp(KZS*T0) * ZS**2
!
!        ZL_Mortality = MZL0 * exp(KZL*T0) * ZL**2
!
!        ZP_Mortality = MZP0 * exp(KZP*T0) * ZP**2
!
!
!        rtemp1 = 1.0 - exp(LamdaS*(PSZS - PS))
!        PS_Grazing_by_ZS = GRmaxS * max (0.0d0, rtemp1 ) * exp (KGS*T0) * ZS
!
!        rtemp1 = 1.0 - exp(LamdaL*(PSZL - PS))
!        PS_Grazing_by_ZL = GRmaxLPS * max (0.0d0, rtemp1 ) * exp (KGL*T0) * ZL
!
!        rtemp1 = 1.0 - exp(LamdaL*(PLZL - PL))
!        PL_Grazing_by_ZL = GRmaxLPL * max (0.0d0, rtemp1 ) * exp (KGL*T0) * ZL
!
!        rtemp1 = 1.0 - exp(LamdaL*(ZSZL - ZS))
!        ZS_Predating_by_ZL = GRmaxLZS *
!     &       max (0.0d0, rtemp1 ) *
!     &       exp (KGL*T0) * ZL
!
!        rtemp1 = 1.0 - exp(LamdaP*(PLZP - PL))
!        PL_Grazing_by_ZP = GRmaxPPL *
!     &       max (0.0d0, rtemp1 ) *
!     &       expPsi( PusaiPL, ZS+ZL)*exp(KGP*T0)*ZP
!
!        rtemp1 = 1.0 - exp(LamdaP*(ZSZP - ZS))
!        ZS_Predating_by_ZP = GRmaxPZS *
!     &       max (0.0d0, rtemp1 ) *
!     &       expPsi( PusaiZS, ZL)*exp(KGP*T0)*ZP
!
!        rtemp1 = 1.0 - exp(LamdaP*(ZLZP - ZL))
!        ZL_Predating_by_ZP = GRmaxPZL *
!     &       max (0.0d0, rtemp1 ) *
!     &       exp (KGP*T0) * ZP
!
!        ZS_Excretion = (GammaZS - BettaZS)* PS_Grazing_by_ZS
!        
!        ZL_Excretion = (GammaZL - BettaZL)* (PS_Grazing_by_ZL +
!     &                 PL_Grazing_by_ZL + ZS_Predating_by_ZL )
!
!        
!        ZP_Excretion = (GammaZP - BettaZP)* (PL_Grazing_by_ZP +
!     &                 ZS_Predating_by_ZP + ZL_Predating_by_ZP )
!
!        ZS_Egestion = (1.0d0 - GammaZS) * PS_Grazing_by_ZS
!
!        
!        ZL_Egestion = (1.0d0 - GammaZL) * (PS_Grazing_by_ZL +
!     &                 PL_Grazing_by_ZL + ZS_Predating_by_ZL)
!
!        ZP_Egestion = (1.0d0 - GammaZP) * (PL_Grazing_by_ZP +
!     &                 ZS_Predating_by_ZP + ZL_Predating_by_ZP)
!
!        
!        POM_Remineralization =  VPA0*exp(KPA*T0)*POM
!
!        POM_Decomposition_to_DOM = VPD0*exp(KPD*T0)*POM
!
!        DOM_Remineralization = VDA0*exp(KDA*T0)*DOM
!
!        Opal_Decomposition = VOpal*exp(KOpal*T0)*Opal
!
!        CaCO3_Decomposition =VCaCO3*exp(kCaCO3*T0)*CaCO3
!        
!        
!
!        Nitrification = NNit0 * exp(KNit*T0)*NH4
!        
!        
!        Si_PL_Shell_Formation = (PL_Photosynthesis - PL_Respiration
!     &            - PL_Extracellular_Excretion)*RSiN
!
!
!        Si_PL_Mortality = PL_Mortality * RSiN
!        
!        Si_ZL_Egestion = PL_Grazing_by_ZL * RSiN
!
!        Si_ZP_Egestion = PL_Grazing_by_ZP * RSiN
!
!
!        Ca_PS_Shell_Formation =( PS_Photosynthesis - PS_Respiration
!     &       - PS_Extracellular_Excretion) * RCN * Rcoco * RCco  
!        
!        
!        Ca_ZS_Shell_Formation = (PS_Grazing_by_ZS) * BettaZS * RCN 
!     &       * Rfora * RCfo
!         
!        
!
!        Ca_PS_Mortality = (PS_Mortality) * RCN * Rcoco * RCco 
!
!        Ca_ZS_Mortality = (ZS_Mortality) * RCN * Rfora * RCfo 
!        
!        Ca_ZS_Egestion = (PS_Grazing_by_ZS) * RCN *  Rcoco * RCco 
!
!        Ca_ZL_Egestion = (PS_Grazing_by_ZL) * RCN * Rcoco * RCco 
!     &                 + (ZS_Predating_by_ZL) * RCN * Rfora * RCfo 
!
!        Ca_ZP_Egestion = (ZS_Predating_by_ZP) * RCN * Rfora * RCfo                
!        PS_f = PS_Photosynthesis - PS_Respiration - 
!     &         PS_Extracellular_Excretion - PS_Mortality -
!     &         PS_Grazing_by_ZS  - PS_Grazing_by_ZL
!
!
!        PL_f = PL_Photosynthesis - PL_Respiration -
!     &         PL_Extracellular_Excretion - PL_Mortality - 
!     &         PL_Grazing_by_ZL - PL_Grazing_by_ZP
!
!        ZS_f = PS_Grazing_by_ZS - ZS_Excretion -
!     &         ZS_Egestion - ZS_Mortality -
!     &         ZS_Predating_by_ZL - ZS_Predating_by_ZP
!
!
!        ZL_f = PS_Grazing_by_ZL + PL_Grazing_by_ZL +
!     &         ZS_Predating_by_ZL - ZL_Excretion -
!     &         ZL_Egestion - ZL_Mortality -
!     &         ZL_Predating_by_ZP
!        
!        ZP_f = PL_Grazing_by_ZP + ZS_Predating_by_ZP + 
!     &         ZL_Predating_by_ZP - ZP_Excretion -
!     &         ZP_Egestion - ZP_Mortality
!
!        NO3_f = Nitrification  - 
!     &          (PS_Photosynthesis - PS_Respiration)*RnewS -
!     &          (PL_Photosynthesis - PL_Respiration)*RnewL
!        NH4_f =  (ZS_Excretion + ZL_Excretion + ZP_Excretion) *1.0    +
!     &           (DOM_Remineralization + POM_Remineralization)*1.0 - 
!     &  Nitrification - 
!     & (PS_Photosynthesis - PS_Respiration)*(1.0d0 - RnewS) -
!     & (PL_Photosynthesis - PL_Respiration)*(1.0d0 - RnewL)
!                
!                
!        POM_f = PS_Mortality + PL_Mortality + ZS_Mortality +
!     &          ZL_Mortality + ZP_Mortality + ZS_Egestion +
!     &          ZL_Egestion + ZP_Egestion - POM_Remineralization - 
!     &          POM_Decomposition_to_DOM
!
!        DOM_f = PS_Extracellular_Excretion + PL_Extracellular_Excretion +
!     &          POM_Decomposition_to_DOM - DOM_Remineralization 
!        
!        SIOH4_f = Opal_Decomposition - Si_PL_Shell_Formation
!
!
!        Opal_f = Si_ZL_Egestion + Si_ZP_Egestion + Si_PL_Mortality -
!     &           Opal_Decomposition
!
!
!
!        
!        
!        Ca_f = CaCO3_Decomposition - Ca_PS_Shell_Formation - 
!     &  Ca_ZS_Shell_Formation        
!
!
!        CaCO3_f = Ca_ZS_Egestion + Ca_ZL_Egestion + Ca_ZP_Egestion +
!     &  Ca_PS_Mortality + Ca_ZS_Mortality - CaCO3_Decomposition 
!        
!
!         CO2_Air_Sea_Gas_Exchange = 0.0
!        TCO2_f = ((NO3_f + NH4_f) * RCN) +  Ca_f
!     &          + CO2_Air_Sea_Gas_Exchange 
!       
!       
!        TALK_f = (2.0 * Ca_f) - NO3_f + NH4_f    
!        
!
!
!        
!        !EXCHANGE AND EXPORT TERMS ie vertical advection.
!        
!        
!        !        shapiro filter to get rid of numerical instability due to the
!        !        advection problem
!!        goto 105
!        do kk = k-1, k+1  
!
!        if (mod(loop,100).eq.0)then 
!
!        if ( k.ge.4 .and. k.le.kmm ) then                               ! k gt 4 or 3 dbt????? 4 only we can
!                                                                        !choose
!        pomin_sm(kk) =(1.0/16.0) * ((-1*pomin(kk-2))+(4*pomin(kk-1))
!     &          +(10*pomin(kk))+(4*pomin(kk+1))+(-1*pomin(kk+2)))
!        opalin_sm(kk)=(1.0/16.0) * ((-1*opalin(kk-2))+(4*opalin(kk-1))
!     &        +(10*opalin(kk))+(4*opalin(kk+1))+(-1*opalin(kk+2)))
!
!        caco3in_sm(kk)=(1.0/16.0)*((-1*caco3in(kk-2))+(4*caco3in(kk-1))
!     &     +(10*caco3in(kk))+(4*caco3in(kk+1))+(-1*caco3in(kk+2)))
!
!        no3in_sm(kk) =(1.0/16.0) *((-1*no3in(kk-2))+(4*no3in(kk-1))
!     &     +(10*no3in(kk))+(4*no3in(kk+1))+(-1*no3in(kk+2)))
!
!        sioh4in_sm(kk) =(1.0/16.0)*((-1*sioh4in(kk-2))+(4*sioh4in(kk-1))
!     &    +(10*sioh4in(kk))+(4*sioh4in(kk+1))+(-1*sioh4in(kk+2))) 
!        
!        tempin_sm(kk)=(1.0/16.0) * ((-1*tempin(kk-2))+(4*tempin(kk-1))
!     &        +(10*tempin(kk))+(4*tempin(kk+1))+(-1*tempin(kk+2)))
!        
!        saltin_sm(kk)=(1.0/16.0) * ((-1*saltin(kk-2))+(4*saltin(kk-1))
!     &        +(10*saltin(kk))+(4*saltin(kk+1))+(-1*saltin(kk+2)))
!        
!        uvelin_sm(kk)=(1.0/16.0) * ((-1*uvelin(kk-2))+(4*uvelin(kk-1))
!     &        +(10*uvelin(kk))+(4*uvelin(kk+1))+(-1*uvelin(kk+2)))
!        
!        vvelin_sm(kk)=(1.0/16.0) * ((-1*vvelin(kk-2))+(4*vvelin(kk-1))
!     &        +(10*vvelin(kk))+(4*vvelin(kk+1))+(-1*vvelin(kk+2)))
!
!                
!!        pomin_sm(kk) = pomin after smoothing by shapiro filter , equations got
!!        it from vinu sr. ask for the reference.
!
!        pomin(kk)=pomin_sm(kk)
!        opalin(kk)=opalin_sm(kk)
!        caco3in(kk)=caco3in_sm(kk)
!        no3in(kk) =no3in_sm(kk) 
!        sioh4in(kk)= sioh4in_sm(kk)
!        tempin(kk) = tempin_sm(kk) 
!        saltin(kk) = saltin_sm(kk) 
!        uvelin(kk) = uvelin_sm(kk) 
!        vvelin(kk) = vvelin_sm(kk) 
!
!
!        end if
!        end if
!        end do
!                
!        if (mod(loop,100).eq.0)then 
!        if ( k.ge.4 .and. k.le.kmm ) then 
!
!        POM = pomin_sm(k)
!        Opal= opalin_sm(k)
!        CaCO3 = caco3in_sm(k)
!        NO3 = no3in_sm(k)
!        SIOH4 = sioh4in_sm(k)       
!        T0    = tempin_sm(k) 
!        S0    = saltin_sm(k)
!        UVEL0 = uvelin_sm(k)
!        VVEL0 = vvelin_sm(k)
! 
!        end if
!        end if
!!105    continue
!        
!        if (k.eq.1) then
!
!       ! w_diff = wvelocity(k)-wvelocity(k+1)
!       ! wb = wvelocity(k)+(w_diff)              !wb -virtual layer
!                                              !velocity
!       ! w1 = ((wvelocity(k))+(wb))/2
!        
!                !***************GRAD CALCULATION ****************
!
!        w1 = (wvelocity(k)+wvelocity(k))/2
!        w2 = (wvelocity(k)+wvelocity(k+1))/2
!        
!       !w1 = min(1e-6,w1)
!       !w2 = min(1e-6,w1)
!        !w1 = 1e-5
!        !w2 = 1e-5
!
!        no31 = (no3in(k)+no3in(k))/2
!        no32 = (no3in(k)+no3in(k+1))/2
!        sioh41 = (sioh4in(k)+sioh4in(k))/2
!        sioh42 = (sioh4in(k)+sioh4in(k+1))/2
!        t01    = (tempin(k)+tempin(k))/2 
!        t02    = (tempin(k)+tempin(k+1))/2
!        s01    = (saltin(k)+saltin(k))/2                
!        s02    = (saltin(k)+saltin(k+1))/2              
!!        uvel01 = (uvelin(k)+uvelin(k))/2               
!!        uvel02 = (uvelin(k)+uvelin(k+1))/2             
!!        vvel01 = (vvelin(k)+vvelin(k))/2               
!!        vvel02 = (vvelin(k)+vvelin(k+1))/2           
! 
!        pom_grad(k) = (((pomin(k)+pomin(k))/2) 
!     &       -((pomin(k)+pomin(k+1))/2)) /(dz)
!                                                !dz direction is positive
!        
!        opal_grad(k)=(((opalin(k)+opalin(k))/2)
!     &       -((opalin(k)+opalin(k+1))/2))/(dz)
!        
!        caco3_grad(k)=(((caco3in(k)+caco3in(k))/2)
!     &     -((caco3in(k)+caco3in(k+1))/2))/(dz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!        no3_grad(k) = ((w1*no31)- (w2*no32))/(dz)
!!        sioh4_grad(k) = ((w1*sioh41) - (w2*sioh42))/(dz)
!!        t0_grad(k)    = ((w1*t01) - (w2*t02))/(dz)
!!        s0_grad(k)    = ((w1*s01) - (w2*s02))/(dz)
!!        uvel0_grad(k) = ((w1*uvel01) - (w2*uvel02))/(dz)
!!        vvel0_grad(k) = ((w1*vvel01) - (w2*vvel02))/(dz) 
!        
!        no3_grad(k) =  ((w1+w2)/2) * (no31 - no32)/(dz)
!        sioh4_grad(k) = ((w1+w2)/2) * (sioh41 - sioh42)/(dz)
!        t0_grad(k)    = ((w1+w2)/2) * (t01 - t02)/(dz)         
!        s0_grad(k)    = ((w1+w2)/2) * (s01 - s02)/(dz) 
!!        uvel0_grad(k) = ((w1+w2)/2) * (uvel01 - uvel02)/(dz)
!!        vvel0_grad(k) = ((w1+w2)/2) * (vvel01 - vvel02)/(dz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                
!
!        else if ((k.ge.2) .and. (k.lt.km)) then
!        w1 = (wvelocity(k)+wvelocity(k-1))/2
!        w2 = (wvelocity(k)+wvelocity(k+1))/2
!
!        !w1 = min(1e-6,w1)
!        !w2 = min(1e-6,w1)
!
!        !w1=1e-5 !delete
!        !w2=1e-5 !delete
!
!        no31 = (no3in(k)+no3in(k))/2
!        no31 = (no3in(k)+no3in(k-1))/2
!        no32 = (no3in(k)+no3in(k+1))/2
!        sioh41 = (sioh4in(k)+sioh4in(k-1))/2
!        sioh42 = (sioh4in(k)+sioh4in(k+1))/2
!        t01    = (tempin(k)+tempin(k-1))/2      
!        t02    = (tempin(k)+tempin(k+1))/2      
!        s01    = (saltin(k)+saltin(k-1))/2     
!        s02    = (saltin(k)+saltin(k+1))/2     
!!        uvel01 = (uvelin(k)+uvelin(k-1))/2    
!!        uvel02 = (uvelin(k)+uvelin(k+1))/2    
!!        vvel01 = (vvelin(k)+vvelin(k-1))/2    
!!        vvel02 = (vvelin(k)+vvelin(k+1))/2     
!
!
!        pom_grad(k) = (((pomin(k) + pomin(k-1))/2) 
!     &      - ((pomin(k)+pomin(k+1))/2))/(dz)
!
!        opal_grad(k) =(((opalin(k) + opalin(k-1))/2)
!     &     - ((opalin(k)+opalin(k+1))/2))/(dz)
!!
!        caco3_grad(k) = (((caco3in(k) + caco3in(k-1))/2)
!     &        -((caco3in(k)+caco3in(k+1))/2))/(dz)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!        no3_grad(k) = ((w1*no31)- (w2*no32))/(dz)
!!        sioh4_grad(k) = ((w1*sioh41) - (w2*sioh42))/(dz)
!!        t0_grad(k)    = ((w1*t01) - (w2*t02))/(dz) 
!!        s0_grad(k)    = ((w1*s01) - (w2*s02))/(dz)
!!        uvel0_grad(k) = ((w1*uvel01) - (w2*uvel02))/(dz)
!!        vvel0_grad(k) = ((w1*vvel01) - (w2*vvel02))/(dz) 
!
!        no3_grad(k) =  ((w1+w2)/2) * (no31 - no32)/(dz)
!        sioh4_grad(k) = ((w1+w2)/2) * (sioh41 - sioh42)/(dz)
!        t0_grad(k)    = ((w1+w2)/2)  * (t01 - t02)/(dz)         
!        s0_grad(k)    = ((w1+w2)/2) * (s01 - s02)/(dz)
!!        uvel0_grad(k) = ((w1+w2)/2) * (uvel01 - uvel02)/(dz)
!!        vvel0_grad(k) = ((w1+w2)/2) * (vvel01 - vvel02)/(dz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        else
!        w1 = (wvelocity(k-1)+wvelocity(k))/2
!        w2 = (wvelocity(k)+wvelocity(k))/2
!        
!        !w1 = min(1e-6,w1)
!        !w2 = min(1e-6,w1)
!
!        !w1 = 1e-5
!        !w2 = 1e-5
!
!        no31 = (no3in(k)+no3in(k))/2
!        no31 = (no3in(k-1)+no3in(k))/2
!        no32 = (no3in(k)+no3in(k))/2
!        sioh41 = (sioh4in(k-1)+sioh4in(k))/2
!        sioh42 = (sioh4in(k)+sioh4in(k))/2
!        t01   = (tempin(k-1)+tempin(k))/2               
!        t02   = (tempin(k)+tempin(k))/2                
!        s01   = (saltin(k-1)+saltin(k))/2              
!        s02   = (saltin(k)+saltin(k))/2                
!!        uvel01= (uvelin(k-1)+uvelin(k))/2              
!!        uvel02= (uvelin(k)+uvelin(k))/2           
!!        vvel01= (vvelin(k-1)+vvelin(k))/2              
!!        vvel02= (vvelin(k)+vvelin(k))/2              
!
!
!        pom_grad(k) = (((pomin(k-1)+pomin(k))/2)
!     &      -((pomin(k)+ pomin(k))/2))/(dz)
!!
!        opal_grad(k) = (((opalin(k-1)+opalin(k))/2)
!     &         - ((opalin(k)+opalin(k))/2))/(dz)
!!
!        caco3_grad(k) = (((caco3in(k-1)+caco3in(k))/2)
!     &    - ((caco3in(k) +caco3in(k))/2) )/(dz)
!
!!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!!        no3_grad(k) = ((w1*no31)- (w2*no32))/(dz)
!!        sioh4_grad(k) = ((w1*sioh41) - (w2*sioh42))/(dz)
!!        t0_grad(k)    = ((w1*t01) - (w2*t02))/(dz) 
!!        s0_grad(k)    = ((w1*s01) - (w2*s02))/(dz)
!!        uvel0_grad(k) = ((w1*uvel01) - (w2*uvel02))/(dz)
!!        vvel0_grad(k) = ((w1*vvel01) - (w2*vvel02))/(dz) 
!
!        no3_grad(k) =  ((w1+w2)/2) * (no31 - no32)/(dz)
!        sioh4_grad(k) = ((w1+w2)/2) * (sioh41 - sioh42)/(dz)
!        t0_grad(k)    = ((w1+w2)/2)  * (t01 - t02)/(dz)         
!        s0_grad(k)    = ((w1+w2)/2) * (s01 - s02)/(dz)
!!        uvel0_grad(k) = ((w1+w2)/2) * (uvel01 - uvel02)/(dz)
!!        vvel0_grad(k) = ((w1+w2)/2) * (vvel01 - vvel02)/(dz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        end if
!!**************instead of this a profile has to give for export****    
!        ExpPON   = SPOM  * pom_grad(k)                          !DOUBT, y ekman_vel is not here
!        ExpOpal  = SOpal * opal_grad(k)                         !DOUBT, " "
!        ExpCaCO3 = SCaCO3  * caco3_grad(k)                      !DOUBT, " "
!!*****************************************************************************
!        
!        ExcNO3   = no3_grad(k)                                  !(ExcvNO3) * no3_grad(k)!instead of
!                                                                !Exchange velocity ekman velocity is given here. no3_grad(k) is not
!                                                                !simply a gradient term it is gradient * wvelocity term.same in
!                                                                !sioh4_grad(k)
!        ExcSiOH4 = sioh4_grad(k)!(ExcvSIOH4) * sioh4_grad(k)
!        ExcT0    = t0_grad(k)
!        ExcS0    = s0_grad(k)
!!        ExcUVEL0 = uvel0_grad(k)
!!        ExcVVEL0 = vvel0_grad(k)
!
!        
!!        ExpPON   =   1.0*SPOM  /100.0 * POM
!!         ExpOpal  =   1.0*SOpal / 100.0 * Opal
!!         ExpCaCO3 =   1.0*SCaCO3 /100.0 * CaCO3 !added by anju. check sigh
!!**********************************************************************************************************
!       !Giving exchange in the model 
!!         ExcNO3   =  ExcTime * ( TNO3d   - NO3   )
!!         ExcSiOH4 = ExcTime * ( TSiOH4d - SIOH4 )
!!!**********************************************************************************************************
!!       go to 101 !PROFILE METHOD
!         Z_c = 55.0  !meter
!         kzc = 51    !11         
!         kzcp = kzc + 1
!        
!         if(k .le. kzc) then
!          !POM_ave = (pomin(k) + pomin(k+1)) /2 !NOW
!         grad_POM = pomin(k)/dz!(pomin(k) - pomin(k+1))/dz
!         grad_Opal = opalin(k)/dz!(opalin(k) - opalin(k+1))/dz
!         grad_CaCO3 = caco3in(k)/dz!(caco3in(k) - caco3in(k+1))/dz
!         
!        
!         
!         POM_f = POM_f - (grad_POM * SPOM)
!         Opal_f = Opal_f - (grad_Opal * SOpal )
!!         CaCO3_f = CaCO3_f -  (grad_CaCO3 * SCaCO3)
!
!
!
!       !POM_f = POM_f - (POM/dz * SPOM)
!        !POM_f = POM_f - (POM_ave * SPOM)
!
!
!! if (POM_f .gt. 0.0)POM_f = POM_f !POM_f * 0.67 !DIS NOW
!
!         POM_F_depth = POM_F_depth  + max(0.0, POM) * dz !(max(0.0,POM_f*(1-0.67)) * dz)
!         Opal_F_depth = Opal_F_depth + max(0.0,Opal) * dz
!
!!         CaCO3_F_depth = CaCO3_F_depth + (CaCO3_f * dz)
!        ! POM_F_depth,Opal_F_depth,CaCO3_F_depth - flux of
!        ! POM,Opal,CaCO3 at depth z=75
!         POM_F_Z(k)   = 0.0
!         Opal_F_Z(k)  = 0.0
!!         CaCO3_F_Z(k) = 0.0
!         endif
!
!         POM_F_depth = max(0.0, POM_F_depth)
!         Opal_F_depth = max(0.0, Opal_F_depth)
!        
!         if((k.ge.kzc).and.(k.le.km)) then
!         POM_F_Z(k) = POM_F_depth *  ((k*dz)/Z_c)**(-0.9) !a=0.9 ref the papr for
!
!!this equation  
!         Opal_F_Z(k) = Opal_F_depth *  ((k*dz)/Z_c)**(-0.9) !a=0.9 ref the papr for
!!         Opal_F_Z(k) =  r_Si * Opal_F_depth * exp(-((k*dz)-(Z_c))/d) 
!!         CaCO3_F_Z(k) = r_Cp * CaCO3_F_depth * exp(-((k*dz)-(Z_c))/d) 
!                                                                        !R-rain ratio
!                                                                        !r_CP - C to Phosphorus Redfield ratio.
!                                                                        !d = scale depth
!         !POM_F_Z,Opal_F_Z,CaCO3_F_Z = flux of pom opal and caco3 at a
!                                                                        !each depth
!         end if
!!         if(k.eq.kzcp) then
!!         POM_f = POM_f - ((POM_F_Z(k-1)-POM_F_Z(k))/dz)/(86400*10)
!!         Opal_f = Opal_f - ((Opal_F_Z(k-1)-Opal_F_Z(k))/dz)/(86400*10)
!!         CaCO3_f = CaCO3_f - ((CaCo3_F_Z(k)-CaCO3_F_Z(k))/dz) 
!!         end if
!         
!         if(k.ge.kzcp) then
!         POM_f = POM_f + ((POM_F_Z(k-1)-POM_F_Z(k))/dz)
!         Opal_f = Opal_f + ((Opal_F_Z(k-1)-Opal_F_Z(k))/dz)
!!         Opal_f = Opal_f - ((Opal_F_Z(k-1)-Opal_F_Z(k))/dz)/(86400*10)
!        
!!         CaCO3_f = CaCO3_f - ((CaCO3_F_Z(k)-CaCO3_F_Z(k-1))/dz)
!         end if
!
!!101     continue
!
!       !exchange 
!      ! ExpPON = -1.0 * POM *(40/86400)   !NOW
!         if (k.le.8) then
!         relaxation_NO3 =(NO3 - no3_obss)/86400.0/15.0             !180.0         
!         relaxation_SIOH4=(SIOH4 - sioh4_obss)/86400.0/15.0        !180.0   
!         end if
!
!!!!!CHECKKKKKKKKKK
!!ABOVEEEEEEEEEEEEEEEEEEEEEEEE***********************************
!         if ((k.gt.8).and.(k.le.16)) then
!         relaxation_NO3 = 5*(NO3 - no3_obss)/86400.0/90.0        !90.0
!         relaxation_SIOH4 = 5*(SIOH4 - sioh4_obss)/86400.0/90.0  !90.0
!         end if
!         if (k.gt.16) then
!         relaxation_NO3 = 5*(NO3 - no3_obss)/86400.0/30.0        !30.0
!         relaxation_SIOH4 =5*(SIOH4 - sioh4_obss)/86400.0/30.0  !30.0
!         end if 
!
!
!         NO3_f   = NO3_f - ExcNO3  - relaxation_NO3*1.0 !1e-7*(no3in(k-1)- no3in(k+1))/ 2/dz ! 1e-5*(no3in(k-1)- no3in(k+1))/ 2/dz !+ SPOM * (NO3/dz)! 1.0 !NOW
!         SIOH4_f = SIOH4_f - ExcSiOH4 -  relaxation_SIOH4*1.0
!!*****************************************************
!         POM_f   = POM_f   - ExpPON *0.0 !here -ve came because SPOM is negative
!         Opal_f  = Opal_f  - ExpOpal *0.0
!         CaCO3_f = CaCO3_f - ExpCaCO3 *0.0!1.0
!        
!        !print*,loop,k,ExcNO3,ExpPON
!        PS = max(0.0, PS + PS_f*dt)
!        PL = max(0.0, PL + PL_f*dt)
!        ZS = max(0.0, ZS + ZS_f*dt)
!        ZL = max(0.0, ZL + ZL_f*dt)
!        ZP = max(0.0, ZP + ZP_f*dt)
!        NO3 = max (0.0 , NO3 + NO3_f*dt)
!        NH4 = max (0.0, NH4 + NH4_f*dt)
!        POM = max (0.0, POM + POM_f*dt)
!        DOM = max(0.0, DOM + DOM_f*dt)
!        SIOH4 = max(0.0, SIOH4 + SIOH4_f*dt)
!        Opal = max(0.0, Opal + Opal_f*dt)
!        Ca = Ca + Ca_f*dt  
!        CaCO3 = CaCO3 + CaCO3_f * dt
!        TCO2 = TCO2 + TCO2_f * dt               
!        TALK =TALK + TALK_f * dt                
!        T0 = T0 - (ExcT0 * dt)                 !doubt here.
!        S0 = S0 - (ExcS0 * dt)                 !doubt here
!
!!        UVEL0 = UVEL0 - (ExcUVEL0 * dt)         
!!        VVEL0 = VVEL0 - (ExcVVEL0 * dt)         
!
!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,DELETE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!        T0 = T0 + (ExcT0 * dt)                         !doubt here.
!!        S0 = S0 + (ExcS0 * dt)                         !doubt here
!!        UVEL0 = UVEL0 + (ExcUVEL0 * dt)                !doubt here
!!        VVEL0 = VVEL0 + (ExcVVEL0 * dt)                !doubt here
!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!        
!
!
!
!        return
!        end subroutine
!
!
!        function Monod (a, b)
!        real(8) a, b, Monod
!        Monod = a/(a+b)
!        end function
!
!        function expPsi (a, b)
!        real(8) expPsi, a, b
!        expPsi = exp(-1.0*a*b)
!        end function

        end module nemuro_mod 

