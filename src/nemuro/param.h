cccccccccccccccccccccccccccc include param.h ccccccccccccccccccccccccccccc
c
c
c	parameter definitions for offline model
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	double precision rx, Ah, rintegral_C, rvolume_integral, rxb
	double precision adv_ve, adv_vw, adv_vn, adv_vs,adv_vt,adv_vb
	double precision adv_vtb,adv_vbb
	double precision adv_fe, adv_fw, adv_fn, adv_fs,adv_ft,adv_fb
	double precision adv_feb, adv_fwb, adv_fnb, adv_fsb,adv_ftb
	double precision adv_fbb
	double precision dif_fe, dif_fw, dif_fn, dif_fs, dif_ft, dif_fb
	double precision ube, ubw, vbn, vbs
	double precision ADV_U, ADV_V, ADV_W , shape_function
	double precision ADV_U_Baro, ADV_V_Baro, ADV_W_Baro
	double precision ADV_U_iso, ADV_V_iso, ADV_W_iso
	double precision DIFF_U, DIFF_V, DIFF_W , rVELO_mask
	double precision DIFF_U_bh, DIFF_V_bh
	double precision diff_1_for_w, diff_2_for_w
	double precision rtemp_doubleprecision, reps

	double precision rVolInt, deg2rad
	double precision co2_spatial_function, pi
	double precision frac
	double precision shear_critical

	common /param/ ADV_U, ADV_V, ADV_W, rx, reps, rxb
	common /param/ adv_ve, adv_vw, adv_vn, adv_vs,adv_vt,adv_vb
	common /param/ adv_fe, adv_fw, adv_fn, adv_fs,adv_ft,adv_fb
	common /param/ adv_feb, adv_fwb, adv_fnb, adv_fsb,adv_ftb,adv_fbb
	common /param/ ube, ubw, vbn, vbs,  adv_vtb, adv_vbb
	common /param/ dif_fe, dif_fw, dif_fn, dif_fs
	common /param/ ADV_U_iso, ADV_V_iso, ADV_W_iso
	common /param/ ADV_U_Baro, ADV_V_Baro, ADV_W_Baro
	common /param/ DIFF_U, DIFF_V, DIFF_W
	common /param/ DIFF_U_bh, DIFF_V_bh
	common /param/ Ah, Ah_back, Ah_max, rintegral_C
	common /param/ rvolume_integral, diff_1_for_w, diff_2_for_w
	common /param/ rintegral_DEN_taup, rintegral_DEN_taum
	common /param/ shape_function_gmeddy, rAITH, shape_function

	common /param/ deg2rad
	common /param/ pi, rVolInt

	common /param/ i,j,k
        common /param/ taum, taun, taup
        common /param/ dt, dtts, loop, loop_total, loop_day, run_days
        common /param/ split_dt, split_dtts, loop_split
        common /param/ day2sec


        common /param/ days, loop_start
        common /param/ number_of_snap, month, month_start, iyear_start
        common /param/ month_start_snap, month_co2
        common /param/ lpm, month_offline, iday_month, iday_start
        common /param/ lp6hourly, iday_offline
        common /param/ offline_switch, offline_switch_daily

        common /param/ A1, A2, B1, B2, C
        common /param/ gamma1, gamma2
        common /param/ alpha
        common /param/ van_kar
        common /param/ rlzero
        common /param/ betaT, betaS
        common /param/ gravity
        common /param/ Rf_background
        common /param/ rho0, cp
        common /param/ tx, ty, shf, swf
        common /param/ fcor, omega
        common /param/ alpha_asselin
        common /param/ i_closure_level
        common /param/ itf_source
        common /param/ p_region_value

	common /param/ rintegral_B


c
c====================== polar.h =====================================
c polar co-ordinate specification for the model
c
        common /param/ dphi, dthe, d2phi, d2the, dphi2, dthe2, Re
	real variable_dyd (jmt)
c
c====================== matric.h ===================================
c
	common /param/ deg2meter
        common /param/ dxd, dyd, reflon, reflat, reflat_start
	common /param/ dx, dy, d2x, d2y, dx2, dy2, dx4, dy4, ymatric
	common /param/ nt
c
c=================== include co2.h ===================================
c
	common /param/ anthropogenic_co2, stf
	common /param/ co2_flux_region_sum, area_of_region
	common /param/ jcfc_count, jcfc_base, rpCFC
	common /param/ pCFC_11_South_Tropic, pCFC_11_North_Tropic
	logical aimpl_vmix_NaN_check, asmooth_NaN_check
	common /param/ rimpl_factor
	common /param/ rdelta
	common /param/ relaxation_time
	common /param/ rSg, rDICg, rVirtualFlux, rVirtSaltError
	common /param/ rPminusE, ff_save, rsc, rALKg, rPMEg, rPminusE_Alk
	common /param/ rCO2exchange, rCO2, rdCO2, rdPCO2, rOceanPCO2
	common /param/ rsalt, wx, rkw, rmodified_rbgc_Sb
	common /param/ rkwco2, rkwo2, izone(19), izone_opt(19)
	real stf
	integer izone
c
c==================== end co2.h ======================================
c
	
c
c====================== file.h ======================================
c for file index definitions
c
        integer ifile_unit
        common /param/ ifile_unit(5)
        character*5 ifile_name
        character*82 ifile_name_yearly
        character*84 ifile_name_daily_u10m
        character*84 ifile_name_daily_v10m
        dimension ifile_name(5)
	data ifile_name(1)/'u.trm'/
	data ifile_name(2)/'v.trm'/
	data ifile_name(3)/'w.trm'/
	data ifile_name(4)/'t.trm'/
	data ifile_name(5)/'s.trm'/

c
c===================== snaps.h ===========================================
c
        parameter (itime=1)
        integer ncid,vlonid,vlatid,vdepthid,vtimeid
        integer start4,start3(3),count4(4),count3(3), count33(3)
        integer startp(2)
        integer startinv(2), countinv(2)
        integer loop_ind, start2(2),count2(2)
        integer vardim(4), vtid, vsid, vuid, vvid, vwid, vardim3(3)
        integer vetaid
        integer vardim33(3), vardim2(2)
        integer varpdim(2)
        integer varinvdim(2)
        integer vqnetid, vqsaltid, vtauxid, vtauyid
        integer vsaltid, vdensid, vvortid
	integer vtempid(nn)
        integer vpxid, vpyid
        integer vthetaid
        integer vmaskid
        real rlon(imt),rlat(jmt),rdepth(km),rtime(itime)
	real var3d(imt,km,jmt)
        real undef
        character*21 stamp_day
        character*20 stamp_time
        character*4 Tr
        dimension stamp_day(12)
        dimension stamp_time(12)
	dimension Tr(nn)
        character*21 stamp_day_daily
        character*20 stamp_time_daily
        common stamp_day_daily(365)
        common stamp_time_daily(365)
        dimension start4(4)
        integer number_of_snap
	
	integer vvar3did(nvar3d), ncid_var3d(nvar3d)
	integer vlonid_var3d(nvar3d), vlatid_var3d(nvar3d)
	integer vdepthid_var3d(nvar3d), vtimeid_var3d(nvar3d)
	integer lonid_var3d(nvar3d), latid_var3d(nvar3d)
	integer idepthid_var3d(nvar3d), itimeid_var3d(nvar3d)
	integer vardim_var3d(4,nvar3d), vardim_var3d_pass(4)
	character*4 var3d_units
	character*4 var3d_char
	character*88 var3d_file_name  
	character*21 ncaptc_var3d
	dimension ncaptc_var3d(nvar3d)
	dimension var3d_file_name(nvar3d)
	dimension var3d_char(nvar3d), var3d_units(nvar3d)
	data var3d_char /"DICT", "FLUX", "PHOS", "ALKT", "DOPT"/ 
	data var3d_units /"m/m3", "m/m2", "m/m3", "m/m3", "m/m3"/ 
	
	data start3/1,1,1/
        data start2/1,1/
        data count2/imt,jmt/
        data count4/imt,km,jmt,1/
        data count3/imt,jmt,itime/
        data count33/imt,km,jmt/

#include "stamp.h"	
c
c==================== end of snaps.h ======================================c
c
c

c
c==================== include isopyc.h ===============================
c
	 risopicSmax = 1.0e-2
c         rAITH = 1.e3   ! horizontal GM bolus velocity coeff.
         rAI = 1.e3     ! horizontal Redi flux coeff.
         rAIvert = 1.e-3 ! vertical Redi flux coeff.
         !!! These coefficients are re-defined inside the
         !!! gm-iso-.F and diffusion.F files with a latitudinal
         !!! factor in order to reduce the mangnitude of coefficient
         !!! in a converging meridion. (3-Mar-2007)
c
c=================== end isopyc.h =====================================
c


	im = i - 1
	jm = j - 1
	kmone = k - 1
	ip = i + 1
	jp = j + 1
	kp = k + 1
	
	A1 = 0.92
        A2 = 0.74
        B1 = 16.6
        B2 = 10.1
        C  = 0.08

        gamma1 = (1.0/3.0) - (2.0*A1/B1)
        gamma2 = (B2 + 6.0*A1)/B1
        gravity = 9.8
        alpha = 0.1
        van_kar = 0.4

        omega = 0.729e-4
        pi = 3.14159265359
        deg2rad = pi/180.0
	Re = 6368.0e3
	
	
	deg2meter = deg2rad*Re
	dxd= 1.0/1.0
	dyd= 1.0/1.0
	reflon = 72.5  ! global 1x1 ORAS4_ASSIM
	reflat_start = -79.25  ! global 2x2
	reflat_start = -81.50  ! java
	reflat = -0.25 - (-reflat_start - (jmt*dyd + reflat_start))/2.0
	dx=dxd*deg2meter
        dy=dyd*deg2meter
        d2x=2.0*dx
        d2y=2.0*dy
        dx2=dx*dx
        dy2=dy*dy
        dx4=dx2*dx2
        dy4=dy2*dy2
        dphi=dxd*deg2rad
c        dthe=dyd*deg2rad
	d2phi=2.0*dphi
c        d2the=2.0*dthe
	

        rho0 = 1022.0
        cp = 4217.4
	
	day2sec = 86400.0
	year2sec = day2sec*365.0
	dt = 86400  ! time step for interpolation
	dtts = 2.0*dt
	split_dt = 7200 ! tracer time step
	split_dtts = 2.0*split_dt
	rimpl_factor = 0.5
	alpha_asselin = 0.05
	relaxation_time = 30*day2sec
	
c	Ah = 1.0e-3  ! most reliable version for six months
c	Ah = 1.0e-6
	Ah_back = 1.0e-4
	Ah_max = 0.01
c	Am = 2.0e4
	
	shf = 0.0
        swf = 0.0001
c        stf = 0.0 ! suppressing initialization here
c
c co2 injection region in lon lat cooridnate
c
	lon_start = 280    ! add reflon in mind
	lon_end = 365      ! add reflon in mind
	lat_start = 65.0   ! region 3 for GFDL_ASSIM
	lat_end = 78.0  ! region 3 for GFDL_ASSIM
	
	lon_start = int ( ((lon_start/dxd) - reflon/dxd))
	lon_end = int ( ((lon_end/dxd) - reflon/dxd))
	lat_start = int ( ((lat_start/dyd) + abs(reflat_start)/dyd))
	lat_end = int ( ((lat_end/dyd + abs(reflat_start)/dyd)))
	
	i_co2_injection_region(1) = lon_start
        i_co2_injection_region(2) = lon_end
        j_co2_injection_region(1) = lat_start
        j_co2_injection_region(2) = lat_end
	
	
	iread_u = 1
	iread_v = 2
	iread_w = 3
	iread_t = 4
	iread_s = 5
