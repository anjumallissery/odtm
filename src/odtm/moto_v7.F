        
        module moto_mod

        use implicit_thomas_mod, only : implicit_vmix

        implicit none

        integer, parameter :: kmax=51

        contains

        SUBROUTINE MT1 (tempin, saltin, tr_in, uvelin, vvelin, tauxin,tauyin, 
     &     qnet, sswin, pmein, ntr, rLATin,dt_rg,SM,SH,we_prof, we_baro, icval,R,
     &     entr_flux, speed, pme_eps, rvr_runoff,
     &     diag1, diag2, diag3, diag4, diag5, diag6,ii,jj) 

        implicit none
        

 
! ONE DIMENSIONAL MIXED LAYER MODEL

      real UN(kmax),VN(kmax),TEN(kmax),SAN(kmax), TRN(kmax,ntr),
     1 U1(0:kmax),V1(0:kmax),TE1(0:kmax),SA1(0:kmax),TR1(0:kmax,ntr),Z(kmax),Z1(kmax),
     2 rQ(kmax),SM(kmax),SH(kmax),R(kmax),YQ(kmax),SM1(kmax),SH1(kmax)
        real tempin(kmax), saltin(kmax), tr_in(kmax,ntr)
        real uvelin(kmax), vvelin(kmax), we_prof(0:kmax+1), we_baro
        real sswin(kmax), zave
        real tauxin, tauyin, qnet, pmein, rLATin, dt_rg
        real :: dq, mf, qf, cp, h, ds, pme_eps
        integer :: nloop
        integer, intent(in) :: ntr
        integer ic(2)
        real :: uvloop(kmax), entr_flux(kmax)
        real :: diag1(kmax), diag2(kmax), diag3(kmax)
        real :: diag4(kmax), diag5(kmax), diag6(kmax)
        integer :: ii, jj, ntrn
        
        real :: a1, a2, b1, b2, c1, day, deg2rad, dt, dtz, dz, dz2, dzz, f0, g1, g2
        real :: ga, r5, ro, y0, z0, z0d, al, bs, e1, e2, e3, e4, e5, e6, e7, f1, f2
        real :: f3, f4, ft, g3, g4, gc, h1, h2, h3, hh, hj1, hj2, twx, twy, tx
        !real :: t
        real :: xk, xk2, afkph, bt, bx, dfkph, diff_bryan_lewis, pi, q1, q2, qx
        real :: rcm2tom2, rtemp, sfkph, sz, txy, ty, UO, uu, uv, uvs, VO, vv
        real :: x1, x2, x3, ya, zfkph, zk, rvr_runoff, speed, UC2, VC2
        integer :: icval, m, naz, naz1, naz2, k1, j

! VARIABLES
! U,V=VELOCITY,TE=TEMP.,SA=SAL.,Z=DEPTH
! WITH N=NEW VALUE AT EACH TIME STEP, 1=OLD VALUE
! Q=TURBULENT KINETIC ENERGY
! SM,SH=MIXING COEFFICIENTS
! R=RICHARDSON NUMBER,YQ=MIXING LENGTH
! NUMBER OF GRIDS
      NAZ=kmax-1
      NAZ1=NAZ+1
      NAZ2=NAZ+2
! GRID SIZE (THICKNESS IN METER)
      DZ=5.0
      DZZ=DZ*DZ
      DZ2=DZ*2.0
! ROUGHNESS
      Z0=0.1
! INITIAL CONDITION
      Y0=1.0
      DO 10 M=1,NAZ1
        UN(M)=uvelin(M)        !0.0
        VN(M)=vvelin(M)  
!    UN(M)=0.0
!      VN(M)=0.0
! STABLE STRATIFICATION
        TEN(M)=tempin(M) !10.0-0.02*(M-1)! temp variation
        SAN(M)=saltin(M) 
        do ntrn = 1, ntr
            TRN(M,ntrn)=tr_in(M,ntrn) 
        end do
      U1(M)=UN(M)
      V1(M)=VN(M)
      TE1(M)=TEN(M)
      SA1(M)=SAN(M)
      do ntrn = 1, ntr
        TR1(M,ntrn)=TRN(M,ntrn)
      end do
   10 CONTINUE

        entr_flux(:) = 0.0

      DO 11 M=1,NAZ1
   11 Z(M)=DZ*(1-M)-Z0
      DO 12 M=1,NAZ1
   12 Z1(M)=DZ*(1-M)

        
! SECONDS IN ONE DAY
      DAY=86400.0
! TIME STEP
      !DT=0.02*DAY !now one day run is done in 6 loops 
      DT=dt_rg/10 ! seconds , same as RG model
        ic(:) = 1
        ic(1) = int(abs(maxval(we_prof)*1e5 )/ 10)
        ic(2) = int(abs(minval(we_prof)*1e5 )/ 10)
        icval = max(1,maxval(ic))
        DT=dt_rg/10/icval
      nloop=dt_rg/DT
      DTZ=DT/DZZ
! WATER DENSITY
      RO=1024.5
      R5=10.0**(-5)
      Z0D=Z0/DZ
! GRAVITY ACCELERATION
      GA=9.8
! CORIOLIS PARAMETER
        deg2rad = 3.14/180.0
        F0= 2.0*(0.729e-4)*sin(rLATin*deg2rad) 
!      F0=1.2*10.0**(-4) 
      FT=F0*DT
! COEFFICIENTS IN THE TURBULENT MODEL
      A1=0.92
      A2=0.74
      B1=16.6
      B2=10.1
      C1=0.08
      G1=(1.0/3.0)-2.0*A1/B1
      G2=(B2+6.0*A1)/B1
      E1=3.0*A2
      E2=G1+G2
      GC=G1/E2
      E3=A1/A2
      E4=B1*(G1-C1)
      E5=E4+6.0*A1+3.0*A2
      E6=B1*G1
      E7=B1*E2-3.0*A1
      F1=2.0*A1*E5
      F2=A1*E4
      F3=A2*E7
      F4=A2*E6
      G3=SQRT(B1*(3.0*A1*(G1-C1))**3) 
      G4=3.0*A2*G1
! SALINE EXPANSION COEFFICIENT
      BS=(-76.0)*10.0**(-5)
! VON KARMAN CONSTANT
      XK=0.4
      XK2=XK*XK
      G3=G3*XK2
      G4=G4*SQRT(B1*3.0*A1*(G1-C1))*XK2*Z0
      XK=(-1.0)*0.5*XK
      AL=0.4 !anju changed from 0.04 to 0.4 08Aug2020

! COEFFICIENTS FOR NEAR SURFACE INTERPOLATION
      H2=ALOG((1.0+Z0D)/(2.0+Z0D))
      H3=ALOG(Z0D/(1.0+Z0D))
      H1=H3-H2
      HH=H2/DZ+1.0/Z0
      HJ1=(H2+1.0/(Z0D+0.5))/(H2-H3)
      HJ2=1.0-HJ1
        
! TIME INTEGRATION STARTS
!     K=0
!   1 K=K+1
      SM1(:) = 0.0
      SH1(:) = 0.0
      we_prof(1) = we_prof(1) + pmein
      K1=0
    2 K1=K1+1
! TIME INCREASING
!      T=T+DT !Prajeesh  
! WIND STRESS
      TWX = tauxin
      TWY = tauyin

        TE1(0) = TE1(1)
        SA1(0) = SA1(1)
        do ntrn = 1, ntr
            TR1(0,ntrn) = TR1(1,ntrn)
        end do
        U1(0)  = U1(1)
        V1(0)  = V1(1)
        TE1(kmax) = TE1(kmax-1)
        SA1(kmax) = SA1(kmax-1)
        do ntrn = 1, ntr
            TR1(kmax,ntrn) = TR1(kmax-1,ntrn)
        end do
        U1(kmax)  = U1(kmax-1)
        V1(kmax)  = V1(kmax-1)
        
! TWX=1.0
! TWY=1.0
! VERTICAL MIXING
      TX=TWX/RO
      TY=TWY/RO
      TXY=SQRT(TX*TX+TY*TY)
      DO 22 M=1,NAZ 
      R(M)=0.0
      SM(M)=0.0
      SH(M)=0.0
   22 CONTINUE

      DO 23 M=1,NAZ
      SZ=0.5*(SA1(M)+SA1(M+1))-30.0
! THERMAL EXPANSION COEFFICIENT
      !BT=((1.3-0.025*SZ)*(TE1(M)+TE1(M+1)) + (3.3+0.325*SZ))*R5  
      BT=((1.3-0.025*SZ)*(TE1(M)+TE1(M+1))*0.5 + (3.3+0.325*SZ))*R5 !0.5
                                                                    !in this equation is added by anju 13Aug2020 
!NOTE FOR CORRECTION, NEED one EXTRA *0.5
      IF(M.EQ.1) GO TO 24
! STRATIFICATION
      BX=BT*(TE1(M-1)-TE1(M+1))+BS*(SA1(M-1)-SA1(M+1))
! VERTICAL SHEAR
      UV=(U1(M-1)-U1(M+1))**2+(V1(M-1)-V1(M+1))**2
      GO TO 25
! NEAR SURFACE STRATIFICATION AND SHEAR
   24 BX=HJ1*(BT*(TE1(1)-TE1(2))+BS*(SA1(1)-SA1(2)))+ 
     & HJ2*(BT*(TE1(2)-TE1(3))+BS*(SA1(2)-SA1(3)))
      UU=HJ1*(U1(1)-U1(2))+HJ2*(U1(2)-U1(3))
      VV=HJ1*(V1(1)-V1(2))+HJ2*(V1(2)-V1(3))
      UV=UU**2+VV**2
        uvloop(M) = UV
   25 IF(UV .LT. 10.0**(-8)) GO TO 23 
! R,SH AND SM CALCULATED
      BX=BX/UV
      BX=BX*GA*DZ2
      X1=F2+F3*BX
      X2=F4*BX
      X3=SQRT(X1**2-2.0*F1*X2)
      R(M)=(X1-X3)/F1 
      SH(M)=E1*(G1-E2*R(M))/(1.0-R(M))
      IF(R(M).GT.GC) SH(M)=0.0
      SM(M)=SH(M)*E3*(E4-E5*R(M))/(E6-E7*R(M))  
      !if (ii .eq. 44 .and. jj .eq. 31) write(*,*)'from moto=',M, SM(M), YQ(M), rQ(M)
   23 CONTINUE
! ITERATIVE CALCULATION OF MIXING LENGTH
      J=0
   30 J=J+1 
      Q1=0.0
      Q2=0.0
      DO 31 M=1,NAZ
      ZK=XK*(Z(M)+Z(M+1))
      YQ(M)=ZK/(1.0+ZK/Y0)
      QX=(U1(M-1)-U1(M+1))**2+(V1(M-1)-V1(M+1))**2
      QX=QX*B1*SM(M)*(1.0-R(M))/(DZZ*4.0)
      QX=SQRT(QX)
      rQ(M)=QX*YQ(M) 
      Q1=Q1+rQ(M)
      Q2=Q2-rQ(M)*0.5*(Z(M)+Z(M+1))
   31 CONTINUE
      IF(Q1.EQ.0.0) GO TO 32
      YA=Y0
      Y0=AL*Q2/Q1
      YA=ABS(Y0-YA)
      IF(J.EQ.10) GO TO 32
      IF(YA.GT.0.01) GO TO 30 
   32 DO 33 M=1,NAZ
      SM(M)=SM(M)*YQ(M)*rQ(M)
      SH(M)=SH(M)*YQ(M)*rQ(M) 
   33 CONTINUE


        DO M=1,NAZ
        diag1(M) = rQ(M)
        diag2(M) =  R(M)
        diag3(M) = YQ(M)
        diag4(M) = E1*(G1-E2*R(M))/(1.0-R(M))
        diag5(M) = diag3(M)*E3*(E4-E5*R(M))/(E6-E7*R(M))
        diag6(M) = 0.0
        ENDDO

! BYAN & LEWIS (1979) MIXING
         afkph = 0.8
         dfkph = 1.05
         sfkph = 4.5e-5
         zfkph = 100.0e2
         rcm2tom2 = 1.0e-4
         pi = 3.14
        DO M=1,NAZ
        rtemp = sfkph*(DZ*M*(100.0) - zfkph )
        diff_bryan_lewis = afkph + dfkph*atan(rtemp)/pi
        SH(M) = SH(M) + diff_bryan_lewis*rcm2tom2 * 0.0
        ENDDO
! SMOOTHING
      DO 34 M=1,NAZ1
      SM1(M)= SM1(M) + SM(M)
      SH1(M)= SH1(M) + SH(M)
   34 CONTINUE
      DO 35 M=2,NAZ
      SM(M)=(SM1(M-1)+2.0*SM1(M)+SM1(M+1))/4.0/K1
      SH(M)=(SH1(M-1)+2.0*SH1(M)+SH1(M+1))/4.0/K1
   35 CONTINUE
! TIME INTEGRATION OF U,V,TE AND SA


        
        
      DO 40 M=1,NAZ
        UN(M) = UN(M) + FT*V1(M) - U1(M)*1.e-3 ! DAMPING
        VN(M) = VN(M) - FT*U1(M) - V1(M)*1.e-3 ! PREVIOUSLY IT WAS UN(M)  
   40 CONTINUE
      DO 60 M=1,NAZ
        TEN(M) = TEN(M) - we_prof(M) * (TE1(M-1) - TE1(M+1))*DT/DZ2
        SAN(M) = SAN(M) - we_prof(M) * (SA1(M-1) - SA1(M+1))*DT/DZ2
        do ntrn = 1, ntr
            TRN(M,ntrn) = TRN(M,ntrn) - we_prof(M) * (TR1(M-1,ntrn) - TR1(M+1,ntrn))*DT/DZ2
        end do
        UN(M) = UN(M) - we_prof(M) * (U1(M-1) - U1(M+1))*DT/DZ2
        VN(M) = VN(M) - we_prof(M) * (V1(M-1) - V1(M+1))*DT/DZ2
        entr_flux(M) = entr_flux(M) + we_prof(M) * (SA1(M-1) -
     &                 SA1(M+1))*DT/DZ2
   60 CONTINUE
      DO 61 M=2,NAZ-1
        TEN(M) = TEN(M) - (we_baro*TE1(M-1) - we_baro*TE1(M+1))*DT/DZ2
        SAN(M) = SAN(M) - (we_baro*SA1(M-1) - we_baro*SA1(M+1))*DT/DZ2
        do ntrn = 1, ntr 
            TRN(M,ntrn) = TRN(M,ntrn) - (we_baro*TR1(M-1,ntrn) - we_baro*TR1(M+1,ntrn))*DT/DZ2
        end do
   61 CONTINUE
!        UN(1) = UN(1) + TX*1.0*DT/DZ
!        UN(1) = UN(1) + TX*1.0*DT/DZ
!        VN(1) = VN(1) + TY*1.0*DT/DZ ! PREVIOUSLY IT WAS UN(M)

!        GOTO 42
! CALCULATION OF THE TOP GRID VALUES
      UO=0.0
      VO=0.0
      UVS=0.01
      UV=SQRT(TXY/G3)
      IF(UV.LT.0.00001) GO TO 41 ! TOO WEAK SHEAR
      UO=TX/G3/UV 
      VO=TY/G3/UV 
      UVS=SQRT(UO**2+VO**2) 
   41 UC2=((UN(2)-UN(3))/DZ-UO/Z0)/HH 
      UN(1)=UN(2)+UC2*H1+UN(2)-UN(3)
      VC2=((VN(2)-VN(3))/DZ-VO/Z0)/HH 
      VN(1)=VN(2)+VC2*H1+VN(2)-VN(3)
! SURFACE BOUNDARY CONDITION
! COOLING 0.1 DEG PER DAY
   42 CONTINUE
        qf = qnet ! read value of heat flux here  
        IF (speed .le. 2.0) qf = qf - 25.0  ! CONVECTIVE LATENT HEAT
        mf = pmein - pme_eps ! read value of freshwater flux in here in m/s
        cp = 4200.0
        h = DZ
        dq= qf*DT/(RO*cp*h)
        ds= mf*35.0*DT/DZ ! Sref= 35.0     


        TEN(1)= TEN(1) + dq
        SAN(1)= SAN(1) + ds
        SAN(2)= SAN(2) - rvr_runoff * 35.0*DT/DZ

        !print *, "What about this ??? Prajeesh  line no: 333 moto_v7.F"

        DO M=2,NAZ
        TEN(M)=TEN(M) + (sswin(M-1)-sswin(M+1-1)) * DT/(RO*cp*DZ)
        ENDDO

        call implicit_vmix (SM, DZ, DT, NAZ1, UN, UN)
        call implicit_vmix (SM, DZ, DT, NAZ1, VN, VN)
        call implicit_vmix (SH, DZ, DT, NAZ1, TEN, TEN)
        call implicit_vmix (SH, DZ, DT, NAZ1, SAN, SAN)
        do ntrn = 1, ntr
            call implicit_vmix (SH, DZ, DT, NAZ1, TRN(:,ntrn), TRN(:,ntrn))
        end do
! REPLACE THE OLD VALUES WITH THE NEW VALUES
   43 DO 52 M=1,NAZ 
      U1(M)=UN(M) 
      V1(M)=VN(M)
      TE1(M)=TEN(M)
      SA1(M)=SAN(M)
      do ntrn = 1, ntr
        TR1(M,ntrn)=TRN(M,ntrn)
      end do
   52 CONTINUE
! TIME STEPS FOR 0.1 DAY
! PRINT OUT THE RESULTS
!      IF(K1.LT.50) GO TO 2 
!c      WRITE(6,101) T/DAY
!  101 FORMAT(1H ,E15.7)
!      WRITE(6,102) (10.0*UN(M),M=1,51,1)
  102 FORMAT(1H ,51F8.2)
!      WRITE(6,102) (10.0*VN(M),M=1,51,1)
!      WRITE(6,102) (TEN(M),M=1,51,1)
!      WRITE(6,102) (100.0*SM(M),M=1,11,1)
!      WRITE(6,102) (100.0*SH(M),M=1,11,1)
! TIME STEPS FOR 0.1 DAY X 50

      IF(K1.LT. nloop) GO TO 2

!      WRITE(6,102) (100.0*SM(1))!,M=1,11,1)
!      WRITE(6,102) (100.0*SH(1))!,M=1,11,1)
      
          DO M=1,NAZ
        tempin(M) = TEN(M)
        saltin(M) = SAN(M)
        do ntrn = 1, ntr
            tr_in(M,ntrn) = TRN(M,ntrn)
        end do
        uvelin(M) = UN(M)
        vvelin(M) = VN(M)
        entr_flux(M) = entr_flux(M) /nloop
        ENDDO

        RETURN
      STOP
      END SUBROUTINE

        end module moto_mod

