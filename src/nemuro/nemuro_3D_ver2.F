ccccc CHECK include old3_nemuro.h_ver3 cccccccccccccccccccccccccc
ccMUSTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT ccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	This module calculates the Nemuro biogeochemistry model
c	as in Yamanaka et al., (2004), J. Oceanography, Vol. 60, 227-241
c
c
c
c
c
c
c	Model coding start on:	04-Dec-2014
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine nemuro_3D (
     I  PS,PL,ZS,ZL,ZP,NO3,NH4,POM,DOM,SIOH4,Opal,Ca,CaCO3,TCO2,TALK,
     I  T0,S0,UVEL0,VVEL0, loop, dt, rlight,k,km,pomin,opalin,caco3in,
     I  no3in,sioh4in,
     I  wvelocity,tempin,saltin,uvelin,vvelin,t0_grad, s0_grad,
     I  uvel0_grad, vvel0_grad,no3_grad,POM_F_depth,Opal_F_depth,
     I  CaCO3_F_depth,no3_obss,sioh4_obss,
     I  RnewS,RnewL,NO3_lim,NH4_lim,N_lim,Si_lim,diff_N_and_Si_lim,
     I  no3_relax_obs,sioh4_relax_obs,sfx)  !PSCChla,PLCChla,mlight)
        
!show error here
#include "nemuro_new.h" !"07_may_2018_nemuro.h_ver3_paper"
        
        
	real dz !real 8 or 4 
	real(8) Monod, expPsi
	real(8) rkappa, rlight, rsum, sfx, day_night
	real(8) rtemp1, rtemp2
        real w1,w2,no31,no32,sioh41,sioh42,t01,t02,s01,s02,uvel01,uvel02
        real vvel01,vvel02 
        real pom_grad(51),opal_grad(51),caco3_grad(51),no3_grad(51)
        real sioh4_grad(51)
        real t0_grad(51),s0_grad(51),uvel0_grad(51),vvel0_grad(51) 
        real pomin(51),opalin(51),caco3in(51),no3in(51),sioh4in(51) 
        real pomin_sm(51),opalin_sm(51),caco3in_sm(51),no3in_sm(51)
        real sioh4in_sm(51) 
        real wvelocity(51)
        real tempin_sm(51),saltin_sm(51),uvelin_sm(51),vvelin_sm(51)
        real tempin(51),saltin(51),uvelin(51),vvelin(51) 
        real POM_F_depth,Opal_F_depth,CaCO3_F_depth,a
        real POM_F_Z(51),Opal_F_Z(51),CaCO3_F_Z(51),R
        real r_CP,d
        real w_diff,wb
        real delta,grad_POM, grad_Opal, grad_CaCO3 
        real POM_ave
        real no3_obss,sioh4_obss
        real relaxation_NO3,relaxation_SIOH4      
        real NO3_lim,NH4_lim,N_lim,Si_lim,diff_N_and_Si_lim
        real no3_relax_obs,sioh4_relax_obs
!        real(8) NLLGR,ter1,PSCCh,PSC,PSCChla(51)
!        real(8) PLCCh,PLC,PLCChla(51),mlight
	integer loop, dt
        integer k,km,kmm,Z_c, kzc, kzcp                                 !depth upto which flux calculating
        
	real(8) PS_f, PL_f, ZS_f, ZL_f, ZP_f, NO3_f, NH4_f, POM_f, DOM_f,
     &             SIOH4_f, Opal_f, Ca_f, CaCO3_f, TCO2_f, TALK_f 
        
        kmm = km -3 
        dz  = 5.0
        a   = 0.9   
       
             
                                     !constant used in POM flux calculation refr more....
        R= 0.08                           !0.07 also possible ref paper
        r_CP = 117.0                      ! 106 is also possible. check ref. paper..
!        r_Si = 2.0 !From NEmuro          
        d = 3500.0   
                                                   ! meter
	PS_Photosynthesis = 
     &  VmaxS *( Monod (NO3, KNO3S)*expPsi(PsiS, NH4) +
     &           Monod (NH4, KNH4S))                 *
     &  exp(KS*T0)*(rlight/IoptS)*exp(1.0d0 - rlight/IoptS) * PS
        

	RnewS = Monod(Monod(NO3, KNO3S)*
     &   expPsi(PsiS, NH4) ,Monod(NH4, KNH4S))
       
	rtemp1 = Monod (NO3, KNO3L)*expPsi(PsiL, NH4) +
     &               Monod (NH4, KNH4L)
	rtemp2 = Monod(SIOH4,KSIL)

	PL_Photosynthesis = VmaxL *
     &  min (rtemp1, rtemp2)* 
     &  exp(KL*T0)*(rlight/IoptL)*exp(1.0 - rlight/IoptL) * PL


!!!DELETE THIS****************NO SILICA LIMITATION CASE
!         PL_Photosynthesis = VmaxL *
!     &  rtemp1 * 
!     &  exp(KL*T0)*(rlight/IoptL)*exp(1.0 - rlight/IoptL) * PL
!*****************************

        
!************************************************************************************************
                ! depth profile of C:Chla ratio for PS
!        NLLGR = VmaxS* exp(KS*T0) * ( Monod (NO3, KNO3S)*expPsi(PsiS,
!     &  NH4) + Monod (NH4, KNH4S))
               
              
!        PSC = PSmax - (KPPS * NLLGR)
        
!        ter1 = (log(sfx)-log(mlight))/4.605
!        PSCCh = (PSC - (PSC - PSmin))* ter1
        !PSCCh = PSC - ((PSC - PSmin)* ter1)


!        PSCChla(k) = PSCCh
        
        
               !depth profile of C:Chla ratio for PL
        
!        NLLGR = VmaxL * exp(KL*T0) * min (rtemp1,rtemp2)
!        PLC  =PLmax - (KPPL * NLLGR)
!        PLCCh = (PLC -(PLC - PLmin))*ter1
!        !PLCCh = PLC - ((PLC - PLmin)* ter1)
!        PLCChla(k) = PLCCh
        
        
!********************************************************************************************************                        

        NO3_lim =  Monod (NO3, KNO3L)*expPsi(PsiL, NH4)         !nitrate lim.
        NH4_lim= Monod (NH4, KNH4L)                             !ammonium limitation
        N_lim= rtemp1                                           !nitrogen limitation
                

        NO3_lim =  Monod (NO3, KNO3L)*expPsi(PsiL, NH4)         !nitrate lim.
        NH4_lim= Monod (NH4, KNH4L)                             !ammonium limitation
        N_lim= rtemp1                                           !nitrogen limitation
        Si_lim = rtemp2                                         !silicate limitation
        

        NO3_lim =  Monod (NO3, KNO3L)*expPsi(PsiL, NH4)!nitrate lim.
        NH4_lim= Monod (NH4, KNH4L)!ammonium limitation
        N_lim= rtemp1 ! nitrogen limitation
        Si_lim = rtemp2 !silicate limitation
        diff_N_and_Si_lim = rtemp1-rtemp2 ! difference b/w N and Si 
        !limitation
        
        
!        if((loop.le.200).and.(k.eq.1))then
!        print*, loop,k,NH4_lim
!        end if
        
	RnewL = Monod(Monod(NO3, KNO3L)*
     & expPsi(PsiL, NH4) ,Monod(NH4, KNH4L))
         
	
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
	PS_Grazing_by_ZS = GRmaxS * 
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGS*T0) * ZS

	rtemp1 = 1.0 - exp(LamdaL*(PSZL - PS))
	PS_Grazing_by_ZL = GRmaxLPS *
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGL*T0) * ZL

	rtemp1 = 1.0 - exp(LamdaL*(PLZL - PL))
	PL_Grazing_by_ZL = GRmaxLPL *
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGL*T0) * ZL

	rtemp1 = 1.0 - exp(LamdaL*(ZSZL - ZS))
	ZS_Predating_by_ZL = GRmaxLZS *
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGL*T0) * ZL

	rtemp1 = 1.0 - exp(LamdaP*(PLZP - PL))
	PL_Grazing_by_ZP = GRmaxPPL *
     &       max (0.0d0, rtemp1 ) *
     &       expPsi( PusaiPL, ZS+ZL)*exp(KGP*T0)*ZP

	rtemp1 = 1.0 - exp(LamdaP*(ZSZP - ZS))
	ZS_Predating_by_ZP = GRmaxPZS *
     &       max (0.0d0, rtemp1 ) *
     &       expPsi( PusaiZS, ZL)*exp(KGP*T0)*ZP

	rtemp1 = 1.0 - exp(LamdaP*(ZLZP - ZL))
	ZL_Predating_by_ZP = GRmaxPZL *
     &       max (0.0d0, rtemp1 ) *
     &       exp (KGP*T0) * ZP

	ZS_Excretion = (GammaZS - BettaZS)* PS_Grazing_by_ZS
	
	ZL_Excretion = (GammaZL - BettaZL)* (PS_Grazing_by_ZL +
     &                 PL_Grazing_by_ZL + ZS_Predating_by_ZL )

    	
	ZP_Excretion = (GammaZP - BettaZP)* (PL_Grazing_by_ZP +
     &                 ZS_Predating_by_ZP + ZL_Predating_by_ZP )

	ZS_Egestion = (1.0d0 - GammaZS) * PS_Grazing_by_ZS

	
	ZL_Egestion = (1.0d0 - GammaZL) * (PS_Grazing_by_ZL +
     &                 PL_Grazing_by_ZL + ZS_Predating_by_ZL)

	ZP_Egestion = (1.0d0 - GammaZP) * (PL_Grazing_by_ZP +
     &                 ZS_Predating_by_ZP + ZL_Predating_by_ZP)

	
	POM_Remineralization =  VPA0*exp(KPA*T0)*POM

	POM_Decomposition_to_DOM = VPD0*exp(KPD*T0)*POM

	DOM_Remineralization = VDA0*exp(KDA*T0)*DOM

	Opal_Decomposition = VOpal*exp(KOpal*T0)*Opal

	CaCO3_Decomposition =VCaCO3*exp(kCaCO3*T0)*CaCO3
	
	

	Nitrification = NNit0 * exp(KNit*T0)*NH4
        
	
	Si_PL_Shell_Formation = (PL_Photosynthesis - PL_Respiration
     &            - PL_Extracellular_Excretion)*RSiN


	Si_PL_Mortality = PL_Mortality * RSiN
	
	Si_ZL_Egestion = PL_Grazing_by_ZL * RSiN

	Si_ZP_Egestion = PL_Grazing_by_ZP * RSiN


        Ca_PS_Shell_Formation =( PS_Photosynthesis - PS_Respiration
     &       - PS_Extracellular_Excretion) * RCN * Rcoco * RCco  
	
	
        Ca_ZS_Shell_Formation = (PS_Grazing_by_ZS) * BettaZS * RCN 
     &       * Rfora * RCfo
         
        

        Ca_PS_Mortality = (PS_Mortality) * RCN * Rcoco * RCco 

        Ca_ZS_Mortality = (ZS_Mortality) * RCN * Rfora * RCfo 
        
        Ca_ZS_Egestion = (PS_Grazing_by_ZS) * RCN *  Rcoco * RCco 

        Ca_ZL_Egestion = (PS_Grazing_by_ZL) * RCN * Rcoco * RCco 
     &                 + (ZS_Predating_by_ZL) * RCN * Rfora * RCfo 

        Ca_ZP_Egestion = (ZS_Predating_by_ZP) * RCN * Rfora * RCfo                
	PS_f = PS_Photosynthesis - PS_Respiration - 
     &         PS_Extracellular_Excretion - PS_Mortality -
     &         PS_Grazing_by_ZS  - PS_Grazing_by_ZL


	PL_f = PL_Photosynthesis - PL_Respiration -
     &         PL_Extracellular_Excretion - PL_Mortality - 
     &         PL_Grazing_by_ZL - PL_Grazing_by_ZP

	ZS_f = PS_Grazing_by_ZS - ZS_Excretion -
     &         ZS_Egestion - ZS_Mortality -
     &         ZS_Predating_by_ZL - ZS_Predating_by_ZP


	ZL_f = PS_Grazing_by_ZL + PL_Grazing_by_ZL +
     &         ZS_Predating_by_ZL - ZL_Excretion -
     &         ZL_Egestion - ZL_Mortality -
     &         ZL_Predating_by_ZP
	
	ZP_f = PL_Grazing_by_ZP + ZS_Predating_by_ZP + 
     &         ZL_Predating_by_ZP - ZP_Excretion -
     &         ZP_Egestion - ZP_Mortality

	NO3_f = Nitrification  - 
     &          (PS_Photosynthesis - PS_Respiration)*RnewS -
     &          (PL_Photosynthesis - PL_Respiration)*RnewL
	NH4_f =  (ZS_Excretion + ZL_Excretion + ZP_Excretion) *1.0    +
     &           (DOM_Remineralization + POM_Remineralization)*1.0 - 
     &  Nitrification - 
     & (PS_Photosynthesis - PS_Respiration)*(1.0d0 - RnewS) -
     & (PL_Photosynthesis - PL_Respiration)*(1.0d0 - RnewL)
        	
		
	POM_f = PS_Mortality + PL_Mortality + ZS_Mortality +
     &          ZL_Mortality + ZP_Mortality + ZS_Egestion +
     &          ZL_Egestion + ZP_Egestion - POM_Remineralization - 
     &          POM_Decomposition_to_DOM

	DOM_f = PS_Extracellular_Excretion + PL_Extracellular_Excretion +
     &          POM_Decomposition_to_DOM - DOM_Remineralization 
	
	SIOH4_f = Opal_Decomposition - Si_PL_Shell_Formation


	Opal_f = Si_ZL_Egestion + Si_ZP_Egestion + Si_PL_Mortality -
     &           Opal_Decomposition



	
	
        Ca_f = CaCO3_Decomposition - Ca_PS_Shell_Formation - 
     &  Ca_ZS_Shell_Formation        


        CaCO3_f = Ca_ZS_Egestion + Ca_ZL_Egestion + Ca_ZP_Egestion +
     &  Ca_PS_Mortality + Ca_ZS_Mortality - CaCO3_Decomposition 
        

	 CO2_Air_Sea_Gas_Exchange = 0.0
        TCO2_f = ((NO3_f + NH4_f) * RCN) +  Ca_f
     &          + CO2_Air_Sea_Gas_Exchange 
       
       
        TALK_f = (2.0 * Ca_f) - NO3_f + NH4_f    
	


        
        !EXCHANGE AND EXPORT TERMS ie vertical advection.
        
        
        !        shapiro filter to get rid of numerical instability due to the
        !        advection problem
!        goto 105
        do kk = k-1, k+1  

        if (mod(loop,100).eq.0)then 

        if ( k.ge.4 .and. k.le.kmm ) then                               ! k gt 4 or 3 dbt????? 4 only we can
                                                                        !choose
        pomin_sm(kk) =(1.0/16.0) * ((-1*pomin(kk-2))+(4*pomin(kk-1))
     &          +(10*pomin(kk))+(4*pomin(kk+1))+(-1*pomin(kk+2)))
        opalin_sm(kk)=(1.0/16.0) * ((-1*opalin(kk-2))+(4*opalin(kk-1))
     &        +(10*opalin(kk))+(4*opalin(kk+1))+(-1*opalin(kk+2)))

        caco3in_sm(kk)=(1.0/16.0)*((-1*caco3in(kk-2))+(4*caco3in(kk-1))
     &     +(10*caco3in(kk))+(4*caco3in(kk+1))+(-1*caco3in(kk+2)))

        no3in_sm(kk) =(1.0/16.0) *((-1*no3in(kk-2))+(4*no3in(kk-1))
     &     +(10*no3in(kk))+(4*no3in(kk+1))+(-1*no3in(kk+2)))

        sioh4in_sm(kk) =(1.0/16.0)*((-1*sioh4in(kk-2))+(4*sioh4in(kk-1))
     &    +(10*sioh4in(kk))+(4*sioh4in(kk+1))+(-1*sioh4in(kk+2))) 
        
        tempin_sm(kk)=(1.0/16.0) * ((-1*tempin(kk-2))+(4*tempin(kk-1))
     &        +(10*tempin(kk))+(4*tempin(kk+1))+(-1*tempin(kk+2)))
        
        saltin_sm(kk)=(1.0/16.0) * ((-1*saltin(kk-2))+(4*saltin(kk-1))
     &        +(10*saltin(kk))+(4*saltin(kk+1))+(-1*saltin(kk+2)))
        
        uvelin_sm(kk)=(1.0/16.0) * ((-1*uvelin(kk-2))+(4*uvelin(kk-1))
     &        +(10*uvelin(kk))+(4*uvelin(kk+1))+(-1*uvelin(kk+2)))
        
        vvelin_sm(kk)=(1.0/16.0) * ((-1*vvelin(kk-2))+(4*vvelin(kk-1))
     &        +(10*vvelin(kk))+(4*vvelin(kk+1))+(-1*vvelin(kk+2)))

                
!        pomin_sm(kk) = pomin after smoothing by shapiro filter , equations got
!        it from vinu sr. ask for the reference.

        pomin(kk)=pomin_sm(kk)
        opalin(kk)=opalin_sm(kk)
        caco3in(kk)=caco3in_sm(kk)
        no3in(kk) =no3in_sm(kk) 
        sioh4in(kk)= sioh4in_sm(kk)
        tempin(kk) = tempin_sm(kk) 
        saltin(kk) = saltin_sm(kk) 
        uvelin(kk) = uvelin_sm(kk) 
        vvelin(kk) = vvelin_sm(kk) 


        end if
        end if
        end do
                
        if (mod(loop,100).eq.0)then 
        if ( k.ge.4 .and. k.le.kmm ) then 

        POM = pomin_sm(k)
        Opal= opalin_sm(k)
        CaCO3 = caco3in_sm(k)
        NO3 = no3in_sm(k)
        SIOH4 = sioh4in_sm(k)       
        T0    = tempin_sm(k) 
        S0    = saltin_sm(k)
        UVEL0 = uvelin_sm(k)
        VVEL0 = vvelin_sm(k)
 
        end if
        end if
!105    continue
        
        if (k.eq.1) then

       ! w_diff = wvelocity(k)-wvelocity(k+1)
       ! wb = wvelocity(k)+(w_diff)              !wb -virtual layer
                                              !velocity
       ! w1 = ((wvelocity(k))+(wb))/2
        
                !***************GRAD CALCULATION ****************

        w1 = (wvelocity(k)+wvelocity(k))/2
        w2 = (wvelocity(k)+wvelocity(k+1))/2
        
       !w1 = min(1e-6,w1)
       !w2 = min(1e-6,w1)
        !w1 = 1e-5
        !w2 = 1e-5

        no31 = (no3in(k)+no3in(k))/2
        no32 = (no3in(k)+no3in(k+1))/2
        sioh41 = (sioh4in(k)+sioh4in(k))/2
        sioh42 = (sioh4in(k)+sioh4in(k+1))/2
        t01    = (tempin(k)+tempin(k))/2 
        t02    = (tempin(k)+tempin(k+1))/2
        s01    = (saltin(k)+saltin(k))/2                
        s02    = (saltin(k)+saltin(k+1))/2              
!        uvel01 = (uvelin(k)+uvelin(k))/2               
!        uvel02 = (uvelin(k)+uvelin(k+1))/2             
!        vvel01 = (vvelin(k)+vvelin(k))/2               
!        vvel02 = (vvelin(k)+vvelin(k+1))/2           
 
        pom_grad(k) = (((pomin(k)+pomin(k))/2) 
     &       -((pomin(k)+pomin(k+1))/2)) /(dz)
                                                !dz direction is positive
        
        opal_grad(k)=(((opalin(k)+opalin(k))/2)
     &       -((opalin(k)+opalin(k+1))/2))/(dz)
        
        caco3_grad(k)=(((caco3in(k)+caco3in(k))/2)
     &     -((caco3in(k)+caco3in(k+1))/2))/(dz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!        no3_grad(k) = ((w1*no31)- (w2*no32))/(dz)
!        sioh4_grad(k) = ((w1*sioh41) - (w2*sioh42))/(dz)
!        t0_grad(k)    = ((w1*t01) - (w2*t02))/(dz)
!        s0_grad(k)    = ((w1*s01) - (w2*s02))/(dz)
!        uvel0_grad(k) = ((w1*uvel01) - (w2*uvel02))/(dz)
!        vvel0_grad(k) = ((w1*vvel01) - (w2*vvel02))/(dz) 
        
        no3_grad(k) =  ((w1+w2)/2) * (no31 - no32)/(dz)
        sioh4_grad(k) = ((w1+w2)/2) * (sioh41 - sioh42)/(dz)
        t0_grad(k)    = ((w1+w2)/2) * (t01 - t02)/(dz)         
        s0_grad(k)    = ((w1+w2)/2) * (s01 - s02)/(dz) 
!        uvel0_grad(k) = ((w1+w2)/2) * (uvel01 - uvel02)/(dz)
!        vvel0_grad(k) = ((w1+w2)/2) * (vvel01 - vvel02)/(dz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                

        else if ((k.ge.2) .and. (k.lt.km)) then
        w1 = (wvelocity(k)+wvelocity(k-1))/2
        w2 = (wvelocity(k)+wvelocity(k+1))/2

        !w1 = min(1e-6,w1)
        !w2 = min(1e-6,w1)

        !w1=1e-5 !delete
        !w2=1e-5 !delete

        no31 = (no3in(k)+no3in(k))/2
        no31 = (no3in(k)+no3in(k-1))/2
        no32 = (no3in(k)+no3in(k+1))/2
        sioh41 = (sioh4in(k)+sioh4in(k-1))/2
        sioh42 = (sioh4in(k)+sioh4in(k+1))/2
        t01    = (tempin(k)+tempin(k-1))/2      
        t02    = (tempin(k)+tempin(k+1))/2      
        s01    = (saltin(k)+saltin(k-1))/2     
        s02    = (saltin(k)+saltin(k+1))/2     
!        uvel01 = (uvelin(k)+uvelin(k-1))/2    
!        uvel02 = (uvelin(k)+uvelin(k+1))/2    
!        vvel01 = (vvelin(k)+vvelin(k-1))/2    
!        vvel02 = (vvelin(k)+vvelin(k+1))/2     


        pom_grad(k) = (((pomin(k) + pomin(k-1))/2) 
     &      - ((pomin(k)+pomin(k+1))/2))/(dz)

        opal_grad(k) =(((opalin(k) + opalin(k-1))/2)
     &     - ((opalin(k)+opalin(k+1))/2))/(dz)
!
        caco3_grad(k) = (((caco3in(k) + caco3in(k-1))/2)
     &        -((caco3in(k)+caco3in(k+1))/2))/(dz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!        no3_grad(k) = ((w1*no31)- (w2*no32))/(dz)
!        sioh4_grad(k) = ((w1*sioh41) - (w2*sioh42))/(dz)
!        t0_grad(k)    = ((w1*t01) - (w2*t02))/(dz) 
!        s0_grad(k)    = ((w1*s01) - (w2*s02))/(dz)
!        uvel0_grad(k) = ((w1*uvel01) - (w2*uvel02))/(dz)
!        vvel0_grad(k) = ((w1*vvel01) - (w2*vvel02))/(dz) 

        no3_grad(k) =  ((w1+w2)/2) * (no31 - no32)/(dz)
        sioh4_grad(k) = ((w1+w2)/2) * (sioh41 - sioh42)/(dz)
        t0_grad(k)    = ((w1+w2)/2)  * (t01 - t02)/(dz)         
        s0_grad(k)    = ((w1+w2)/2) * (s01 - s02)/(dz)
!        uvel0_grad(k) = ((w1+w2)/2) * (uvel01 - uvel02)/(dz)
!        vvel0_grad(k) = ((w1+w2)/2) * (vvel01 - vvel02)/(dz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        else
        w1 = (wvelocity(k-1)+wvelocity(k))/2
        w2 = (wvelocity(k)+wvelocity(k))/2
        
        !w1 = min(1e-6,w1)
        !w2 = min(1e-6,w1)

        !w1 = 1e-5
        !w2 = 1e-5

        no31 = (no3in(k)+no3in(k))/2
        no31 = (no3in(k-1)+no3in(k))/2
        no32 = (no3in(k)+no3in(k))/2
        sioh41 = (sioh4in(k-1)+sioh4in(k))/2
        sioh42 = (sioh4in(k)+sioh4in(k))/2
        t01   = (tempin(k-1)+tempin(k))/2               
        t02   = (tempin(k)+tempin(k))/2                
        s01   = (saltin(k-1)+saltin(k))/2              
        s02   = (saltin(k)+saltin(k))/2                
!        uvel01= (uvelin(k-1)+uvelin(k))/2              
!        uvel02= (uvelin(k)+uvelin(k))/2           
!        vvel01= (vvelin(k-1)+vvelin(k))/2              
!        vvel02= (vvelin(k)+vvelin(k))/2              


        pom_grad(k) = (((pomin(k-1)+pomin(k))/2)
     &      -((pomin(k)+ pomin(k))/2))/(dz)
!
        opal_grad(k) = (((opalin(k-1)+opalin(k))/2)
     &         - ((opalin(k)+opalin(k))/2))/(dz)
!
        caco3_grad(k) = (((caco3in(k-1)+caco3in(k))/2)
     &    - ((caco3in(k) +caco3in(k))/2) )/(dz)

!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
!        no3_grad(k) = ((w1*no31)- (w2*no32))/(dz)
!        sioh4_grad(k) = ((w1*sioh41) - (w2*sioh42))/(dz)
!        t0_grad(k)    = ((w1*t01) - (w2*t02))/(dz) 
!        s0_grad(k)    = ((w1*s01) - (w2*s02))/(dz)
!        uvel0_grad(k) = ((w1*uvel01) - (w2*uvel02))/(dz)
!        vvel0_grad(k) = ((w1*vvel01) - (w2*vvel02))/(dz) 

        no3_grad(k) =  ((w1+w2)/2) * (no31 - no32)/(dz)
        sioh4_grad(k) = ((w1+w2)/2) * (sioh41 - sioh42)/(dz)
        t0_grad(k)    = ((w1+w2)/2)  * (t01 - t02)/(dz)         
        s0_grad(k)    = ((w1+w2)/2) * (s01 - s02)/(dz)
!        uvel0_grad(k) = ((w1+w2)/2) * (uvel01 - uvel02)/(dz)
!        vvel0_grad(k) = ((w1+w2)/2) * (vvel01 - vvel02)/(dz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end if
!**************instead of this a profile has to give for export****    
        ExpPON   = SPOM  * pom_grad(k)                          !DOUBT, y ekman_vel is not here
        ExpOpal  = SOpal * opal_grad(k)                         !DOUBT, " "
        ExpCaCO3 = SCaCO3  * caco3_grad(k)                      !DOUBT, " "
!*****************************************************************************
        
        ExcNO3   = no3_grad(k)                                  !(ExcvNO3) * no3_grad(k)!instead of
                                                                !Exchange velocity ekman velocity is given here. no3_grad(k) is not
                                                                !simply a gradient term it is gradient * wvelocity term.same in
                                                                !sioh4_grad(k)
        ExcSiOH4 = sioh4_grad(k)!(ExcvSIOH4) * sioh4_grad(k)
        ExcT0    = t0_grad(k)
        ExcS0    = s0_grad(k)
!        ExcUVEL0 = uvel0_grad(k)
!        ExcVVEL0 = vvel0_grad(k)

        
! 	 ExpPON   =   1.0*SPOM  /100.0 * POM
!         ExpOpal  =   1.0*SOpal / 100.0 * Opal
!         ExpCaCO3 =   1.0*SCaCO3 /100.0 * CaCO3 !added by anju. check sigh
!**********************************************************************************************************
       !Giving exchange in the model 
!         ExcNO3   =  ExcTime * ( TNO3d   - NO3   )
!         ExcSiOH4 = ExcTime * ( TSiOH4d - SIOH4 )
!!**********************************************************************************************************
!       go to 101 !PROFILE METHOD
         Z_c = 55.0  !meter
         kzc = 51    !11         
         kzcp = kzc + 1
        
         if(k .le. kzc) then
          !POM_ave = (pomin(k) + pomin(k+1)) /2 !NOW
         grad_POM = pomin(k)/dz!(pomin(k) - pomin(k+1))/dz
         grad_Opal = opalin(k)/dz!(opalin(k) - opalin(k+1))/dz
         grad_CaCO3 = caco3in(k)/dz!(caco3in(k) - caco3in(k+1))/dz
         
        
         
         POM_f = POM_f - (grad_POM * SPOM)
         Opal_f = Opal_f - (grad_Opal * SOpal )
!         CaCO3_f = CaCO3_f -  (grad_CaCO3 * SCaCO3)



       !POM_f = POM_f - (POM/dz * SPOM)
        !POM_f = POM_f - (POM_ave * SPOM)


! if (POM_f .gt. 0.0)POM_f = POM_f !POM_f * 0.67 !DIS NOW

         POM_F_depth = POM_F_depth  + max(0.0, POM) * dz !(max(0.0,POM_f*(1-0.67)) * dz)
         Opal_F_depth = Opal_F_depth + max(0.0,Opal) * dz

!         CaCO3_F_depth = CaCO3_F_depth + (CaCO3_f * dz)
        ! POM_F_depth,Opal_F_depth,CaCO3_F_depth - flux of
        ! POM,Opal,CaCO3 at depth z=75
         POM_F_Z(k)   = 0.0
         Opal_F_Z(k)  = 0.0
!         CaCO3_F_Z(k) = 0.0
         endif

         POM_F_depth = max(0.0, POM_F_depth)
         Opal_F_depth = max(0.0, Opal_F_depth)
        
         if((k.ge.kzc).and.(k.le.km)) then
         POM_F_Z(k) = POM_F_depth *  ((k*dz)/Z_c)**(-0.9) !a=0.9 ref the papr for

!this equation  
         Opal_F_Z(k) = Opal_F_depth *  ((k*dz)/Z_c)**(-0.9) !a=0.9 ref the papr for
!         Opal_F_Z(k) =  r_Si * Opal_F_depth * exp(-((k*dz)-(Z_c))/d) 
c         CaCO3_F_Z(k) = r_Cp * CaCO3_F_depth * exp(-((k*dz)-(Z_c))/d) 
                                                                        !R-rain ratio
                                                                        !r_CP - C to Phosphorus Redfield ratio.
                                                                        !d = scale depth
         !POM_F_Z,Opal_F_Z,CaCO3_F_Z = flux of pom opal and caco3 at a
                                                                        !each depth
         end if
!         if(k.eq.kzcp) then
!         POM_f = POM_f - ((POM_F_Z(k-1)-POM_F_Z(k))/dz)/(86400*10)
!         Opal_f = Opal_f - ((Opal_F_Z(k-1)-Opal_F_Z(k))/dz)/(86400*10)
c         CaCO3_f = CaCO3_f - ((CaCo3_F_Z(k)-CaCO3_F_Z(k))/dz) 
!         end if
         
         if(k.ge.kzcp) then
         POM_f = POM_f + ((POM_F_Z(k-1)-POM_F_Z(k))/dz)
         Opal_f = Opal_f + ((Opal_F_Z(k-1)-Opal_F_Z(k))/dz)
!         Opal_f = Opal_f - ((Opal_F_Z(k-1)-Opal_F_Z(k))/dz)/(86400*10)
        
c         CaCO3_f = CaCO3_f - ((CaCO3_F_Z(k)-CaCO3_F_Z(k-1))/dz)
         end if

!101     continue

       !exchange 
      ! ExpPON = -1.0 * POM *(40/86400)   !NOW
         if (k.le.8) then
         relaxation_NO3 =(NO3 - no3_obss)/86400.0/15.0             !180.0         
         relaxation_SIOH4=(SIOH4 - sioh4_obss)/86400.0/15.0        !180.0   
         end if

!!!!CHECKKKKKKKKKK
!ABOVEEEEEEEEEEEEEEEEEEEEEEEE***********************************
         if ((k.gt.8).and.(k.le.16)) then
         relaxation_NO3 = 5*(NO3 - no3_obss)/86400.0/90.0        !90.0
         relaxation_SIOH4 = 5*(SIOH4 - sioh4_obss)/86400.0/90.0  !90.0
         end if
         if (k.gt.16) then
         relaxation_NO3 = 5*(NO3 - no3_obss)/86400.0/30.0        !30.0
         relaxation_SIOH4 =5*(SIOH4 - sioh4_obss)/86400.0/30.0  !30.0
         end if 


         NO3_f   = NO3_f - ExcNO3  - relaxation_NO3*1.0 !1e-7*(no3in(k-1)- no3in(k+1))/ 2/dz ! 1e-5*(no3in(k-1)- no3in(k+1))/ 2/dz !+ SPOM * (NO3/dz)! 1.0 !NOW
         SIOH4_f = SIOH4_f - ExcSiOH4 -  relaxation_SIOH4*1.0
!*****************************************************
         POM_f   = POM_f   - ExpPON *0.0 !here -ve came because SPOM is negative
         Opal_f  = Opal_f  - ExpOpal *0.0
	 CaCO3_f = CaCO3_f - ExpCaCO3 *0.0!1.0
        
        !print*,loop,k,ExcNO3,ExpPON
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
        T0 = T0 - (ExcT0 * dt)                 !doubt here.
        S0 = S0 - (ExcS0 * dt)                 !doubt here

!        UVEL0 = UVEL0 - (ExcUVEL0 * dt)         
!        VVEL0 = VVEL0 - (ExcVVEL0 * dt)         

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,DELETE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!        T0 = T0 + (ExcT0 * dt)                         !doubt here.
!        S0 = S0 + (ExcS0 * dt)                         !doubt here
!        UVEL0 = UVEL0 + (ExcUVEL0 * dt)                !doubt here
!        VVEL0 = VVEL0 + (ExcVVEL0 * dt)                !doubt here
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	



	return
	end


	function Monod (a, b)
	real(8) a, b, Monod
	Monod = a/(a+b)
	end

	function expPsi (a, b)
	real(8) expPsi, a, b
	expPsi = exp(-1.0*a*b)
	end


