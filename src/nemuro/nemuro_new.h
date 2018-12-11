ccccccccccccccccccccccccccccc nemuro.h cccccccccccccccccccccccccccccccccccccccc
c
c	Header informations for nemuro_3D.F
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Define local variables for each Tracer	c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) PS, PL, ZS, ZL, ZP
	real(8) NO3, NH4, POM, DOM, SIOH4, Opal, CaCO3
	real(8) Ca, TCO2, TALK    
	real(8) T0, S0, UVEL0, VVEL0

	real(8) ExpPON, ExpOpal, ExcNO3, ExcSiOH4
	real(8) ExpCaCO3,ExcT0,ExcS0,ExcUVEL0,ExcVVEL0

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Small Phyto Plankton, Tracer = 3	c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) PS_Photosynthesis
	real(8) PS_Respiration
	real(8) PS_Extracellular_Excretion
	real(8) PS_Mortality
	real(8) PS_Grazing_by_ZS
	real(8) PS_Grazing_by_ZL

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Large Phyto Plankton, Tracer = 4	c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) PL_Photosynthesis
	real(8) PL_Respiration
	real(8) PL_Extracellular_Excretion
	real(8) PL_Mortality
	real(8) PL_Grazing_by_ZL
	real(8) PL_Grazing_by_ZP

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Small Zoo Plankton, Tracer = 5		c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) ZS_Excretion
	real(8) ZS_Egestion
	real(8) ZS_Mortality
	real(8) ZS_Predating_by_ZL
	real(8) ZS_Predating_by_ZP

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Large Zoo Plankton, Tracer = 6		c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) ZL_Excretion
	real(8) ZL_Egestion
	real(8) ZL_Mortality
	real(8) ZL_Predating_by_ZP

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Predatory Zoo Pankton, Tracer = 7	c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) ZP_Excretion
	real(8) ZP_Egestion
	real(8) ZP_Mortality

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Nitrate (NO3), Tracer = 8		c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) Nitrification
	real(8) RnewS
	real(8) RnewL

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Amonia (NH4), Tracer = 9		c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) DOM_Remineralization
	real(8) POM_Remineralization

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Purt. Org. Matt (POM), Tracer = 10	c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) POM_Decomposition_to_DOM
	real(8) Sinking_POM 

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Diss. Org. Matt. (DOM), Tracer = 11	c
ccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Silicate, (Si(OH)4), Tracer = 12	c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) Opal_Decomposition
	real(8)Si_PL_Shell_Formation

ccccccccccccccccccccccccccccccccccccccccccccccccc
c	Opal, Tracer = 13			c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) Si_ZL_Egestion
        real(8) Si_ZP_Egestion
	real(8) Si_PL_Mortality
	real(8) Sinking_Opal

ccccccccccccccccccccccccccccccccccccccccccccccccc  
c        Ca, Tracer=14                          c
ccccccccccccccccccccccccccccccccccccccccccccccccc
        real(8) CaCO3_Decomposition
	real(8) Ca_PS_Shell_Formation
	real(8) Ca_ZS_Shell_Formation


ccccccccccccccccccccccccccccccccccccccccccccccccc 
c	CaCO3, Tracer=15                        c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) Ca_ZS_Egestion 
	real(8) Ca_ZL_Egestion
	real(8) Ca_ZP_Egestion
	real(8) Ca_PS_Mortality
	real(8) Ca_ZS_Mortality
	real(8) Sinking_CaCO3


ccccccccccccccccccccccccccccccccccccccccccccccccc !added by ANJU/22-06-2016 -it is not complete.. 
c      CO2, Tracer=16                           c
ccccccccccccccccccccccccccccccccccccccccccccccccc
	real(8) CO2_Air_Sea_Gas_Exchange


cccccccccccccccccccccccccccccccccccccccccccccccccc !added by ANJU/22-06-2016 -it is not complete.. check in the paper and complete it.
c      Alkalinity, Tracer=17                     c
cccccccccccccccccccccccccccccccccccccccccccccccccc
	



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Biological Parameters.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	real(8), parameter :: D2S	=	86400.0		! days in seconds                                !refer murthus paper!!!
	real(8), parameter :: PSmax	=	200.0		!  Max C:Chla ratio for PS, [g:g]             
	real(8), parameter :: PLmax	=	120.0		! Max C:Chla ratio for PL, [g:g]            
	real(8), parameter :: KPPS	=	95.0*D2S		! Photoassimilation coefficient for PS [(g:g)day] 
	real(8), parameter :: KPPL 	=	70.0*D2S		! photoassimilation coefficient for PL [(g:g)day] !day or secnd refer above paper
	real(8), parameter :: PSmin 	=	30.0		! min C:Chla ratio for PS, [g:g]
	real(8), parameter :: PLmin 	=	15.0		! min C:Chla ratio for PL, [g:g]

	real(8), parameter :: VmaxS	=	0.3/D2S   !0.4/D2S         ! PS Maximum Photosynthetic rate @0degC   [/s] !delete
	real(8), parameter :: KNO3S	=	1.0e-6    !1.0e-6	        ! PS Half satuation constant for Nitrate  [molN/l] 
	real(8), parameter :: KNH4S	=	0.1e-6    !0.1e-6		! PS Half satuation constant for Ammonium [molN/l] !
	real(8), parameter :: PsiS	=	1.5/1e-6	! PS Ammonium Inhibition Coefficient      [l/molN]         !
	real(8), parameter :: KS	=	0.0693		! PS Temp. Coeff. for Photosynthetic Rate [/degC]  !
	real(8), parameter :: IoptS	=	104.7 		! PS Optimum Light Intensity              [W/m2]    !
	real(8), parameter :: MPS0	=       0.0585/1e-6/D2S	! PS Mortality Rate @0degC                [l/mol/s]
	real(8), parameter :: KMPS	=	0.0693		! PS Temp. Coeff. for Mortality           [/degC]
	real(8), parameter :: RPS0	=	0.03/D2S	! PS Respiration Rate at @0degC           [/s]
	real(8), parameter :: KRS	=	0.0519		! PS Temp. Coeff. for Respiration         [/degC]
	real(8), parameter :: GammaS	=	0.135		! PS Ratio of Extracell. Excret. to Photo.[(nodim)]

	real(8), parameter :: VmaxL	=	2.0/D2S   !0.8/D2S	        ! PL Maximum Photosynthetic rate @0degC   [/s] !delet
	real(8), parameter :: KNO3L	=	3.0e-6    !3.0e-6		! PL Half satuation constant for Nitrate  [molN/l]
	real(8), parameter :: KNH4L	=	0.3e-6    !0.3e-6		! PL Half satuation constant for Ammonium [molN/l]
	real(8), parameter :: KSIL	=	6.0e-6    !6.0e-6		! PL Half satuation constant for Silicate [molSi/l]
	real(8), parameter :: PsiL	=	1.5/1e-6	! PL Ammonium Inhibition Coefficient      [l/molN]
	real(8), parameter :: KL	=	0.0693		! PL Temp. Coeff. for Photosynthetic Rate [/degC]
	real(8), parameter :: IoptL	=	104.7 		! PL Optimum Light Intensity              [W/m2]
	real(8), parameter :: MPL0	=	0.029/1e-6/D2S  ! PL Mortality Rate @0degC                [l/mol/s]
	real(8), parameter :: KMPL	=	0.0693		! PL Temp. Coeff. for Mortality           [/degC]
	real(8), parameter :: RPL0	=	0.0300/D2S	! PL Respiration Rate at @0degC           [/s]
	real(8), parameter :: KRL	=	0.0519		! PL Temp. Coeff. for Respiration         [/degC]
	real(8), parameter :: GammaL	=	0.135		! PL Ratio of Extracell. Excret. to Photo.[(nodim)]
	
	real(8), parameter :: GRmaxS	=	0.4/D2S		! ZS Maximum Rate of Grazing PS @0degC    [/s]
	real(8), parameter :: KGS	=	0.0693		! ZS Temp. Coeff. for Grazing             [/degC]
	real(8), parameter :: LamdaS	=	1.4/1e-6	! ZS Ivlev constant                       [l/molN]
	real(8), parameter :: PSZS	=	0.043e-6	! ZS Threshold Value for Grazing PS       [molN/l]
	real(8), parameter :: GammaZS	=	0.7		! ZS Assimilation Efficiency              [(nodim)]
	real(8), parameter :: BettaZS	=	0.3		! ZS Growth Efficiency                    [(nodim)]
	real(8), parameter :: MZS0	=	0.0585/1e-6/D2S ! ZS Mortality Rate @0degC                [l/mol/s]
	real(8), parameter :: KZS	=	0.0693		! ZS Temp. Coeff. for Mortality           [/degC]

	real(8), parameter :: GRmaxLPS	=	0.1/D2S		! ZL Maximum Rate of Grazing PS @0degC    [/s]
	real(8), parameter :: GRmaxLPL	=	0.4/D2S		! ZL Maximum Rate of Grazing PL @0degC    [/s]
	real(8), parameter :: GRmaxLZS	=	0.4/D2S		! ZL Maximum Rate of Grazing ZS @0degC    [/s]
	real(8), parameter :: KGL	=	0.0693		! ZL Temp. Coeff. for Grazing             [/degC]
	real(8), parameter :: LamdaL	=	1.4/1e-6	! ZL Ivlev constant                       [l/molN]
	real(8), parameter :: PSZL	=	0.04e-6		! ZL Threshold Value for Grazing PS       [molN/l]
	real(8), parameter :: PLZL	=	0.04e-6		! ZL Threshold Value for Grazing PL       [molN/l]
	real(8), parameter :: ZSZL	=	0.04e-6		! ZL Threshold Value for Grazing ZS       [molN/l]
	real(8), parameter :: GammaZL	=	0.7		! ZL Assimilation Efficiency              [(nodim)]
	real(8), parameter :: BettaZL	=	0.3		! ZL Growth Efficiency                    [(nodim)]
	real(8), parameter :: MZL0	=	0.0585/1e-6/D2S ! ZL Mortality Rate @0degC                [l/mol/s]
	real(8), parameter :: KZL	=	0.0693		! ZL Temp. Coeff. for Mortality		  [/degC]
	
	real(8), parameter :: GRmaxPPL	=	0.4/D2S !was 0.1	! ZP Maximum Rate of Grazing PS @0degC, !"Grazing PL"(in kishi study.)[/s] !it is not frm the paper. from A7station or PAPA station
	real(8), parameter :: GRmaxPZS	=	0.4/D2S	!was 0.2 ZP Maximum Rate of Grazing PL @0degC !"Predating ZS" (in kishi) [/s] !""
	real(8), parameter :: GRmaxPZL	=	0.4/D2S	!was 0.2 ZP Maximum Rate of Grazing ZS @0degC !"predatitng ZL"(in kishi) [/s] !""
	real(8), parameter :: KGP	=	0.0693	! ZP Temp. Coeff. for Grazing !/predation   [/degC]!""
	real(8), parameter :: LamdaP	=	1.5/1e-6! ZP Ivlev constant                       [l/molN]!""
	real(8), parameter :: PLZP	=	0.04e-6	! ZP Threshold Value for Grazing PL       [molN/l]!"", but it was PSZP, we dont have that why???
	real(8), parameter :: ZSZP	=	0.04e-6	! ZP Threshold Value for Grazing ZS       [molN/l]
	real(8), parameter :: ZLZP	=	0.04e-6	! ZP Threshold Value for Grazing ZL       [molN/l]
	real(8), parameter :: GammaZP	=	0.7	! ZP Assimilation Efficiency              [(nodim)]
	real(8), parameter :: BettaZP	=	0.3	! ZP Growth Efficiency                    [(nodim)]
	real(8), parameter :: MZP0	=	0.0500/1e-6/D2S ! ZP Mortality Rate @0degC                [/s] !l/mol/s unit was error
	real(8), parameter :: KZP	=	0.0693		! ZP Temp. Coeff. for Mortality		  [/degC]
	real(8), parameter :: PusaiPL   =       4.605e6    	! ZP Preference Coeff. for PL             [l/molN]
	real(8), parameter :: PusaiZS   =       3.010e6 	! ZP Preference Coeff. for ZS             [l/molN]

	real(8), parameter :: Alpha1	=	0.04		! Light Dissip Coeff. of sea water	  [/m]
	real(8), parameter :: Alpha2	=	0.04/1e-6	! Self Shading Coeff.			  [1/mol m]
	
	real(8), parameter :: NNit0	=	0.03/D2S	! NH4 Nitrification Rate @0degC           [/s]
	real(8), parameter :: KNit	=	0.0693		! NH4 Temp. coefficient for Nitrification [/degC]

	real(8), parameter :: SPOM	=	40.0/D2S	! Part. Org. Matt. Sinking velocity	  [m/sec]
	real(8), parameter :: VPA0	=	0.10/D2S	! PON Decomp. Rate to Ammonium @0degC     [/s]         !
	real(8), parameter :: KPA	=	0.0693		! PON Temp. Coeff. for Decomp. to Ammon.  [/degC]
	real(8), parameter :: VPD0	=	0.10/D2S	! PON Decomp. Rate to DON @0degC          [/s]
	real(8), parameter :: KPD	=	0.0693		! PON Temp. Coeff. for Decomp. to DON     [/degC]
	real(8), parameter :: VDA0	=	0.2/D2S      !0.02/D2S	! DON Decomp. Rate to Ammonium @0degC     [/s] !
	real(8), parameter :: KDA	=	0.0693		! DON Temp. Coeff. for Decomp. to Ammon.  [/degC]
	real(8), parameter :: SOpal	=	40.0/D2S!(100m/day in org. nemuro)	! Sinking Velocity of Opal		  [m/sec]
	real(8), parameter :: VOpal	=	0.10/D2S	! Opal Decomp. Rate to Silicate @0degC    [/s]
	real(8), parameter :: KOpal	=	0.0693		! Opal Temp. Coeff. for Decomp.to Silicate[/degC]
	real(8), parameter :: ExcvNO3   =       1e-6            ! exchange velocity of no3                [m/s] !by anju doubt in it
	real(8), parameter :: ExcvSIOH4 =       1e-6            ! exchange velocity of sioh4              [m/s] !by anju doubt in it

	real(8), parameter :: RSiN	=	2.0		! Si/N ratio				  [Nodim]

	real(8), parameter :: ExcTime   = 	3.3d0 / (100.0d0*D2S) ! Exch. Coeff. between Sur-Deep 	  [/s]
	real(8), parameter :: TNO3d	= 	10.0d-6    	! Nitrate Concentraion in the Deep Layer  [molN/l]
	real(8), parameter :: TSiOH4d	=	10.0d-6    	! Silicate Concentraion in the Deep Layer [molSi/l]
	
	real(8), parameter :: SCaCO3    =       40.0/D2S           ! Calcium Carbonate Sinking Velocity   [m/s]   added by ANJU/22-06-2016
        real(8), parameter :: VCaCO3    =       0.050/D2S           ! Decomposition rate of CaCO3 at 0 degC [/s]     ! added by ANJU/22-06-2016
	real(8), parameter :: kCaCO3    =       0.0693          !Tem. Coefficient for CaCO3 decomposition [/degC]  ! added by ANJU/22-06-2016

	real(8), parameter :: RCN       =       6.625           !Stoichiomery of carbon to Nitrogen       [Nodim]  ! added by ANJU/22-06-2016
	
	real(8), parameter :: Rcoco     =       0.1  !Ratio of Cocolithophoroids to small phytoplankton   [Nodim]  ! added by ANJU/22-06-2016
	real(8), parameter :: Rfora     =       0.1  !Ratio of foraminifera to small phytoplankton        [Nodim]  ! added by ANJU/22-06-2016
	real(8), parameter :: RCco      =       0.5  !Ratio of Inorganic C to total C in Cocolithophoroids[Nodim]  ! added by ANJU/22-06-2016
	real(8), parameter :: RCfo      =       0.5  !Ratio of Inorganic C to total C in Foraminifera     [Nodim]  ! added by ANJU/22-06-2016

