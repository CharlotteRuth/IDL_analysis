Function EdotInstant,T,rho,abunds,R_Cool_Metal,R_Heat_Metal,Zmetal 
NOCOMPTON = 1
MOLECULARH = 0
NOMETALCOOLING = 1
SHIELDHI = 1
DOLOWTCOOL = 0
CL_B_gm = (6.022e23*(938.7830/931.494))
zIn = 0
CL_k_Boltzmann = 1.38066e-16
CL_eV_erg = 1.60219e-12
CL_eV_per_K = (CL_k_Boltzmann/CL_eV_erg)
CL_eH2 = (4.476*CL_eV_erg) ;/*Energy lost during H2 dissociation, Shapira & Kang (1987) through Abel 1997 */
;/*#define CL_eH2     (14.5*CL_eV_erg) Dissasociation energy for H2(?)CC  Gould and Salpeter 1963*/
CL_eHI = (13.60*CL_eV_erg)
CL_eHeI = (24.59*CL_eV_erg)
CL_eHeII = (54.42*CL_eV_erg)
CL_E2HeII = (3.0*13.6*CL_eV_erg)
CL_Ccomp0 = 0.565e-9 
CL_Tcmb0 = 2.735
CL_Ccomp = (CL_Ccomp0*CL_Tcmb0)


Y_HI = REFORM(abunds[0,*])
Y_H2 = REFORM(abunds[1,*]) 
Y_HII = REFORM(abunds[2,*]) 
Y_HeI = REFORM(abunds[3,*])
Y_HeII = REFORM(abunds[4,*])
Y_HeIII = REFORM(abunds[5,*]) 
Y_e = REFORM(abunds[6,*]) 
s_dust = REFORM(abunds[7,*])
s_self = REFORM(abunds[8,*])

en_B = rho*CL_B_gm
n_e = en_B*Y_e

R_Radr_HII = clRateRadrHII(T) 
R_Phot_HI = 9.79e-14 ;1d-100
R_Coll_HI = clRateCollHI(T)
R_TotrForm_H2 = clRateTotrFormH2(ZMetal)
R_Coll_e_H2 = clRateColl_e_H2(T)
R_Coll_H_H2 = clRateColl_H_H2(T)
R_Coll_H2_H2 = clRateColl_H2_H2(T)
R_Phot_H2 = 1.5917e-11

R_Phot_HeI = 2.850000e-14 ;1d-100
R_Radr_HeII = clRateRadrHeII(T) 
R_Diel_HeII = clRateDielHeII(T)
R_Phot_HeII = 2.940000e-16 ;1d-100
R_Radr_HeIII =clRateRadrHeII(T)
R_Coll_HeI = clRateCollHeI(T) 
R_Coll_HeII = clRateCollHeII(T) 

R_Cool_Coll_HI = CL_eHI*CL_B_gm; 
R_Cool_Coll_HeI = CL_eHeI*CL_B_gm;
R_Cool_Coll_HeII = CL_eHeII*CL_B_gm;
R_Cool_Diel_HeII = (CL_E2HeII+CL_eHeI)*CL_B_gm;
R_Cool_Coll_H2 = CL_eH2*CL_B_gm;
R_Cool_Comp = ((1+zIn)*CL_Ccomp)^4.0*CL_B_gm 
R_Tcmb = CL_Tcmb0*(1+zIn)
R_Heat_Phot_HI = 3134276950407.8179
R_Heat_Phot_H2 = 9660222337904.6992
R_Heat_Phot_HeI = 5306780052799.2021
R_Heat_Phot_HeII = 24606650403217.199

IF (NOCOMPTON)  THEN  Edot_cool = 0 ELSE Edot_cool =  Y_e*R_Cool_Comp*(T - R_Tcmb) ; 
Edot_cool = Edot_cool + n_e*(clCoolBrem1(T)*(Y_HII + Y_HeII) + $
                             clCoolBrem2(T)*Y_HeIII + $	  
                             R_Cool_Diel_HeII*Y_HeII*R_Diel_HeII + $
                             clCoolRadrHII(T)*Y_HII*R_Radr_HII  + $ 
                             clCoolRadrHeII(T)*Y_HeII*R_Radr_HeII + $
                             clCoolRadrHeIII(T)*Y_HeIII*R_Radr_HeIII + $
                             clCoolLineHI(T)*Y_HI + $
                             clCoolLineHeI(T)*Y_HeI + $
                             clCoolLineHeII(T)*Y_HeII)

IF (MOLECULARH) THEN Edot_cool =  Edot_cool + n_e*(clCoolLineH2(T, Y_HI + Y_HII + 2.0*Y_H2)*Y_H2 +  $
                                                    R_Cool_Coll_H2*Y_H2*R_Coll_e_H2*s_dust*s_self)   

Edot_cool = Edot_cool + n_e*(R_Cool_Coll_HI*Y_HI*R_Coll_HI + $
                             R_Cool_Coll_HeI*Y_HeI*R_Coll_HeI +  $
                             R_Cool_Coll_HeII*Y_HeII*R_Coll_HeII)

IF (MOLECULARH) THEN Edot_cool = Edot_cool + R_Cool_Coll_H2*Y_H2*R_Coll_H_H2*Y_HI*en_B*s_dust*s_self $ 
                                           + R_Cool_Coll_H2*Y_H2*R_Coll_H2_H2*Y_H2*en_B*s_dust*s_self

IF (DOLOWTCOOL)         THEN Edot_cool = Edot_cool  + LowTCool
IF NOT (NOMETALCOOLING) THEN Edot_cool = Edot_cool + R_Cool_Metal

Edot_heat = 0
IF NOT (NOMETALCOOLING) THEN Edot_heat = Edot_heat + R_Heat_Metal
IF (MOLECULARH)         THEN Edot_heat = Edot_heat + Y_H2*R_Heat_Phot_H2*R_Phot_H2*s_dust*s_self
IF (SHIELDHI)           THEN Edot_heat = Edot_heat + Y_HI*R_Heat_Phot_HI*R_Phot_HI*s_dust $
                        ELSE Edot_heat = Edot_heat + Y_HI * R_Heat_Phot_HI * R_Phot_HI 

Edot_heat = Edot_heat + Y_HeI*R_Heat_Phot_HeI*R_Phot_HeI + Y_HeII*R_Heat_Phot_HeII*R_Phot_HeII ;

Edot = Edot_heat - Edot_cool

Edot_phot = 0
IF (MOLECULARH) THEN Edot_phot = Edot_phot + Y_H2*R_Heat_Phot_H2*R_Phot_H2*s_dust*s_self
IF (SHIELDHI)   THEN Edot_phot = Edot_phot + Y_HI*R_Heat_Phot_HI*R_Phot_HI*s_dust $
                ELSE Edot_phot = Edot_phot + Y_HI*R_Heat_Phot_HI*R_Phot_HI 
Edot_phot = Edot_phot + Y_HeI*R_Heat_Phot_HeI*R_Phot_HeI + Y_HeII*R_Heat_Phot_HeII*R_Phot_HeII ;
Edot_metal = R_Heat_Metal - R_Cool_Metal
Edot_LineHI = n_e*clCoolLineHI(T)*Y_HI
Edot_CollHI = n_e*R_Cool_Coll_HI*Y_HI*R_Coll_HI
Edot_Brem = n_e*(clCoolBrem1(T)*(Y_HII + Y_HeII) + clCoolBrem2(T)*Y_HeIII) 
Edot_Rad =  n_e*(clCoolRadrHII(T)*Y_HII*R_Radr_HII + clCoolRadrHeII(T)*Y_HeII*R_Radr_HeII + clCoolRadrHeIII(T)*Y_HeIII*R_Radr_HeIII)
;IF (3 lt alog10(T) AND  alog10(T)lt 3.5) then stop
RETURN,[EDOT,Edot_metal,Edot_phot,Edot_heat,Edot_cool,Edot_LineHI,Edot_CollHI,Edot_Brem,Edot_Rad,clCoolBrem1(T)]
END


FUNCTION clCoolBrem1,T
CL_Cbremss1 = 1.426d-27
CL_al = 0.79464
CL_bl = 0.1243
CL_ar = 2.13164
CL_br = -0.1240
CL_B_gm = (6.022d23*(938.7830/931.494))

Tlog10 = alog10(T)              ;
Tsq = sqrt(T)                   ;
if (T le 3.2e5) THEN return,Tsq*CL_Cbremss1*(CL_al+CL_bl*Tlog10)*CL_B_gm $
                ELSE return,Tsq*CL_Cbremss1*(CL_ar+CL_br*Tlog10)*CL_B_gm ;
END

FUNCTION clCoolBrem2,T
CL_al = 0.79464
CL_bl = 0.1243
CL_ar = 2.13164
CL_br = -0.1240
CL_alog4  = 0.602059991
CL_alII = (4.0*(CL_al-CL_bl*CL_alog4))
CL_blII = (4.0*CL_bl)
CL_arII = (4.0*(CL_ar-CL_br*CL_alog4))
CL_brII = (4.0*CL_br)
CL_Cbremss1 = 1.426d-27
CL_B_gm = (6.022e23*(938.7830/931.494))

Tlog10 = alog10(T)               ;
Tsq = sqrt(T)                   ;
if (T le 12.8e5) THEN return,Tsq*CL_Cbremss1*(CL_alII+CL_blII*Tlog10)*CL_B_gm $
                 ELSE return,Tsq*CL_Cbremss1*(CL_arII+CL_brII*Tlog10)*CL_B_gm;
END

FUNCTION clCoolRadrHII,T
CL_aHII =  0.0215964
CL_b = 0.270251
CL_k_Boltzmann = 1.38066e-16
CL_B_gm = (6.022e23*(938.7830/931.494))

Tpow = T^CL_b                ;
;  /* return CL_B_gm*(CL_eHI+exp(-CL_aHII*Tpow)*CL_k_Boltzmann*T); */
;  /* Though 13.6eV is lost to the Gas as radiation, calculating the
;   * Energy using u = 3/2 k T requires we don't subtract it here.
;   */
return,CL_B_gm*(exp(-1.0*CL_aHII*Tpow)*CL_k_Boltzmann*T) ;
END
 
FUNCTION clCoolRadrHeII,T
CL_b = 0.270251
CL_B_gm = (6.022e23*(938.7830/931.494))
CL_eV_erg = 1.60219e-12
CL_eHI = (13.60*CL_eV_erg)
CL_aHII =  0.0215964
CL_k_Boltzmann = 1.38066e-16

Tpow = T^CL_b                ;
;  /* return CL_B_gm*(CL_eHeI+exp(-(CL_aHII*pow(13.6/24.59,CL_b))*Tpow)*CL_k_Boltzmann*T); */
return,CL_B_gm*(exp(-1.0*(CL_aHII*(13.6/24.59)^CL_b)*Tpow)*CL_k_Boltzmann*T) ;
END

FUNCTION clCoolRadrHeIII,T 
CL_B_gm = (6.022e23*(938.7830/931.494))
CL_aHII =  0.0215964
CL_b = 0.270251
CL_k_Boltzmann = 1.38066e-16
    
Tpow=T^CL_b                     ;
;  /* return CL_B_gm*(CL_eHeII+exp(-(CL_aHII*pow(13.6/54.42,CL_b))*Tpow)*CL_k_Boltzmann*T); */
  return,CL_B_gm*(exp(-1.0*(CL_aHII*(13.6/54.42)^CL_b)*Tpow)*CL_k_Boltzmann*T);
END

;/*-----------------------------------------------------------------
; *     Line Cooling
; *-----------------------------------------------------------------*/
;/*      CEN (1992, Ap.J.Suppl 78,341) ADVOCATES MULTIPLYING EACH OF 
; *      THESE RATES BY Cen_correctn - HE CLAIMS THIS GIVES THE RIGHT
; *      HIGH T LIMIT FOR PROCESSES INVOLVING A FREE EL INTERACTING 
; *      WITH AN ORBITAL ELECTRON ?? */

;FUNCTION clCoolLineHI,T
;CL_B_gm = (6.022e23*(938.7830/931.494))
;CL_aHI = -1275.29
;CL_bHI = 374.535
;CL_cHI = -42.8662
;CL_dHI = 2.19186
;CL_fHI = -0.0422634
;CL_bHI_old = 1.18348e05
;CL_MAX_NEG_EXP_ARG = -500.0

;arg = -CL_bHI_old/T             ;
;lnT = alog(T)                    ; 
;if (arg lt CL_MAX_NEG_EXP_ARG) return,0 ;
;ln_cool = CL_aHI +  CL_bHI*lnT +  CL_cHI * lnT^2.0 +  CL_dHI * lnT^3.0 +  CL_fHI*lnT^4.0 
;return,CL_B_gm*exp(ln_cool)     ;
;END

FUNCTION clCoolLineHI,T
CL_B_gm = (6.022e23*(938.7830/931.494))
CL_aHI = 7.5e-19
CL_bHI = 1.18348e05
Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5)) ;
CL_MAX_NEG_EXP_ARG = -500.0

T_inv=1.0/T                     ;
arg = -1.0*CL_bHI*T_inv             ;
if (arg lt CL_MAX_NEG_EXP_ARG) THEN return,0 ELSE return,CL_B_gm*CL_aHI*exp( arg )*Cen_correctn ;
END

;FUNCTION clCoolLineHeI,T
;CL_B_gm = (6.022e23*(938.7830/931.494))
;CL_aHeI = -3133.67
;CL_bHeI = 980.557
;CL_cHeI = -117.344
;CL_dHeI = 6.26898
;CL_fHeI = -0.125976
;CL_bHeI_old = 1.3179e04
;CL_MAX_NEG_EXP_ARG = -500.0

;arg = -CL_bHeI_old/T            ;
;lnT = alog(T)                    ; 
;if (arg lt CL_MAX_NEG_EXP_ARG) return,0 ;
;ln_cool = CL_aHeI +  CL_bHeI*lnT +  CL_cHeI * lnT^2.0 +  CL_dHeI * lnT^3.0 +  CL_fHeI*lnT^4.0 ; 
;return,CL_B_gm*exp(ln_cool)     ;
;END

FUNCTION clCoolLineHeI,T
CL_B_gm = (6.022e23*(938.7830/931.494))
CL_aHeI = 9.10e-27
CL_bHeI = 1.3179e04
CL_p_HeI = 0.1687
CL_MAX_NEG_EXP_ARG = -500.0

Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));
T_inv=1.0/T                     ;
arg = -1.0*CL_bHeI*T_inv            ;
if (arg lt CL_MAX_NEG_EXP_ARG) THEN return,0 ELSE return,CL_B_gm*CL_aHeI*exp(-1.0*CL_bHeI*T_inv)*T_inv^CL_p_HeI*Cen_correctn ;
END

;#ifdef CLOUDY 
;#define CL_aHeII  -4002.63
;#define CL_bHeII  1197.95 
;#define CL_cHeII  -136.703 
;#define CL_dHeII  6.96923 
;#define CL_fHeII  -0.133835
;#define CL_bHeII_old  4.73638e05
;FUNCTION clCoolLineHeII,T
;  double lnT, ln_cool, arg;
;  arg = -CL_bHeII_old/T;
;  lnT = log(T); 
;  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
;  ln_cool = CL_aHeII +  CL_bHeII*lnT +  CL_cHeII * pow(lnT, 2.0) +  CL_dHeII * pow(lnT, 3.0) +  CL_fHeII*pow(lnT, 4.0); 
;  return CL_B_gm*exp(ln_cool);
;END    

FUNCTION clCoolLineHeII,T
CL_B_gm = (6.022e23*(938.7830/931.494))
CL_aHeII =  5.54e-17
CL_bHeII =  4.73638e05
CL_p_HeII = 0.397
CL_MAX_NEG_EXP_ARG = -500.0

Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

T_inv=1.0/T                     ;
arg = -1.0*CL_bHeII*T_inv           ;
if (arg lt CL_MAX_NEG_EXP_ARG) THEN return,0 ELSE return,CL_B_gm*CL_aHeII*exp(-1.0*CL_bHeII*T_inv)*T_inv^CL_p_HeII*Cen_correctn ;
END

;//Cooling from rot-vib transitions or H2 CC -- Martin, Schwarz & Mandy, 1996
FUNCTION clCoolLineH2,T,YH
if (T gt 45000 OR YH lt 1e-3) THEN return,0 ;

alpha1 =-1.058656e3
alpha2 = 6.412920e2
alpha3 =-1.330667e2
alpha4 = 9.285717
alpha5 = 9.160106e1
alpha6 = 2.680075e3
alpha7 = 9.500433e3
alpha8 =-1.253746e3
alpha9 = 7.792906e2
alpha10=-1.628687e2
alpha11= 1.145088e1
alpha12= 1.057438e2
alpha13= 2.889762e3
alpha14= 1.382858e4
alpha15=-1.175497e2
alpha16= 7.886144e1
alpha17=-1.777082e1
alpha18= 1.338843
alpha19= 3.406424e3

logt = alog10(T)                   ;
logyh = alog10(YH)                 ;
  
a = -1.0*logyh + alpha1 + alpha2*logt + alpha3*logt*logt + alpha4*logt*logt*logt + alpha5*alog10(1.0 + alpha6/T) + alpha7/T ;
b = alpha8 + alpha9*logt + alpha10*logt*logt + alpha11*logt*logt*logt + alpha12*alog10(1.0 + alpha13/T) + alpha14/T ;
f = logyh - b                   ; 
w = alpha15 + alpha16*logt + alpha17*logt*logt + alpha18*logt*logt*logt + alpha19/T ;
lograd = a + 0.5*(f - sqrt(f*f+2.*w*w)) ;
return,10.0^lograd         ;
END
