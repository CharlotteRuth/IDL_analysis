FUNCTION LTE_abund_extended,T,rho,zmetal,h
;LTE_abund,temperature,density,metallicity,length radiation travels over in cm ( = c/gas_shear)
;print,'LTE Calc: ',T,rho,Y_HI,Y_H2
CL_B_gm  = 6.022e23*(938.7830/931.494)
en_B = CL_B_gm*rho     
EPS = 1e-5
mu_metal = 17.600299999999997

;print,T,rho

if (ZMetal le 0.1) then  yHe = (0.236 + 2.1*ZMetal)/4.0 else yHe = (-0.446*(ZMetal - 0.1)/0.9 + 0.446)/4.0
yH = 1.0 - yHe*4.0 - ZMetal;
yH2 = 0                      
yHI = yH                      
yHII = 0                        
yHeI = 0                        
yHeII = 0                   
yHeIII = 0       
YeMax = yH + yHe*2               

R_Radr_HII = clRateRadrHII(T)  ;radiative dissociation
R_Phot_HI = -9.300920e-13; photon dissociation
R_Coll_HI = clRateCollHI(T) ; collisional 
R_TotrForm_H2 = clRateTotrFormH2(ZMetal)
R_Coll_e_H2 = clRateColl_e_H2(T)
R_Coll_H_H2 = clRateColl_H_H2(T)
R_Coll_H2_H2 = clRateColl_H2_H2(T)
R_Phot_H2 = -1.591700e-11;1d-100;1.5917e-11

R_Phot_HeI = 2.850000e-14 ;1d-100
R_Totr_HeII =clRateRadrHeII(T) + clRateDielHeII(T)
R_Phot_HeII = 2.940000e-16 ;1d-100
R_Radr_HeIII =clRateRadrHeII(T)
R_Coll_HeI = clRateCollHeI(T) 
R_Coll_HeII = clRateCollHeII(T) 

rcirrHeI  = (R_Coll_HeI)/(R_Totr_HeII) ;
rcirrHeII = (R_Coll_HeII)/(R_Radr_HeIII) ;
rpirrHeI  = (R_Phot_HeI)/(R_Totr_HeII * en_B) ;
rpirrHeII = (R_Phot_HeII)/(R_Radr_HeIII * en_B) ;

S_H2 = S_H2(yH2*rho*5.9753790e+23/yH,h)
S_d  = S_d(yHI*rho*5.9753790e+23,yH2*rho*5.9753790e+23,h)

;IF( t lt 6e3 AND t gt 5.95e3) then print,'Rates: ',R_Radr_HII,R_Coll_HI,R_Coll_e_H2,R_Coll_H_H2,R_Coll_H2_H2
;if (R_Coll_e_H2 lt 1e-40) then stop
for  i=0, 20 do begin
    yHI_old   = yHI             
    yHeII_old = yHeII           
    yH2_old   = yH2             
   
    ye = (yeMax-(yHI + 2 * yH2 + 2 * yHeI + yHeII)) ; //Free electrons
    if (ye le 0) then begin
        y_e = 0                 ;
        yHII = 0    
        yHeI = yHe             
        yHeII = 0              
        yHeIII = 0       

        fHI = 2.0*(R_TotrForm_H2*en_B*yHI_old)/
                 ((R_Coll_H_H2*en_B*yHI_old  + 
                   R_Coll_H2_H2*en_B*yH2_old +  
                   R_Phot_H2)*S_d*S_H2)
        yHI = yH/(1 + fHI)
        yH2 = (yH - yHI)/2.0       ; /*Now H2, so what do I add CC*/
    endif else begin
        rye = 1/ye              
    
        IF (s_d le 1d-150) THEN BEGIN
            yHI = 0
            yH2 = yH/2.0
            yHII = 0
        ENDIF ELSE BEGIN
            fHII = (R_Radr_HII*en_B*ye)/(R_Phot_HI*S_d + R_Coll_HI*en_B*ye) ;       //rcirrH2 + rpirrH2 * rye;
            fHI = 2.0*(R_TotrForm_H2*en_B*yHI_old)/
                     ((R_Coll_e_H2*en_B*ye       + 
                       R_Coll_H_H2*en_B*yHI_old  +  
                       R_Coll_H2_H2*en_B*yH2_old +  
                       R_Phot_H2)*S_d*S_H2)
            rfH  =  1 / ( 1 + fHII * (1 + fHI) )  
            yHII =  yH * rfH        
            yHI  =  yH * rfH * fHII 
            yH2  = (yH - yHI - yHII)/2.0 
        ENDELSE

        fHeI  = rcirrHeI + rpirrHeI * rye ;/* HeI->HeII/HeII->HeI */
        fHeII = rcirrHeII+rpirrHeII * rye ; /* HeII->HeIII/HeIII->HeII */ 
        rfHe  = 1 / ( 1 + fHeI * (1 + fHeII) ) ;
        yHeI  = yHe * rfHe      ;
        yHeII = yHe * fHeI * rfHe ;
        yHeIII = yHe / ((1.0/fHeI+1.0)/fHeII+1.0) ;
    endelse
    IF(~FINITE(YHI) OR ~FINITE(YH2)) THEN BEGIN
        print,rho*5.9753790e+23,T
        print,yHI, yH2
        print,yHI_old,yH2_old
        print,s_d,s_H2
        yHI = yHI_old
        yH2 = yH2_old
    ENDIF
    if ( ABS(yHeII_old-yHeII) lt EPS * yHeII AND abs(yHI_old-yHI) lt EPS * yHI ) then  break
;    S_H2 = S_H2(yH2,h)
;    S_d  = S_d(yHI,yH2,h)
;    S_H2 = S_H2(yH2*rho*5.9753790e+23/yH/2.0,h)
;    S_d  = S_d(rho*5.9753790e+23*(yH - YHII),0,h)
    S_H2 = S_H2(yH2*rho*5.9753790e+23/yH,h)
    S_d  = S_d(yHI*rho*5.9753790e+23,yH2*rho*5.9753790e+23,h)
    if (YH2 lt 0) then begin
        YH2 = 0                 ;
        YHII = yH - yHI
    ENDIF
    if (YHI lt 0) then yHI = 0;
    if (YHII lt 0) then begin
        YHII = 0 
        YH2 = (YH - yHI)/2;
    ENDIF
    if (YH2 gt yH/2.0) then YH2 = yH/2.0 ;
    if (YHI gt yH) then YHI = yH ;
    if (YHII gt yH) then YHII = yH ;
;    if (rho*5.9753790e+23 gt 10 AND rho*5.9753790e+23 lt 20) then begin
;        print,s_H2,s_d
;        print,YHI,YH2,YHII
;       stop
;    endif
endfor
Y_e = ye                        ;
;rye = 1/ye
;fHeI  = rcirrHeI + rpirrHeI * rye

Y_HI = yHI                      ;
Y_HII = yHII                    ;
Y_H2 = yH2                      ;
Y_HeI = yHeI                    ;
Y_HeII = yHeII                  ;
Y_HeIII = yHeIII
;fHeII = rcirrHeII + rpirrHeII*rye ;
;Y_HeIII = yHe / ((1.0/fHeI+1.0)/fHeII+1.0) ;
;if (Y_HeI lt 1e-20) then Y_HeI = 0 ;
;if (Y_HeII lt 1e-20) then Y_HeII = 0 ;
;if (Y_HeIII lt 1e-20) then Y_HeIII = 0 ;
;if (Y_HeI gt yHe) then Y_HeI = yHe ;
;if (Y_HeII gt yHe) then Y_HeII = yHe ;
;if (Y_HeIII gt yHe) then Y_HeIII = yHe ;
;IF(  t lt 6e3 AND t gt 5.95e3) then begin
;print,'LTE Calc: T: ',T,',rho: ',rho,',Zmetal: ',Zmetal,',Y_HI: ',Y_HI,', Y_H2: ',Y_H2
;print,'InitAbund',T,',',rho,',',Y_HI,',',Y_H2,',',en_B
;print,'en_B: ',en_B,R_Radr_HeIII
;print,'FinalAbund: ',yH,Y_HII,Y_HI,Y_H2,Y_HeI,Y_HeII,Y_HeIII
;print,'Rates: ',rcirrHeI,rpirrHeI,rcirrHeII,rpirrHeII
;print,'Fractions',fHI,fHII,rfH
;stop
;endif
IF(~FINITE(Y_HI) OR ~FINITE(Y_H2)) THEN stop
print,"HI: ",Y_HI,"; H2: ",Y_H2,"; HII: ",Y_HII,";"
RETURN,[Y_HI,Y_H2,Y_HII,Y_HeI,Y_HeII,Y_HeIII,Y_e,s_d,s_H2]
;return,number fraction of HI, number fraction of H2, number fraction of HeI,  number fraction of HeII,  number fraction of HeIII,  number fraction of electrons, dust shielding parameter, self-shielding parameter
;By number fraction, I mean (the number of HI)/(total number of atoms)  
;This is the same format as the XXXXX.HI files outputted from gasoline are in

end


;CL_Rgascode      =   8.2494e7
;CL_Eerg_gm_degK    =    CL_Rgascode
;CL_ev_degK         =    1.0/1.1604e4
;CL_Eerg_gm_ev      =    CL_Eerg_gm_degK/CL_ev_degK
;CL_Eerg_gm_degK3_2   =  1.5*CL_Eerg_gm_degK
;CL_MAX_NEG_EXP_ARG  =   -500.
;CL_B_gm        =    (6.022e23*(938.7830/931.494))  ;/*Avegadro's Number * Mass_Hydrogen/Energy_AMU */
;CL_k_Boltzmann   =  1.38066e-16
;CL_eV_erg       =   1.60219e-12
;CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

; /*-----------------------------------------------------------------
; *     Collisional Ionization rates
; *-----------------------------------------------------------------*/
;/*     H + e- -> H+ + 2e-  Janev et al. 1987 (Abel 1996) */
function clRateCollHI,T
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066d-16
CL_eV_erg       =   1.60219d-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

TL = alog(T*CL_eV_per_K)        
arg = -32.713967867 + TL*(13.536556     + TL*(-5.73932875 + $
                                              TL*(1.56315498 + $
                                                  TL*(-0.2877056     + TL*(3.48255977e-2 + TL*(-2.63197617e-3 + $
                                                                                               TL*(1.11954395e-4  + TL*(-2.03914985e-6)))))))) 
if (arg lt CL_MAX_NEG_EXP_ARG) then return, 0 
return, DOUBLE(exp ( arg ))     
end
  
;/*     He + e- -> He+ + 2e-  Janev et al. 1987 (Abel 1996) */
function clRateCollHeI,T 
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066d-16
CL_eV_erg       =   1.60219d-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

TL = alog(T*CL_eV_per_K)        
arg = -44.09864886  + TL*(23.91596563   + TL*(-10.7532302 +  $
                      TL*(3.05803875 +  $
                      TL*(-0.56851189    + TL*(6.79539123e-2 + TL*(-5.00905610e-3 + $
                      TL*(2.06723616e-4  + TL*(-3.64916141e-6)))))))) 
if (arg lt CL_MAX_NEG_EXP_ARG) then return, 0 
return, DOUBLE(exp( arg ))
end

;/*     He+ + e- -> He++ + 2e- Aladdin Database 1989 (Abel 1996) */
function clRateCollHeII,T  
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

TL = alog(T*CL_eV_per_K)        
arg = -68.71040990  + TL*(43.93347633 + TL*(-18.4806699 + $
                      TL*(4.70162649  + TL*(-0.76924663 + $
                      TL*(8.113042e-2 + TL*(-5.32402063e-3 +  $
                      TL*(1.97570531e-4  + TL*(-3.16558106e-6))))))))
if (arg lt CL_MAX_NEG_EXP_ARG) then return, 0
return, DOUBLE(exp( arg ))
end

;//Proton Collision
;//Lepp & Shull, 1983
function clRateColl_H2_H2,T
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

if (T lt 7291) then ratecollH2 = 5.22e-14*exp(-3.22e4/T) else ratecollH2 = 3.17e-15*exp(-4060./T - (7500./T)*(7500./T))
return, DOUBLE(ratecollH2)
end

;//Electron Collision
;//Donahue & Shull, Abel 1997
function clRateColl_e_H2,T
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

ratecollH2 = 5.6e-11*SQRT(CL_eV_per_K*T)*exp(-8.8/(CL_eV_per_K*T)) ; //Is 8.8eV the right number?  Abel has 102124
return, DOUBLE(ratecollH2)
end

;//Neutral H
;//Dove and Mandy 1986, Abel 1997
function clRateColl_H_H2,T
;  //Donahue & Shull  ratecollH2 = 6.11e-14*exp(-4.48*CL_eV_per_K/T);
CL_MAX_NEG_EXP_ARG  =-500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

ratecollH2 = 1.067e-10*(CL_eV_per_K*T)^2.012*exp(-1.0*(4.463/T/CL_eV_per_K)*(1+0.2472*CL_eV_per_K*T)^3.512) 
return, DOUBLE(ratecollH2)      
end

;'/*-----------------------------------------------------------------
; *     Radiative Recombination rates
; *-----------------------------------------------------------------*/
;/*     H+ + e- -> H + gam  Verner & Ferland 1996 */
function clRateRadrHII,T 
Tsq = sqrt(T)                   
return, DOUBLE(7.982e-11/( Tsq*0.563615 * (1+Tsq*0.563615)^0.252 * (1+Tsq*1.192167e-3)^1.748)) 
end

;/*     He+ + e- -> He + gam  radiative  Verner & Ferland 1996 */
function clRateRadrHeII,T 
lgT = alog10(T)                  
lg_rec = -10.3217 + 0.749031* lgT - 1.66439*lgT^2 + 1.13114 * lgT^3 -0.459420 * lgT^4.0 + 0.114337 * lgT^5.0 -0.0174612 * lgT^6.0 + 0.00158680 * lgT^7.0 -7.86098e-05 * lgT^8 +  1.63398e-06 * lgT^9.0 
return, DOUBLE(10.0^lg_rec)      
end

function clRateTotrFormH2,z
clump = 10.0   ; /*Ranges from 2-10 to 30-100, Gendin et al 2008 CC*/ 
Rate_dust = 3.5e-17*z/0.0177*clump ; /*Formation rate coefficient of molecular hydrogen on dust, Gnedin et al 2008, Wolfire 2008, unit of cc per s CC*/  
return, DOUBLE(Rate_dust)       
end

function clRateGasFormH2,T 
  return, 0;
end

function clRateDielHeII,T
CL_MAX_NEG_EXP_ARG  =   -500.
T_inv = 1.0/T
arg = -4.7e5*T_inv              
if (arg lt CL_MAX_NEG_EXP_ARG)THEN return, 0 
return, DOUBLE(1.9e-3*T^(-1.5)*exp(arg)*(1+0.3*exp(-9.4e4*T_inv))) ;
end

function clRateChtrHeII,T
T_4 = T/1e4                      
return, 0.                       
end

;/*     He++ + e- -> He+ + gam  Verner & Ferland 1996 */
function clRateRadrHeIII,T
Tsq = sqrt(T)                   ;
return, DOUBLE(1.891e-10/( Tsq*0.326686 * (1+Tsq*0.326686)^0.2476 * (1+Tsq*6.004084e-4)^1.7524)) ;
end

function S_H2,yH2,h
omega_H2 = 0.2
x = yH2*h/5d14
return, (1 - omega_H2)/(1 + x)^2 + omega_H2/SQRT(1 + x)*exp(-0.00085*SQRT(1 + x))
end

function S_d,yHI,yH2,h
sigma_d = 2d-21
return, exp(-1.0*sigma_d*(yHI*h + 2.0*yH2*h))
end

;Length: 1.332104e+21, Rho HI: 1.067018e-02, Rho H2: 1.670782e-11
