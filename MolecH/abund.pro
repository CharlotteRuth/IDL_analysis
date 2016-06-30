pro abund,dt = dt, last=last, start = start

IF (NOT KEYWORD_SET(last)) THEN last = 3000
IF (NOT KEYWORD_SET(dt)) THEN dt = 20
IF (NOT KEYWORD_SET(start))THEN start = 1;10/dt + 1

base = "/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_H2_UV_2grad_soft/MW_disk" ;abund,dt = 1,start = 1, last = 1
;base = "/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/MW_disk"  ;abund,dt = 1,start = 10, last = 10
;base = "/astro/net/scratch1/christensen/MolecH/12M/CoolingCloud/CoolingCloud"
;base = "/astro/net/scratch1/christensen/MolecH/12M/H2_cooling/12M"

!Y.STYLE = 1
!X.STYLE = 1

cm_per_kpc = 3.0857d21
gm_per_msol = 1.989d33
amu_per_gm = 6.022142d23
molec_weight = (0.76*1 + 0.24*4.0)

msol_per_sysmass = 1.36e17
kpc_per_syslength = 1e5

lengthunit =  kpc_per_syslength*cm_per_kpc;system length unit in cm (=1kpc)
massunit   = msol_per_sysmass*gm_per_msol;system mass unit in gm
dens_convert =  msol_per_sysmass * gm_per_msol*amu_per_gm/molec_weight/kpc_per_syslength^3/ cm_per_kpc^3 ;massunit*amu_per_gm/molec_weight/lengthunit^3

FOR i = start, last/dt DO BEGIN 
    if (i*dt lt 10) THEN step = '0000'+STRTRIM(i*dt,2) ELSE BEGIN
        if (i*dt lt 100) THEN step = '000'+STRTRIM(i*dt,2) ELSE step = '00'+STRTRIM(dt*i,2)
    ENDELSE
    readcol,base+"."+step+".HI",HI_prime,/silent
    readcol,base+"."+step+".H2",H2_prime,/silent
    rtipsy,base+"."+step,h,g,d,s
;    h.ngas=1000
    HI = HI_prime[1:h.ngas]
    H2 = H2_prime[1:h.ngas]
    t = g.tempg
    rho = g.dens*dens_convert
    zmetal = g.zmetal

    ind = where(H2 ne 0)
;    if ind[0] ne -1 then begin
    HI = HI[ind]
    H2 = H2[ind]         
 ;   ENDIF
    ind = SORT(H2)
    IF((where(H2 lt 1e-7))[0] ne -1) THEN  BEGIN
        HI_noH2 = HI[where(H2 lt 1e-7)] 
        H2_noH2 = H2[where(H2 lt 1e-7)]
        t_noH2 = t[where(H2 lt 1e-7)]
        rho_noH2 = rho[where(H2 lt 1e-7)]
        zmetal_noH2 = zmetal[where(H2 lt 1e-7)]
        ind_noH2 = SORT(HI_noH2)
    ENDIF
    IF((where(H2 ge 1e-7))[0] ne -1) THEN BEGIN
        HI_yesH2 = HI[where(H2 ge 1e-7)] 
        H2_yesH2 = H2[where(H2 ge 1e-7)]
        t_yesH2 = t[where(H2 ge 1e-7)]
        rho_yesH2 = rho[where(H2 ge 1e-7)]
        zmetal_yesH2 = zmetal[where(H2 ge 1e-7)]        
        ind_yesH2 = SORT(H2_yesH2)
        IF((where(H2 lt 1e-7))[0] ne -1) THEN  BEGIN
            y1 = [HI_noH2[ind_noH2],HI_yesH2[ind_yesH2]]
            y2 = 2.0*[H2_noH2[ind_noH2],H2_yesH2[ind_yesH2]]
            t = [t_noH2[ind_noH2],t_yesH2[ind_yesH2]]
            rho = [rho_noH2[ind_noH2],rho_yesH2[ind_yesH2]]
            zmetal = [zmetal_noH2[ind_noH2],zmetal_yesH2[ind_yesH2]]
        ENDIF ELSE BEGIN
            y1 = HI_yesH2[ind_yesH2]
            y2 = 2.0*H2_yesH2[ind_yesH2]
            t =  t_yesH2[ind_yesH2]
            rho = rho_yesH2[ind_yesH2]  
            zmetal = zmetal_yesH2[ind_yesH2]       
        ENDELSE
    ENDIF ELSE BEGIN
        y1 = [HI_noH2[ind_noH2]]
        y2 = 2.0*[H2_noH2[ind_noH2]]
        t = [t_noH2[ind_noH2]]
        rho = [rho_noH2[ind_noH2]]
        zmetal = [zmetal_noH2[ind_noH2]]
    ENDELSE

    y3 = y1 + y2
    x = FINDGEN(N_ELEMENTS(y1))
;    window,0
;    plot,x,y1,/ylog,yrange=[1e-7,1],xrange=[0,N_ELEMENTS(y1)]
;    oplot,x,y2,linestyle = 2
;    print,i*dt,MINMAX(y2)
;    stop

    Y_HI = HI;dblarr(N_ELEMENTS(t))
    Y_H2 = H2;Y_HI
    print,'Computed Range'
    print,'HI: ',MINMAX(HI)
    print,'H2: ',MINMAX(H2)
    all = fltarr(2,N_ELEMENTS(t))
    for i=0LL,LONG64(N_ELEMENTS(t)) - 1 do  BEGIN
        all[*,i] = LTE_abund(T[i],rho[i],zmetal[i],Y_HI[i],Y_H2[i])
    ENDFOR
    Y_HI = all[0,*]
    Y_H2 = all[1,*]    
    print,''
    print,'LTE Range'
    print,'HI: ',MINMAX(Y_HI)
    print,'H2: ',MINMAX(Y_H2)
    print,'HI: ',aLOG10(MINMAX(Y_HI))
    print,'H2: ',ALOG10(MINMAX(Y_H2))
    Y_HI = HI_prime
    Y_H2 = H2_prime
    Y_HI[1:h.ngas] = all[0,*]
    Y_H2[1:h.ngas] = all[1,*]
    openw,1,base+"."+step+".HI_LTE"
    openw,2,base+"."+step+".H2_LTE"
    for i=0LL,LONG64(N_ELEMENTS(Y_HI)) - 1 do  BEGIN
        printf,1,Y_HI[i]
        printf,2,Y_H2[i]
    ENDFOR
    close,1
    close,2
    stop
;    wait,0.5
endfor
stop
end


FUNCTION LTE_abund,T,rho,zmetal,Y_HI,Y_H2
;print,'LTE Calc: ',T,rho,Y_HI,Y_H2
CL_B_gm  = 6.022e23*(938.7830/931.494)
en_B = CL_B_gm*rho     
;print,rho
;stop                               ;
EPS = 1e-5
;ZMetal = 0.01
mu_metal = 17.6003
;T = 1e4;

if (ZMetal le 0.1) then  yHe = (0.236 + 2.1*ZMetal)/4.0 else yHe = (-0.446*(ZMetal - 0.1)/0.9 + 0.446)/4.0;

yH = 1.0 - yHe*4.0 - ZMetal; 

;yH = 0.764                      ; //Total H atoms
yH2 = Y_H2                      ;0; //H2 molecules
yHI = Y_HI                      ;0.765; //HI atoms
yHII = 0                        ; //HII atoms

;yHe = 0.235                     ; //Total He atoms
yHeI = yHe                        ;.235;
yHeII = 0                   ;
yHeIII = 0                      ;

R_Radr_HII = clRateRadrHII(T) 
R_Phot_HI = 1e-100
R_Coll_HI = clRateCollHI(T)
R_TotrForm_H2 = clRateTotrFormH2(ZMetal)
R_Coll_e_H2 = clRateColl_e_H2(T)
R_Coll_H_H2 = clRateColl_H_H2(T)
R_Coll_H2_H2 = clRateColl_H2_H2(T)
R_Phot_H = 1e-100
R_Phot_HeI = 1e-100
R_Phot_H2 = 1e-100
R_Totr_HeII =clRateRadrHeII(T) + clRateDielHeII(T)
R_Phot_HeII = 1e-100
R_Radr_HeIII =clRateRadrHeII(T)
R_Coll_HeI = clRateCollHeI(T) 
R_Coll_HeII = clRateCollHeII(T) 

rcirrHeI  = (R_Coll_HeI)/(R_Totr_HeII) ;
rcirrHeII = (R_Coll_HeII)/(R_Radr_HeIII) ;
rpirrHeI  = (R_Phot_HeI)/(R_Totr_HeII * en_B) ;
rpirrHeII = (R_Phot_HeII)/(R_Radr_HeIII * en_B) ;

;print,'Rates: ',R_Radr_HII,R_Phot_HI,R_Coll_HI,R_TotrForm_H2,R_Coll_e_H2,R_Coll_H_H2,R_Coll_H2_H2,R_Phot_H2;
;if (R_Coll_e_H2 lt 1e-40) then stop
for  i=0, 20 do begin
    yHI_old   = yHI             ;
    yHeII_old = yHeII           ;
    yH2_old   = yH2             ;
   
    ye = (1.235-(yHI + 2 * yHeI + yHeII)) ; //Free electrons
    if (ye le 0) then begin
        y_e = 0                 ;
        y_HI = yH - yH2*2       ; /*Now H2, so what do I add CC*/
        Y_HII = 0               ;
        Y_HeI = yHe             ;
        Y_HeII = 0              ;
        Y_HeIII = 0             ;
        Y_Total = yH + yHe + ZMetal/17.6003 - yH2 ; /*Total particles, subtract off yH2 to avoid double counting H atoms CC*/
        break
    endif
    rye = 1/ye                  ;
    
    fHII = (R_Radr_HII*en_B*ye)/(R_Phot_HI + R_Coll_HI*en_B*ye) ;       //rcirrH2 + rpirrH2 * rye;
    fHI = 2.0*(R_TotrForm_H2*en_B*yHI_old)/(R_Coll_e_H2*en_B*ye + R_Coll_H_H2*en_B*yHI_old +  R_Coll_H2_H2*en_B*yH2_old +  R_Phot_H2)
    
    rfH  =  1 / ( 1 + fHII * (1 + fHI) ) ; 
    yHII =  yH * rfH            ;
    yHI  =  yH * rfH * fHII     ;
    yH2  = (yH - yHI - yHII)/2.0 ;
;    print,"Itir: ",yH,rfH,fHII,yHI,yH2
    
    fHeI  = rcirrHeI + rpirrHeI * rye ;/* HeI->HeII/HeII->HeI */
    fHeII = rcirrHeII+rpirrHeII * rye ; /* HeII->HeIII/HeIII->HeII */ 
    rfHe  = 1 / ( 1 + fHeI * (1 + fHeII) ) ;
    yHeI  = yHe * rfHe          ;
    yHeII = yHe * fHeI * rfHe   ;
    if ( ABS(yHeII_old-yHeII) gt EPS * yHeII && abs(yHI_old-yHI) lt EPS * yHI ) then  break
endfor
Y_e = ye                        ;
rye = 1/ye
fHeI  = rcirrHeI + rpirrHeI * rye

Y_HI = yHI                      ;
Y_HII = yHII                    ;
Y_H2 = yH2                      ;
if (Y_H2 lt 1e-20) then Y_H2 = 0 ;
if (Y_HI lt 1e-20) then Y_HI = 0 ;
if (Y_HII lt 1e-20) then Y_HII = 0 ;
Y_HeI = yHeI                    ;
Y_HeII = yHeII                  ;
fHeII = rcirrHeII + rpirrHeII*rye ;
Y_HeIII = yHe / ((1.0/fHeI+1.0)/fHeII+1.0) ;
if (Y_HeI lt 1e-20) then Y_HeI = 0 ;
if (Y_HeII lt 1e-20) then Y_HeII = 0 ;
if (Y_HeIII lt 1e-20) then Y_HeIII = 0 ;
;print,'LTE Calc: ',T,rho,Y_HI,Y_H2
;print,'InitAbund',T,rho,Y_HI,Y_H2,en_B
;print,'Fractions',fHI,fHII,rfH
RETURN,[Y_HI,Y_H2]
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

  TL = alog(T*CL_eV_per_K);
  arg = -32.713967867 + TL*(13.536556     + TL*(-5.73932875 + $
      TL*(1.56315498 + $
      TL*(-0.2877056     + TL*(3.48255977e-2 + TL*(-2.63197617e-3 + $
      TL*(1.11954395e-4  + TL*(-2.03914985e-6))))))));
  if (arg lt CL_MAX_NEG_EXP_ARG) then return, 0;
  return, DOUBLE(exp ( arg ));
end
  
;/*     He + e- -> He+ + 2e-  Janev et al. 1987 (Abel 1996) */
function clRateCollHeI,T 
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066d-16
CL_eV_erg       =   1.60219d-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

  TL = alog(T*CL_eV_per_K);
  arg = -44.09864886  + TL*(23.91596563   + TL*(-10.7532302 +  $
      TL*(3.05803875 +  $
      TL*(-0.56851189    + TL*(6.79539123e-2 + TL*(-5.00905610e-3 + $
      TL*(2.06723616e-4  + TL*(-3.64916141e-6))))))));
  if (arg lt CL_MAX_NEG_EXP_ARG) then return, 0;
  return, DOUBLE(exp( arg ));
end

;/*     He+ + e- -> He++ + 2e- Aladdin Database 1989 (Abel 1996) */
function clRateCollHeII,T  
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

  TL = alog(T*CL_eV_per_K);
  arg = -68.71040990  + TL*(43.93347633   + TL*(-18.4806699 + $
      TL*(4.70162649 + $
      TL*(-0.76924663    + TL*(8.113042e-2   + TL*(-5.32402063e-3 +  $
      TL*(1.97570531e-4  + TL*(-3.16558106e-6))))))));
  if (arg lt CL_MAX_NEG_EXP_ARG) then return, 0;
  return, DOUBLE(exp( arg ));
end

;//Proton Collision
;//Lepp & Shull, 1983
function clRateColl_H2_H2,T
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

  if (T lt 7291) then ratecollH2 = 5.22e-14*exp(-3.22e4/T) else ratecollH2 = 3.17e-15*exp(-4060./T - (7500./T)*(7500./T));
  return, DOUBLE(ratecollH2);
end

;//Electron Collision
;//Donahue & Shull, Abel 1997
function clRateColl_e_H2,T
CL_MAX_NEG_EXP_ARG = -500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

  ratecollH2 = 5.6e-11*SQRT(CL_eV_per_K*T)*exp(-8.8/(CL_eV_per_K*T)); //Is 8.8eV the right number?  Abel has 102124
  return, DOUBLE(ratecollH2);
end

;//Neutral H
;//Dove and Mandy 1986, Abel 1997
function clRateColl_H_H2,T
;  //Donahue & Shull  ratecollH2 = 6.11e-14*exp(-4.48*CL_eV_per_K/T);
CL_MAX_NEG_EXP_ARG  =-500.
CL_k_Boltzmann   =  1.38066e-16
CL_eV_erg       =   1.60219e-12
CL_eV_per_K     =   (CL_k_Boltzmann/CL_eV_erg)

  ratecollH2 = 1.067e-10*(CL_eV_per_K*T)^2.012*exp(-1.0*(4.463/T/CL_eV_per_K)*(1+0.2472*CL_eV_per_K*T)^3.512);
  return, DOUBLE(ratecollH2);
end

;'/*-----------------------------------------------------------------
; *     Radiative Recombination rates
; *-----------------------------------------------------------------*/
;/*     H+ + e- -> H + gam  Verner & Ferland 1996 */
function clRateRadrHII,T 
   Tsq = sqrt(T)                   ; 
   return, DOUBLE(7.982e-11/( Tsq*0.563615 * (1+Tsq*0.563615)^0.252 * (1+Tsq*1.192167e-3)^1.748));  
end

;/*     He+ + e- -> He + gam  radiative  Verner & Ferland 1996 */
function clRateRadrHeII,T 
  lgT = alog10(T)                  ; 
  lg_rec = -10.3217 + 0.749031* lgT - 1.66439*lgT^2 + 1.13114 * lgT^3 -0.459420 * lgT^4.0 + 0.114337 * lgT^5.0 -0.0174612 * lgT^6.0 + 0.00158680 * lgT^7.0 -7.86098e-05 * lgT^8 +  1.63398e-06 * lgT^9.0; 
  return, DOUBLE(10.0^lg_rec); 
end

function clRateTotrFormH2,z
  clump = 50.0; /*Ranges from 2-10 to 30-100, Gendin et al 2008 CC*/ 
  Rate_dust = 3.5e-17*z*clump; /*Formation rate coefficient of molecular hydrogen on dust, Gnedin et al 2008, Wolfire 2008, unit of cc per s CC*/  
  return, DOUBLE(Rate_dust);
end

function clRateGasFormH2,T 
  return, 0;
end

function clRateDielHeII,T
CL_MAX_NEG_EXP_ARG  =   -500.
  T_inv = 1.0/T

  arg = -4.7e5*T_inv;
  if (arg lt CL_MAX_NEG_EXP_ARG)THEN return, 0;
  return, DOUBLE(1.9e-3*T^(-1.5)*exp(arg)*(1+0.3*exp(-9.4e4*T_inv)));
end

function clRateChtrHeII,T
  T_4 = T/1e4; 
  return, 0. ; 
end

;/*     He++ + e- -> He+ + gam  Verner & Ferland 1996 */
function clRateRadrHeIII,T
  Tsq = sqrt(T);

  return, DOUBLE(1.891e-10/( Tsq*0.326686 * (1+Tsq*0.326686)^0.2476 * (1+Tsq*6.004084e-4)^1.7524));
end
