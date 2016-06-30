 ;coolingCurve,filename ="/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/glass_highd/cooldebug.txt" 
;coolingCurve,filename = "/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/glass_highd_nometal/cooldebug.txt"
;coolingCurve,filename ="/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/glass_highd_z0/cooldebug.txt"
;coolingCurve,filename ="/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/glass_highd_noH2_metal/cooldebug.txt"

pro coolingCurveDetail, filename = filename,outplot = outplot
outfilebase = '~/Scratch2/MolecH/results/'
outfilename = '_z0'
outfilename = '_noH2_metal'
CL_B_gm = (6.022d23*(938.7830/931.494))

nel = 400.
min = 1.
max = 5.
tarray = 10.0^(findgen(nel)/nel*(max - min) + min)
CoolLineH2_H = dblarr(nel)
CoolLineH2_H2 = dblarr(nel)
CoolLineH2_He = dblarr(nel)
CoolLineH2_e = dblarr(nel)
CoolLineH2_HII = dblarr(nel)

FOR i = 0, nel - 1 DO BEGIN
    CoolLineH2_H[i] = clCoolLineH2_H(Tarray[i])
    CoolLineH2_H2[i] = clCoolLineH2_H2(Tarray[i])
    CoolLineH2_He[i] = clCoolLineH2_He(Tarray[i])
    CoolLineH2_e[i] = clCoolLineH2_e(Tarray[i])
    CoolLineH2_HII[i] = clCoolLineH2_HII(Tarray[i])
ENDFOR

readcol,filename,T,     en_B,  n_e, H2, HI,HII,HeI,HeII,HeIII,Sheild, Edot, Comp,Brems, Diel, Radr, Line,LineH2_H,LineH2_H2, LineH2_He,  LineH2_e, LineH2_HII, Coll,CollH2, LowT,MCool, Cool, Phot, PhotH2, MHeat
;                 %#2f, %#2f, %#2e, %g, %g, %g, %g,  %g,   %g,  %#2e,
;                 %#2e, %#2e, %#2e, %#2e, %#2e, %#2e,  %#2e, %#2e,
;                 %#2e, %#2e, %#2e, %#2e, %#2e,  %#2e, %#2e
nel = N_ELEMENTS(en_B)
x = indgen(nel)
rho = en_B
rho_gm_cc = en_B/CL_B_gm
ye = n_e/CL_B_gm
LowT = LowT[1:nel - 1]
MCool = MCool[1:nel - 1]
Cool = Cool[1:nel - 1]
Phot = Phot[1:nel - 1]
PhotH2 = PhotH2[1:nel - 1]
MHeat = MHeat[1:nel - 1]
T = T[1:nel - 1]
en_B = en_B[1:nel - 1]
n_e = n_e[1:nel - 1]
H2 = H2[1:nel - 1]
HI = HI[1:nel - 1]
HII = HII[1:nel - 1]
HeI = HeI[1:nel - 1]
HeII = HeII[1:nel - 1]
HeIII = HeIII[1:nel - 1]
Sheild = Sheild[1:nel - 1]
Edot = Edot[1:nel - 1]
Comp = Comp[1:nel - 1]
Brems = Brems[1:nel - 1]
Diel = Diel[1:nel - 1]
Radr = Radr[1:nel - 1]
Line = Line[1:nel - 1]
LineH2_H = LineH2_H[1:nel - 1] 
LineH2_H2 = LineH2_H2[1:nel - 1] 
LineH2_He = LineH2_He[1:nel - 1] 
LineH2_e = LineH2_e[1:nel - 1] 
LineH2_HII = LineH2_HII[1:nel - 1] 
Coll = Coll[1:nel - 1]
CollH2 = CollH2[1:nel - 1]
YH = (2.0*H2 + HI + HII)

CL_Rgascode        = 8.2494e7
CL_Eerg_gm_degK    = CL_Rgascode
CL_ev_degK         = 1.0/1.1604e4
CL_Eerg_gm_ev      = CL_Eerg_gm_degK/CL_ev_degK
CL_Eerg_gm_degK3_2 = 1.5*CL_Eerg_gm_degK
Y_Total =  n_e/rho + H2 + HI + HII + HeI + HeII + HEII 
ergs = Y_Total*CL_Eerg_gm_degK3_2*T
set_plot,'x'
loadct,39
window,0
plot,n_e/rho + 2.0*H2 + HI + HeI*2.0 + HeII

;stop
plot,T,-1.0/en_B/CL_B_gm*cool/yH/yH,ytitle = "Cooling [ergs s^-1 g^-1]",xtitle = "Output",/ylog,xrange=[10,1e8],/xlog,yrange=[1e-27,1e-21],ystyle = 1, xstyle = 1
oplot,T,-1.0/en_B/CL_B_gm*Brems/yH/yH,linestyle = 1,color = 50 ;Minescule
;oplot,x,-1.0*n_e * (clCoolBrem1(T) * (HII + HeII )+ clCoolBrem2(T) *HeIII),linestyle = 2,color = 20 
oplot,T,-1.0/en_B/CL_B_gm*Radr/yH/yH,linestyle = 3,color = 50
oplot,T,-1.0/en_B/CL_B_gm*Line/yH/yH,linestyle = 2,color = 100
;oplot,x,-1.0*n_e * (clCoolLineHI(T)*HI + clCoolLineHeI(T)*HeI + clCoolLineHeII(T)*HeII),linestyle = 2
oplot,T,-1.0/en_B/CL_B_gm*MCool/yH/yH,linestyle = 2,color = 180
oplot,T,-1.0/en_B/CL_B_gm*(LineH2_H + LineH2_H2 + LineH2_He + LineH2_e + LineH2_HII)/yH/yH,linestyle = 2, color = 120
oplot,T,-1.0/en_B/CL_B_gm*Coll/yH/yH,linestyle = 4,color = 240
oplot,T,-1.0/en_B/CL_B_gm*CollH2/yH/yH,linestyle = 4,color = 200
legend,['Brems','Radr','Line (H, He)','Line (Metal)','Line (H2)','Coll','Coll H2','Total'],linestyle = [1,3,2,2,2,4,4,0],color = [20,50,100,180,120,240,200,0],/top,/left

window,1
;plot,T,-1.0*(CollH2 + LineH2_H +  LineH2_H2 + LineH2_He + LineH2_e + LineH2_HII),xtitle = 'T [K]',ytitle = 'Cooling erg s^-1 g^-1',/xlog,/ylog,xrange = [1e1, 1e4]

;stop
;plot, T,-1.0/en_B/CL_B_gm*(LineH2_H/HI/H2 +  LineH2_H2/H2/H2 + LineH2_He/HeI/H2 + LineH2_e/ye/H2 + LineH2_HII/HII/H2),xtitle = 'T [K]',ytitle = 'Cooling erg cm^3 s^-1 ',/xlog,/ylog,yrange=[1e-28,1e-21],xrange=[100,1e4],ystyle = 1
;oplot,T,-1.0/en_B/CL_B_gm*(LineH2_H/HI/H2 +  LineH2_H2/H2/H2 + LineH2_He/HeI/H2 + LineH2_e/ye/H2 + LineH2_HII/HII/H2),linestyle = 0,color = 240
;oplot,Tarray,CoolLineH2_H + CoolLineH2_H2 + CoolLineH2_He + CoolLineH2_e + CoolLineH2_HII,linestyle = 0,color = 180

plot,T,-1.0/en_B/CL_B_gm*LineH2_H/H2/HI,xtitle = 'T [K]',ytitle = 'Cooling erg cm^3 s^-1 ',/xlog,/ylog,yrange=[1e-28,1e-21],xrange=[100,1e4],ystyle = 1,linestyle = 1
oplot,Tarray,CoolLineH2_H,linestyle = 1,color = 70

oplot,T,-1.0/en_B/CL_B_gm*LineH2_H2/H2/H2,linestyle = 2
oplot,Tarray,CoolLineH2_H2,linestyle = 2,color = 70

oplot,T,-1.0/en_B/CL_B_gm*LineH2_He/H2/HeI,linestyle = 3
oplot,Tarray,CoolLineH2_He,linestyle = 3,color = 70

oplot,T,-1.0/CL_B_gm*LineH2_e/H2/n_e,linestyle = 4
oplot,Tarray,CoolLineH2_e,linestyle = 4,color = 70

oplot,T,-1.0/en_B/CL_B_gm*LineH2_HII/H2/HII,linestyle = 5
oplot,Tarray,CoolLineH2_HII,linestyle = 5,color = 70
legend,['H','H2','He','e','p'],linestyle = [1,2,3,4,5],/top,/left

window,2
plot,t,hI,/xlog,/ylog,linestyle = 1,yrange = [1e-8,1e0],xrange = [10,1e8],xstyle = 1
oplot,t,hI,color = 80,linestyle = 1
oplot,t,HII,color = 120,linestyle = 1
oplot,t,h2*2.0,color = 50,linestyle = 1
oplot,t,HeI,color = 180,linestyle = 2
oplot,t,HeII,color = 200,linestyle = 2
oplot,t,HeIII,color = 240,linestyle = 2
oplot,t,n_e,linestyle = 0
legend,['H2','HI','HII','HeI','HeII','HeIII','N_e'],color = [30,80,120,180,200,240,0],linestyle = [1,1,1,2,2,2,0],/right,/top

window,3
xHe = 0.0825
xH2 = 0.001
xHII = 1e-4
xe = 1e-4
xHI = 1 - xH2*2 - xHII 
plot, T,-1.0/en_B/CL_B_gm*(LineH2_H*xHI/HI +  LineH2_H2*xH2/H2 + LineH2_He*xHe/HeI + LineH2_e*xe/n_e + LineH2_HII*xHII/HII)/H2,xtitle = 'T [K]',ytitle = 'Cooling erg cm^3 s^-1 ',/ylog,yrange=[1e-28,1e-23],xrange=[0,2000],ystyle = 1
;oplot,T,-1.0/en_B/CL_B_gm*(LineH2_H/HI/H2 +  LineH2_H2/H2/H2 + LineH2_He/HeI/H2 + LineH2_e/ye/H2 + LineH2_HII/HII/H2),linestyle = 0,color = 240
oplot,Tarray,CoolLineH2_H*xHI + CoolLineH2_H2*xH2 + CoolLineH2_He*xHe + CoolLineH2_e*xe + CoolLineH2_HII*xHII,color = 100,linestyle = 2

;plot,T,-1.0*en_B/CL_B_gm*Edot,xtitle = 'T [K]',ytitle = 'Cooling erg cm^3 s^-1 ',/xlog,/ylog
;stop

IF KEYWORD_SET(outpolot) THEN BEGIN
    set_plot,'ps'
    
    device,filename = outfilebase + 'coolingCurve'+outfilename+'.ps',/color,bits_per_pixel=8,/times
    plot,T,-1.0/en_B/CL_B_gm*cool/yH/yH,ytitle = "Cooling [ergs s^-1 g^-1]",xtitle = "Output",/ylog,xrange=[10,1e8],/xlog,yrange=[1e-27,1e-21],ystyle = 1, xstyle = 1
    oplot,T,-1.0/en_B/CL_B_gm*Brems/yH/yH,linestyle = 1,color = 20 ;Minescule
;oplot,x,-1.0*n_e * (clCoolBrem1(T) * (HII + HeII )+ clCoolBrem2(T) *HeIII),linestyle = 2,color = 20 
    oplot,T,-1.0/en_B/CL_B_gm*Radr/yH/yH,linestyle = 3,color = 50
    oplot,T,-1.0/en_B/CL_B_gm*Line/yH/yH,linestyle = 2,color = 100
;oplot,x,-1.0*n_e * (clCoolLineHI(T)*HI + clCoolLineHeI(T)*HeI + clCoolLineHeII(T)*HeII),linestyle = 2
    oplot,T,-1.0/en_B/CL_B_gm*MCool/yH/yH,linestyle = 2,color = 180
    oplot,T,-1.0/en_B/CL_B_gm*(LineH2_H + LineH2_H2 + LineH2_He + LineH2_e + LineH2_HII)/yH/yH,linestyle = 2, color = 120
    oplot,T,-1.0/en_B/CL_B_gm*Coll/yH/yH,linestyle = 4,color = 240
    oplot,T,-1.0/en_B/CL_B_gm*CollH2/yH/yH,linestyle = 4,color = 200
    legend,['Brems','Radr','Line (H, He)','Line (Metal)','Line (H2)','Coll','Coll H2','Total'],linestyle = [1,3,2,2,2,4,4,0],color = [20,50,100,180,120,240,200,0],/top,/left
    device,/close
    
    device,filename = outfilebase + 'coolingH2'+outfilename+'.ps',/color,bits_per_pixel=8,/times
    plot,T,-1.0/en_B/CL_B_gm*LineH2_H/H2/HI,xtitle = 'T [K]',ytitle = 'Cooling erg cm^3 s^-1 ',/xlog,/ylog,yrange=[1e-28,1e-21],xrange=[100,1e4],ystyle = 1,linestyle = 1
    oplot,Tarray,CoolLineH2_H,linestyle = 1,color = 70
    
    oplot,T,-1.0/en_B/CL_B_gm*LineH2_H2/H2/H2,linestyle = 2
    oplot,Tarray,CoolLineH2_H2,linestyle = 2,color = 70
    
    oplot,T,-1.0/en_B/CL_B_gm*LineH2_He/H2/HeI,linestyle = 3
    oplot,Tarray,CoolLineH2_He,linestyle = 3,color = 70
    
    oplot,T,-1.0/CL_B_gm*LineH2_e/H2/n_e,linestyle = 4
    oplot,Tarray,CoolLineH2_e,linestyle = 4,color = 70
    
    oplot,T,-1.0/en_B/CL_B_gm*LineH2_HII/H2/HII,linestyle = 5
    oplot,Tarray,CoolLineH2_HII,linestyle = 5,color = 70
    legend,['H','H2','He','e','p'],linestyle = [1,2,3,4,5],/top,/left
    device,/close

    device,filename = outfilebase + 'fracabund'+outfilename+'.ps',/color,bits_per_pixel=8,/times
    plot,t,hI,/xlog,/ylog,linestyle = 1,yrange = [1e-8,1e0],xrange = [10,1e8],xstyle = 1
    oplot,t,hI,color = 80,linestyle = 1
    oplot,t,HII,color = 120,linestyle = 1
    oplot,t,h2*2.0,color = 30,linestyle = 1
    oplot,t,HeI,color = 180,linestyle = 2
    oplot,t,HeII,color = 200,linestyle = 2
    oplot,t,HeIII,color = 240,linestyle = 2
    oplot,t,n_e,linestyle = 0
    legend,['H2','HI','HII','HeI','HeII','HeIII','N_e'],color = [30,80,120,180,200,240,0],linestyle = [1,1,1,2,2,2,0],/right,/top

    device,/close

    device,filename = outfilebase + 'coolingH2_total'+outfilename+'.ps',/color,bits_per_pixel=8,/times
    plot, T,-1.0/en_B/CL_B_gm*(LineH2_H*xHI/HI +  LineH2_H2*xH2/H2 + LineH2_He*xHe/HeI + LineH2_e*xe/n_e + LineH2_HII*xHII/HII)/H2,xtitle = 'T [K]',ytitle = 'Cooling erg cm^3 s^-1 ',/ylog,yrange=[1e-28,1e-23],xrange=[0,2000],ystyle = 1
;oplot,T,-1.0/en_B/CL_B_gm*(LineH2_H/HI/H2 +  LineH2_H2/H2/H2 + LineH2_He/HeI/H2 + LineH2_e/ye/H2 + LineH2_HII/HII/H2),linestyle = 0,color = 240
    oplot,Tarray,CoolLineH2_H*xHI + CoolLineH2_H2*xH2 + CoolLineH2_He*xHe + CoolLineH2_e*xe + CoolLineH2_HII*xHII,color = 100,linestyle = 2
    device,/close
ENDIF
stop
END

PRO coolingFits
nel = 400.
min = 1.
max = 5.
t = 10.0^(findgen(nel)/nel*(max - min) + min)
CoolLineH2_H = dblarr(nel)
CoolLineH2_H2 = dblarr(nel)
CoolLineH2_He = dblarr(nel)
CoolLineH2_e = dblarr(nel)
CoolLineH2_HII = dblarr(nel)

FOR i = 0, nel - 1 DO BEGIN
    CoolLineH2_H[i] = clCoolLineH2_H(T[i])
    CoolLineH2_H2[i] = clCoolLineH2_H2(T[i])
    CoolLineH2_He[i] = clCoolLineH2_He(T[i])
    CoolLineH2_e[i] = clCoolLineH2_e(T[i])
    CoolLineH2_HII[i] = clCoolLineH2_HII(T[i])
ENDFOR
plot,t,CoolLineH2_H,/xlog,/ylog
oplot,t,CoolLineH2_H2,linestyle = 1
oplot,t,CoolLineH2_He,linestyle = 2
oplot,t,CoolLineH2_e,linestyle = 3
oplot,t,CoolLineH2_HII,linestyle = 4
oplot,t,clCoolBrem1(t),linestyle = 5
oplot,t,clCoolBrem2(t),linestyle = 6
oplot,t,clCoolLineHI(t),linestyle = 1
oplot,t,clCoolLineHeI(t),linestyle = 2
legend,["H2-H","H2-H2","H2-He","H2-e","H2-HII"],linestyle =[0,1,2,3,4]
stop

END

FUNCTION clCoolLineH2_H,T     ;{ /* Cooling based on radiating out of a H2-H collisionally-induced excited state, Glover & Abel 08 */
    a00 = -16.818342
    a10 = 37.383713
    a20 = 58.145166
    a30 = 48.656103
    a40 = 20.159831
    a50 = 3.8479610;
  a01 = -24.311209
    a11 = 3.5692468
    a21 = -11.332860
    a31 = -27.850082
    a41 = -21.328264
    a51 = -4.2519023
  a02 = -24.311209
    a12 = 4.6450521
    a22 = -3.7209846
    a32 = 5.9369081
    a42 = -5.5108047
    a52 = 1.5538288

  xint = 6000.

 ; print,'Slope of H2-H Curve at T = 6000.0:       ',
;  slope = a12*1.0/1000.0 + $
;    a22*alog10(xint/1000.0)*2.0/1000.0  + $
;    a32*alog10(xint/1000.0)*alog10(xint/1000.0)*3.0/1000.0 + $ 
;    a42*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*4.0/1000.0 +  $ 
;    a52*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*5.0/1000.0 
  slope = 2.10095

  ;print,'Y-intercept of H2-H Curve at T = 6000.0: ', 
;  yint = 10.0^(a02 +  $ 
;               a12*alog10(xint/1000.0) +  $ 
;               a22*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;               a32*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;               a42*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;               a52*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0))
  yint = 1.86368e-22

  if (T le 100) THEN BEGIN return, 10.0^(a00 + $
                                 a10*alog10(T/1000.0) + $ 
                                 a20*alog10(T/1000.0)*alog10(T/1000.0) + $ 
                                 a30*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +  $ 
                                 a40*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +  $ 
                                 a50*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0));
  ENDIF ELSE BEGIN
      if (T le 1000) THEN BEGIN return, 10.0^(a01 +  $ 
                                 a11*alog10(T/1000.0) +  $ 
                                 a21*alog10(T/1000.0)*alog10(T/1000.0) +  $ 
                                 a31*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +  $ 
                                 a41*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +  $ 
                                 a51*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0));
      ENDIF else BEGIN  
          if (T le xint) THEN BEGIN return, 10.0^(a02 +  $ 
                        a12*alog10(T/1000.0) +  $ 
                        a22*alog10(T/1000.0)*alog10(T/1000.0) +  $ 
                        a32*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +  $ 
                        a42*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +  $ 
                        a52*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)) ;
          ENDIF ELSE return, 10.0^(slope*(alog10(T/1000.0) - alog10(xint/1000.0)) + alog10(yint))
      ENDELSE
  ENDELSE
END

;-------------------------------------------------------------------------------------------------------------------------------
FUNCTION clCoolLineH2_H2,T ;/* Cooling based on radiating out of a H2-H collisionally-induced excited state, Glover & Abel 08 */
  a0 = -23.962112
  a1 = 2.09433740
  a2 = -0.77151436
  a3 = 0.43693353
  a4 = -0.14913216
  a5 = -0.033638326             ;
  
  xint = 6000.0

;  print,'Slope of H2-H2 Curve at T = 6000.0:       ', a1*1.0 + $
;    a2*alog10(xint/1000.0)*2.0  + $
;    a3*alog10(xint/1000.0)*alog10(xint/1000.0)*3.0 + $ 
;    a4*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*4.0 +  $ 
;    a5*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*5.0 
  slope = 1.34460

;  print,'Y-intercept of H2-H2 Curve at T = xint.0: ', 10.0^(a0 +  $ 
;                                                            a1*alog10(xint/1000.0) +  $ 
;                                                            a2*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a3*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a4*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a5*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0))
  yint = 2.19802e-23
  
  IF (T le xint) THEN BEGIN
      return, 10.0^(a0 +   $
	     a1*alog10(T/1000.0) +   $
	     a2*alog10(T/1000.0)*alog10(T/1000.0) +   $
	     a3*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
	     a4*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
	     a5*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0));
  ENDIF ELSE return, 10.0^(slope*(alog10(T/1000.0) - alog10(xint/1000.0)) + alog10(yint))
END

;-----------------------------------------------------------------------------------------------------------------------------------
FUNCTION clCoolLineH2_He,T ;{ /* Cooling based on radiating out of a H2-H collisionally-induced excited state, Glover & Abel 08 */
  a0 = -23.689237
    a1 = 2.1892372
    a2 = -0.81520438
    a3 = 0.29036281
    a4 = -0.16596184
    a5 = 0.19191375
  
  xint = 6000.0

;  print,'Slope of H2-He Curve at T = 6000.0:       ', a1*1.0 + $
;    a2*alog10(xint/1000.0)*2.0  + $
;    a3*alog10(xint/1000.0)*alog10(xint/1000.0)*3.0 + $ 
;    a4*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*4.0 +  $ 
;    a5*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*5.0 
  slope = 1.48703

;  print,'Y-intercept of H2-He Curve at T = xint.0: ', 10.0^(a0 +  $ 
;                                                            a1*alog10(xint/1000.0) +  $ 
;                                                            a2*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a3*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a4*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a5*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0))
  yint = 4.48145e-23

  IF (T le xint) THEN BEGIN
      return, 10.0^(a0 +   $
                    a1*alog10(T/1000.0) +   $
                    a2*alog10(T/1000.0)*alog10(T/1000.0) +   $
                    a3*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
                    a4*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
                    a5*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)) ;
  ENDIF ELSE return, 10.0^(slope*(alog10(T/1000.0) - alog10(xint/1000.0)) + alog10(yint))
END

;-----------------------------------------------------------------------------------------------------------------------------------
FUNCTION clCoolLineH2_HII,T ;{ /* Cooling based on radiating out of a H2-H collisionally-induced excited state, Glover & Abel 08 */
  a0 = -21.716699
    a1 = 1.3865783
    a2 = -0.37915285
    a3 = 0.11453688
    a4 = -0.23214154
    a5 = 0.058538864

  xint = 10000.0

;  print,'Slope of H2-HII Curve at T = 10000.0:       ', a1*1.0 + $
;    a2*alog10(xint/1000.0)*2.0  + $
;    a3*alog10(xint/1000.0)*alog10(xint/1000.0)*3.0 + $ 
;    a4*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*4.0 +  $ 
;    a5*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*5.0 
  slope = 0.336011

;  print,'Y-intercept of H2-HII Curve at T = xint.0: ', 10.0^(a0 +  $ 
;                                                            a1*alog10(xint/1000.0) +  $ 
;                                                            a2*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a3*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a4*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a5*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0))
  yint = 1.70474e-21

  IF (T le xint) THEN BEGIN  
      return, 10.0^( a0 +   $
                     a1*alog10(T/1000.0) +   $
                     a2*alog10(T/1000.0)*alog10(T/1000.0) +   $
                     a3*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
                     a4*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
                     a5*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)) ;
  ENDIF ELSE return, 10.0^(slope*(alog10(T/1000.0) - alog10(xint/1000.0)) + alog10(yint))      
END

;-------------------------------------------------------------------------------------------------------------------------------
FUNCTION clCoolLineH2_e,T;{ /* Cooling based on radiating out of a H2-H collisionally-induced excited state, Glover & Abel 08 */
  a00 = -34.286155
    a10 = -48.537163
    a20 = -77.121176
    a30 = -51.352459
    a40 = -15.169160
    a50 = -0.98120322;
  a01 = -22.190316
    a11 = 1.5728955
    a21 = -0.21335100
    a31 = 0.96149759
    a41 = -0.91023495
    a51 = 0.13749749

  xint = 10000.0

;  print,'Slope of H2-e Curve at T = 10000.0:       ', a11*1.0 + $
;    a21*alog10(xint/1000.0)*2.0  + $
;    a31*alog10(xint/1000.0)*alog10(xint/1000.0)*3.0 + $ 
;    a41*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*4.0 +  $ 
;    a51*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*5.0 
  slope = 1.07723

;  print,'Y-intercept of H2-e Curve at T = xint.0: ', 10.0^(a01 +  $ 
;                                                            a11*alog10(xint/1000.0) +  $ 
;                                                            a21*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a31*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a41*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0) +  $ 
;                                                            a51*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0)*alog10(xint/1000.0))
  yint = 2.28029e-21

  
  if (T le 200) THEN BEGIN return, 10.0^( a00 +   $
                                 a10*alog10(T/1000.0) +   $
                                 a20*alog10(T/1000.0)*alog10(T/1000.0) +   $
                                 a30*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
                                 a40*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
                                 a50*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0));
      ENDIF else BEGIN    
          IF (T le xint) THEN BEGIN 
              return, 10.0^(a01 +   $
                            a11*alog10(T/1000.0) +   $
                            a21*alog10(T/1000.0)*alog10(T/1000.0) +   $
                            a31*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
                            a41*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0) +   $
                            a51*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)*alog10(T/1000.0)) ;
          ENDIF ELSE return, 10.0^(slope*(alog10(T/1000.0) - alog10(xint/1000.0)) + alog10(yint)) 
      ENDELSE
END

function clCoolBrem1,T 
CL_B_gm    =(6.022e23*(938.7830/931.494))
CL_Cbremss1= 1.426e-27
CL_al      = 0.79464
CL_bl      = 0.1243
CL_ar      = 2.13164
CL_br      = (-0.1240)

  Tlog10 = alog10(T);
  Tsq = sqrt(T);
     return,Tsq*CL_Cbremss1*(CL_al+CL_bl*Tlog10)*CL_B_gm;
END

function clCoolBrem2,T 
CL_B_gm      =(6.022e23*(938.7830/931.494))
CL_Cbremss1 =1.426e-27
CL_al       =0.79464
CL_bl      = 0.1243
CL_ar      = 2.13164
CL_br      = (-0.1240)
CL_alog4   =0.602059991
CL_alII    =(4.0*(CL_al-CL_bl*CL_alog4))
CL_blII    =(4.0*CL_bl)
CL_arII    =(4.0*(CL_ar-CL_br*CL_alog4))
CL_brII    =(4.0*CL_br)

  Tlog10 = alog10(T);
  Tsq = sqrt(T);

   return,Tsq*CL_Cbremss1*(CL_alII+CL_blII*Tlog10)*CL_B_gm;
END


FUNCTION clCoolLineHeI,T
CL_B_gm      =(6.022e23*(938.7830/931.494))
CL_aHeI =  9.10e-27
CL_bHeI =  1.3179e04
CL_p_HeI = 0.1687
Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5)) ;
CL_MAX_NEG_EXP_ARG = -500

T_inv=1.0/T                     ;
arg = -CL_bHeI*T_inv            ;
;if (arg lt CL_MAX_NEG_EXP_ARG)then  return,0 else
 return, CL_B_gm*CL_aHeI*exp(-CL_bHeI*T_inv)*T_inv^(CL_p_HeI)*Cen_correctn ;
END

FUNCTION clCoolLineHI,T 
CL_B_gm      =(6.022e23*(938.7830/931.494))
CL_aHI  =7.5e-19
CL_bHI = 1.18348e05
  Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));
CL_MAX_NEG_EXP_ARG = -500

  T_inv=1.0/T;
  arg = -CL_bHI*T_inv;
;  if (arg lt CL_MAX_NEG_EXP_ARG) THEN return,0 ELSE 
return,CL_B_gm*CL_aHI*exp( arg )*Cen_correctn;
END

FUNCTION clCoolLineHeII,T 
CL_B_gm      =(6.022e23*(938.7830/931.494))
CL_aHeII  = 5.54e-17
CL_bHeII  = 4.73638e05
CL_p_HeII = 0.397
  Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));
CL_MAX_NEG_EXP_ARG = -500

  T_inv=1.0/T;
  arg = -CL_bHeII*T_inv;
;  if (arg lt CL_MAX_NEG_EXP_ARG) THEN return,0 ELSE 
return,CL_B_gm*CL_aHeII*exp(-CL_bHeII*T_inv)*T_inv^CL_p_HeII*Cen_correctn;
END
