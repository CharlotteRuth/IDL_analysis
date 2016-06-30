PRO LTEcoolingcurve

!X.STYLE = 1
!Y.STYLE = 1
M_H = 1.6719999999999997e-24
Y_H = 0.68649999884516
cm_per_kpc = 3.0857d21
kpc_per_syslength = 1e5
h = 2.60000e-07*cm_per_kpc*kpc_per_syslength
loadct,0

;********************* Read in Metal Cooling Tables *********************

filename_c = '/astro/net/scratch2/christensen/MolecH/cooltable.txt'
filename_h = '/astro/net/scratch2/christensen/MolecH/heattable.txt'
nnH = 121
nHminlog = -9
nHmaxlog = 3
dnH = 0.1
nt = 142
tminlog = 2
tmaxlog = 9.05
dt = 0.05
h_array = findgen(nnH)*dnH + nHminlog
temp_array = findgen(nt)*dt + tminlog
h_table = fltarr(nnH,nt)
for i = 0, N_ELEMENTS(temp_array)-1 do h_table[*,i] = 10^h_array
cooltable = fltarr(nt,nnH)
heattable = fltarr(nt,nnH)

openr,1,filename_c
readf,1,cooltable
close,1
openr,1,filename_h
readf,1,heattable
close,1
cooltable = rotate(cooltable,4)
heattable = rotate(heattable,4)
cooltablefinal = exp(cooltable)*Y_H/M_H*h_table
heattablefinal = exp(heattable)*Y_H/M_H*h_table

;***************************** Fit Phase Diagram **************************

result = fit_curve()
nel = 100.0
minlogrho = -12
maxlogrho = 3 
logrho = findgen(nel)/nel * (maxlogrho - minlogrho) + minlogrho
logT = phasecurve(logrho,result)
set_plot,'x'
window,3
plot,logrho,logT,xtitle='Log (nH)',ytitle = 'Log T'
;******************************* Calculate Abundances ***********************

all = fltarr(9,N_ELEMENTS(logT))
zmetal = fltarr(N_ELEMENTS(logT)) + 0.025
smooth = fltarr(N_ELEMENTS(logT)) + 8.0228202e19 
for i=0LL,LONG64(N_ELEMENTS(logT)) - 1 do  BEGIN
    all[*,i] = LTE_abund(10.0^logT[i],10.0^logrho[i]/5.9753790e+23,zmetal[i],smooth[i])
ENDFOR
Y_HI = REFORM(all[0,*])
Y_H2 = REFORM(all[1,*]) 
Y_HII = REFORM(all[2,*]) 
Y_HeI = REFORM(all[3,*])
Y_HeII = REFORM(all[4,*])
Y_HeIII = REFORM(all[5,*]) 
Y_e = REFORM(all[6,*]) 

;lowT_ind = where(logT lt alog10(9000.0))
;lowT_Y_HII = Y_HII[lowT_ind]
;lowT_Y_e = Y_e[lowT_ind]
;lowT_H = Y_HI[lowT_ind] + 2.0*Y_H2[lowT_ind] + Y_HII[lowT_ind]
;maxfracion = 1
;lowT_Y_HII[where(lowT_Y_e/lowT_H gt maxfracion)] = maxfracion
;lowT_Y_e[where(lowT_Y_e/lowT_H gt maxfracion)] = maxfracion

;Y_e[lowT_ind] = lowT_Y_e
;Y_HII[lowT_ind] = lowT_Y_HII
;all[6,lowT_ind] = lowT_Y_e
;all[2,lowT_ind] = lowT_Y_HII
;all = [[Y_HI],[Y_H2],[Y_HII],[Y_HeI],[Y_HeII],[Y_HeIII],[Y_e],REFORM(all[7,*]),REFORM(all[8,*])]

window,0
;set_plot,'ps'
;device,filename = 'ions.ps',/color,bits_per_pixel= 8,/times
plot,logT,Y_HI,yrange=[1e-4,1],/ylog,xtitle = 'Log(T) K'
oplot,logT,Y_H2,linestyle = 2
oplot,logT,Y_HII,linestyle = 3
oplot,logT,Y_HeI,linestyle = 0,color = 100
oplot,logT,Y_HeII,linestyle = 2,color = 100
oplot,logT,Y_HeIII,linestyle = 3,color = 100
oplot,logT,Y_e,linestyle = 1,color = 150
legend,['Y_HI','Y_H2','Y_HII','Y_HeI','Y_HeII','Y_HeIII','e'],color=[255,255,255,100,100,100,150], $
  linestyle = [0,2,3,0,2,3,1],/bottom,/right
;device,/close

;********************************** Calculate clEdot *********************
;edot_all = fltarr(10,N_ELEMENTS(logT))
;for i=0LL,LONG64(N_ELEMENTS(logT)) - 1 do  BEGIN
;    temp = MIN(abs(h_array - logrho[i]),irho)
;    temp = MIN(abs(temp_array - logT[i]),it)
;    edot_all[*,i] = EdotInstant(10.0^logT[i],10.0^logrho[i],all[*,i],cooltablefinal[it,irho],heattablefinal[it,irho],Zmetal[i])
;ENDFOR
;edot = REFORM(edot_all[0,*])
;edot_metal = REFORM(edot_all[1,*])
;edot_phot = REFORM(edot_all[2,*])
;edot_heat = REFORM(edot_all[3,*])
;edot_cool = REFORM(edot_all[4,*])
;edot_LineHI = REFORM(edot_all[5,*])
;edot_CollHI = REFORM(edot_all[6,*])
;edot_Brem = REFORM(edot_all[7,*])
;edot_Rad = REFORM(edot_all[8,*])
;edot_BremS = REFORM(edot_all[9,*])

;window,1
;set_plot,'ps'
;device,filename = 'cooling.ps',/color,bits_per_pixel= 8,/times
;plot,logT,edot_cool/10.0^logrho,xtitle='Log(T) K',/ylog,ytitle = 'Change in Energy'
;oplot,logT,edot_LineHI/10.0^logrho,linestyle = 1
;oplot,logT,edot_CollHI/10.0^logrho,linestyle = 2
;oplot,logT,edot_Brem/10.0^logrho,linestyle = 3
;oplot,logT,edot_Rad/10.0^logrho,linestyle = 4
;legend,['Total Cooling','HI Line', 'HI Collisional','Brems','Radiative'],linestyle = [0,1,2,3,4],/top,/left
;device,/close
stop
END

FUNCTION fit_curve,filename = filename
IF NOT KEYWORD_SET(filename) THEN filename = "/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_UV_solar_soft/MW_disk.00005"
rtipsy,filename,h,g,d,s
msol_per_sysmass = 1.36e17
kpc_per_syslength = 1e5
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3
total_H = fltarr(N_ELEMENTS(g))
total_He = fltarr(N_ELEMENTS(g))
total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
IF(comple[0] ne -1) THEN total_He[comple] = (-0.446*(g2[comple].zmetal - 0.1)/0.9 + 0.446)/4.0 ;
total_H = 1.0 - total_He*4.0 - g.zmetal ; 
rho_particle = alog10(g.dens * dens_convert * total_H)

x_data = alog10(g[where(g.tempg gt 100)].dens*dens_convert*total_H)
y_data = alog10(g[where(g.tempg gt 100)].tempg)

RETURN,poly_fit(x_data,y_data,10)
END

FUNCTION phasecurve,x_fit,result
curve = result[0] + x_fit*result[1] + x_fit*x_fit*result[2] + x_fit*x_fit*x_fit*result[3] + x_fit*x_fit*x_fit*x_fit*result[4] + x_fit*x_fit*x_fit*x_fit*x_fit*result[5] + x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[6]+ x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[7] + x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[8]+ x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[9]+ x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[10]
highrho = where(x_fit ge -5.5)
lowrho = where(x_fit lt -5.5)
slope = (curve[highrho[0]] - curve[highrho[1]])/(x_fit[highrho[0]] - x_fit[highrho[1]])
curve[lowrho] = slope*(x_fit[lowrho] - x_fit[highrho[0]]) + curve[highrho[0]]
RETURN,curve
END

PRO makephaseD
set_plot,'x'
file = '/astro/net/scratch1/abrooks/FABIO/MW1.1024g1bwK/MW1.1024g1bwK.00512/MW1.1024g1bwK.00512'
file = '/astro/net/scratch2/christensen/MolecH/cosmoMW.std'
cm_per_kpc = 3.0857d21
gm_per_msol = 1.989d33
amu_per_gm = 6.022142d23
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)

msol_per_sysmass = 3.171e15
kpc_per_syslength = 2.85714e4
lengthunit =  kpc_per_syslength*cm_per_kpc;system length unit in cm (=1kpc)
massunit   = msol_per_sysmass*gm_per_msol;system mass unit in gm
dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/ kpc_per_syslength^3/ cm_per_kpc^3 ;massunit*amu_per_gm/molec_weight/lengthunit^3 ;Converts to grams/cm^3

rtipsy,file,h,g,d,s
;g = g[0:10000]
;plot,g.dens*dens_convert,g.tempg,psym = 3,/xlog,/ylog,xtitle = "Density",ytitle = 'Temperature"

contour_plus,alog10(g.dens*dens_convert),alog10(g.tempg),$
  xtitle = "Log(Density)",ytitle = 'Log(Temperature)',xmin = -11,ymin = 2,xmax = 3, ymax = 7, $
  xbinsize = 0.1,ybinsize = 0.05,nlevels = 20,threshold = 10000
END
