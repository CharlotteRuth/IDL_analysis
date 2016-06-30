;cd,'/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_repl/z1_lw3e7_cp10'
;file = 'MW_disk.00010'

;H2surface_den_column,dir,file
PRO H2frac_column,dir,file,outplot = outplot
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790d+23
molec_weight = (0.76*1 + 0.24*4.0)
loadct,39
f_H = 0.764
zsol = 0.0177                   ;0.025

dKpcUnit        = 1e5
dMsolUnit       = 1.36e17
;dMsolUnit       = 2.310e15
;dKpcUnit        = 25000.
dens_convert =  dMsolUnit*gm_per_msol*amu_per_gm/dKpcUnit^3/cm_per_kpc^3

cd,dir
rtipsy,file,h,g,d,s
readcol,file+'.H2',h2,/silent
readcol,file+'.HI',hI,/silent 
readcol,file+'.shear',mach,/silent
HI = HI[1:h.ngas]
H2 = H2[1:h.ngas] 
column = 1.0/mach[1:h.ngas]*dKpcUnit*cm_per_kpc  
readarr,file + '.smoothlength',h,smooth,/ascii, part = 'gas'
rho = alog10(g.dens * dens_convert * (HI + 2.0*H2))
zmetal = MEAN(g.zmetal)
print,'Mean Metallicity: ',zmetal/zsol,' Zsol'

IF keyword_set(outplot) THEN BEGIN
    set_plot,'ps'
    !P.CHARSIZE = 1.5             ;3
    !P.THICK=1.5                ;4
    !P.CHARTHICK=1.5            ;4
    !P.font=0
    !X.THICK=1.5                ;4
    !Y.THICK=1.5                ;4
    device,filename = file + '_logSD.eps',/color,bits_per_pixel = 8,/times,xsize = 6,ysize = 4.5,/inch
ENDIF ELSE BEGIN
    set_plot,'x'
    window,0
ENDELSE

plot,alog10(g.dens*dens_convert*(HI + 2.0*H2)*column),2.0*H2/(HI + 2.0*H2),psym = 3,/ylog,xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.0],xrange = [17,24],xstyle=1,ystyle=1,symsize =0.5

sprad = sqrt(g.x*g.x + g.y*g.y + g.z*g.z)
ind = sprad lt 0.000044
radii = SQRT(g.x*g.x + g.y*g.y)*dKpcUnit
nel = 100.
maxr = 5.
dr = maxr/nel
rarray = findgen(nel)/nel*(maxr - dr)
marray = dblarr(nel)
hmarray = dblarr(nel)
h2marray = dblarr(nel)
aarray = !PI*((rarray + dr)*(rarray + dr) - rarray*rarray)*cm_per_kpc*cm_per_kpc
darray = dblarr(nel)
FOR i = 0, nel - 1 DO BEGIN
    ind = where(radii lt rarray[i]+dr AND radii ge rarray[i])
    IF (ind[0] NE -1) THEN BEGIN
        marray[i] = TOTAL(g[ind].mass)*amu_per_gm*gm_per_msol*dMsolUnit 
        hmarray[i] = TOTAL((HI[ind]+2.0*H2[ind])*g[ind].mass)*amu_per_gm*gm_per_msol*dMsolUnit
        h2marray[i] = TOTAL(2.0*H2[ind]*g[ind].mass)*amu_per_gm*gm_per_msol*dMsolUnit
        darray[i] = MEAN(g[ind].dens)*dens_convert
    ENDIF ELSE BEGIN
        marray[i] = 0
        hmarray[i] = 0
    ENDELSE
ENDFOR
sden = marray/aarray
shden = hmarray/aarray
oplot,alog10(shden),h2marray/hmarray,psym = 4,color = 100

IF 1 THEN BEGIN
    zmetal = 0.025
    readcol,'/astro/users/christensen/code/MolecH/Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
    readcol,'/astro/users/christensen/code/MolecH/Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'
    oplot,nhw,10^logfH2w,psym = 1,color = 240
    oplot,nhg,10^logfH2g,psym = 4,color = 240
    legend,['Simulated Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [3,4,1],color=[0,240,240],/right,/bottom,charsize = 1.2
ENDIF

IF keyword_set(outplot) THEN BEGIN
    device,/close
    device,filename = file + '_krumholtz.eps',/color,bits_per_pixel = 8,/times,xsize = 6,ysize = 4.5,/inch
ENDIF ELSE window,1
surf_den_convert = 1.258075e20
plot,(g.dens*dens_convert*(HI + 2.0*H2)*column)/surf_den_convert,2.0*H2/(HI + 2.0*H2),psym = 3,xtitle = textoidl("N_{HI} + 2N_{H_2} (M_{\odot} pc^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [0,1.0],xrange = [0.1,1e4],xstyle=1,ystyle=1,symsize =0.5,/xlog
;oplot,shden/surf_den_convert,h2marray/hmarray,psym = 4,color = 100
krumholtz,MEAN(g.zmetal),x,y
oplot,x,y,linestyle = 2,color = 50
;krumholz,0.0177*10.0,x,y
;oplot,x,y,color = 210,linestyle = 3 
;krumholz,0.0177/10.0,x,y
;oplot,x,y,color = 210,linestyle = 2
;krumholz,0.0177/100.0,x,y
;oplot,x,y,color = 210,linestyle = 1

IF 0 THEN BEGIN
    oplot,10^nhw/surf_den_convert,10^logfH2w,psym = 1,color = 240
    oplot,10^nhg/surf_den_convert,10^logfH2g,psym = 4,color = 240
    legend,['Simulated Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [3,4,1],color=[0,240,240],/right,/bottom,charsize = 1.2
ENDIF


omega_H2 = 0.2                ;0.2
sigma_d = 2d-21     ;2d-21                 ;4d-21
zmetal = MEAN(g.zmetal)
x = (g.dens*dens_convert*(HI + 2.0*H2)*column)/5d14
h_shield = 1.0 - (1.0 - omega_H2)/(1.0 + x)^2.0 + omega_H2/SQRT(1.0 + x)*exp(-0.00085*SQRT(1.0 + x))
d_shield = 1.0 - exp(-1.0*(zmetal/zsol)*sigma_d*(g.dens*dens_convert*(HI + 2.0*H2)*column))
;oplot,(g.dens*dens_convert*(HI + 2.0*H2)*column)/surf_den_convert,d_shield*h_shield,psym = 3;,color = 150
;oplot,(g.dens*dens_convert*(HI +
;2.0*H2)*column)/surf_den_convert,2.0*H2/(HI + 2.0*H2),psym = 3,color
;= 180
;d_shield = 1.0 - exp(-1.0*sigma_d*(g.dens*dens_convert*(HI + 2.0*H2)*column))
;oplot,(g.dens*dens_convert*(HI + 2.0*H2)*column)/surf_den_convert,d_shield*h_shield,psym = 3,color = 120
;stop

IF keyword_set(outplot) THEN BEGIN
    device,/close
    device,filename = file + '_phaseD.eps',/color,bits_per_pixel = 8,/times,xsize = 6,ysize = 4.5,/inch
ENDIF ELSE window,2
plot,g.dens*dens_convert,g.tempg,psym = 3,/xlog,/ylog,xrange = [1e-8,1e4],yrange = [10,1e7],xtitle = 'Density [amu cm^-3]',ytitle = 'Temperature [K]'
colors = 2.0*H2/(HI + 2.0*H2)/0.75*200 + 20
for i = 0L,N_ELEMENTS(H2) -1 do oplot,[g[i].dens*dens_convert,g[i].dens*dens_convert],[g[i].tempg,g[i].tempg,g[i].tempg,g[i].tempg],color = colors[i]
IF keyword_set(outplot) THEN device,/close
;stop
set_plot,'x'

END

; rphot = 1.5917e-9
;IDL> rdust = 3.5e-17
;IDL> en_B = fH2*rphot/(0.758*rdust*(1 - fH2)

;en_B = fH2*rphot/(0.758*rdust*(1 - fH2)
;IDL> en_B = fH2*rphot/(0.758*rdust*(1 - fH2))
;IDL> plot,en_B,fH2,/xlog,/ylog
;IDL> en_B2 = fH2*rphot*1000.0/(0.758*rdust*(1 - fH2))
;IDL> oplot,en_B2,fH2,linestyle = 2
;IDL> en_B3 = fH2*rphot*0.01/(0.758*rdust*(1 - fH2))
;IDL> oplot,en_B3,fH2,linestyle = 1

