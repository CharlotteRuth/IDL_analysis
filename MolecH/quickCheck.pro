;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g7HBWK_00480'
;file ='steps.noFB/h516.cosmo25cmb.1536g8HBWK_noFB.00408.00014.dir/h516.cosmo25cmb.1536g8HBWK_noFB.00408.00014.halo.1'

;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g8HBWK'
;file ='steps/h516.cosmo25cmb.1536g9HBWK.00300.dir/h516.cosmo25cmb.1536g9HBWK.00300.halo.1'


;setsphere 2 -1.484926e-02 5.009335e-02 7.029962e-03 5e-4
;y
;abox 2
;angleup 2 gas
;rotate up 90
;viewgas logrho 0 7

pro quickCheck,file,starlog = starlog,outplot = outplot,pfile = pfile
loadct,39
formatplot,outplot = outplot
cm_in_kpc = 3.08568025e21
ZSOLAR    = 0.0130215
;starlog = 'h516.cosmo25cmb.1536g8HBWK.JeansSF.newShear.00408.starlog'
;starlog = 'h516.cosmo25cmb.1536g8HBWK.JeansSF.shield.00408.starlog'
;starlog ='h516.cosmo25cmb.1536g8HBWK.JeansSF.T3000.shield.newShear.00408.starlog'
;starlog = 'h516.cosmo25cmb.1536g8HBWK.starlog'
;starlog = 'h516.cosmo25cmb.1536g9HBWK.starlog'
;starlog = 'h516.cosmo25cmb.1536g8HBWK_noFB.00408.starlog'

IF NOT KEYWORD_SET(pfile) THEN BEGIN
    spawn,'ls *.param',pfilelist
    pfile = pfilelist[0]
ENDIF
units = tipsyunits(pfile)
rtipsy,file,h,g,d,s
g.dens = g.dens*units.rhounit/h.time/h.time/h.time
readarr,file+'.HI',h,HI,/ascii,part = 'gas'
readarr,file+'.H2',h,H2,/ascii,part = 'gas'
readarr,file+'.lw',h,lw,/ascii,part = 'gas'
readarr,file+'.smoothlength',h,hl,/ascii,part = 'gas'
hl = hl*units.lengthunit*cm_in_kpc
IF (FILE_TEST(file+".correL")) THEN BEGIN
    readarr,file+'.correL',h,correL,/ascii,part = 'gas'
    correL = correL*units.lengthunit*cm_in_kpc
    length = correL
ENDIF ELSE correL = hl
totalH = HI + 2.0*H2
length = hl

minlw = 15 ;20.0
maxlw = 32.0
colors_lw = (alog10(lw) - minlw)*240/(maxlw - minlw)
ind = where(alog10(lw) lt minlw)
if (ind[0] ne -1) THEN colors_lw[ind] = 0
ind = where(alog10(lw) gt maxlw)
if (ind[0] ne -1) THEN colors_lw[ind] = 240

minz = 0
maxz = 0.03 ;0.0015
colors_z = (g.zmetal - minz)*240/(maxz - minz)
ind = where(g.zmetal lt minz)
if (ind[0] ne -1) THEN colors_z[ind] = 0
ind = where(g.zmetal gt maxz)
if (ind[0] ne -1) THEN colors_z[ind] = 240

;minL = 18.0
;maxL = 22.0
;colors_L = (alog10(correL) - minL)/(maxz - minz)*240
;ind = where(alog10(correL) lt minL)
;if (ind[0] ne -1) THEN colors_L[ind] = 0
;ind = where(alog10(correL) gt maxL)
;if (ind[0] ne -1) THEN colors_L[ind] = 240


minT = 1
maxT = 4
colors_T = (alog10(g.tempg) - minT)*240/(maxT - minT)
ind = where(alog10(g.tempg) lt minT)
if (ind[0] ne -1) THEN colors_T[ind] = 0
ind = where(alog10(g.tempg) gt maxT)
if (ind[0] ne -1) THEN colors_T[ind] = 0 ;240

minH = -6
maxH = 0
colors_H = (alog10(H2/totalH*2) - minH)*240/(maxH - minH)
ind = where(alog10(H2/totalH*2) lt minH)
if (ind[0] ne -1) THEN colors_H[ind] = 0
ind = where(alog10(H2/totalH*2) gt maxH)
if (ind[0] ne -1) THEN colors_H[ind] = 240

set_plot,'x'
window,0
!p.multi = [0,3,1]
plot,g.dens,H2/totalH*2,psym = 3,yrange = [1e-6,1],/xlog,xrange = [1e-1,1e4],xtitle = 'Density [amu/cc]',ytitle = 'H_2/H',title = 'Lyman Werner'
FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens],[H2[i],H2[i]]/totalH[i]*2,psym = 3,color = colors_lw[i]
plot,g.dens*length,H2/totalH*2,psym = 3,yrange = [1e-6,1],/xlog,/ylog,xrange = [1e18,1e24],xtitle = 'Surface Density [amu/cm^2]',ytitle = 'H_2/H',title = 'Lyman Werner'
FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens]*length[i],[H2[i],H2[i]]/totalH[i]*2,psym = 3,color = colors_lw[i]
histogramp,alog10(lw),min = minlw, max = maxlw, xtitle = 'Log(LW)',nbins = 100
stop

!p.multi = [0,3,1]
plot,g.dens,H2/totalH*2,psym = 3,yrange = [1e-6,1],/xlog,/ylog,xrange = [1e-1,1e4],xtitle = 'Density [amu/cc]',ytitle = 'H_2/H',title = 'Metallicity'
FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens],[H2[i],H2[i]]/totalH[i]*2,psym = 3,color = colors_z[i]
plot,g.dens*length,H2/totalH*2,psym = 3,yrange = [1e-6,1],/xlog,/ylog,xrange = [1e18,1e24],xtitle = 'Surface Density [amu/cm^2]',ytitle = 'H_2/H',title = 'Metallicity'
FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens]*length[i],[H2[i],H2[i]]/totalH[i]*2,psym = 3,color = colors_z[i]
histogramp,alog10(g.zmetal/ZSOLAR),min = -3, max = alog10(maxz/ZSOLAR), xtitle = 'Metallicity',nbins = 100,/ylog
stop

!p.multi = [0,3,1]
plot,g.dens,H2/totalH*2,psym = 3,yrange = [1e-6,1],xrange = [1e-1,1e4],xtitle = 'Density [amu/cc]',ytitle = 'H_2/H',title = 'Temperature',/xlog;,/ylog,
FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens],[H2[i],H2[i]]/totalH[i]*2,psym = 3,color = colors_T[i]
plot,g.dens*length,H2/totalH*2,psym = 3,yrange = [1e-6,1],xrange = [1e18,1e24],xtitle = 'Surface Density [amu/cm^2]',ytitle = 'H_2/H',title = 'Temperature',/xlog;,/ylog
FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens]*length[i],[H2[i],H2[i]]/totalH[i]*2,psym = 3,color = colors_T[i]
histogramp,alog10(g.tempg),min = mint, max = maxt, xtitle = 'Temperature',nbins = 100
stop

;!p.multi = [0,3,1]
;contour_plus,alog10(g.dens),alog10(H2/totalH*2),yrange = [-6,0],xrange = [-1,4],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),threshold = 10,levels = [5,10,20,30,40,50]
;plot,g.dens,H2/totalH*2,psym = 3,yrange = [1e-6,1],xrange = [1e-1,1e4],xtitle = 'Density [amu/cc]',ytitle = 'H_2/H',title = 'Temperature',/xlog;,/ylog,
;FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens],[H2[i],H2[i]]/totalH[i]*2,psym = 3,color = colors_T[i]
;plot,g.dens*length,H2/totalH*2,psym = 3,yrange = [1e-6,1],xrange = [1e18,1e24],xtitle = 'Surface Density [amu/cm^2]',ytitle = 'H_2/H',title = 'Temperature',/xlog;,/ylog
;FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens]*length[i],[H2[i],H2[i]]/totalH[i]*2,psym = 3,color = colors_T[i]
;stop

!p.multi = [0,2,1]
plot,alog10(g.dens),H2/totalH*2,psym = 3,yrange = [0,1],xrange = [-1,4],xtitle = 'Log(Density) [amu/cc]',ytitle =  textoidl('Log(H_2/H)')
;contour_plus,alog10(g.dens),H2/totalH*2,yrange = [0,1],xrange = [-1,4],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),threshold = 40,levels = [40,80,120,160],/nofill
plot,alog10(g.dens*length),alog10(H2/totalH*2),yrange = [-6,0],xrange = [18,24],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),psym = 3
;contour_plus,alog10(g.dens*length),alog10(H2/totalH*2),yrange = [-6,0],xrange = [18,24],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),threshold = 20,levels = [20,100,200,400,600,800,1000],/nofill
stop

!p.multi = 0
plot,g.dens,H2/totalH*2,psym = 3,yrange = [1e-6,1],/xlog,/ylog,xrange = [1e-2,1e3],xtitle = 'Density [amu/cc]',ytitle = 'H_2/H',title = 'Correlation Length'
;FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens],[H2[i],H2[i]]/totalH[i]*2,psym = 3,color = colors_L[i]

window,1
!p.multi = [0,2,1]
plot,g.dens,correL,psym = 3,/xlog,/ylog,xrange = [1e-2,1e3],xtitle = 'Dens [amu/cc]',ytitle = 'Length [cm]',yrange = [1e19,1e22]
oplot,g.dens,hl,psym = 3, color = 240
oplot,[0.01,1000],[1e21/0.01,1e21/1000]

plot,g.dens,g.dens*correL,psym = 3,/xlog,/ylog,xrange = [1e-2,1e3],xtitle = 'Dens [amu/cc]',ytitle = 'Surface Density [amu/cm^2]',yrange = [1e18,1e23]
oplot,g.dens,g.dens*hl,psym = 3, color = 240
oplot,[0.01,1000],[1e21,1e21]
oplot,[1,1],[1e18,1e23],linestyle = 2
oplot,[10,10],[1e18,1e23],linestyle = 2
oplot,[100,100],[1e18,1e23],linestyle = 2

window,2
!p.multi = 0
;device,filename = file+'phase.eps',/color,bits_per_pixel= 8,/times
plot,g.dens,g.tempg,/xlog,/ylog,psym = 3,xrange = [1e-6,1e6],yrange = [10,1e7],xstyle = 1,xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]'
FOR i = 0L, N_ELEMENTS(g.dens) - 1 DO oplot,[g[i].dens,g[i].dens],[g[i].tempg,g[i].tempg],psym = 3,color = colors_H[i]
IF keyword_set(starlog) THEN BEGIN
    starform = rstarlog(starlog,/MOLECULARH)
    oplot,starform.rhoform*units.rhounit,starform.tempform,psym = 2,color = 240
    window,3
    sfr,s,timeunit  = units.timeunit,massunit = units.massunit,binsize = 1e7
    stop
ENDIF

;window,3
;!p.multi = 0
;y_H2 = histogram(alog10(H2/totalH*2),locations = x,max = 1, min = -6,nbins = 100)
;plot,x,FLOAT(y_H2)/N_ELEMENTS(g),psym = 10,title = file
;stop
end

