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

pro quickCheckPaper,file,outplot = outplot
loadct,39
IF KEYWORD_SET(outfile) THEN BEGIN
    fgcolor = 0
    l_charsize = 0.75
ENDIF ELSE BEGIN
    l_charsize = 1.0
    fgcolor = 255
ENDELSE

formatplot,outplot = outplot
cm_in_kpc = 3.08568025e21
zsolar = 0.0130215
;starlog = 'h516.cosmo25cmb.1536g8HBWK.JeansSF.newShear.00408.starlog'
;starlog = 'h516.cosmo25cmb.1536g8HBWK.JeansSF.shield.00408.starlog'
;starlog ='h516.cosmo25cmb.1536g8HBWK.JeansSF.T3000.shield.newShear.00408.starlog'
;starlog = 'h516.cosmo25cmb.1536g8HBWK.starlog'
;starlog = 'h516.cosmo25cmb.1536g9HBWK.starlog'
;starlog = 'h516.cosmo25cmb.1536g8HBWK_noFB.00408.starlog'

spawn,'ls *.param',pfilelist
pfile = pfilelist[0]
units = tipsyunits(pfile)
rtipsy,file,h,g,d,s
g.dens = g.dens*units.rhounit/h.time/h.time/h.time
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,filename=outplot+'_phase.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2
ENDIF ELSE window,3
!p.multi = 0
;plot,g.dens,g.tempg,/xlog,/ylog,psym = 3,xrange = [1e-6,1e6],yrange = [10,1e7],xstyle = 1,xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]'
contour_plus,alog10(g.dens),alog10(g.tempg),xrange = [-6,4],yrange = [1,7],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(Temperature) [K]'),threshold = 100,levels = [100,200,400,800,1600,3200,6400],/nofill


ind = where(g.dens gt 1e-3 and g.tempg lt 2e4)
g = g[ind]
readarr,file+'.HI',h,HI,/ascii,part = 'gas'
readarr,file+'.H2',h,H2,/ascii,part = 'gas'
readarr,file+'.lw',h,lw,/ascii,part = 'gas'
readarr,file+'.smoothlength',h,hl,/ascii,part = 'gas'
hl = hl*units.lengthunit*cm_in_kpc
IF (FILE_TEST(file+".correL")) THEN BEGIN
    readarr,file+'.correL',h,correL,/ascii,part = 'gas'
    correL = correL*units.lengthunit*cm_in_kpc
    length = correL
    correL = correL[ind]
;    print,'here'
ENDIF
totalH = HI + 2.0*H2
length = hl
HI = HI[ind]
H2 = H2[ind]
lw = lw[ind]
hl = hl[ind]
totalH = totalH[ind]
length = length[ind]


IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename=outplot+'_smoothL.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  
    ENDIF ELSE window,0
!p.multi = [0,2,1]
plot,alog10(g.dens),H2/totalH*2,psym = 3,yrange = [0,1],xrange = [-1,4],xtitle = 'Log(Density) [amu/cc]',ytitle =  textoidl('Log(H_2/H)')
;contour_plus,alog10(g.dens),H2/totalH*2,yrange = [0,1],xrange = [-1,4],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),threshold = 40,levels = [40,80,120,160],/nofill
plot,alog10(g.dens*length),alog10(H2/totalH*2),yrange = [-6,0],xrange = [18,24],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),psym = 3
;contour_plus,alog10(g.dens*length),alog10(H2/totalH*2),yrange = [-6,0],xrange = [18,24],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),threshold = 20,levels = [20,100,200,400,600,800,1000],/nofill

IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename=outplot+'_phase.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2
ENDIF ELSE window,4
y = histogram(g.zmetal,locations = x,min = 0, max = 0.028,nbins = 500)
plot,x,y,/ylog,yrange = [1,1e4],xtitle = 'Metallicity',ytitle = 'Distribution'

IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename=outplot+'_correL.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 
ENDIF ELSE window,1
!p.multi = [0,2,1]
plot,alog10(g.dens),H2/totalH*2,psym = 3,yrange = [0,1],xrange = [-1,4],xtitle = 'Log(Density) [amu/cc]',ytitle =  textoidl('Log(H_2/H)')
;contour_plus,alog10(g.dens),H2/totalH*2,yrange = [0,1],xrange = [-1,4],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),threshold = 40,levels = [40,80,120,160],/nofill
plot,alog10(g.dens*correL),alog10(H2/totalH*2),yrange = [-6,0],xrange = [18,24],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),psym = 3
;contour_plus,alog10(g.dens*length),alog10(H2/totalH*2),yrange = [-6,0],xrange = [18,24],xtitle = 'Log(Density) [amu/cc]',ytitle = textoidl('Log(H_2/H)'),threshold = 20,levels = [20,100,200,400,600,800,1000],/nofill

IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename=outplot+'_lengths.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 
ENDIF ELSE window,2
!p.multi = [0,2,1]
plot,g.dens,correL,psym = 3,/xlog,/ylog,xrange = [1e-2,1e3],xtitle = 'Dens [amu/cc]',ytitle = 'Length [cm]',yrange = [1e19,1e22]
oplot,g.dens,hl,psym = 3, color = 240
oplot,[0.01,1000],[1e21/0.01,1e21/1000]
legend,['Correlation Length','Smoothing Length'],linestyle = [1,1],color = [fgcolor,240],charsize = l_charsize,/right

plot,g.dens,g.dens*correL,psym = 3,/xlog,/ylog,xrange = [1e-2,1e3],xtitle = 'Dens [amu/cc]',ytitle = 'Surface Density [amu/cm^2]',yrange = [1e18,1e23]
oplot,g.dens,g.dens*hl,psym = 3, color = 240
oplot,[0.01,1000],[1e21,1e21]
oplot,[1,1],[1e18,1e23],linestyle = 2
oplot,[10,10],[1e18,1e23],linestyle = 2
oplot,[100,100],[1e18,1e23],linestyle = 2



IF (KEYWORD_SET(outplot)) THEN  device,/close ELSE stop
end

