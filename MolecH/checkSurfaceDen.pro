pro checkSurfaceDen,file,pfile = pfile,LTEarr = LTEarr,outplot = outplot
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790d+23
molec_weight = (0.76*1 + 0.24*4.0)
zsolar = 0.0130215
f_H = 0.764
;amucm2smpc = 7.94938572e-21
amucm2smpc = 1/gm_per_msol/amu_per_gm*cm_per_kpc*cm_per_kpc/1000/1000

loadct,39
formatplot,outplot = outplot
IF KEYWORD_SET(outplot) THEN BEGIN
    fgcolor = 0
    l_charsize = 0.75
ENDIF ELSE BEGIN
    l_charsize = 1.0
    fgcolor = 255
ENDELSE
set_plot,'x'
cm_in_kpc = 3.08568025e21
;starlog = 'h516.cosmo25cmb.1536g8HBWK.JeansSF.newShear.00408.starlog'
;starlog = 'h516.cosmo25cmb.1536g8HBWK.JeansSF.shield.00408.starlog'
;starlog ='h516.cosmo25cmb.1536g8HBWK.JeansSF.T3000.shield.newShear.00408.starlog'
;starlog = 'h516.cosmo25cmb.1536g8HBWK.starlog'
;starlog = 'h516.cosmo25cmb.1536g9HBWK.starlog'
;starlog = 'h516.cosmo25cmb.1536g8HBWK_noFB.00408.starlog'

IF (NOT KEYWORD_SET(pfile)) THEN BEGIN
    spawn,'ls *.param',pfilelist
    pfile = pfilelist[0]
ENDIF
units = tipsyunits(pfile)
rtipsy,file,h,g,d,s
readarr,file+'.HI',h,HI,/ascii,part = 'gas'
readarr,file+'.H2',h,H2,/ascii,part = 'gas'
IF KEYWORD_SET(LTEarr) THEN BEGIN
    readarr,file+'.HI_LTEc10',h,HI_LTE,/ascii,part = 'gas'
    readarr,file+'.H2_LTEc10',h,H2_LTE,/ascii,part = 'gas'
ENDIF
readarr,file+'.lw',h,lw,/ascii,part = 'gas'
readarr,file+'.smoothlength',h,hl,/ascii,part = 'gas'
;readarr,file+'.correL',h,correL,/ascii,part = 'gas'
correL = hl
g.dens = g.dens*units.rhounit/h.time/h.time/h.time
correL = correL*units.lengthunit*cm_in_kpc
hl = hl*units.lengthunit*cm_in_kpc
totalH = (HI + 2.0*H2)*0 + 0.761192
marray = g.mass*units.massunit;amu_per_gm*gm_per_msol
hmarray = (HI+2.0*H2)*g.mass*units.massunit
hImarray = (HI)*g.mass*units.massunit
h2marray = (2.0*H2*g.mass)*units.massunit
IF KEYWORD_SET(LTEarr) THEN BEGIN
    hImarray_LTE = (HI_LTE*g.mass)*units.massunit
    h2marray_LTE = (2.0*H2_LTE*g.mass)*units.massunit
ENDIF
ind = where(g.tempg lt 1e4)
zmetal = TOTAL(g[ind].zmetal*g[ind].mass)/TOTAL(g[ind].mass)
print,'Total H2 Mass: ',TOTAL((2.0*H2)*g.mass*units.massunit)*f_H,' Solar Masses'
print,'Total HI + H2 Mass: ',TOTAL((HI + 2.0*H2)*g.mass*units.massunit)*f_H,' Solar Masses'
print,'Mean Metallicity: ',zmetal/0.0177 ,' Zsol'



!p.multi = [0,2,1]
window,1,xsize = 712,ysize = 392
plot,g.dens,correL/(0.75/!PI*g.mass*units.massunit*gm_per_msol*amu_per_gm/g.dens)^(1.0/3.0),psym = 3,/ylog,/xlog,xtitle = 'Density [amu/cc]',ytitle = 'Ratio',xrange = [1e-1,1e3];,xrange = [1e21,4e22],yrange = [1e21,4e22]
oplot,g.dens,hl/(0.75/!PI*g.mass*units.massunit*gm_per_msol*amu_per_gm/g.dens)^(1.0/3.0),psym = 3,color = 240
;oplot,[1e21,4e22],[1e21,4e22]
length = hl

xsigmalow = findgen(300)/100 - 1.
xsigma = 10.0^(findgen(400)/100 + 1.)
ysigma=2.5e-4*xsigma^1.4
ysigma1=2.5e-4*xsigma^(1.4 + 0.15)
ysigma2=2.5e-4*xsigma^(1.4 - 0.15)
ysigmalow = xsigmalow*2.4 - 5.0
range = 5.0;10.0
delta = 0.0400 ;0.750 ;0.750;/dKpcUnit

g.x = g.x *h.time*units.lengthunit
g.y = g.y *h.time*units.lengthunit
g.z = g.z *h.time*units.lengthunit
rg = sqrt(g.x*g.x + g.y*g.y)
s.x = s.x *h.time*units.lengthunit
s.y = s.y *h.time*units.lengthunit
s.z = s.z *h.time*units.lengthunit
xmin  = -1.0*range
xmax  =  1.0*range
ymin  = -1.0*range 
ymax  =  1.0*range
nx = (xmax - xmin)/delta
ny = (ymax - ymin)/delta    
grid_sd = fltarr(nx,ny)
grid_sd_ave = fltarr(nx,ny)
grid_sd_msol = fltarr(nx,ny)
grid_sd_ave_msol = fltarr(nx,ny)
grid_frac = grid_sd
IF KEYWORD_SET(LTEarr) THEN grid_frac_LTE = grid_sd
grid_z = grid_sd
grid_sfr = grid_sd
xarray = findgen(nx + 1)*delta + xmin
yarray = findgen(ny + 1)*delta + ymin
deltat = 100.e6 ;Appearently, in Bigiel, this is based on the FUV flux, so OB stars with have a lifetime of 10^8 years

timeunit = 3.872d3*SQRT((units.lengthunit*3.0856805d21)^3/(units.massunit*1.98892d33))/31556926.0
tform=s.tform*(timeunit)[0]
tcurrent = max(tform)           ;1e10
tclip = tcurrent - deltat

minr = 0
maxr = range
nbins = 100
r_bin = findgen(nbins)*range/(nbins)
y_H2 = weighted_HISTOGRAM(rg,input=H2marray,binsize=range/nbins);,min=minr,max=maxr);*amu_per_gm*gm_per_msol
y_HI = weighted_HISTOGRAM(rg,input=HImarray,binsize=range/nbins);,min=minr,max=maxr);*amu_per_gm
IF KEYWORD_SET(LTEarr) THEN BEGIN
    y_H2_LTE = weighted_HISTOGRAM(rg,input=H2marray_LTE,binsize=range/nbins,min=minr,max=maxr) ;*amu_per_gm
    y_HI_LTE = weighted_HISTOGRAM(rg,input=HImarray_LTE,binsize=range/nbins,min=minr,max=maxr) ;*amu_per_gm
ENDIF
dr = range/(nbins)
area = ((r_bin + dr)*(r_bin + dr) - r_bin*r_bin)*!PI*cm_per_kpc*cm_per_kpc
plot,r_bin,y_HI*amu_per_gm*gm_per_msol/area,psym=10,xtitle="Radius [kpc]",/ylog;,ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' pc^2')+"]" ,title='Gas Surface Density Profile -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years',/ylog,yrange = [1e2,3e7]
oplot,r_bin,y_H2*amu_per_gm*gm_per_msol/area,linestyle = 2,psym = 10
IF KEYWORD_SET(LTEarr) THEN BEGIN
    oplot,r_bin,y_HI_LTE*amu_per_gm*gm_per_msol/area,linestyle = 2,psym = 10,color = 240
    oplot,r_bin,y_H2_LTE*amu_per_gm*gm_per_msol/area,linestyle = 2,psym = 10,color = 240
ENDIF
stop
!p.multi = 0
FOR ix = 0, nx - 1 DO BEGIN
    FOR iy = 0, ny - 1 DO BEGIN
        ind = where(g.x gt xarray[ix] AND g.x lt xarray[ix + 1] AND g.y gt yarray[iy] AND g.y lt yarray[iy + 1])
        inds = where(s.x gt xarray[ix] AND s.x lt xarray[ix + 1] AND s.y gt yarray[iy] AND s.y lt yarray[iy + 1] AND tform gt tclip) 
        inds_new = where(tform gt tclip) 
        window,2,xsize = 392, ysize = 392
        plot,g.x,g.y,psym = 3,xrange = [-1.0*range,range],yrange = [-1*range,range], xtitle = 'X [kpc]', ytitle = 'Y [kpc]',xstyle = 1, ystyle = 1
        oplot,s[inds_new].x,s[inds_new].y,psym = 2, color = 50,symsize = 0.5
        if N_ELEMENTS(ind) gt 1 THEN  oplot,g[ind].x,g[ind].y,psym = 3,color = 240
        if N_ELEMENTS(inds) gt 1 THEN oplot,s[inds].x,s[inds].y,psym = 2,color = 100,symsize = 0.5       
        if ind[0] NE -1 THEN BEGIN
            grid_frac[ix,iy] = TOTAL(h2marray[ind])/TOTAL(hmarray[ind])
            IF KEYWORD_SET (LTEarr) THEN grid_frac_LTE[ix,iy] = TOTAL(h2marray_LTE[ind])/TOTAL(hmarray[ind])
            grid_sd[ix,iy] = TOTAL(hmarray[ind])*amu_per_gm*gm_per_msol/delta/delta/cm_per_kpc/cm_per_kpc;/units.lengthunit/units.lengthunit/delta[i]/delta[i]  
            grid_sd_msol[ix,iy] = TOTAL(hmarray[ind])/delta/delta/1000.0/1000.0;/units.lengthunit/units.lengthunit/delta[i]/delta[i] 
            grid_z[ix,iy] = TOTAL(g[ind].zmetal*g[ind].mass)/TOTAL(g[ind].mass)

            window,0,xsize = 712,ysize = 392
;            plot,g[ind].dens*length[ind],H2[ind]/totalH[ind]*2,psym = 2,yrange = [1e-6,1],/xlog,/ylog,xtitle = 'Density [Msol/pc^2]',ytitle = 'H_2/H',xrange = [1e-3,100];amucm2smpc
            plot,g[ind].dens*length[ind],H2[ind]/totalH[ind]*2,psym = 2,yrange = [1e-6,1],/xlog,/ylog,xtitle = 'Surface Density [amu/cm^2]',ytitle = 'H_2/H',xrange = [1e19,1e23];amucm2smpc
            loadct,0
            IF KEYWORD_SET (LTEarr) THEN oplot,g[ind].dens*correL[ind],H2_LTE[ind]/totalH[ind]*2,psym = 2,color = 100
            loadct,39
            ave_sd = MEAN(g[ind].dens*length[ind])
            ave_sd_msol = MEAN(g[ind].dens*length[ind])*amucm2smpc
            grid_sd_ave[ix,iy] = ave_sd
            grid_sd_ave_msol[ix,iy] = ave_sd_msol
;            oplot,grid_sd_ave,grid_frac,psym = 2,color = 70
            oplot,grid_sd,grid_frac,psym = 2,color = 240
;            oplot,grid_sd_ave,grid_frac_LTE,psym = 2,color = 120
            IF KEYWORD_SET (LTEarr) THEN oplot,grid_sd,grid_frac_LTE,psym = 2,color = 200
;            oplot,[ave_sd,ave_sd],[grid_frac[ix,iy],grid_frac[ix,iy]],psym = 5,color = 100
            oplot,[ave_sd,ave_sd],[1e-6,1],linestyle = 2
            oplot,[grid_sd[ix,iy],grid_sd[ix,iy]],[1e-6,1]
            oplot,[1e19,1e23],[grid_frac[ix,iy],grid_frac[ix,iy]]

            window,5,xsize = 712,ysize = 392
            plot,g[ind].dens*length[ind]*amucm2smpc,H2[ind]/totalH[ind]*2,psym = 2,yrange = [1e-6,1],/xlog,/ylog,xtitle = 'Density [Msol/pc^2]',ytitle = 'H_2/H',xrange = [1e-3,100];amucm2smpc
;            plot,g[ind].dens*length[ind],H2[ind]/totalH[ind]*2,psym = 2,yrange = [1e-6,1],/xlog,/ylog,xtitle = 'Surface Density [amu/cm^2]',ytitle = 'H_2/H',xrange = [1e19,1e23];amucm2smpc
            loadct,0
            IF KEYWORD_SET (LTEarr) THEN oplot,g[ind].dens*correL[ind]*amucm2smpc,H2_LTE[ind]/totalH[ind]*2,psym = 2,color = 100
            loadct,39
;            oplot,grid_sd_ave,grid_frac,psym = 2,color = 70
            oplot,grid_sd_msol,grid_frac,psym = 2,color = 240
;            oplot,grid_sd_ave,grid_frac_LTE,psym = 2,color = 120
            IF KEYWORD_SET (LTEarr) THEN oplot,grid_sd_msol,grid_frac_LTE,psym = 2,color = 200
;            oplot,[ave_sd,ave_sd],[grid_frac[ix,iy],grid_frac[ix,iy]],psym = 5,color = 100
            oplot,[ave_sd_msol,ave_sd_msol],[1e-6,1],linestyle = 2
            oplot,[grid_sd_msol[ix,iy],grid_sd_msol[ix,iy]],[1e-6,1]
            oplot,[1e-4,1e2],[grid_frac[ix,iy],grid_frac[ix,iy]]

        ENDIF ELSE BEGIN
            grid_frac[ix,iy] = 0
            grid_sd[ix,iy] = 0
            grid_z[ix,iy] = 0
        ENDELSE
        if inds[0] NE -1 THEN BEGIN
            grid_sfr[ix,iy] = TOTAL(s[inds].mass*units.massunit)/delta/delta/deltat
        ENDIF ELSE BEGIN
            grid_sfr[ix,iy] = 0
        ENDELSE
        window,3,xsize = 392, ysize = 392
        plot,[-5,-5],[-5,-5],xrange = [-1,5],yrange = [-5,3],psym = 2,xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]",xstyle = 1,ystyle = 1
        oplot,[1,1],[-5,3],linestyle = 1
        oplot,xsigmalow,ysigmalow
        oplot,alog10(xsigma),alog10(ysigma)
        oplot,alog10(grid_sd),alog10(grid_sfr),psym = 2
        oplot,[alog10(grid_sd_msol[ix,iy]),alog10(grid_sd_msol[ix,iy])],[alog10(grid_sfr[ix,iy]),alog10(grid_sfr[ix,iy])],psym = 2,color = 240
;        wait, 2
;        stop
    ENDFOR
ENDFOR    
plot,grid_sd,grid_frac,/ylog,/xlog,psym = 1,xrange = [1e19,1e23],yrange = [1e-6,1],xtitle = 'Density [Msol/pc^2]',ytitle = textoidl('H_2/H')
    readcol,'/astro/users/christensen/code/MolecH/Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
    readcol,'/astro/users/christensen/code/MolecH/Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'
    oplot,nhw,10^logfH2w,psym = 4
    oplot,nhg,10^logfH2g,psym = 6
    legend,['Simulated Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [1,4,6],/right,/bottom,charsize = l_charsize
end
