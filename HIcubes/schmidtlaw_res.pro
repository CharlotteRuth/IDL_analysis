;filename = 'h603.cosmo50cmb.2304g5bwK.00512'
;cd,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG'
pro schmidtlaw_res,filename,kpcunit,massunit,massform,res = res,range = range,obsH = obsH,i = i,color = color,psym = psym,overplot = overplot,starlogfile = starlogfile
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
timeunit=SQRT((kpcunit*3.086d21)^3/(6.67d-8*massunit*1.99d33))/(3600.*24.*365.24)
rhounit = massunit*gm_per_msol*H_per_gm/kpcunit^3/cm_per_kpc^3
massform = massform*massunit
loadct,39

IF NOT keyword_set(deltat) THEN deltat=100.e6 ; in yrs
IF NOT keyword_set(psym)   THEN psym = 2
IF NOT keyword_set(i)      THEN i = 90
IF NOT keyword_set(res)    THEN res = 0.75 ;kpc

IF keyword_set(obsH) THEN BEGIN
    rtipsy,filename + '.std',h,g,d,s,/justhead
    rhounit = rhounit/(h.time^3)
;tcurrent = 13.7e9 - wmap3_lbfunc(h.time)
    z = 1/h.time - 1
    tcurrent=13.6958e9-wmap3_lookback(z) ;tcurrent=13.73d9
    tclip=tcurrent-deltat  ; only use starformation events since tclip

    cube_tform = read_cube_fits(filename+'.cube.tform.fits',header_tform)
    mom0_HI = mrdfits(filename + '.mom0.fits',0,header_HI)/amu_per_gm/gm_per_msol*cm_per_kpc*cm_per_kpc/1d6
    IF NOT (where(~finite(cube_tform)))[0] EQ -1 THEN cube_tform[where(~finite(cube_tform))] = 0
    IF NOT (where(~finite(mom0_HI)))[0]    EQ -1 THEN mom0_HI[where(~finite(mom0_HI))] = 0

    spatial_axes = findgen(header_tform.NAXIS1)*header_tform.CDELT1 + header_tform.CRVAL1
    time_axes    = findgen(header_tform.NAXIS3)*header_tform.CDELT3 + header_tform.CRVAL3
    mintime_axes = time_axes - header_tform.CDELT3/2.0
    maxtime_axes = time_axes + header_tform.CDELT3/2.0
    lastti = (where(maxtime_axes - tcurrent GE 0))[0]
    IF (lastti EQ -1) THEN BEGIN
        lastti = N_ELEMENTS(maxtime_axes) - 1 
        print,'Check time range: ',maxtime_axes[lastti],'>=',tcurrent
    ENDIF
    foo = MIN(abs(time_axes - tclip),startti)
    totalt = tcurrent - mintime_axes[startti]

    dcell = ROUND(res/header_tform.CDELT1)
    dxcell = dcell*header_tform.CDELT1
    ncell = ROUND(header_tform.NAXIS1/dcell)
    cell_axes = findgen(ncell)*dxcell + header_tform.CRVAL1 - header_tform.CDELT1/2.0
    areapc = 1e6*dxcell*dxcell
    areakpc = dxcell*dxcell

    HIrebin = fltarr(ncell,ncell)
    tformrebin = fltarr(ncell,ncell)
    schmidt = fltarr(ncell,ncell)

    FOR ix = 0, ncell - 1 DO BEGIN
        FOR iy = 0,ncell - 1 DO BEGIN 
            HIrebin[ix,iy] = MEAN(mom0_HI[ix*dcell : ix*dcell + dcell - 1, iy*dcell :iy*dcell + dcell - 1])
            tformrebin[ix,iy] = TOTAL(cube_tform[ix*dcell : ix*dcell + dcell - 1, iy*dcell :iy*dcell + dcell - 1,startti : lastti])/areakpc/totalt
        ENDFOR
    ENDFOR
ENDIF ELSE BEGIN
    rtipsy,filename + '.std',h,g,d,s
    base = strmid(filename,0,strlen(filename) - 13)
    IF NOT keyword_set(starlogfile) THEN starlogfile = '../../'+base + '.starlog'
    starlog = rstarlog(starlogfile,/molecularH)
    rhounit = rhounit/(h.time^3)
    s.x = s.x*kpcunit*h.time
    s.y = s.y*kpcunit*h.time
    s.z = s.z*kpcunit*h.time
    g.x = g.x*kpcunit*h.time
    g.y = g.y*kpcunit*h.time
    g.z = g.z*kpcunit*h.time
    s.tform = s.tform*timeunit
    tcurrent= max(s.tform)
    tclip=tcurrent-deltat  ; only use starformation events since tclip
    g.mass = g.mass*massunit
    readarr,filename + '.HI',h,HI,/ascii,part = 'gas'
    IF file_test(filename + '.H2') THEN $
      readarr,filename + '.H2',h,H2,/ascii,part = 'gas' ELSE $
      H2 = HI*0
    readarr,filename + '.coolontime',h,coolon,/ascii,part = 'gas'
    readarr,filename + '.iord',h,iordg,/ascii,part = 'gas',type = 'long'
    readarr,filename + '.iord',h,iords,/ascii,part = 'star',type = 'long'
    coolonclean = coolon
    coolonclean = coolonclean - tcurrent/timeunit
    plot,alog10(g.dens*rhounit),alog10(H2*2/(HI + 2.0*H2)),xrange = [0,3],psym = 3
    plot_colorscale,alog10(g.dens*rhounit),alog10(H2*2/(HI + 2.0*H2)),coolonclean,xrange = [0,3],/overplot

    IF NOT keyword_set(range) THEN range = 5
    minrange = range/(-2.0) 
    ncell = round(range/res)
    res = 1.0*range/ncell
    tformrebin = fltarr(ncell,ncell)
    HIrebin = fltarr(ncell,ncell)
    H2rebin = fltarr(ncell,ncell)
    H2rebincool = fltarr(ncell,ncell)
    gasrebin = fltarr(ncell,ncell)
    FOR ix = 0, ncell - 1 DO BEGIN
        FOR iy = 0,ncell - 1 DO BEGIN 
            inds = where((s.x GT (res*ix + minrange)) AND (s.x LE (res*(ix+1) + minrange)) AND (s.y GT (res*iy + minrange)) AND (s.y LT (res*(iy+1) + minrange)) AND s.tform GT tclip)
            IF inds[0] NE -1 THEN tformrebin[ix,iy] = TOTAL(s[inds].mass*massunit)/res/res/deltat $
            ELSE tformrebin[ix,iy] = 0
            indg = where((g.x gt (res*ix + minrange)) AND (g.x le (res*(ix+1) + minrange)) AND (g.y gt (res*iy + minrange)) AND (g.y lt (res*(iy+1) + minrange)))
            IF indg[0] NE -1 THEN HIrebin[ix,iy]  = TOTAL(g[indg].mass* HI[indg])                /res/res/1000.0/1000.0 ELSE HIrebin[ix,iy]  = 0
            IF indg[0] NE -1 THEN H2rebin[ix,iy]  = TOTAL(g[indg].mass*            2.0*H2[indg]) /res/res/1000.0/1000.0 ELSE H2rebin[ix,iy]  = 0
            IF indg[0] NE -1 THEN gasrebin[ix,iy] = TOTAL(g[indg].mass*(HI[indg] + 2.0*H2[indg]))/res/res/1000.0/1000.0 ELSE gasrebin[ix,iy] = 0

            IF tformrebin[ix,iy] NE 0 THEN BEGIN
                print, H2rebin[ix,iy], H2rebin[ix,iy] + HIrebin[ix,iy], tformrebin[ix,iy]
;                plot,alog10(g[indg].dens*rhounit),alog10(H2[indg]*2/(HI[indg] + 2.0*H2[indg])),xrange = [0,3],psym = 3,yrange = [-6,0]
;                plot_colorscale,alog10(g[indg].dens*rhounit),alog10(H2[indg]*2/(HI[indg] + 2.0*H2[indg])),coolonclean[indg],xrange = [0,3],/overplot,min = -0.33                
;                stop
;                plot,alog10(g[indg].dens*rhounit),alog10(g[indg].tempg),psym = 3
;                plot_colorscale,alog10(g[indg].dens*rhounit),alog10(g[indg].tempg),coolonclean[indg],min = 1e-6,/overplot
                
;Set the gas particles that just had their cooling turned off to the
;H2 fraction of the particles that just formed stars when they were
;forming the stars
                match,iords[inds],starlog.iorderstar,temp,inds_sl
                meanH2frac = mean(starlog[inds_sl].H2form)
                
                indcoolon = where(coolonclean[indg] GE 0) 
                H2cool = H2
                HIcool = HI
                IF (indcoolon)[0] NE -1 THEN BEGIN
                   H2cool[indg[indcoolon]] = meanH2frac/2*MAX(HI)
                   HIcool[indg[indcoolon]] = MAX(HI) - meanH2frac/2*MAX(HI)
                ENDIF
                FOR is = 0, n_elements(inds_sl) - 1 DO BEGIN
                   iordg_sl = where(iordg EQ starlog[inds_sl[is]].iordergas)
                   IF iordg_sl NE -1 THEN H2cool[iordg_sl] = starlog[inds_sl[is]].H2form
                ENDFOR
                H2rebincool[ix,iy]  = TOTAL(g[indg].mass*2.0*H2cool[indg]) /res/res/1000.0/1000.0
;                stop
            ENDIF
        ENDFOR
    ENDFOR
ENDELSE

ind = where(tformrebin GT 0)
xsigmalow = findgen(300)/100 - 1.
xsigma = 10.0^(findgen(400)/100 + 1.) ;xsigmalow
ysigma=2.5e-4*xsigma^1.4
ysigma1=2.5e-4*xsigma^(1.4 + 0.15)
ysigma2=2.5e-4*xsigma^(1.4 - 0.15)
ysigmalow = xsigmalow*2.4 - 5.0

IF ind[0] NE -1 THEN BEGIN
   IF keyword_set(overplot) THEN $
      oplot,alog10(sin(i)*HIrebin[ind] + H2rebin[ind]),alog10(sin(i)*tformrebin[ind]), psym = psym, color = color ELSE $
         plot, alog10(sin(i)*HIrebin[ind] + H2rebin[ind]),alog10(sin(i)*tformrebin[ind]), psym = psym, color = color, xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1, _EXTRA=_extra, pos=[.15,.15,.95,.95], xrange = [-1,5], yrange = [-5,3]
ENDIF
;oplot,alog10(xsigma),alog10(ysigma)
;oplot,xsigmalow,ysigmalow 

IF ind[0] NE -1 THEN BEGIN
   IF keyword_set(overplot) THEN $
      oplot,alog10(sin(i)*H2rebin[ind]),alog10(sin(i)*tformrebin[ind]), psym = psym, color = color ELSE $
         plot, alog10(sin(i)*H2rebin[ind]),alog10(sin(i)*tformrebin[ind]), psym = psym, color = color, xtitle=textoidl('Log \Sigma')+"!H2!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1, _EXTRA=_extra, pos=[.15,.15,.95,.95], xrange = [-5,3], yrange = [-5,-0.5]
   oplot,alog10(sin(i)*H2rebincool[ind]),alog10(sin(i)*tformrebin[ind]), psym = psym, color = 50
   fit = linfit(alog10(sin(i)*H2rebincool[ind]*1e6),alog10(sin(i)*tformrebin[ind]))
   oplot,alog10(sin(i)*H2rebincool[ind]),fit[0] + fit[1]*alog10(sin(i)*H2rebincool[ind]*1e6)
ENDIF
histogramp,alog10(H2rebincool[ind]*1e6/tformrebin[ind]),nbins = 20

stop
END


PRO schmidtlaw_ave,filename,kpcunit,massunit,massform,colors,psyms, outplot = outplot,overplot = overplot,guo = guo
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
;massunit = 1.84793e16
;kpcunit = 50000.
massform = massform*massunit
;htime = 1.0
diskheight = 2.5/2.0
rtipsy,filename + '.halo.1.std',h,g,d,s
rdfloat, filename + '.halo.1.HI', hI, skipline=1 
hI = hI[0:h.ngas-1]
z = (h.time)-1.
tcurrent=13.7e9-wmap3_lookback(-1.0*z) ;tcurrent=13.73d9
IF (keyword_set(deltat) EQ 0) THEN deltat=100.e6 ; in yrs
tclip=tcurrent-deltat      ; only use starformation events since tclip
timeunit = 3.872d3*SQRT((kpcunit*3.0856805d21)^3/(massunit*1.98892d33))/31556926.0
tform=s.tform*timeunit
tcurrent = max(tform)
tclip = tcurrent - deltat
rform=SQRT(s.x^2+s.y^2)*kpcunit
rgas=SQRT(g.x^2+g.y^2)*kpcunit

radius = MAX(rform[where(tform GT tclip AND ABS(rform) LT 50 AND ABS(s.z)*kpcunit LT diskheight)])
radius = 5
indform=WHERE(rform LT radius AND tform GT tclip AND ABS(s.z)*kpcunit LT diskheight ,nindform)
indgas=WHERE(rgas LT radius AND ABS(g.z)*kpcunit LT diskheight ,nindgas)
g_surfacerho = TOTAL(g[indgas].mass*hI[indgas]*massunit)/(!PI*radius^2*1.e6)
s_surfacerho = TOTAL(massform*nindform)/(!PI*radius^2*deltat)

IF keyword_set(overplot) THEN oplot,alog10([g_surfacerho,g_surfacerho]),alog10([s_surfacerho,s_surfacerho]),psym = psyms,color = colors
;MT = 'X,X,X,X,X,F,X,F,F,F,X,X,X,X,X,X,X,X,X,X,X'
;eadcol,filename + '.amiga.stat',F=FMT,vir,gas,star,dark
;uo = [(TOTAL(d.mass))*massunit,TOTAL(s.mass)*massunit]

END

;schmidtlaw_res_master,datafile = 'schmidtlaw_files_z0.txt'
PRO schmidtlaw_res_master,outplot = outplot,datafile = datafile,colors = colors,res = res
loadct,39
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN base = '/home/christensen/Code/' ELSE base = '/astro/users/christensen/code/'
readcol,base + 'Datafiles/HIcubes/'+datafile,dir,file,massunit,kpcunit,massform,inclin,color,psym,key,format = '(A,A,F,F,F,F,F,D,A)'
;IF keyword_set(outplot) THEN BEGIN
;    squareplot,filename=outbase+'schmidt.ps'
;   set_plot, 'ps'
;    device, filename=outbase+'schmidt.ps',/COLOR,bits_per_pixel= 8,/times
;ENDIF ELSE set_plot,'x'
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 12
    ysize = 12
    formatplot,/outplot
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 400
    ysize = 400
    formatplot
 ENDELSE
n = n_elements(file)
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF colors[0] EQ 1 THEN  colors = (findgen(n) + 1)*240/n ELSE colors = colors
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(psyms) THEN psyms = fltarr(n) + 4 ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n)
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(psyms) THEN psyms = findgen(n) + 4;(findgen(n) + 2)*2   
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
ENDELSE

IF keyword_set(outplot) THEN BEGIN
   device,filename=outplot+'_resks_tipsy.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2,/encapsul
ENDIF ELSE BEGIN
   window,0,xsize = xsize,ysize = ysize ;392
ENDELSE

guo_mass = fltarr(2,N_ELEMENTS(dir))
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    IF (i EQ 0) THEN $
       schmidtlaw_res,file[i],kpcunit[i],massunit[i],massform[i], psym = psym[i],i = inclin[i],res = res
    schmidtlaw_res,   file[i],kpcunit[i],massunit[i],massform[i],color = color[i],psym = psym[i],i = inclin[i],/overplot,res = res
;    schmidtlaw_ave,  file[i],kpcunit[i],massunit[i],massform[i],color[i],psym[i] + 3, /overplot;,guo = guo
;    FMT = 'X,X,X,X,X,F,X,F,F,F,X,X,X,X,X,X,X,X,X,X,A'
;    readcol,file[i] + '.amiga.stat',F=FMT,vir,gas,star,dark,contam,/silent                                                          
;    guo_mass[*,i] = [dark[0],star[0]]
ENDFOR

xsigmalow = findgen(200)/100 - 1.
xsigma = findgen(600)/100 - 1.
xsigmahigh = findgen(400)/100 + 1.
;xsigma = 10.0^(findgen(400)/100 + 1.) ;xsigmalow
;ysigma=2.5e-4*xsigma^1.4
;ysigma1=2.5e-4*xsigma^(1.4 + 0.15)
;ysigma2=2.5e-4*xsigma^(1.4 - 0.15)
ysigma = xsigma
ysigmahigh = xsigmahigh - 3.35
ysigmalow = xsigmalow*2.6 - 5.0
;oplot,alog10(xsigma),alog10(ysigma)
oplot,xsigmahigh,ysigmahigh
oplot,xsigmalow,ysigmalow 
oplot,xsigma,ysigma - 3.9,linestyle = 1
oplot,xsigma,ysigma - 2.9,linestyle = 1
oplot,xsigma,ysigma - 1.9,linestyle = 1
legend,key,psym = psym,color = color
IF keyword_set(outplot) THEN device,/close ELSE stop

IF keyword_set(outplot) THEN device, filename=outplot+'guo.ps',/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch ELSE window,1

x = findgen(700)/100 + 9
y = alog10(0.1257*((10^x/10^(11.36))^(-0.9147) + (10^x/10^(11.36))^(0.2485))^(-2.574)*10^x);guo 09, (3)
plot,x,y,xstyle = 1, ystyle = 1,xrange = [9.5,16],yrange = [6,12],xtitle = textoidl('log(M_{halo}[M')+sunsymbol()+'])',ytitle = textoidl('log(M_{*}[M')+sunsymbol()+'])'
FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
   loadct,ctables[i]
    oplot,[alog10(guo_mass[0,i]),alog10(guo_mass[0,i])],[alog10(guo_mass[1,i]),alog10(guo_mass[1,i])],psym = psym[i]+3,color = color[i]
ENDFOR
legend,key,psym = psym+3,color = color,/bottom,/right
IF keyword_set(outplot) THEN device,/close


END
