PRO toomreRadius,files,dKpcUnit,dMsolunit,outplot = outplot,maxdistance = maxdistance,color = color,keys = keys,thicks = thicks,symbols = symbols,linestyle = linestyle
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_amu = 1/amu_per_gm
boltzmann = 1.380658e-16 ;cgs
sec_per_year = 31556926
grav = 6.67e-8 ;cgs

n = N_ELEMENTS(files)
IF KEYWORD_SET(outplot) THEN     fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT KEYWORD_SET(symbols) THEN symbols = (fltarr(n)+4)
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT KEYWORD_SET(thicks) THEN thicks = (findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = REVERSE(findgen(n)*2)   
    IF NOT KEYWORD_SET(symbols) THEN symbols = (findgen(n)+2)*2
ENDELSE
yrange = [0,70]

IF NOT KEYWORD_SET(maxdistance) THEN maxdistance = 4.0
nbins = 100.0
bind = maxdistance/nbins
vc_x = (findgen(nbins)+1)*bind*cm_per_kpc
tmass = fltarr(nbins) 

nbins1 = 20.0
bind1 = maxdistance/nbins1*cm_per_kpc
vc_x1 = (findgen(nbins1)+1)*bind1
toomre = fltarr(nbins1,n)
sd = fltarr(nbins1,n)
sigma = fltarr(nbins1,n)
kapa = fltarr(nbins1,n)

n = N_ELEMENTS(files)
;!p.multi = [0,1,n]
;multiplot,[1,n];,mxtitle = 'Radius [kpc]',mytitle = 'Density [amu/cc]';,/doxaxis,/doyaxis;gap=0.04,
IF KEYWORD_SET(outplot) THEN BEGIN
    bgcolor = 255
    fgcolor = 0 
ENDIF ELSE BEGIN
    bgcolor = 0
    fgcolor = 255
ENDELSE
FOR i = 0, n - 1 DO BEGIN
    rtipsy,files[i] + '.halo.1',h,g,d,s
    a = h.time
    IF (N_ELEMENTS(dKpcUnit) GT 1) THEN kpc_per_syslength = dKpcUnit[i] ELSE kpc_per_syslength = (dKpcUnit)[0]*a
    IF (N_ELEMENTS(dMsolUnit) GT 1) THEN msol_per_sysmass = dMsolUnit[i] ELSE msol_per_sysmass = dMsolUnit[0]
    timeunit = SQRT((cm_per_kpc*kpc_per_syslength)^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year
    dens_convert =  dMsolunit[i] * gm_per_msol * 5.9753790e+23/dKpcUnit[i]^3/cm_per_kpc^3/h.time^3

;----------------- Vcirc -------
    massg = g.mass*msol_per_sysmass
    IF N_ELEMENTS(s) gt 0 THEN masss = s.mass*msol_per_sysmass
    massd = d.mass*msol_per_sysmass
    IF N_ELEMENTS(s) gt 0 THEN mass = [massg,masss,massd] ELSE mass = [massg,massd]

    distancesg = SQRT(g.x^2.0 + g.y^2.0 + g.z^2.0)*kpc_per_syslength*cm_per_kpc
    IF N_ELEMENTS(s) gt 0 THEN distancess = SQRT(s.x^2.0 + s.y^2.0 + s.z^2.0)*kpc_per_syslength*cm_per_kpc
    distancesd = SQRT(d.x^2.0 + d.y^2.0 + d.z^2.0)*kpc_per_syslength*cm_per_kpc
    IF N_ELEMENTS(s) gt 0 THEN distances = [distancesg,distancess,distancesd] ELSE distances = [distancesg,distancesd]    

    FOR j = 0, nbins-1 DO BEGIN
        ind = where(distances LE vc_x[j])
        if (ind[0] ne -1)then tmass[j] = TOTAL(mass[ind])else tmass[j] = 0
    ENDFOR
;    vc_x = vc_x*cm_per_kpc ;now in cm
;    bind = bind*cm_per_kpc
    vc_y = SQRT(tmass*gm_per_msol*grav/vc_x);cm/s
    ;IF i eq 0 THEN 
 ;   plot,alog10(vc_x),alog10(vc_y),linestyle = linestyle[i],xtitle = 'Radius [cm]',ytitle = textoidl('V_{circ} [cm/s]'),yrange = [5,7] ;,xra ge = alog10([0,maxdistance])
 ;   oplot,alog10(vc_x),alog10(vc_y),linestyle = linestyle[i],color = colors[i]
    fitpar = poly_fit(alog10(vc_x),alog10(vc_y),3)
    fit = fitpar[0] + alog10(vc_x)*fitpar[1] + alog10(vc_x)^2*fitpar[2] + alog10(vc_x)^3*fitpar[3]; + alog10(vc_x)^4*fitpar[4] + alog10(vc_x)^5*fitpar[5]
 ;   oplot,alog10(vc_x),fit,linestyle = linestyle[i],color = 50

;----------------- Toomre -------
    disk = where(g.tempg lt 5e4)
    g1 = g[disk]
    g1mass = g1.mass*msol_per_sysmass*gm_per_msol
    radius = distancesg[disk] ;now in cm
    FOR j = 0, nbins1-1 DO BEGIN
        ind = where(radius LE vc_x1[j] AND radius GE vc_x1[j] - bind1)
        radmid = vc_x1[j] + bind1/2
        area = !PI*(vc_x1[j]^2 - (vc_x1[j] - bind1)^2)
        mass = TOTAL(g1mass[ind])
        sd[j,i] = mass/area
        avetemp = TOTAL(g1[ind].tempg*g1mass[ind])/mass
        sigma[j,i] = SQRT(3.0*boltzmann*avetemp/gm_per_amu)
;        stop
        vcirc = fitpar[0] + alog10(radmid)*fitpar[1] + alog10(radmid)^2*fitpar[2] + alog10(radmid)^3*fitpar[3] ;+ alog10(radmid)^4*fitpar[4] + alog10(radmid)^5*fitpar[5]
        beta = fitpar[1] + 2.0*alog10(radmid)*fitpar[2] + 3.0*alog10(radmid)^2*fitpar[3] ;+ 4.0*alog10(radmid)^3*fitpar[4] + alog10(radmid)^5*fitpar[5]
        kapa[j,i] = 1.41*10^vcirc/radmid*SQRT(1 + beta)

        toomre[j,i] = sigma[j,i]*kapa[j,i]/(!PI*grav*sd[j,i])
 ;       stop
    ENDFOR 
ENDFOR
IF KEYWORD_SET(outplot) THEN device,filename=outplot + 'toomre.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 ELSE window,0
FOR i = 0, n - 1 DO BEGIN
    if i eq 0 THEN plot,(vc_x1 - nbins1/2.0)/cm_per_kpc,toomre[*,i],psym = -1.0*symbols[i],ytitle = 'Q',xtitle = 'Radius [kpc]',yrange = [0.5,20],/ylog, xrange = [0,maxdistance],linestyle = linestyle[i],thick = thicks[i]
    oplot,(vc_x1 - nbins1/2.0)/cm_per_kpc,toomre[*,i],psym = -1.0*symbols[i],color = colors[i],linestyle = linestyle[i],thick = thicks[i]
;    multiplot;,/doyaxis 
ENDFOR
IF KEYWORD_SET(outplot) THEN device,/close else stop

;----------- Toomre obs ---------
angle = !pi/2
;nbins = 60
;r25 = opticalRadii(filename,verbose = 0)
;maxr = r25*4
;radii_bins = findgen(nbins1 + 1)*maxdistance/nbins1
;***************** Gas Cubes **********************

;sdr = fltarr(nbins)
;sigmar = fltarr(nbins)
;sigmar_norm = sigmar
;radii_x = findgen(nbins)*maxr/nbins + maxr/nbins/2.0

toomre_obs = fltarr(nbins1,n)
sd_obs = fltarr(nbins1,n)
sigma_obs = fltarr(nbins1,n)
sigma_norm_obs = fltarr(nbins1,n)
sigma_all_obs = fltarr(nbins1,n)
kapa_obs = fltarr(nbins1,n)
r25 = fltarr(n)

FOR i = 0, n - 1 DO BEGIN
    r25[i] = opticalRadii(filename = files[i])
    filebase = files[i] + '.halo.1.90.smoothed'
    cubeHI_surface_den = read_dencube_fits(filebase + '.mom0.fits',headerHI,/noscale)/amu_per_gm;*4787.23257/amu_per_gm
    cubeHI_sigma = read_dencube_fits(filebase + '.mom2.fits',headerHI,/noscale)*1e5
    xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1
    yaxes = (findgen(headerHI.NAXIS2) - headerHI.NAXIS1/2.0)*headerHI.CDELT2
    radius_obs = fltarr(N_ELEMENTS(xaxes),N_ELEMENTS(yaxes))
    FOR ix = 0,N_ELEMENTS(xaxes) - 1 DO $
      FOR iy = 0,N_ELEMENTS(yaxes) - 1 DO $
      radius_obs[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy]/SIN(angle)/SIN(angle))*cm_per_kpc
    for j = 0 , nbins1 - 1 DO BEGIN 
        ind = where(radius_obs LT vc_x1[j] and radius_obs GE vc_x1[j] - bind1)
        radmid = vc_x1[j] + bind1/2
        sigma_obs[j,i] = MEAN(cubeHI_sigma[ind])
        sigma_norm_obs[j,i] = TOTAL(cubeHI_sigma[ind]*cubeHI_surface_den[ind])/TOTAL(cubeHI_surface_den[ind])
        sigma_all_obs[j,i]= sqrt(sigma_norm_obs[j,i]^2 + sigma[j,i]^2)
        sd_obs[j,i] = MEAN(cubeHI_surface_den[ind]) 
        vcirc = fitpar[0] + alog10(radmid)*fitpar[1] + alog10(radmid)^2*fitpar[2] + alog10(radmid)^3*fitpar[3] ;+ alog10(radmid)^4*fitpar[4] + alog10(radmid)^5*fitpar[5]
        beta = fitpar[1] + 2.0*alog10(radmid)*fitpar[2] + 3.0*alog10(radmid)^2*fitpar[3] ;+ 4.0*alog10(radmid)^3*fitpar[4] + alog10(radmid)^5*fitpar[5]
        kapa_obs[j,i] = 1.41*10^vcirc/radmid*SQRT(1 + beta)

        toomre_obs[j,i] = sigma_all_obs[j,i]*kapa[j,i]/(!PI*grav*sd_obs[j,i])        
;    stop
    ENDFOR
ENDFOR
IF KEYWORD_SET(outplot) THEN device,filename=outplot + 'toomre_obs.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 ELSE window,1
FOR i = 0, n - 1 DO BEGIN
    if i eq 0 THEN plot,(vc_x1 - nbins1/2.0)/cm_per_kpc,toomre_obs[*,i],psym = -1*symbols[i],ytitle = 'Q',xtitle = 'Radius [kpc]',yrange = [0.5,20],/ylog, xrange = [0,maxdistance],linestyle = linestyle[i],thick = thicks[i]
    oplot,(vc_x1 - nbins1/2.0)/cm_per_kpc,toomre_obs[*,i],psym = -1*symbols[i],color = colors[i],linestyle = linestyle[i],thick = thicks[i]
    oplot,[r25[i],r25[i]],[1e-4,1e4],linestyle = linestyle[i],color = colors[i],thick = thicks[i]
;    multiplot;,/doyaxis  
ENDFOR
if KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,psym = -1*symbols,/right,/bottom
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop

END

PRO toomreRadius_Master,outplot = outplot
prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
tfiles = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
         'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
dKpcUnit = [25000.,25000.]
dMsolUnit = [ 2.310e15, 2.310e15]
formatplot,outplot = outplot

toomreRadius,prefix + tfiles,dKpcUnit,dMsolunit,xrange = [0,8],outplot = outplot

END
