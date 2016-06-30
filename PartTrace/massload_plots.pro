;Charlotte Christensen
;5/31/12
;This program, very much like outflow_plots.pro, plots the outflowing mass

PRO eject_plots, dirs, halo = halo, unscale = unscale, outplot = outplot, keys = keys, colors = colors, thicks = thicks, label = label, ctables = ctables, yrange_fbcum = yrange_fbcum, massloading = massloading, molecularH = molecularH, formatthick = formatthick, redshift = redshift, linestyles = linestyles,ratio = ratio

formatplot,outplot = outplot,thick = formatthick
IF keyword_set(redshift) THEN $
   IF keyword_set(outplot) THEN !Y.MARGIN = [8,4] ELSE !Y.MARGIN = [6,4]

n = n_elements(dirs)
IF NOT keyword_set(halo) THEN halo = strarr(n) + '1'

IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 800
    ysize = 500
ENDELSE
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2)  
ENDELSE
IF NOT keyword_set (yrange_fbcum) THEN yrange_fbcum = [0,250]

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
maxtime = wmap3_lookback(1000)

ageUniverse = 13.7346*1e9       ;wmap3_lookback(100)
IF keyword_set(redshift) THEN BEGIN
    z = reverse(findgen(100)*10.0/100.0)
    t = ageUniverse - wmap3_lookback(z)
    tickred_in_t = reverse(findgen(5))
    ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e9
ENDIF

IF keyword_set(massloading) THEN outplot_ext = '_massloading2.eps' ELSE BEGIN
    IF keyword_set(unscale) THEN outplot_ext = '_fbcum.eps' ELSE outplot_ext = '_fbcum_scale.eps'
    IF keyword_set(ratio)   THEN outplot_ext = '_fbratio.eps'
ENDELSE

IF keyword_set(outplot) THEN  device,filename = outplot + outplot_ext,/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
;stop
FOR i = 0, N_ELEMENTS(dirs) - 1 DO BEGIN
    loadct, ctables[i]
    spawn,'ls ' + dirs[i] + '*512/*512.coolontime',file_coolon
;    spawn,'ls ' + dirs[i] + '*512/*512.halo.1.coolontime',file_coolonhalo
    spawn,'ls ' + dirs[i] + '*512/*512.iord',file_iord
    spawn,'ls ' + dirs[i] + '*512/*cosmo**512',file
;    spawn,'ls ' + dirs[i] + '*512/*cosmo**512.halo.1',file_halo
    spawn,'ls ' + dirs[i] + 'h*param',pfile
    rtipsy,file,h,g,d,s,/justhead
    units = tipsyunits(pfile[0])
    IF keyword_set(massloading) THEN BEGIN
;       spawn,'ls ' + dirs[i] + '*.starlog',file_starlog 
;       IF strlen(file_starlog) GT 1 THEN BEGIN
;           starlog = rstarlog(file_starlog[0],molecularH = molecularH)
;           single = rem_dup(starlog.iorderstar) ;remove any duplicate stars
;           starlog = starlog[single]
;           slmassform = starlog.massform
;           sltimeform = starlog.timeform
;       ENDIF ELSE BEGIN
;           rtipsy,file + '.halo.1.std',h,g,d,s
           starlog = mrdfits(dirs[i] + 'starlog.' + halo[i] + '.fits',1)
           sltimeform = starlog.timeform;*units.timeunit
           slmassform = starlog.massform*units.massunit
;           sltimeform = s.tform
;       ENDELSE
;           stop
    ENDIF

    zend = 1.0/h.time - 1.0
    timeend = 13.734377

    reejecti = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_iord.fits',0)
    reejectm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reeject.fits',0)
    reejectz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_z.fits',0)
    reejectz[where(reejectz LT 0)] = 0
    reejectt = z_to_t(reejectz)
    reejecti_uniq = reejecti[uniq(reejecti,sort(reejecti))]
    reejectmcum  = weighted_histogram(reejectt, weight =  reejectm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)
    
    reexpelli = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_iord.fits',0)
    reexpellm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reexpell.fits',0) 
    reexpellz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_z.fits',0)
    reexpellz[where(reexpellz LT 0)] = 0
    reexpellt = z_to_t(reexpellz)
    reexpelli_uniq = reexpelli[uniq(reexpelli,sort(reexpelli))]
    reexpellmcum  = weighted_histogram(reexpellt, weight =  reexpellm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)
    reexpellm1 = fltarr(n_elements(reexpellm))
    reexpellz1 = fltarr(n_elements(reexpellz))
    exi = 0L
    FOR eji = 0, n_elements(reejecti) -1 DO BEGIN
       IF exi GE n_elements(reexpelli) THEN BREAK
       IF reejecti[eji] EQ reexpelli[exi] THEN BEGIN
          reexpellm1[exi] = reejectm[eji]
          reexpellz1[exi] = reejectz[eji]
          exi = exi + 1
;          IF exi MOD 100 EQ 0 THEN print,exi
       ENDIF
    ENDFOR
    reexpellt1 = z_to_t(reexpellz1)

;------------------------- Bin the Data -----------------------    
    reejectmasscum  = weighted_histogram(reejectt, weight =  reejectm,   binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    reexpellmasscum = weighted_histogram(reexpellt,weight =  reexpellm,  binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    reexpellmasscum1= weighted_histogram(reexpellt1,weight =  reexpellm1,  binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)

    binsize_large = 1;0.5
    reejectmasshist  = weighted_histogram(reejectt, weight =  reejectm,    binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)
    reexpellmasshist = weighted_histogram(reexpellt,weight =  reexpellm,   binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)
    reexpellmasshist1 = weighted_histogram(reexpellt1,weight =  reexpellm1,   binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)

    IF keyword_set(massloading) THEN BEGIN
       sfrhist = weighted_histogram(sltimeform*units.timeunit/1e9,weight = slmassform,binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)
    ENDIF; ELSE sfrhist = fltarr(n_elements(ejectmasshist)) + 1.e8
;    stop

    yrange = [0,20]             ;[0,14]
    IF i EQ 0 THEN $
       plot, timearrhist, reejectmasshist/  sfrhist, xtitle='Time [Gyr]', ytitle=textoidl('Outflow / SFR'), title=label, psym=10, xrange=[0,14], /xstyle, /nodata,yrange = yrange ;,yrange = yrange ;h516 [0,30]; h986 [0,6] ;h516.H2 
    oplot,   timearrhist, reejectmasshist/  sfrhist, thick=thicks[i], color=colors[i],psym = 10, linestyle = 0 ;linestyles[i]; + 1
    oplot,[0,14],[total(reejectm)/total(slmassform),total(reejectm)/total(slmassform)],thick=thicks[i], color=colors[i],psym = 10, linestyle = 2
    stop
 ENDFOR
IF keyword_set(massloading) THEN legend,keys,linestyle = fltarr(n),thick = thicks, color = colors,ctables = ctables,/left,/top,box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop
END

