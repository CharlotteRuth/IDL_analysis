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
    spliti = strsplit(dirs[i],'/')
    filename = strmid(dirs[i],spliti[n_elements(spliti) -1],strlen(dirs[i])-spliti[n_elements(spliti) -1] - 1)
    stop
    spawn,'ls ' + dirs[i] + filename + '.00512/*512.coolontime',file_coolon
;    spawn,'ls ' + dirs[i] + '*512/*512.halo.1.coolontime',file_coolonhalo
    spawn,'ls ' + dirs[i] + filename + '.00512/*512.iord',file_iord
    spawn,'ls ' + dirs[i] + filename + '.00512/*cosmo**512',file
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

;Reading in files
;    outflowz     = mrdfits(dirs[i] + 'grp' + halo[i] + '.rvir.outflow_z.fits',0)
;    outflowiord  = mrdfits(dirs[i] + 'grp' + halo[i]  + '.rvir.outflow_iord.fits',0)
;    outflowm     = mrdfits(dirs[i] + 'grp' + halo[i] + '.rvir.mass_at_outflow.fits',0)
    ejectmass    = mrdfits(dirs[i] + 'grp' + halo[i] + '.mass_at_eject.fits',0)
    ejectz       = mrdfits(dirs[i] + 'grp' + halo[i] + '.eject_z.fits',0)
;    ejectiord    = mrdfits(dirs[i] + 'grp' + halo[i] + '.eject_iord.fits',0)
    expellmass   = mrdfits(dirs[i] + 'grp' + halo[i] + '.mass_at_expell.fits',0)
    expellz      = mrdfits(dirs[i] + 'grp' + halo[i] + '.expell_z.fits',0)

;    ind99 = where(outflowz EQ 99, comp=n99)
;    outflowiord2 = outflowiord[n99]
;Use the coolon file from the last output to select the particles that
;have been heated by SN.
;Should be unnecessary
;    readarr,file_coolon,h,coolon,part = 'gas',/ascii
;    readarr,file_iord,  h,iord,  part = 'gas',/ascii,type = 'long'
;    match,iord,outflowiord2,ind_iord,ind_outflowiord
;    coolonoutind = where(coolon[ind_iord] NE 0)
;    outflowz = outflowz[n99[coolonoutind]]
;    outflowlb = wmap3_lookback(outflowz)
;    outflowm = outflowm[n99[coolonoutind]]
;    stop

;    rtipsy,file_halo,h_halo,g_halo,d_halo,s_halo
;   readarr,file_coolonhalo,h_halo,coolon_halo,part = 'gas',/ascii

;    alliord = mrdfits(dirs[i] + 'grp' + halo[i] + '.allgas.iord.fits',0)
;    rtipsy,file,h,g,d,s
;    match,iord,alliord,ind_iord,ind_alliord
;    gall = g[ind_iord]
;    coolonall = coolon[ind_iord]
;    mcoolon = total(g[where(coolonall NE 0)].mass)*units.massunit

    ind99 = where(expellz EQ 99, comp=n99)
    expellmass    = expellmass[n99]
    expellz       = expellz[n99]
    expellmass1   = ejectmass[n99] ;When the particle that will be lost from the halo first left this disk
    expellz1      = ejectz[n99]
    ind99 = where(ejectz  EQ 99, comp=n99)
    ejectmass     = ejectmass[n99]
    ejectz        = ejectz[n99]

    uniqz = expellz[uniq(expellz,sort(expellz))]
    uniqt = (ageUniverse - wmap3_lookback(uniqz))/1e9
    
    ejectt  = z_to_t(ejectz)
    expellt = z_to_t(expellz)
    expellt1= z_to_t(expellz1)
;    ejectlb  = wmap3_lookback(ejectz)
;    expelllb = wmap3_lookback(expellz)

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
;    histogramp,timeend - outflowlb/1e9,weight = outflowm, binsize = 0.2, min = 1e-6, max = timeend,/cum,yrange = [0,1e10]
;    histogramp,ejectt, weight =  ejectmass, binsize = 0.2, min = 1e-6, max = timeend,/cum,/overplot,color = 240
;    histogramp,expellt,weight =  expellmass,binsize = 0.2, min = 1e-6, max = timeend,/cum,/overplot,color = 60

;    outflowcum    = weighted_histogram(timeend - outflowlb/1e9,weight = outflowm, binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    ejectmasscum    = weighted_histogram(ejectt,   weight =  ejectmass,  binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    expellmasscum   = weighted_histogram(expellt,  weight =  expellmass, binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    expellmasscum1  = weighted_histogram(expellt1, weight =  expellmass1,binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    reejectmasscum  = weighted_histogram(reejectt, weight =  reejectm,   binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    reexpellmasscum = weighted_histogram(reexpellt,weight =  reexpellm,  binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    reexpellmasscum1= weighted_histogram(reexpellt1,weight =  reexpellm1,  binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)

    binsize_large = 1;0.5
;    outflowhist    = weighted_histogram(timeend - outflowlb/1e9,weight = outflowm, binsize = 1, min = 1e-6, max = timeend, locations = timearrhist)
    ejectmasshist    = weighted_histogram(ejectt,   weight =  ejectmass,  binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)
    expellmasshist   = weighted_histogram(expellt,  weight =  expellmass, binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)
    expellmasshist1  = weighted_histogram(expellt1, weight =  expellmass1,binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)
    reejectmasshist  = weighted_histogram(reejectt, weight =  reejectm,    binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)
    reexpellmasshist = weighted_histogram(reexpellt,weight =  reexpellm,   binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)
    reexpellmasshist1 = weighted_histogram(reexpellt1,weight =  reexpellm1,   binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)

    IF keyword_set(massloading) THEN BEGIN
       sfrhist = weighted_histogram(sltimeform*units.timeunit/1e9,weight = slmassform,binsize = binsize_large, min = 1e-6, max = timeend, locations = timearrhist)
    ENDIF; ELSE sfrhist = fltarr(n_elements(ejectmasshist)) + 1.e8
;    stop
    IF keyword_set(massloading) THEN BEGIN
          IF keyword_set(redshift) THEN BEGIN
             yrange = [0,5]
              IF i EQ 0 THEN BEGIN
                 plot, timearrhist, expellmasshist/sfrhist, xtitle='Time [Gyr]', ytitle=textoidl('Outflow / SFR'), title=label, psym=10, xrange=[0,14], xstyle = 9, ystyle = 1,/nodata,yrange = yrange;[0,17] ;h516 [0,30]; h986 [0,6] ;h516.H2 
                 axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
              ENDIF
          ENDIF ELSE BEGIN
             yrange = [0,20];[0,14]
             IF i EQ 0 THEN $
                plot, timearrhist, reejectmasshist/  sfrhist, xtitle='Time [Gyr]', ytitle=textoidl('Outflow / SFR'), title=label, psym=10, xrange=[0,14], /xstyle, /nodata,yrange = yrange;,yrange = yrange ;h516 [0,30]; h986 [0,6] ;h516.H2 
          ENDELSE
;          oplot,   timearrhist, expellmasshist/  sfrhist, thick=thicks[i], color=colors[i],psym = 10, linestyle = 2 ;linestyles[i]
;          oplot,   timearrhist, reexpellmasshist/ sfrhist, thick=thicks[i], color=colors[i],psym = 10, linestyle = linestyles[i] + 2
;          oplot,   timearrhist, expellmasshist1/ sfrhist, thick=thicks[i], color=colors[i] + 100,psym = 10, linestyle = linestyles[i]
;          oplot,   timearrhist, reexpellmasshist1/sfrhist, thick=thicks[i], color=colors[i],psym = 10, linestyle = linestyles[i] + 2
;          oplot,   timearrhist, ejectmasshist/    sfrhist, thick=thicks[i], color=colors[i],psym = 10, linestyle = linestyles[i]
          oplot,   timearrhist, reejectmasshist/  sfrhist, thick=thicks[i], color=colors[i],psym = 10, linestyle = 0;linestyles[i]; + 1
          oplot,[0,14],[total(reejectm)/total(slmassform),total(reejectm)/total(slmassform)],thick=thicks[i], color=colors[i],psym = 10, linestyle = 2
;          stop
    ENDIF ELSE BEGIN ;No Mass Loading
       IF keyword_set(ratio) THEN BEGIN
          IF keyword_set(redshift) THEN BEGIN
              IF i EQ 0 THEN BEGIN
                 plot, timearr, expellmasscum/ejectmasscum, xtitle='Time [Gyr]',ytitle=textoidl('Cumulative Expelled/Ejected'), title=label, psym=10, xrange=[0,14], xstyle = 9, ystyle = 1,/nodata, yrange = [0,1] ;h799[0,0.4], ;h516 [0,1] ;h986 [0,10] ;h516.H2
                 axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
              ENDIF
          ENDIF ELSE BEGIN
                plot, timearr, expellmasscum/ejectmasscum, xtitle='Time [Gyr]',ytitle=textoidl('Cumulative Expelled/Ejected'), title=label, psym=10, xrange=[0,14], /xstyle, ystyle = 1,/nodata, yrange = [0,1] ;h799[0,0.4], ;h516 [0,1] ;h986 [0,10] ;h516.H2
          ENDELSE
;          oplot,   timearr, expellmasscum/ejectmasscum,     thick=thicks[i], color=colors[i], linestyle = 0
;          oplot,   timearr, reexpellmasscum/reejectmasscum, thick=thicks[i], color=colors[i], linestyle = 1
 ;         oplot,   timearr, expellmasscum1/ejectmasscum,    thick=thicks[i], color=colors[i] + 100, linestyle = 0
          oplot,   timearr, reexpellmasscum1/reejectmasscum,thick=thicks[i], color=colors[i];, linestyle = 1
       ENDIF ELSE BEGIN
       IF keyword_set(unscale) THEN BEGIN
          yrange = [0,50]
          IF keyword_set(redshift) THEN BEGIN
              IF i EQ 0 THEN BEGIN
                 plot, timearr, ejectmasscum/ 1.e9, xtitle='Time [Gyr]',ytitle=textoidl('Cumulative Outflow / 10^{9} [M')+sunsymbol()+']', title=label, psym=10, xrange=[0,14], xstyle = 9, ystyle = 1,/nodata, yrange = yrange ;h799[0,0.4], ;h516 [0,1] ;h986 [0,10] ;h516.H2
                 axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
              ENDIF
          ENDIF ELSE BEGIN
                IF i EQ 0 THEN plot, timearr, reejectmasscum/ 1.e9, xtitle='Time [Gyr]',ytitle=textoidl('Cumulative Outflow / 10^{9} [M')+sunsymbol()+']', title=label, psym=10, xrange=[0,14], /xstyle, ystyle = 1,/nodata, yrange = [0,0.55] ;h799[0,0.55], ;h516 [0,1] ;h986 [0,10] ;h516.H2
          ENDELSE
;          oplot,   timearr, expellmasscum/1.e9,   thick=thicks[i], color=colors[i], linestyle = 2
          oplot,   timearr, reexpellmasscum/1.e9, thick=thicks[i], color=colors[i], linestyle = linestyles[i] + 2
;          oplot,   timearr, expellmasscum1/1.e9,  thick=thicks[i], color=colors[i] + 100, linestyle = linestyles[i]
;          oplot,   timearr, reexpellmasscum1/1.e9,thick=thicks[i], color=colors[i], linestyle = linestyles[i] + 3
;          oplot,   timearr, ejectmasscum/ 1.e9,   thick=thicks[i], color=colors[i], linestyle = linestyles[i] ;linestyle = 0
          oplot,   timearr, reejectmasscum/1.e9,  thick=thicks[i], color=colors[i], linestyle = linestyles[i]; + 1 ;linestyle = 0
          print,max(ejectmasscum/1e9),max(expellmasscum/1e9)

;          IF i EQ 0 THEN $
;            plot,timearr,expellmasscum/ejectmasscum, xtitle='Time [Gyr]', ytitle='Expell/Eject]', title=label,yrange = [0,1]
;          oplot,timearr,expellmasscum/ejectmasscum,thick =thicks[i], color = colors[i]
 
;           IF i EQ 0 THEN  plot, timearr, ejectmasscum/mcoolon, xtitle='Time [Gyr]', ytitle=textoidl('Cumulative Outflow / 10^{9} [M')+sunsymbol()+']', title=label, psym=10, xrange=[0,14], /xstyle, /nodata, yrange = [0,1] ;h516 [0,1] ;h986 [0,10] ;h516.H2
;          oplot,   timearr, ejectmasscum/mcoolon, thick=thicks[i], color=colors[i]
;          oplot,   timearr, expellmasscum/mcoolon, thick=thicks[i], color=colors[i], linestyle = 2          
;          stop
       ENDIF ELSE BEGIN         
;Calculate the total amount of mass accreted
          massaccr  = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdisk.fits',0);mass_at_accr
          zaccr     = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0);accr_z
          massaccr  = massaccr[where(zaccr GE zend)]
          igmass    = median(massaccr)
          early     = mrdfits(dirs[i] + '/early.iord.fits',0)
          massearly = n_elements(early)*igmass
          totalmass = total(massaccr) + massearly

          yrange = [0,0.75]
          IF keyword_set(redshift) THEN BEGIN          
              IF i EQ 0 THEN BEGIN
                 plot, timearr, ejectmasscum/totalmass, xtitle='Time [Gyr]', ytitle=textoidl('Cumulative Outflow / Total Gas Mass'), title=label, psym=10, xrange=[0,14], xstyle = 9, ystyle = 1,/nodata, yrange = yrange ;h986[0,0.35];h603[0,0.2]
                 axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
              ENDIF
          ENDIF ELSE BEGIN
              IF i EQ 0 THEN $
                plot, timearr, ejectmasscum/totalmass, xtitle='Time [Gyr]', ytitle=textoidl('Cumulative Outflow / Total Gas Mass'), title=label, psym=10, xrange=[0,14], /xstyle,/nodata, yrange = yrange ;h986[0,0.35];h603[0,0.2]              
          ENDELSE
;          oplot,   timearr, expellmasscum/totalmass, thick=thicks[i], color=colors[i], linestyle = linestyles[i]
          oplot,   timearr, reexpellmasscum/totalmass,thick=thicks[i], color=colors[i], linestyle = linestyles[i] + 2
;          oplot,   timearr, expellmasscum1/totalmass, thick=thicks[i], color=colors[i] + 100, linestyle = linestyles[i]
;          oplot,   timearr, reexpellmasscum1/totalmass,thick=thicks[i], color=colors[i], linestyle = linestyles[i] + 3
;          oplot,   timearr, ejectmasscum/ totalmass, thick=thicks[i], color=colors[i],   linestyle = linestyles[i] ;linestyle = 0
          oplot,   timearr, reejectmasscum/totalmass,thick=thicks[i], color=colors[i],   linestyle = linestyles[i];linestyle = 0

;         oplot,   timearr, outflowcum/totalmass,thick=thicks[i], color=colors[i], linestyle = 1
          print,file,totalmass,max(ejectmasscum),max(ejectmasscum/totalmass),max(expellmasscum),max(expellmasscum/totalmass),max(reejectmasscum),max(reexpellmasscum)
       ENDELSE
;       IF keyword_set(keys) AND i lt N_ELEMENTS(keys) THEN BEGIN
;          l_keys = lb_keys
;          l_keys[i] = keys[i]
;          l_thicks = lb_thicks
;          l_thicks[i] = thicks[i]
;          l_color = lb_color
;          l_color[i] = colors[i]
;          legend,l_keys,linestyle = fltarr(n_elements(l_keys)),color = l_color,thick = l_thicks,/top,/left,box = 0
;       ENDIF
    ENDELSE
    ENDELSE
;    stop
ENDFOR
IF keyword_set(massloading) OR keyword_set(ratio) THEN BEGIN
    IF keyword_set(ratio) AND keyword_set(keys) THEN legend,keys,linestyle = linestyles,thick = thicks, color = colors,ctables = ctables,/left,/top,box = 0
    IF keyword_set(massloading) OR keyword_set(keys) THEN legend,keys,linestyle = fltarr(n),thick = thicks, color = colors,ctables = ctables,/left,/top,box = 0
ENDIF ELSE BEGIN
    IF keyword_set(keys) THEN legend,keys,linestyle = fltarr(n),thick = thicks, color = colors,ctables = ctables,/left,/top,box = 0
;    legend,keys,linestyle = fltarr(n),    thick = thicks, color = colors,ctables = ctables,/left,/top,box = 0  
    legend,['Ejected','Expelled'],linestyle = [0,2],thick = fltarr(2) + thicks[0],/right,box = 0    
ENDELSE
IF keyword_set(outplot) THEN device, /close ELSE stop
END

PRO eject_plots_master,outplot = outplot
  formatplot,outplot = outplot  ;,/thick
  spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN prefix = '/nobackupp8/crchrist/MolecH/' $
ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

  steps  = ['00492']
  dirtop = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/']

  steps  = ['00512']
  dirtop = [prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/']
  massloading = 0
  ratio = 0
  unscale = 0
  redshift = 1

  eject_plots,dirtop,/colors,massloading = massloading, ratio = ratio,unscale = unscale, redshift = redshift
END
