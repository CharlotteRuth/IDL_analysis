;Charlotte Christensen
;5/31/12
;This program, very much like outflow_plots.pro, plots the outflowing mass

PRO eject_plots, dirs, unscale = unscale, outplot = outplot, keys = keys, colors = colors,thicks = thicks,label = label,ctables = ctables,yrange_fbcum = yrange_fbcum,massloading = massloading, molecularH = molecularH,formatthick = formatthick
formatplot,outplot = outplot,thick = formatthick
n = N_ELEMENTS(dirs)

IF KEYWORD_SET(outplot) THEN BEGIN
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
IF KEYWORD_SET(colors) THEN BEGIN
    loadct,39
    IF NOT keyword_set(ctables) THEN ctables = [39,39,39]
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN ctables = [0,0,0]
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT KEYWORD_SET(thicks) THEN thicks = (findgen(n) + 1)*6/n - 1
ENDELSE
IF NOT KEYWORD_SET (yrange_fbcum) THEN yrange_fbcum = [0,250]

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
maxtime = wmap3_lookback(1000)

IF KEYWORD_SET(outplot) THEN  device,filename = outplot + '_fbcum.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize

FOR i = 0, N_ELEMENTS(dirs) - 1 DO BEGIN
    loadct, ctables[i]
    spawn,'ls ' + dirs[i] + '*512/*512.coolontime',file_coolon
    spawn,'ls ' + dirs[i] + '*512/*512.iord',file_iord
    spawn,'ls ' + dirs[i] + '*512/*cosmo**512',file
    spawn,'ls ' + dirs[i] + 'h*param',pfile
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
           starlog = mrdfits('starlog.1.fits',1)
           sltimeform = starlog.timeform
           slmassform = starlog.massform*units.massunit
;           sltimeform = s.tform
;       ENDELSE
    ENDIF
    rtipsy,file,h,g,d,s,/justhead
    readarr,file_coolon,h,coolon,part = 'gas',/ascii
    readarr,file_iord,  h,iord,  part = 'gas',/ascii,type = 'long'
    zend = 1.0/h.time - 1.0
    timeend = 13.734377

;Reading in files
;    outflowz     = mrdfits(dirs[i] + 'grp1.rvir.outflow_z.fits',0)
;    outflowiord  = mrdfits(dirs[i] + 'grp1.rvir.outflow_iord.fits',0)
;    outflowm     = mrdfits(dirs[i] + 'grp1.rvir.mass_at_outflow.fits',0)
    ejectmass    = mrdfits(dirs[i] + 'grp1.mass_at_eject.fits',0)
    ejectz       = mrdfits(dirs[i] + 'grp1.eject_z.fits',0)
    expellmass   = mrdfits(dirs[i] + 'grp1.mass_at_expell.fits',0)
    expellz      = mrdfits(dirs[i] + 'grp1.expell_z.fits',0)

;    ind99 = where(outflowz EQ 99, comp=n99)
;    outflowiord2 = outflowiord[n99]
;Use the coolon file from the last output to select the particles that
;have been heated by SN.
;Should be unnecessary
;    match,iord,outflowiord2,ind_iord,ind_outflowiord
;    coolonoutind = where(coolon[ind_iord] NE 0)
;    outflowz = outflowz[n99[coolonoutind]]
;    outflowlb = wmap3_lookback(outflowz)
;    outflowm = outflowm[n99[coolonoutind]]
;    stop

    ind99 = where(ejectz  EQ 99, comp=n99)
    ejectmass    = ejectmass[n99]
    ejectz       = ejectz[n99]
    ind99 = where(expellz EQ 99, comp=n99)
    expellmass   = expellmass[n99]
    expellz      = expellz[n99]

    ejectlb  = wmap3_lookback(ejectz)
    expelllb = wmap3_lookback(expellz)
    
;    histogramp,timeend - outflowlb/1e9,weight = outflowm, binsize = 0.2, min = 1e-6, max = timeend,/cum,yrange = [0,1e10]
;    histogramp,timeend - ejectlb/1e9, weight =  ejectmass, binsize = 0.2, min = 1e-6, max = timeend,/cum,/overplot,color = 240
;    histogramp,timeend - expelllb/1e9,weight =  expellmass,binsize = 0.2, min = 1e-6, max = timeend,/cum,/overplot,color = 60

;    outflowcum    = weighted_histogram(timeend - outflowlb/1e9,weight = outflowm, binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    ejectmasscum  = weighted_histogram(timeend - ejectlb/1e9, weight =  ejectmass, binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)
    expellmasscum = weighted_histogram(timeend - expelllb/1e9,weight =  expellmass,binsize = 0.2, min = 1e-6, max = timeend, locations = timearr,/cum)

;    outflowhist    = weighted_histogram(timeend - outflowlb/1e9,weight = outflowm, binsize = 1, min = 1e-6, max = timeend, locations = timearrhist)
    ejectmasshist  = weighted_histogram(timeend - ejectlb/1e9, weight =  ejectmass, binsize = 1, min = 1e-6, max = timeend, locations = timearrhist)
    expellmasshist = weighted_histogram(timeend - expelllb/1e9,weight =  expellmass,binsize = 1, min = 1e-6, max = timeend, locations = timearrhist)

    IF keyword_set(massloading) THEN BEGIN
       sfrhist = weighted_histogram(sltimeform,weight = slmassform,binsize = 1, min = 1e-6, max = timeend, locations = timearrhist)
    ENDIF; ELSE sfrhist = fltarr(n_elements(ejectmasshist)) + 1.e8

;    stop

    IF keyword_set(massloading) THEN BEGIN
          IF i EQ 0 THEN $
             plot, timearrhist, ejectmasshist/sfrhist, xtitle='Time [Gyr]', ytitle=textoidl('Outflow / SFR'), title=label, psym=10, xrange=[0,14], /xstyle, /nodata,yrange = [0,6] ;h516.H2
          oplot,   timearrhist, ejectmasshist/sfrhist, thick=thicks[i], color=colors[i],psym = 10
          oplot,   timearrhist, expellmasshist/sfrhist, thick=thicks[i], color=colors[i], linestyle = 2,psym = 10
;          oplot,   timearrhist, outflowhist/sfrhist, thick=thicks[i], color=colors[i], linestyle = 1,psym = 10
    ENDIF ELSE BEGIN
       IF keyword_set(unscale) THEN BEGIN
          IF i EQ 0 THEN $
             plot, timearr, ejectmasscum/ 1.e9, xtitle='Time [Gyr]', ytitle=textoidl('Cumulative Outflow / 10^{9} [M')+sunsymbol()+']', title=label, psym=10, xrange=[0,14], /xstyle, /nodata, yrange = [0,10] ;h516.H2
          oplot,   timearr, ejectmasscum/ 1.e9, thick=thicks[i], color=colors[i]
          oplot,   timearr, expellmasscum/1.e9, thick=thicks[i], color=colors[i], linestyle = 2
       ENDIF ELSE BEGIN
          
;Calculate the total amount of mass accreted
          massaccr  = mrdfits(dirs[i] + '/grp1.mass_at_accr.fits',0)
          zaccr     = mrdfits(dirs[i] + '/grp1.accrz.fits',0)
          massaccr  = massaccr[where(zaccr GE zend)]
          igmass    = median(massaccr)
          early     = mrdfits(dirs[i] + '/early.iord.fits',0)
          massearly = n_elements(early)*igmass
          totalmass = total(massaccr) + massearly
          
          IF i EQ 0 THEN $
             plot, timearr, ejectmasscum/totalmass, xtitle='Time [Gyr]', ytitle=textoidl('Cumulative Outflow / Total Gas Mass'), title=label, psym=10, xrange=[0,14], /xstyle, yrange=[0,0.35], /nodata ;, ymargin=[7,3]
          oplot,   timearr, ejectmasscum/totalmass, thick=thicks[i], color=colors[i]
          oplot,   timearr, expellmasscum /totalmass, thick=thicks[i], color=colors[i], linestyle = 2
;         oplot,   timearr, outflowcum/totalmass,thick=thicks[i], color=colors[i], linestyle = 1
          print,file,totalmass,max(ejectmasscum),max(ejectmasscum/totalmass),max(expellmasscum),max(expellmasscum/totalmass)
       ENDELSE
    ENDELSE

    IF KEYWORD_SET(keys) AND i lt N_ELEMENTS(keys) THEN BEGIN
        l_keys = lb_keys
        l_keys[i] = keys[i]
        l_thicks = lb_thicks
        l_thicks[i] = thicks[i]
        l_color = lb_color
        l_color[i] = colors[i]
        legend,l_keys,color = l_color,thick = l_thicks,/top,/left,box = 0
     ENDIF
;    stop

ENDFOR

stop

IF KEYWORD_SET(outplot) THEN device, /close
END

PRO eject_plots_master,outplot = outplot
  formatplot,outplot = outplot  ;,/thick
  spawn,'hostname',hostname
  IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

  steps  = ['00492']
  dirtop = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/']
  eject_plots,dirtop,/colors
END
