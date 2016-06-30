PRO outflow_plots, dirs, unscale = unscale, outplot = outplot, keys = keys, colors = colors,thicks = thicks, linestyles = linestyles,label = label,ctables = ctables,yrange_fbcum = yrange_fbcum
formatplot,outplot = outplot
n = N_ELEMENTS(dirs)

IF KEYWORD_SET(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
   IF KEYWORD_SET(multiframe) THEN BEGIN
        xsize = 10*n
        ysize = 12
        mxTitSize = 1.5
        mxTitOffset = 2
    ENDIF ELSE BEGIN
        xsize = 18
        ysize = 12
    ENDELSE
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    IF KEYWORD_SET(multiframe) THEN BEGIN
        xsize = 400*n
        ysize = 475
        mxTitSize = 1.5
        mxTitOffset = 1
    ENDIF ELSE BEGIN
        xsize = 400
        ysize = 266
    ENDELSE
ENDELSE
IF KEYWORD_SET(colors) THEN BEGIN
    loadct,39
    if NOT keyword_set(ctables) then ctables = [39,39,39]
    if colors[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
ENDIF ELSE BEGIN
    loadct,0    
    if NOT keyword_set(ctables) then ctables = [0,0,0]
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT KEYWORD_SET(thicks) THEN thicks = (findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = REVERSE(findgen(n)*2)   
ENDELSE
IF NOT KEYWORD_SET (yrange_fbcum) THEN yrange_fbcume = [0,250]

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
lb_linestyle = intarr(n)
maxtime = wmap3_lookback(1000)

IF KEYWORD_SET(outplot) THEN  device,filename = outplot + '_fbcum.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,0
FOR i = 0, N_ELEMENTS(dirs) - 1 DO BEGIN
    loadct, ctables[i]
    spawn,'ls ' + dirs[i] + '*512/*512.coolontime',file_coolon
    spawn,'ls ' + dirs[i] + '*512/*512.iord',file_iord
    spawn,'ls ' + dirs[i] + '*512/*cosmo**512',file
    rtipsy,file,h,g,d,s,/justhead
    readarr,file_coolon,h,coolon,part = 'gas',/ascii
    readarr,file_iord,h,iord,part = 'gas',/ascii

;Ejection defined by Amiga
;    massout = mrdfits(dirs[i]+'/grp1.mass_at_outflow.fits',0)
;    zout = mrdfits(dirs[i]+'/grp1.outflow_z.fits',0)
;    iordout = mrdfits(dirs[i]+'/grp1.outflow_iord.fits',0)
;    ind99 = where(zout EQ 99, comp=n99)
;    massout = massout[n99]
;    zout = zout[n99]
;    iordout = iordout[n99]
;    lb = wmap3_lookback(zout)
;    match,iord,iordout,ind_iord,ind_iordout
;    coolonout = where(coolon[ind_iord] ne 0)

;    massout = massout[ind_iordout[coolonout]]
;    zout = zout[ind_iordout[coolonout]]
;    iordout = iordout[ind_iordout[coolonout]]
;    lb = lb[ind_iordout[coolonout]]

 ;   histogramp,coolon[ind_iord[coolonout]]*units.timeunit,nbins = 100
 ;   histogramp,maxtime - lb[ind_iordout[coolonout]],nbins = 100,/overplot

 ;   histogramp,coolon[ind_iord[coolonout]]*units.timeunit - (maxtime - lb[ind_iordout[coolonout]]),nbins = 100

;    hlb = histogram((13.7-lb/1.e9), binsize=0.5, min=0, reverse_indices=ri)
;    mass = fltarr(n_elements(hlb))
;    masscum = fltarr(n_elements(hlb))
;    FOR j=0,n_elements(hlb)-1 DO IF ri[j+1] GT ri[j] THEN mass[j] = total(massout[ri[ri[j]:ri[j+1]-1]])
;    FOR j=0,n_elements(hlb)-1 DO masscum[j] = total(mass[0:j])

;Ejection defined by virial radius
    massout_rvir = mrdfits(dirs[i]+'/grp1.rvir.mass_at_outflow.fits',0)
    zout_rvir = mrdfits(dirs[i]+'/grp1.rvir.outflow_z.fits',0)
    iordout_rvir = mrdfits(dirs[i]+'/grp1.rvir.outflow_iord.fits',0)
    ind99_rvir = where(zout_rvir EQ 99, comp=n99_rvir) 
    massout_rvir = massout_rvir[n99_rvir]
    zout_rvir = zout_rvir[n99_rvir]
    iordout_rvir = iordout_rvir[n99_rvir]
;    lb_rvir = wmap3_lookback(zout_rvir[n99_rvir])
    lb_rvir = wmap3_lookback(zout_rvir)
    match,iord,iordout_rvir,ind_iord_rvir,ind_iordout_rvir
    coolonout_rvir = where(coolon[ind_iord_rvir] ne 0)

    massout_rvir = massout_rvir[ind_iordout_rvir[coolonout_rvir]]
    zout_rvir = zout_rvir[ind_iordout_rvir[coolonout_rvir]]
    iordout_rvir= iordout_rvir[ind_iordout_rvir[coolonout_rvir]]
    lb_rvir = lb_rvir[ind_iordout_rvir[coolonout_rvir]]

    hlb = histogram((13.7-lb_rvir/1.e9), binsize=0.5, min=0, reverse_indices=ri)
    mass_rvir = fltarr(n_elements(hlb))
    masscum_rvir = fltarr(n_elements(hlb))
    FOR j=0,n_elements(hlb)-1 DO IF ri[j+1] GT ri[j] THEN mass_rvir[j] = total(massout_rvir[ri[ri[j]:ri[j+1]-1]])
    FOR j=0,n_elements(hlb)-1 DO masscum_rvir[j] = total(mass_rvir[0:j])

    IF KEYWORD_SET(unscale) THEN BEGIN
        IF i EQ 0 THEN plot, indgen(28)*0.5, masscum/1.e8, xtitle='Time [Gyr]', ytitle=textoidl('Cumulative Outflow / 10^{8} [M')+sunsymbol()+']', title=label, psym=10, xrange=[0,14], /xstyle, yrange=yrange_fbcum, /nodata
        oplot, indgen(28)*0.5, masscum/1.e8, thick=thicks[i], color=colors[i],linestyle = linestyles[i] ;, psym=10
    ENDIF ELSE BEGIN
        massaccr = mrdfits(dirs[0]+'/grp1.mass_at_accr.fits',0)
        igmass = max(massaccr)
        early = mrdfits(dirs[0]+'/early.iord.fits',0)
        massearly = N_ELEMENTS(early)*igmass
        totalmass = total(massaccr) + massearly
        IF i EQ 0 THEN plot, indgen(28)*0.5, masscum_rvir/totalmass, xtitle='Time [Gyr]', ytitle=textoidl('Cumulative Outflow / Total Gas Mass'), title=label, psym=10, xrange=[0,14], /xstyle, yrange=yrange_fbcum, /nodata, ymargin=[7,3]
;        oplot, indgen(28)*0.5, masscum/totalmass, thick=thicks[i], color=colors[i],linestyle = linestyles[i] ;, psym=10
        oplot, indgen(28)*0.5, masscum_rvir/totalmass, thick=thicks[i], color=colors[i],linestyle = linestyles[i];, psym = -4 ;, psym=10        
        print,file,totalmass,max(masscum_rvir),max(masscum_rvir/totalmass)
    ENDELSE


    IF KEYWORD_SET(keys) AND i lt N_ELEMENTS(keys) THEN BEGIN
        l_keys = lb_keys
        l_keys[i] = keys[i]
        l_thicks = lb_thicks
        l_thicks[i] = thicks[i]
        l_color = lb_color
        l_color[i] = colors[i]
        l_linestyle = lb_linestyle
        l_linestyle[i] = linestyles[i]
        legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/top,/left,box = 0
    ENDIF
ENDFOR

IF KEYWORD_SET(outplot) THEN device, /close
END

PRO outflow_plots_master,outplot = outplot
  formatplot,outplot = outplot  ;,/thick
  IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

  steps  = ['00492']
  dirtop = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.' + steps[0] + '.dir']
  outflow_plots,dirtop
END
