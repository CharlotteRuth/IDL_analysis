;dir = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir']
;dir =['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir','/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir']


pro fft_result,dir,outplot = outplot,keys = keys,color = color,symbols = symbols,ctables = ctables,formatthick = formatthick,label = label,symsizes = symsizes
rad = 12
diameter = rad*2
n = N_ELEMENTS(dir)
formatplot,outplot = outplot,thick = formatthick
names = strarr(n)
IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(color) THEN BEGIN
    obscolor = fgcolor
    if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
    if NOT keyword_set(ctables) then ctables = [39,39,39]
;    IF NOT KEYWORD_SET(THICKS) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(symbols) THEN symbols = fltarr(n) + 4
ENDIF ELSE BEGIN
    obscolor = fgcolor ;150
    colors = fltarr(n) + fgcolor;  fltarr(n) + 5
    if NOT keyword_set(ctables) then ctables = [0,0,0]
;    IF NOT KEYWORD_SET(THICKS) THEN thicks = fltarr(n) + 2 ;(findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(symbols) THEN symbols = (findgen(n)+2)*2
ENDELSE
symthick = 5
IF (KEYWORD_SET(outplot)) THEN device,filename=outplot+'_fft_s.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 392,ysize = 392

lb_keys = strarr(n)
;lb_thicks = intarr(n)
lb_color = intarr(n)
lb_linestyle = intarr(n)
lb_symbols = intarr(n)

FOR id = 0, N_ELEMENTS(dir) - 1 DO BEGIN
    cd,dir[id]
    loadct,ctables[id]    
    real = mrdfits("FFTW_R.fits")
    imag = mrdfits("FFTW_I.fits")
    mag = real*real + imag*imag
    nx = (size(real))[1]
    ny = (size(real))[2]
    ny_orig = 2*(ny - 1)
    all = fltarr(nx,ny_orig)
    
    all[*,0:ny - 1] = mag
    offset = ny - (ny_orig - ny)
    all[*,ny : ny_orig - 1] = (rotate(mag,2))[*,offset : ny - 1]
    
    xaxes = findgen(nx) - (nx - 1)/2.0
    yaxes = findgen(ny_orig) - (ny_orig - 1)/2.0
    radius = fltarr(N_ELEMENTS(xaxes),N_ELEMENTS(yaxes))
    FOR ix = 0,N_ELEMENTS(xaxes) - 1 DO $
      FOR iy = 0,N_ELEMENTS(yaxes) - 1 DO $
      radius[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy])
    
;    loadct,39
;    contour,alog10(all),xaxes,yaxes,nlevels = 20,/fill,/xstyle,/ystyle
    maxr = 512
    dr = 4
    nbins = maxr/dr
    radii = findgen(nbins)*dr
    avepower = fltarr(nbins)
    stdevpower = fltarr(nbins)
    
    FOR ir = 0, nbins - 1 DO BEGIN
        ind = where(radius ge radii[ir] AND radius lt (radii[ir] + dr))
        avepower[ir] = MEAN(mag[ind])
        stdevpower[ir] = stdev(mag[ind])
    ENDFOR
    if id eq 0 then plot,radii,avepower,ytitle = 'Ave Power',xtitle = 'Frequency',/xstyle,/xlog,/ylog,xrange = [1,700],yrange = [1e5,2e9],/nodata,title = label ;[1e4,2e9]
    oplot,radii,avepower,psym = symcat(symbols[id],thick = symthick),color = colors[id],symsize = symsizes[id]
    ind = where(radii gt 5 AND radii le 150)
    xfitp = alog10(radii[ind])
    yfitp = alog10(avepower[ind]) 

    fit = poly_fit(xfitp,yfitp,1)
    oplot,10^(xfitp),10^(xfitp * fit[1] + fit[0]),color = obscolor,thick = 2
;lambda[ind])
    names[id] = 'Powerlaw Fit, k = ' + strtrim(fit[1],2)
;    IF KEYWORD_SET(keys) THEN BEGIN
;        l_keys = lb_keys
;        l_keys[id] = keys[id]
;        l_thicks = lb_thicks
;        l_thicks[id] = thicks[id]
;        l_color = lb_color
;        l_color[id] = colors[id]
;;        l_linestyle = lb_linestyle
;;        l_linestyle[id] = l_linestyle[id]
;        l_symbols = lb_symbols
;        l_symbols[id] = symbols[id]
;        legend,l_keys,color = l_color,thick = l_thicks,psym = l_symbols,/right,box = 0
;    ENDIF
ENDFOR
legend,names,linestyle = fltarr(n),/bottom,/left,ctables = ctables,psym = symbols,color = colors,thick = symthick;,thick = thicks
IF (KEYWORD_SET(outplot)) THEN device,/close ELSE stop


IF (KEYWORD_SET(outplot)) THEN device,filename=outplot+'_fft_k.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 392,ysize = 392

FOR id = 0, N_ELEMENTS(dir) - 1 DO BEGIN
    cd,dir[id]
    loadct,ctables[id]
    real = mrdfits("FFTW_R.fits")
    imag = mrdfits("FFTW_I.fits")
    mag = real*real + imag*imag
    nx = (size(real))[1]
    ny = (size(real))[2]
    ny_orig = 2*(ny - 1)
    all = fltarr(nx,ny_orig)
    
    all[*,0:ny - 1] = mag
    offset = ny - (ny_orig - ny)
    all[*,ny : ny_orig - 1] = (rotate(mag,2))[*,offset : ny - 1]
    
    xaxes = findgen(nx) - (nx - 1)/2.0
    yaxes = findgen(ny_orig) - (ny_orig - 1)/2.0
    radius = fltarr(N_ELEMENTS(xaxes),N_ELEMENTS(yaxes))
    FOR ix = 0,N_ELEMENTS(xaxes) - 1 DO $
      FOR iy = 0,N_ELEMENTS(yaxes) - 1 DO $
      radius[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy])
    
;    loadct,39
;    contour,alog10(all),xaxes,yaxes,nlevels = 20,/fill,/xstyle,/ystyle
    maxr = 512
    dr = 4
    nbins = maxr/dr
    radii = findgen(nbins)*dr
    avepower = fltarr(nbins)
    stdevpower = fltarr(nbins)
    
    FOR ir = 0, nbins - 1 DO BEGIN
        ind = where(radius ge radii[ir] AND radius lt (radii[ir] + dr))
        avepower[ir] = MEAN(mag[ind])
        stdevpower[ir] = stdev(mag[ind])
    ENDFOR

    lambda = diameter/(radii + dr/2)*1000 ;/!PI)*1000
    min_lambda = 170
    ind = where(lambda ge min_lambda AND lambda le 6000)
    if id eq 0 then plot,lambda[ind],avepower[ind],psym = symcat(symbols[id]),ytitle = 'Average Power ' + textoidl('[n_{HI}^2 cm^{-4}]'),xtitle = 'Spatial Scale [pc]',xrange = [10000,100],yrange = [1e5,5e9],/xstyle,/xlog,/ylog,/nodata,title = label;,xmargin = [20,3],[1e5,2e9
    oplot,lambda[ind],avepower[ind],psym = symcat(symbols[id],thick = symthick),color = colors[id],symsize = symsizes[id]
    xfitp = alog10(lambda[ind])
    yfitp = alog10(avepower[ind])

    fit = poly_fit(xfitp,yfitp,1)
    oplot,10^(xfitp),10^(xfitp * fit[1] + fit[0]),color = obscolor,thick = 2
;lambda[ind])
    names[id] = 'Powerlaw Fit, k = ' + strtrim(fit[1],2)
;    IF KEYWORD_SET(keys) THEN BEGIN
;        l_keys = lb_keys
;        l_keys[id] = keys[id]
;        l_thicks = lb_thicks
;        l_thicks[id] = thicks[id]
;        l_color = lb_color
;        l_color[id] = colors[id]
;        l_symbols = lb_symbols
;        l_symbols[id] = symbols[id]
;        legend,l_keys,color = l_color,thick = l_thicks,psym = l_symbols,/right,box = 0
;    ENDIF
ENDFOR
IF keyword_set(keys) THEN legend,keys,linestyle = fltarr(n),/bottom,psym = symbols,color = colors,/left,ctables = ctables,box = 0,thick = symthick;,thick = fltarr(n) + thicks[0],box = 0
print,names
IF (KEYWORD_SET(outplot)) THEN device,/close ELSE stop


end
