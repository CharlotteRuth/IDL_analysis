;Charlotte Christensen
;6/5/12
;This program plots the density and temperature of the gas particles
;that formed stars

;dir = '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.3072g1bwK'
;filebase = 'h603.cosmo50cmb.3072g1bwK'

;dir = '/astro/store/nbody2/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK'
;filebase = 'h603.cosmo50cmb.3072gs1MbwK'

;dir = '/astro/store/nbody3/christensen/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK'
;filebase = 'h603.cosmo50cmb.3072g14HBWK'

;sfgas_contour,dir,filebase

PRO sfgas_contour,dirs,filebases,molecularH = molecularH,zdivide = zdivide,outplot = outplot, keys = keys, colors = colors, thicks = thicks, linestyles = linestyles,label = label,ctables = ctables,min1 = min1, max1 = max1, min2 = min2, max2 = max2,nlevels = nlevels,maxlevel = maxlevel,formatthick = formatthick
n = n_elements(filebases)
formatplot,outplot = outplot,thick = formatthick
IF KEYWORD_SET(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 18
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 600
    ysize = 600
ENDELSE
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF NOT keyword_set(ctables) then ctables = [39,39,39]
    IF colors[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) then ctables = [0,0,0]
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = reverse(findgen(n)*2)   
ENDELSE

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
lb_linestyle = intarr(n)

IF NOT keyword_set(maxlevel) THEN maxlevel = 2.5e9
IF NOT keyword_set(nlevels) THEN nlevels = 5 ;50
IF NOT keyword_set(min1) THEN min1 = 1 ;-1
IF NOT keyword_set(max1) THEN max1 = 3.5 ;4.5
IF NOT keyword_set(min2) THEN min2 = 1.5 ;1
IF NOT keyword_set(max2) THEN max2 = 4.0 ;4.5
binsize1 = 0.1
binsize2 = 0.1
min1 = min1 - binsize1
max1 = max1 + binsize1
min2 = min2 - binsize2
max2 = max2 + binsize2

IF keyword_set(outplot) THEN device,filename=outplot + '_sfgascont.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,0,xsize = xsize, ysize = ysize

FOR i = 0, n - 1 DO BEGIN
    cd,dirs[i]
    starlog = rstarlog(filebases[i] + '.starlog',molecularH = molecularH[i])
    units = tipsyunits(filebases[i] + '.param')
    starlog.timeform = starlog.timeform*units.timeunit
    starlog.massform = starlog.massform*units.massunit
    starlog.rhoform  = starlog.rhoform *units.rhounit
    starlogzform = z_from_time(starlog.timeform)

    phasehist = better_hist2d(alog10(starlog.rhoform),alog10(starlog.tempform),starlog.massform, min1 = min1, max1 = max1, min2 = min2, max2 = max2, binsize1 = binsize1,binsize2 = binsize2)
    xarray = (findgen((size(phasehist))[1]) + 0.5)*binsize1 + min1
    yarray = (findgen((size(phasehist))[2]) + 0.5)*binsize2 + min2

    IF keyword_set(zdivide) THEN BEGIN
        indz4 = where(starlogzform GE 4)
        phasehistz4 = better_hist2d(alog10(starlog[indz4].rhoform),alog10(starlog[indz4].tempform),starlog[indz4].massform, min1 = min1, max1 = max1, min2 = min2, max2 = max2, binsize1 = binsize1,binsize2 = binsize2)
        
;        indz2 = where(starlogzform GE 2 AND starlogzform LT 4)
;        phasehistz2 = better_hist2d(alog10(starlog[indz2].rhoform),alog10(starlog[indz2].tempform),starlog[indz2].massform, min1 = min1, max1 = max1, min2 = min2, max2 = max2, binsize1 = binsize1,binsize2 = binsize2)
        
        indz1 = where(starlogzform GE 1 AND starlogzform LT 4)
        phasehistz1 = better_hist2d(alog10(starlog[indz1].rhoform),alog10(starlog[indz1].tempform),starlog[indz1].massform, min1 = min1, max1 = max1, min2 = min2, max2 = max2, binsize1 = binsize1,binsize2 = binsize2)
        
        indz0 = where(starlogzform LT 1)
        phasehistz0 = better_hist2d(alog10(starlog[indz0].rhoform),alog10(starlog[indz0].tempform),starlog[indz0].massform, min1 = min1, max1 = max1, min2 = min2, max2 = max2, binsize1 = binsize1,binsize2 = binsize2)        
    ENDIF

;contour,alog10(phasehist),xarray,yarray,nlevels = 254,min = 0,/cell_fill,xrange = [min1,max1],yrange = [min2,max2],xstyle = 1,ystyle = 1
;contour,alog10(phasehist)  ,xarray,yarray,nlevels = 10;,/overplot
;contour,alog10(phasehistz4),xarray,yarray,nlevels = 10,/overplot,color = 254;,linestyle = 1
;contour,alog10(phasehistz2),xarray,yarray,nlevels = 10,/overplot,color = 150
;contour,alog10(phasehistz1),xarray,yarray,nlevels = 10,/overplot,color = 80
;contour,alog10(phasehistz0),xarray,yarray,nlevels = 10,/overplot,color = 30

    loadct,ctables[i]
    levels = [4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5]
    IF i EQ 0 THEN contour,alog10(phasehist),xarray,yarray,levels = levels,xrange = [min1 + binsize1,max1 - binsize1],yrange = [min2 + binsize2,max2 - binsize2],xstyle = 1,ystyle = 1,/nodata,xtitle = textoidl('log(\rho)') +' [amu/cc]',ytitle = 'log(T) [K]',title = label,/closed;,min = 1,max = maxlevel,nlevels = nlevel
    
    IF keyword_set(zdivide) THEN BEGIN
        loadct,39
        c4 = 254
        c1 = 170
        c0 = 85
;        contour,phasehistz2,xarray,yarray,nlevels = nlevels,min = 1,max = maxlevel,/overplot,c_colors = 150,c_thick = thicks[i]
        contour,phasehistz1,xarray,yarray,nlevels = nlevels,min = 1,max = maxlevel,/overplot,c_colors = c1,c_thick = thicks[i],/closed
        contour,phasehistz0,xarray,yarray,nlevels = nlevels,min = 1,max = maxlevel,/overplot,c_colors = c0,c_thick = thicks[i],/closed
        contour,phasehistz4,xarray,yarray,nlevels = nlevels,min = 1,max = maxlevel,/overplot,c_colors = c4,c_thick = thicks[i],/closed
        legend,['4 < z','1 < z < 4 ','z < 1'],color = [c4,c1,c0],thick = [thicks[i],thicks[i],thicks[i]],/right,/top,box = 0,linestyle = [0,0,0]
    ENDIF ELSE  contour,alog10(phasehist),xarray,yarray,levels = levels,c_colors = colors[i],c_linestyle = linestyles[i],c_thick = thicks[i],/overplot,/closed;,min = 1, max = maxlevel,nlevels = nlevel
;    stop
;    IF keyword_set(keys) AND i LT n_elements(keys) THEN BEGIN
;        l_keys = lb_keys
;        l_keys[i] = keys[i]
;        l_thicks = lb_thicks
;        l_thicks[i] = thicks[i]
;        l_color = lb_color
;        l_color[i] = colors[i]
;        l_linestyle = lb_linestyle
;        l_linestyle[i] = linestyles[i]
;        legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/right,/top,box = 0
;    ENDIF
ENDFOR
IF keyword_set(keys) THEN legend,keys,thick = thicks, color = colors, linestyle = linestyles,ctables = ctables,/left,/bottom,box = 0
IF 1 THEN BEGIN
    i = 0
    cd,dirs[i]
    starlog = rstarlog(filebases[i] + '.starlog',molecularH = molecularH[i])
    units = tipsyunits(filebases[i] + '.param')
    starlog.timeform = starlog.timeform*units.timeunit
    starlog.massform = starlog.massform*units.massunit
    starlog.rhoform  = starlog.rhoform *units.rhounit
    starlogzform = z_from_time(starlog.timeform)

    phasehist = better_hist2d(alog10(starlog.rhoform),alog10(starlog.tempform),starlog.massform, min1 = min1, max1 = max1, min2 = min2, max2 = max2, binsize1 = binsize1,binsize2 = binsize2)
    xarray = (findgen((size(phasehist))[1]) + 0.5)*binsize1 + min1
    yarray = (findgen((size(phasehist))[2]) + 0.5)*binsize2 + min2


    loadct,ctables[i]
    contour,phasehist,xarray,yarray,nlevels = nlevels,min = 1, max = maxlevel,c_colors = colors[i],c_linestyle = linestyles[i],c_thick = thicks[i],/overplot,/closed
ENDIF

IF keyword_set(outplot) THEN device,/close
END
