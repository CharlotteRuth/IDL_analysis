;------------------------- Radius where feedback comes from ---------
PRO reaccr_radius,dir,filenames,halo = halo,datafile = datafile,halfr = halfr,scalem = scalem,scaler = scaler,normalize = normalize,outplot = outplot,keys = keys,colors = colors,linestyles = linestyles,thicks = thicks,label = label,ctables = ctables,formatthick = formatthick,ratio = ratio

formatplot,outplot = outplot,thick = formatthick
n = n_elements(dir)
xrange_fb = [0,8.0]
IF NOT keyword_set(scaler) THEN scaler = fltarr(n) + 1
IF NOT keyword_set(scalem) THEN scalem = fltarr(n) + 1
IF NOT keyword_set(halo) THEN halo = strtrim(fltarr(n) + '1',2)
IF NOT keyword_set(datafile) THEN datafile = 'reaccrdisk_history'
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
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = findgen(n)*2  
ENDELSE
halfr = fltarr(n)

IF keyword_set(outplot) THEN device,filename = outplot + '_eject_reaccr_r.eps',/encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8 


FOR i = 0, n_elements(dir)-1 DO BEGIN
   loadct,ctables[i]
   cd,dir[i]
   splitbase = strsplit(filenames[i],'.')
   base = strmid(filenames[i],0,splitbase[n_elements(splitbase) - 1] - 1)

   totalmass = 1
   history =  mrdfits(dir[i] + 'grp' + halo[i] + '.' + datafile + '.fits',1)
   print,dir[i] + 'grp' + halo[i] + '.' + datafile
 
   histrad    = weighted_histogram(sqrt(history.x*history.x + history.y*history.y)*scaler[i], nbins = 100,max = xrange_fb[1],weight = history.mass*scalem[i])
   histradcum = weighted_histogram(sqrt(history.x*history.x + history.y*history.y)*scaler[i], nbins = 100,max = xrange_fb[1],weight = history.mass*scalem[i],locations = locations,/cum)
   uniqind = uniq(histradcum)
   halfr[i] = spline(histradcum[uniqind]/max(histradcum),locations[uniqind],0.5)

   IF i EQ 0 AND     keyword_set(normalize) THEN histogramp,sqrt(history.x*history.x   + history.y*history.y)*scaler[i],  nbins = 100,max = xrange_fb[1],weight = history.mass*scalem[i], /normalize,xrange = xrange_fb,xtitle = 'Radius [kpc]',ytitle = '1/M dM/dr',title = label,yrange = [0,0.15],   /nodata
   IF i EQ 0 AND NOT keyword_set(normalize) THEN histogramp,sqrt(history.x*history.x   + history.y*history.y)*scaler[i],  nbins = 100,max = xrange_fb[1],weight = history.mass*scalem[i],            xrange = xrange_fb,xtitle = 'Radius [kpc]',ytitle = 'dM/dr'    ,title = label,yrange = [0,1.2e8],/nodata ;,xmargin = [16,3]
   IF                keyword_set(normalize) THEN histogramp,sqrt(history.x*history.x   + history.y*history.y)*scaler[i],  nbins = 100,max = xrange_fb[1],weight = history.mass*scalem[i], /normalize,/overplot,color = colors[i],thick = thicks[i],linestyle = linestyles[i]
   IF NOT            keyword_set(normalize) THEN histogramp,sqrt(history.x*history.x   + history.y*history.y)*scaler[i],  nbins = 100,max = xrange_fb[1],weight = history.mass*scalem[i],            /overplot,color = colors[i],thick = thicks[i],linestyle = linestyles[i]
   oplot,[halfr[i],halfr[i]],[0,1],linestyle = 2,color = colors[i],thick = thicks[i]

   IF keyword_set(keys) THEN legend,keys,color = colors,linestyle = linestyles,thick = thicks,/top,/right,/box
   print,max(weighted_histogram(sqrt(history.x*history.x + history.y*history.y),weight = totalmass/1e8,max = xrange_fb[1],nbins = 100))
   print,dir[i] + 'grp' + halo[i] + '.' + datafile + '.fits'
   undefine,totalmass
   undefine,history
ENDFOR    
IF keyword_set(outplot) THEN device,/close ELSE stop
END
