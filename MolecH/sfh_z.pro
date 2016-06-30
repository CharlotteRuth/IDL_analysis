
pro sfh_z,dir,files,pfiles,keys = keys,colors = colors,linestyle = linestyle,outplot = outplot,thick = thick,binsize = binsize,cumlative = cumlative,yrange =yrange,xrange = xrange,_EXTRA=_extra,label = label,ctables = ctables, redshift = redshift,formatthick = formatthick,halos = halos


;!Y.STYLE = 1
;!X.STYLE = 1
n = N_ELEMENTS(files)
IF NOT keyword_set(halos) THEN halos = string(strtrim(intarr(n) + 1,2))
formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0
    bgcolor = 255
    device,filename=outplot + '_sfhz.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 27,xoffset = 2,yoffset = 0
;,xsize = 18,ysize= 12
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    window,0,xsize = 600,ysize = 900
ENDELSE
;print,'outplot: ',outplot
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF colors[0] EQ 1 THEN  colors = (findgen(n) + 1)*240/n ELSE colors = colors
    IF NOT keyword_set(ctables) THEN ctables = intarr(n) + 39
    IF NOT keyword_set(thick) THEN thick = fltarr(n) + 2
    IF NOT keyword_set(linestyle) THEN linestyle = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n)
    IF NOT keyword_set(thick) THEN thick = findgen(n)*2
    IF NOT keyword_set(linestyle) THEN linestyle = reverse(findgen(n)*2)
ENDELSE
IF NOT keyword_set(binsize) THEN binsize = 5e7
IF keyword_set(xrange) THEN maxt =   xrange[1]*1e9
print,linestyle

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
lb_linestyle = intarr(n)
lb_psym = intarr(n) + 3

ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
z = reverse(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
tickred_in_t = reverse(findgen(5))
ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e9

multiplot,[1,3]
FOR i = 0, n - 1 DO BEGIN
   cd,dir[i] + '/../../'
   units = tipsyunits(pfiles[i],/silent)
   filebase = strmid(pfiles[i],0,strlen(pfiles[i]) - 6)
   readcol,filebase + '.grp' + halos[i] + '.haloid.dat',names,haloid,format='(A,F)'
   nn = n_elements(names)
   zarray = fltarr(nn)
   tarray = fltarr(nn)
   FOR j = 0, nn - 1 DO BEGIN
      metals = mrdfits(names[j] + '.metals.fits',1)
      ind = where(metals.grp EQ haloid[j])
      rtipsy,names[j],h,g,d,s,/justhead
      znow = 1/h.time - 1
      tarray[j] = ageUniverse - wmap3_lookback(znow)
      IF metals[ind].ox_sfr GT 0 THEN zarray[j] = metals[ind].ox_sfr ELSE zarray[j] = metals[ind].ox_inneratom
   ENDFOR
   print,xrange
   IF i EQ 0 THEN  plot,tarray/1e9,zarray,ytitle = '12 + log(O/H)',yrange = [6.8,8.1],/nodata,xstyle = 9,xrange = xrange,ystyle = 1
   oplot,tarray/1e9,zarray,thick = thick[i],linestyle = linestyle[i],color = colors[i]
ENDFOR
axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
multiplot,/doyaxis
FOR i = 0, n - 1 DO BEGIN
   loadct,ctables[i]
   units = tipsyunits(pfiles[i])
   rtipsy,files[i],h,g,d,s
   print,xrange
   IF i eq 0 THEN  sfr,s,massunit = units.massunit,timeunit = units.timeunit,binsize=binsize,thick = thick[i],linestyle = linestyle[i],sarray = sarray,tarray = tarray,maxt = maxt,xrange = xrange,yrange = yrange,title = label,/cumlative,/nodata,/ylog,xstyle = 5
   sfr,s,                massunit = units.massunit,timeunit = units.timeunit,binsize=binsize,thick = thick[i],linestyle = linestyle[i],sarray = sarray,tarray = tarray,maxt = maxt,/cumlative,color = colors[i],/overplot
   ind = where(tarray/1e9 gt 10)
ENDFOR
multiplot,/doxaxis,/doyaxis
FOR i = 0, n - 1 DO BEGIN
   loadct,ctables[i]
   units = tipsyunits(pfiles[i])
   rtipsy,files[i],h,g,d,s
   print,xrange
   IF i eq 0 THEN  sfr,s,massunit = units.massunit,timeunit = units.timeunit,binsize=binsize,thick = thick[i],linestyle = linestyle[i],sarray = sarray,tarray = tarray,maxt = maxt,srange = xrange,xstyle = 1,yrange = [0,0.39],title = label,/nodata,ystyle = 1
   sfr,s,                massunit = units.massunit,timeunit = units.timeunit,binsize=binsize,thick = thick[i],linestyle = linestyle[i],sarray = sarray,tarray = tarray,maxt = maxt,color = colors[i],/overplot
   ind = where(tarray/1e9 gt 10)
   print,'Mean: ',mean(sarray(ind))
   print,'STDEV: ',stdev(sarray(ind))
   IF keyword_set(keys) THEN BEGIN
       l_keys = lb_keys
       l_keys[i] = keys[i]
       l_thicks = lb_thicks
       l_thicks[i] = thick[i]
       l_color = lb_color
       l_color[i] = colors[i]
       l_linestyle = lb_linestyle
       l_linestyle[i] = linestyle[i]
       l_psym = lb_psym
       l_psym[i] = 0
       legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,psym = l_psym,/right,/top,box = 0
   ENDIF
ENDFOR
;IF keyword_set(key) THEN legend,key,color = colors,linestyle = linestyle,thick = thick,/right
if keyword_set(outplot) then device,/close else stop
multiplot,/reset
end
