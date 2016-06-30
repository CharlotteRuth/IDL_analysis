;------------------------- Metallicity of feedback comes from ---------
PRO eject_metallicity,dir,filenames,scale = scale,normalize = normalize,outplot = outplot,keys = keys,colors = colors,linestyles = linestyles,thicks = thicks,label = label,ctables = ctables,formatthick = formatthick
formatplot,outplot = outplot,thick = formatthick
n = n_elements(dir)
zsolar  =  0.0130215
xrange_fb = [-3,0]
IF NOT keyword_set(scale) THEN scale = fltarr(n) + 1
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
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(linestyles) THEN linestyles = findgen(n)*2  
ENDELSE
lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
lb_linestyle = intarr(n)

IF keyword_set(outplot) AND keyword_set(normalize)     THEN device,filename = outplot + '_fbz_norm.eps',/encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8
IF keyword_set(outplot) AND NOT keyword_set(normalize) THEN device,filename = outplot + '_fbz.eps',     /encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8

FOR i = 0, n_elements(dir)-1 DO BEGIN
   loadct,ctables[i]
   cd,dir[i]
   splitbase = strsplit(filenames[i],'.')
   base = strmid(filenames[i],0,splitbase[n_elements(splitbase) - 1] - 1)

   totalmass = 1
;  eject_gas_character,dirhead[i],ejecthistory,expellhistory,totalmass = totalmass
   ejecthistory =  mrdfits(dir[i] + 'grp1.eject_disk.fits',1)
   expellhistory = mrdfits(dir[i] + 'grp1.expell_disk.fits',1)
   IF i EQ 0 AND     keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,xrange = xrange_fb,xtitle = 'Z/Z'+sunsymbol(),ytitle = '1/M dM/dr',title = label,yrange = [0,1],   /nodata
   IF i EQ 0 AND NOT keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            xrange = xrange_fb,xtitle = 'Z/Z'+sunsymbol(),ytitle = 'dM/dr'    ,title = label,yrange = [0,2e9], /nodata ;,xmargin = [16,3]
   IF                keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,/overplot,color = colors[i],thick = thicks[i],linestyle = linestyles[i]
   IF NOT            keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            /overplot,color = colors[i],thick = thicks[i],linestyle = linestyles[i]
;   histogramp,                                              alog10(expellhistory.metallicity/zsolar),nbins = 100,max = xrange_fb[1],weight = expellhistory.mass*scale[i],           /overplot,color = colors[i],thick = thicks[i],linestyle = 2
   IF keyword_set(keys) AND i lt n_elements(keys) THEN BEGIN
      l_keys = lb_keys
      l_keys[i] = keys[i]
      l_thicks = lb_thicks
      l_thicks[i] = thicks[i]
      l_color = lb_color
      l_color[i] = colors[i]
      l_linestyle = lb_linestyle
      l_linestyle[i] = linestyles[i]
      legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/top,/right,box = 0
   ENDIF
   print,max(weighted_histogram(sqrt(ejecthistory.x*ejecthistory.x + ejecthistory.y*ejecthistory.y),weight = totalmass/1e8,max = xrange_fb[1],nbins = 100))

   undefine,totalmass
   undefine,ejecthistory
   undefine,expellhistory
ENDFOR    
IF keyword_set(outplot) THEN device,/close
END
