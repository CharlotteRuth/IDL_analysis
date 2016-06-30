PRO plot_z_m,modez,dirs,files,halo = halo,colors = colors,outplot = outplot,normalize = normalize,symbols = symbols,formatthick = formatthick,absolute = absolute
zsolar  =  0.0130215
IF NOT keyword_set(finalstep) THEN finalstep = '00512'
n = n_elements(dirs)
;formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
;    xsize = 18
;    ysize = 12
;    IF keyword_set(stellarmass) THEN outplot = outplot + '_sm' ELSE outplot = outplot + '_vm'
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
;    xsize = 800
;    ysize = 500
ENDELSE
IF keyword_set(stellarmass) THEN BEGIN
   IF NOT keyword_set(absolute) THEN xtitle = 'Stellar Mass [M' +sunsymbol() + ']'
   xrange = [1e7,1e11]
ENDIF ELSE BEGIN
   IF NOT keyword_set(absolute) THEN xtitle = 'Virial Mass [M' +sunsymbol() + ']'
   xrange = [1e10,1e12]
ENDELSE
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symbols) THEN symbols = [15]
    symbols_boarder = [6]
    IF NOT keyword_set(symsize) THEN symsize = 2
    colore = [254]
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2)  
    IF NOT keyword_set(symbols) THEN symbols = [15]
    symbols_boarder = [6]
    IF NOT keyword_set(symsize) THEN symsize = 2
    colore = [254]
ENDELSE

vmass = fltarr(n)
smass = fltarr(n)
zgas = fltarr(n)
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group eq halo[i]))[0]
   vmass[i] = stat[ind].m_tot
   smass[i] = stat[ind].m_star
   readcol,dirs[i] + '/grp' + halo[i] + '.metals.txt',z,ox,fe,H2,HI,H
   zgas[i] = z[n_elements(z) - 1]/h[n_elements(z) - 1]/zsolar
ENDFOR

IF keyword_set(absolute) THEN BEGIN
;    outplotext = outplotext + '_zm_abs.eps'
    ytitle = 'Log(Z/Z'+sunsymbol()+')'
    yrange = [-1,0.2]
ENDIF ELSE BEGIN
;    outplotext = outplotext + '_zm_scale.eps'
    ytitle = textoidl('Log(Z_{eject}/Z_{gas})')
    yrange = [-0.05,0.6]
ENDELSE

IF keyword_set(stellarmass) THEN xaxis = smass ELSE xaxis = vmass
;IF keyword_set(outplot) THEN  device,filename = outplot + outplotext,/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 0, xsize = xsize, ysize = ysize
plot,xaxis,alog10(modez),/xlog,xrange = xrange,xtitle = xtitle,yrange = yrange,ytitle = ytitle,psym = symcat(symbols[0]),symsize = symsize,/nodata
IF keyword_set(absolute) THEN oplot,xaxis,alog10(zgas),psym = symcat(16),symsize = symsize,color = fgcolor
oplot,xaxis,alog10(modez),psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
IF keyword_set(symbols_boarder) THEN oplot,xaxis,alog10(modez),psym = symcat(symbols_boarder[0]),symsize = symsize,color = fgcolor,thick = 1
IF NOT keyword_set(absolute) THEN oplot,xrange,[0,0],thick = thicks[0]
;legend,['Gas Ejected','Star Formation'],psym = [symbols[1],46],color = [colore[1],fgcolor],/bottom,/right,box =0
;IF keyword_set(outplot) THEN device, /close ELSE stop
END
