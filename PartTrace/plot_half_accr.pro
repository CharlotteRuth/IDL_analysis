PRO plot_half_accr,dirs,files,halo = halo,outplot = outplot,colors = colors,finalstep = finalstep,stellarmass = stellarmass,symsize = symsize, symbols = symbols,formatthick = formatthick,sfr = sfr
IF NOT keyword_set(finalstep) THEN finalstep = '00512'
n = n_elements(dirs)
formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
    IF keyword_set(stellarmass) THEN outext = '_sm_' ELSE outext = '_vm_'
    IF keyword_set(sfr) THEN outext = outext + 'orig_' ELSE outext = outext + 'final_'
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 800
    ysize = 500
 ENDELSE
IF keyword_set(stellarmass) THEN BEGIN
   xtitle = 'Stellar Mass [M' +sunsymbol() + ']'
   xrange = [1e7,1e11]
ENDIF ELSE BEGIN
   xtitle = 'Virial Mass [M' +sunsymbol() + ']'
   xrange = [3e9,1e12]
ENDELSE
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symbols) THEN symbols = [17,15]
    symbols_boarder = [5,6] 
    IF NOT keyword_set(symsize) THEN symsize = 2
    colore = [120,254]
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2)  
    IF NOT keyword_set(symbols) THEN symbols = [14,15]
    symbols_boarder = [4,6] 
    IF NOT keyword_set(symsize) THEN symsize = 2
    colore = [50,254]
ENDELSE

eject_radius,dirs,files,halo = halo,colors = colors,linestyles = linestyles,thicks = thicks,/normalize,halfr = eject_halfr
reaccr_radius,dirs,files,halo = halo,colors = colors,linestyles = linestyles,thicks = thicks,/normalize,halfr = reaccr_halfr

vmass = fltarr(n)
smass = fltarr(n)
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group eq halo[i]))[0]
   vmass[i] = stat[ind].m_tot
   smass[i] = stat[ind].m_star
ENDFOR

formatplot,outplot = outplot
IF keyword_set(stellarmass) THEN xaxis = smass ELSE xaxis = vmass

;--------------------------- Reeject and Reaccretion
IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'ejectr_accr.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xaxis,eject_halfr,/xlog,xrange = xrange,xtitle = xtitle,yrange = [0,5],ytitle = 'Half Mass Radius [kpc]',psym = symcat(symbols[1]),symsize = symsize,/nodata
oplot,xaxis,eject_halfr,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
oplot,xaxis,eject_halfr,psym = symcat(symbols_boarder[1]),symsize = symsize,color = fgcolor,thick = 1
oplot,xaxis,reaccr_halfr,psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
oplot,xaxis,reaccr_halfr,psym = symcat(symbols_boarder[0]),symsize = symsize,color = fgcolor
legend,['Gas Ejected','Gas Reaccreted'],psym = [symbols[1], symbols[0]],color = [colore[1],colore[0]],/bottom,/right,box =0
IF keyword_set(outplot) THEN device, /close ELSE stop

;----------------------------- Reeject/Reaccretion
IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'accr_by_ejectr.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xaxis,reaccr_halfr/eject_halfr,psym = symcat(symbols[1]),/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Reaccretion Half Radius/Eject Half Radius',symsize = symsize,/nodata,yrange = [1,2.6]
;oplot,xrange,[1,1],thick = thicks[0]
oplot,xaxis,reaccr_halfr/eject_halfr,psym = symcat(16),symsize = symsize,color = fgcolor
IF keyword_set(outplot) THEN device, /close ELSE stop

;----------------------------- Multiplot --------------------------
IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'eject_sfr_multi.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*1.8,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize*1.8
multiplot,[1,2]

plot,xaxis,eject_halfr,/xlog,xrange = xrange,yrange = [0.1,5],ytitle = 'Half Mass Radius [kpc]',psym = symcat(symbols[1]),symsize = symsize,/nodata,/ylog
oplot,xaxis,eject_halfr,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
oplot,xaxis,eject_halfr,psym = symcat(symbols_boarder[1]),symsize = symsize,color = fgcolor,thick = 1
oplot,xaxis,reaccr_halfr,psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
oplot,xaxis,reaccr_halfr,psym = symcat(symbols_boarder[0]),symsize = symsize,color = fgcolor
legend,['Gas Ejected','Gas Reaccreted'],psym = [symbols[1], symbols[0]],color = [colore[1],colore[0]],/bottom,/right,box =0
multiplot

plot,xaxis,reaccr_halfr/eject_halfr,psym = symcat(symbols[1]),/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Reaccretion Half Radius/Eject Half Radius',symsize = symsize,/nodata,yrange = [1,2.6]
;oplot,xrange,[1,1],thick = thicks[0]
oplot,xaxis,reaccr_halfr/eject_halfr,psym = symcat(16),symsize = symsize,color = fgcolor
multiplot,/reset
IF keyword_set(outplot) THEN device, /close ELSE stop


END
