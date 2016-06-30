PRO plot_half_eject,dirs,files,halo = halo,normalize = normalize,outplot = outplot,colors = colors,finalstep = finalstep,stellarmass = stellarmass,symsize = symsize, symbols = symbols,formatthick = formatthick,sfr = sfr
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
    IF NOT keyword_set(symbols) THEN symbols = [14,15,17]
    symbols_boarder = [4,6,5] 
    IF NOT keyword_set(symsize) THEN symsize = 2
    colore = [50,254,120]
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2)  
    IF NOT keyword_set(symbols) THEN symbols = [14,15,17]
    symbols_boarder = [4,6,5] 
    IF NOT keyword_set(symsize) THEN symsize = 2
    colore = [50,254,120]
ENDELSE
symbols = [14,15,17]
eject_radius,dirs,files,halo = halo,colors = colors,linestyles = linestyles,thicks = thicks,/normalize,halfr = halfr
;stop
;eject_radius,dirs,files,halo = halo,colors = colors,linestyles = linestyles,thicks = thicks,/normalize,halfr = halfrexp,datafile = 'expell_disk'
reaccr_radius,dirs,files,halo = halo,colors = colors,linestyles = linestyles,thicks = thicks,halfr = halfr_reaccr

halfsfr = fltarr(n)
vmass = fltarr(n)
smass = fltarr(n)
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group eq halo[i]))[0]
   vmass[i] = stat[ind].m_tot
   smass[i] = stat[ind].m_star

   IF sfr THEN sldist,dirs[i],halo[i],files[i],histr,histsf ELSE BEGIN
       rtipsy,dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.halo.' + halo[i],h,g,d,s
       units = tipsyunits(dirs[i] + files[i] + '.param')
       radius = sqrt(s.x*s.x + s.y*s.y)*units.lengthunit
       histsf = weighted_histogram(radius,locations = histr,max = 50,nbins = 500)  
;       plot,histr,histsf,psym = 10,title = files[i] + '.'  + halo[i] ,xrange = [0,5] ;,yrange = [0,1.3e-8]     
   ENDELSE
   nbins = n_elements(histsf)
   cumarr=fltarr(nbins,nbins)
   ii=lindgen(nbins*nbins)
   cumarr(where(ii mod nbins le ii/nbins)) = 1
;    sfh0cum=reform(sfh0 # cumarr)
   histsf=reform(histsf # cumarr)
   ii=0
   cumarr=0
   uniqind = uniq(histsf)
   IF histsf[0]/max(histsf) GT 0.5 THEN halfsfr[i] = histr[0] ELSE $
      halfsfr[i] = spline(histsf[uniqind]/max(histsf),histr[uniqind],0.5)
   
;   plot,histr[uniqind],histsf[uniqind]/max(histsf),title = files[i] + '.' + halo[i],xrange = [0,5],yrange = [0,1]
;   oplot,[halfsfr[i],halfsfr[i]],[0,1],linestyle = 2
;   oplot,[0,50],[0.5,0.5],linestyle = 2
;   stop
ENDFOR

formatplot,outplot = outplot
IF keyword_set(stellarmass) THEN xaxis = smass ELSE xaxis = vmass
IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'ejectr.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xaxis,halfr,/xlog,xrange = xrange,xtitle = xtitle,yrange = [0.04,5],ytitle = 'Half Mass Radius [kpc]',psym = symcat(symbols[1]),symsize = symsize,/nodata,/ylog;[0,5]
oplot,[1e9,1e12],[0.085,0.085],linestyle = 1
oplot,[1e9,1e12],[0.17,0.17],linestyle = 1
oplot,xaxis,halfr_reaccr,psym = symcat(symbols[2]),symsize = symsize,color = colore[2]
IF keyword_set(symbols_boarder) THEN oplot,xaxis,halfr_reaccr,psym = symcat(symbols_boarder[2]),symsize = symsize,color = fgcolor
oplot,xaxis,halfr,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
IF keyword_set(symbols_boarder) THEN oplot,xaxis,halfr,psym = symcat(symbols_boarder[1]),symsize = symsize,color = fgcolor,thick = 1
;oplot,xaxis,halfrexp,psym = symcat(symbols[0]),symsize =symsize,color = colore[0]
oplot,xaxis,halfsfr,psym = symcat(46),symsize = symsize,color = fgcolor
;legend,['Gas Expelled','Gas Ejected'],psym = symbols,color = colore,/bottom,/right
legend,['Gas Ejected','Gas Reaccreted','Star Formation'],psym = [symbols[1],symbols[2],46],color = [colore[1],colore[2],fgcolor],/top,/left,box =0
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'sfr.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xaxis,halfsfr,psym = symcat(46),/xlog,xrange = xrange,xtitle = xtitle,yrange = [0,5],ytitle = 'Half Radius of Star Formation [kpc]',symsize = symsize
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'eject_sfr_r.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
;ind = where(halfsfr LT 100)
;plot,smass[ind],halfr[ind]/halfsfr[ind],psym = 4,/xlog,xrange = xrange,xtitle = xtitle,yrange = [0,1],ytitle = 'Ratio Eject Half Radius to SF Half Radius'
plot,xaxis,halfr/halfsfr,psym = symcat(symbols[1]),/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Eject Half Radius/SF Half Radius',yrange = [0,4],symsize = symsize,/nodata
oplot,xrange,[1,1],thick = thicks[0]
oplot,xaxis,halfr/halfsfr,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
IF keyword_set(symbols_boarder) THEN oplot,xaxis,halfr/halfsfr,psym = symcat(symbols_boarder[1]),symsize = symsize,color = fgcolor,thick = 1
;oplot,xaxis,halfrexp/halfsfr,psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
;legend,['Gas Expelled','Gas Ejected'],psym = symbols,color =colore,/bottom,/right,box =0
IF keyword_set(outplot) THEN device, /close ELSE stop

;----------------------------- Multiplot --------------------------
IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'eject_sfr_multi.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*1.8,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize*1.8
multiplot,[1,2]

plot,xaxis,halfr,/xlog,xrange = xrange,yrange = [0.04,5],ytitle = 'Half Mass Radius [kpc]',psym = symcat(symbols[1]),symsize = symsize,/nodata,/ylog
oplot,[1e9,1e12],[0.085,0.085],linestyle = 1
oplot,[1e9,1e12],[0.17,0.17],linestyle = 1
oplot,xaxis,halfr_reaccr,psym = symcat(symbols[2]),symsize = symsize,color = colore[2]
IF keyword_set(symbols_boarder) THEN oplot,xaxis,halfr_reaccr,psym = symcat(symbols_boarder[2]),symsize = symsize,color = fgcolor
oplot,xaxis,halfr,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
IF keyword_set(symbols_boarder) THEN oplot,xaxis,halfr,psym = symcat(symbols_boarder[1]),symsize = symsize,color = fgcolor,thick = 1
;oplot,xaxis,halfrexp,psym = symcat(symbols[0]),symsize =symsize,color = colore[0]
oplot,xaxis,halfsfr,psym = symcat(46),symsize = symsize,color = fgcolor
;legend,['Gas Expelled','Gas Ejected'],psym = symbols,color = colore,/bottom,/right
legend,['Gas Ejected','Gas Reaccreted','Star Formation'],psym = [symbols[1],symbols[2],46],color = [colore[1],colore[2],fgcolor],/top,/left,box =0
multiplot

plot,xaxis,halfr/halfsfr,psym = symcat(symbols[1]),/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Eject Half Radius/SF Half Radius',yrange = [0.5,4.1],symsize = symsize,/nodata
oplot,xrange,[1,1],thick = thicks[0]
oplot,xaxis,halfr/halfsfr,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
IF keyword_set(symbols_boarder) THEN oplot,xaxis,halfr/halfsfr,psym = symcat(symbols_boarder[1]),symsize = symsize,color = fgcolor,thick = 1
;oplot,xaxis,halfrexp/halfsfr,psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
;oplot,xaxis,halfrexp/halfrexp,psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
;legend,['Gas Expelled','Gas Ejected'],psym = symbols,color=colore,/bottom,/right,box =0
multiplot,/reset
IF keyword_set(outplot) THEN device, /close ELSE stop


END
