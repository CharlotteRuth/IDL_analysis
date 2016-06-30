PRO plot_angmom,dirs,files,halo = halo,outplot = outplot,colors = colors,finalstep = finalstep,stellarmass = stellarmass,symsize = symsize, symbols = symbols,formatthick = formatthick,sfr = sfr
IF NOT keyword_set(finalstep) THEN finalstep = '00512'
n = n_elements(dirs)
formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
    IF keyword_set(stellarmass) THEN outext = '_sm_' ELSE outext = '_vm_'
;    IF keyword_set(sfr) THEN outext = outext + 'orig_' ELSE outext = outext + 'final_'
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 800
    ysize = 500
 ENDELSE
IF keyword_set(stellarmass) THEN BEGIN
   xtitle = 'Stellar Mass [M' + sunsymbol() + ']'
   xrange = [1e7,1e11]
ENDIF ELSE BEGIN
   xtitle = 'Virial Mass [M' + sunsymbol() + ']'
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


vmass = fltarr(n)
smass = fltarr(n)
ejectL = fltarr(n)
eject2L = fltarr(n)
accrL = fltarr(n)
ejectL_med = fltarr(n)
eject2L_med = fltarr(n)
accrL_med = fltarr(n)
ejectL_stdev = fltarr(n)
eject2L_stdev = fltarr(n)
accrL_stdev = fltarr(n)
diffL = fltarr(n)
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   print,files[i] + '.' + finalstep + '.' + halo[i]
   stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group eq halo[i]))[0]
   vmass[i] = stat[ind].m_tot
   smass[i] = stat[ind].m_star

   history_eject =  mrdfits(dirs[i] + 'grp' + halo[i] + '.reeject_disk.fits',1)
   history_accr = mrdfits(dirs[i] + 'grp' + halo[i] + '.reaccrdisk_history.fits',1)
   history_eject2 = eject_accr_match(dirs[i],halo[i])
   IF 0 THEN BEGIN
      ind_history_eject_r = reverse(indgen(n_elements(history_eject)))
      history_eject_r = history_eject[ind_history_eject_r]
      ind_history_accr_r = reverse(indgen(n_elements(history_accr)))
      history_accr_r = history_accr(ind_history_accr_r)

      history_eject2 = history_accr
      accriord = history_accr[uniq(history_accr.iord,sort(history_accr.iord))].iord
      FOR j = 0, n_elements(accriord) - 1 DO BEGIN
         indaccr = where(history_accr.iord EQ accriord[j])
         indeject = where(history_eject.iord EQ accriord[j])
         history_eject2[indaccr] = history_eject[indeject[0:n_elements(indaccr) - 1]]
      ENDFOR
   ENDIF

   r_eject = [[history_eject.x*history_eject.mass],[history_eject.y*history_eject.mass],[history_eject.z*history_eject.mass]]
   r_eject2 = [[history_eject2.x*history_eject2.mass],[history_eject2.y*history_eject2.mass],[history_eject2.z*history_eject2.mass]]
   r_accr  = [[history_accr.x*history_accr.mass],[history_accr.y*history_accr.mass],[history_accr.z*history_accr.mass]]
   v_eject = [[history_eject.vx],[history_eject.vy],[history_eject.vz]]/3e16
   v_eject2 = [[history_eject2.vx],[history_eject2.vy],[history_eject2.vz]]/3e16
   v_accr  = [[history_accr.vx],[history_accr.vy],[history_accr.vz]]/3e16

   j_eject_v = crossp_multi(r_eject, v_eject)
   j_eject2_v = crossp_multi(r_eject2, v_eject2)
   j_accr_v = crossp_multi(r_accr, v_accr)
   j_eject = sqrt(j_eject_v[*,0]^2 + j_eject_v[*,1]^2 + j_eject_v[*,2]^2)
   j_eject2 = sqrt(j_eject2_v[*,0]^2 + j_eject2_v[*,1]^2 + j_eject2_v[*,2]^2)
   j_accr = sqrt(j_accr_v[*,0]^2 + j_accr_v[*,1]^2 + j_accr_v[*,2]^2)
   ejectL[i] = total(j_eject)
   eject2L[i] = total(j_eject2)
   accrL[i] = total(j_accr)
   ejectL_med[i] = median(j_eject)
   eject2L_med[i] = median(j_eject2)
   accrL_med[i] = median(j_accr)
   ejectL_stdev[i] = stdev(j_eject)
   eject2L_stdev[i] = stdev(j_eject2)
   accrL_stdev[i] = stdev(j_accr)
   diffL[i] = total(j_accr) - total(j_eject2)
   y = histogram(j_eject,nbins = 100);,min = 0,max = 1e6)
   IF NOT keyword_set(outplot) THEN window,2
   histogramp,j_eject,nbins = 100,yrange = [0,max(y)]
   oplot,[ejectL_med[i],ejectL_med[i]],[0,1e7],color = fgcolor
   histogramp,j_eject2,nbins = 100,color = 100,/overplot
   oplot,[eject2L_med[i],eject2L_med[i]],[0,1e7],color = 100
   histogramp,j_accr,nbins = 100,color = 200,/overplot
   oplot,[accrL_med[i],accrL_med[i]],[0,1e7],color = 200
   print,mean(j_eject),mean(j_eject2),mean(j_accr),median(j_eject),median(j_eject2),median(j_accr),total(j_eject),total(j_eject2),total(j_accr),total(j_accr) - total(j_eject2),total(j_accr)/total(j_eject2)
   legend,['All Ejecta','Ejecta to be reaccreted','Reaccretion'],color = [fgcolor,100,200],linestyle = [0,0,0],/right,box = 0
ENDFOR
stop
formatplot,outplot = outplot
IF keyword_set(stellarmass) THEN xaxis = smass ELSE xaxis = vmass

;--------------------------- Reeject and Reaccretion
IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'ejectL_accr.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xaxis,ejectL,/xlog,xrange = xrange,xtitle = xtitle,ytitle = textoidl('L [M') + sunsymbol()+textoidl(' kpc^2 s^{-1}]'),psym = symcat(symbols[1]),symsize = symsize,/nodata,/ylog,yrange = [1e-9,1e-2]
oplot,xaxis,ejectL,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
oplot,xaxis,ejectL,psym = symcat(symbols_boarder[1]),symsize = symsize,color = fgcolor,thick = 1
oplot,xaxis,accrL,psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
oplot,xaxis,accrL,psym = symcat(symbols_boarder[0]),symsize = symsize,color = fgcolor
legend,['Origianl','Reaccreted'],psym = [symbols[1], symbols[0]],color = [colore[1],colore[0]],/bottom,/right,box =0
IF keyword_set(outplot) THEN device, /close ELSE stop

;--------------- Fractional Total Angular Momentum change ------------
IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'ejectL_accr_frac.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xaxis,accrL/eject2L,/xlog,xrange = xrange,xtitle = xtitle,ytitle = textoidl('L_{accr}/L_{eject}'),psym = symcat(16),symsize = symsize,yrange = [0.6,2]
IF keyword_set(outplot) THEN device, /close ELSE stop

;----------------------------- Fractional Median Angular Momentum Change --------
IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'ejectr_accr_fracmed.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xaxis,accrL_med/eject2L_med,/xlog,xrange = xrange,xtitle = xtitle,ytitle = textoidl('L_{accr}/L_{eject}'),psym = symcat(16),symsize = symsize,yrange = [1,4]
IF keyword_set(outplot) THEN device, /close ELSE stop

;----------------------------- Multiplot --------------------------
IF keyword_set(outplot) THEN  device,filename = outplot + outext + 'ejectL_multi.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*1.8,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize*1.8
multiplot,[1,2]

plot,xaxis,ejectL_med,/xlog,xrange = xrange,yrange = [1e-13,1e-9],ytitle = textoidl('L [M') + sunsymbol()+textoidl(' kpc^2 s^{-1}]'),psym = symcat(symbols[1]),symsize = symsize,/nodata,/ylog
;oploterror,xaxis,eject2L_med,eject2L_stdev,psym = 3,color =fgcolor,errcolor = colore[1],thick = 4
oplot,xaxis,ejectL_med,psym = symcat(symbols_boarder[1]),symsize = symsize,color = colore[1]
oplot,xaxis,eject2L_med,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
oplot,xaxis,eject2L_med,psym = symcat(symbols_boarder[1]),symsize = symsize,color = fgcolor,thick = 1
;oploterror,xaxis,accrL_med,accrL_stdev,psym = 3,color = fgcolor,errcolor = colore[0],thick = 4
oplot,xaxis,accrL_med,psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
oplot,xaxis,accrL_med,psym = symcat(symbols_boarder[0]),symsize = symsize,color = fgcolor
legend,['Ejected','Ejected (later reaccreted)','Reaccreted'],psym = [symbols_boarder[1],symbols[1], symbols[0]],color = [colore[1],colore[1],colore[0]],/bottom,/right,box =0
multiplot

plot,xaxis,accrL_med/eject2L_med,psym = symcat(symbols[1]),/xlog,xrange = xrange,xtitle = xtitle,ytitle = textoidl('L_{reaccreted}/L_{ejected}'),symsize = symsize,/nodata,yrange = [1,4.5]
;oplot,xrange,[1,1],thick = thicks[0]
oplot,xaxis,accrL_med/ejectL_med,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
oplot,xaxis,accrL_med/ejectL_med,psym = symcat(symbols_boarder[1]),symsize = symsize,color =fgcolor
multiplot,/reset
IF keyword_set(outplot) THEN device, /close ELSE stop


END
 
