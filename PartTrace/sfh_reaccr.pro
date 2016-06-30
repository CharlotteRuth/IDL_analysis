PRO sfh_reaccr,dir,halo = halo,formatthick = formatthick,color = color,outplot = outplot

IF NOT keyword_set(halo) THEN halo = '1'

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
device,decomposed = 0

IF keyword_set(color) THEN BEGIN
   loadct,39
   colors = [fgcolor,50,254]
   linestyles = [0,0,0]
ENDIF ELSE BEGIN
   loadct,0
   colors = [fgcolor,fgcolor,fgcolor]
   linestyles = [0,1,4]
ENDELSE
IF keyword_set(formatthick) THEN thick = 3 ELSE thick = 1
spawn,'ls ' + dir + '/h*param',pfile
units = tipsyunits(pfile[0])
spawn,'ls ' + dir + '/*.starlog',slfile
sl = rstarlog(slfile[0],/molecularH)
sl.timeform = sl.timeform*units.timeunit/1e9
sl.massform = sl.massform*units.massunit
slcut = mrdfits(dir + 'starlog.cut.'+halo + '.fits',1)
match,sl.iorderstar,slcut.iorderstar,ind1,ind2
slgal = sl[ind1]

accr_iord = mrdfits(dir + '/grp' + halo + '.reaccrdisk_iord.fits',0)
accr_z = mrdfits(dir + '/grp' + halo + '.reaccrdisk_z.fits',0)
accr_t = z_to_t(accr_z)
accr_iord_r = reverse(accr_iord)
accr_t_r = reverse(accr_t)
reaccr_iord_r = reverse(accr_iord)
reaccr_t_r = reverse(accr_t)
hist = histogram(accr_iord,min = min(accr_iord),max = max(accr_iord),locations = iord)

IF 0 THEN BEGIN
;reaccr_iord = iord[where(hist GT 1)]
;match2,reaccr_iord,slgal.iordergas,ind1,ind2
;reaccr_iord = reaccr_iord[where(ind1 NE -1)]

FOR i = 0, n_elements(reaccr_iord) - 1 DO BEGIN
   temp = where(slgal.iordergas EQ reaccr_iord[i])
   IF temp[0] NE -1 THEN BEGIN
      slgas = slgal[temp]
      accr_ti = accr_t[where(accr_iord EQ reaccr_iord[i])]
      temp = where(slgas.timeform*units.timeunit/1e9 GT accr_ti[1])
      IF keyword_set(reaccr_sl) THEN reaccr_sl = [reaccr_sl,slgal[temp]] ELSE reaccr_sl = [slgal[temp]]
   ENDIF ELSE print,'Problem'
ENDFOR
ENDIF

initaccr_iord_r = accr_iord_r[uniq(accr_iord_r,sort(accr_iord_r))] 
initaccr_t_r = accr_t_r[uniq(accr_iord_r,sort(accr_iord_r))] 

reaccr_iord_r[uniq(accr_iord_r,sort(accr_iord_r))] = -1
reaccr_t_r[uniq(accr_iord_r,sort(accr_iord_r))] = -1
reaccr_iord_r = reaccr_iord_r[where(reaccr_iord_r NE -1)]
reaccr_t_r = reaccr_t_r[where(reaccr_t_r NE -1)]

hist0 = histogram(accr_t,min = 1e-6,max = 14,nbins = 100,locations = xtime)
hist0cum = total(hist0,/cumulative)
hist1 = histogram(initaccr_t_r,min = 1e-6,max = 14,nbins = 100)
hist1cum = total(hist1,/cumulative)
hist2 = histogram(reaccr_t_r,min = 1e-6,max = 14,nbins = 100)
hist2cum = total(hist2,/cumulative)

IF keyword_set(outplot) THEN device,filename = outplot + '_sfh_accr.eps',/color,bits_per_pixel= 8 ELSE window,0
multiplot,[1,2]
plot,xtime,hist0cum,xrange = [0,14],/nodata
oplot,xtime,hist0cum,color = colors[0],linestyle = linestyles[0],psym = 10
oplot,xtime,hist1cum,color = colors[1],linestyle = linestyles[1],psym = 10
oplot,xtime,hist2cum,color = colors[2],linestyle = linestyles[2],psym = 10
multiplot
plot,xtime,hist2cum/hist0cum,xrange = [0,14],xtitle = 'Time [Gyr]',psym = 10
multiplot,/reset
IF keyword_set(outplot) THEN device,/close ELSE stop

match2,reaccr_iord_r,slgal.iordergas,ind1,ind2
slgal_sub = slgal[where(ind2 NE -1)]

tarr = reaccr_t_r[uniq(reaccr_t_r,sort(reaccr_t_r))]
tarr = [tarr,14]
FOR it = 0, n_elements(tarr) - 2 DO BEGIN
   slgal_sub_tbin = slgal_sub[where(slgal_sub.timeform GE tarr[it] AND slgal_sub.timeform LT tarr[it + 1])]
   ind_reaccr_r_sub = where(reaccr_t_r LE tarr[it])
   match2,slgal_sub_tbin.iordergas,reaccr_iord_r[ind_reaccr_r_sub],ind1,ind2
   IF keyword_set(reaccr_sl) THEN reaccr_sl = [reaccr_sl,slgal_sub_tbin[where(ind1 NE -1)]] ELSE reaccr_sl = slgal_sub_tbin[where(ind1 NE -1)]
ENDFOR

nbins = 100
max = 14.
binsize = nbins/max
hist0_sf = weighted_histogram(slgal.timeform,weight = slgal.massform/1e9/binsize,max = max,min = 1e-6,nbins = nbins,locations = xtime)
hist0cum_sf = total(hist0_sf,/cumulative)*binsize*1e9
hist2_sf = weighted_histogram(reaccr_sl.timeform,weight = reaccr_sl.massform/1e9/binsize,max = max,min = 1e-6,nbins = nbins,locations = xtime)
hist2cum_sf = total(hist2_sf,/cumulative)*binsize*1e9
hist1_sf = hist0_sf - hist2_sf
hist1cum_sf = hist0cum_sf - hist2cum_sf

IF keyword_set(outplot) THEN device,filename = outplot + '_sfh_accr.eps',/color,bits_per_pixel= 8 ELSE window,0
multiplot,[1,2]
plot,xtime,hist0_sf,ytitle = 'Star Formation Rate [M' + sunsymbol() + textoidl('yr^{-1}]'),psym = 10,xrange = [0,14],/nodata
oplot,xtime,hist0_sf,psym = 10,color = colors[0],linestyle = linestyles[0]
oplot,xtime,hist1_sf,psym = 10,color = colors[1],linestyle = linestyles[1]
oplot,xtime,hist2_sf,psym = 10,color = colors[2],linestyle = linestyles[2]
multiplot
plot,xtime,hist2_sf/hist0_sf,ytitle = 'Accreted/Total Star Formation',xtitle = 'Time [Gyr]',xrange = [0,14],psym = 10
multiplot,/reset
IF keyword_set(outplot) THEN device,/close ELSE stop

IF keyword_set(outplot) THEN device,filename = outplot + '_sfcum_accr.eps',/color,bits_per_pixel= 8 ELSE window,0
multiplot,[1,2]
plot,xtime,hist0cum_sf,ytitle = 'Star Formation Rate [M' + sunsymbol() + textoidl('yr^{-1}]'),psym = 10,xrange = [0,14],/nodata
oplot,xtime,hist0cum_sf,psym = 10,color = colors[0],linestyle = linestyles[0]
oplot,xtime,hist1cum_sf,psym = 10,color = colors[1],linestyle = linestyles[1]
oplot,xtime,hist2cum_sf,psym = 10,color = colors[2],linestyle = linestyles[2]
multiplot
plot,xtime,hist2cum_sf/hist0cum_sf,ytitle = 'Accreted/Total Star Formation',xtitle = 'Time [Gyr]',xrange = [0,14],psym = 10
multiplot,/reset
IF keyword_set(outplot) THEN device,/close ELSE stop

END
