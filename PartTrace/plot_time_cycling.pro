;Plots the output from time_cycling

PRO plot_time_cycling,dirs,files,halo = halo,finalstep = finalstep,colors = colors,outplot = outplot,symbols = symbols,expelled = expelled
n = n_elements(dirs)
IF NOT keyword_set(halo) THEN halo = strarr(n) + '1'
IF NOT keyword_set(finalstep) THEN finalstep = '00512'

formatplot,outplot = outplot,thick = formatthick

IF keyword_set(expelled) THEN BEGIN
   readfile = 'expelltime' 
   outfile = 'expelled'
ENDIF ELSE BEGIN 
   readfile = 'outflowtime'
   outfile = 'ejected'
ENDELSE

xmrange = [3e9,1e12]
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
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    colors_all = (findgen(n) + 1)*240/n
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symbols) THEN symbols = fltarr(n) + 14
    symsize = 2
ENDIF ELSE BEGIN
    loadct,0    
    colors = [fgcolor,fgcolor];  fltarr(n) + 5
    colors_all = (findgen(n) + 1)*fgcolor
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2) 
    IF NOT keyword_set(symbols) THEN symbols = fltarr(n) + 14 
    symsize = 2
ENDELSE

smass = fltarr(n)
vmass = fltarr(n)
FOR i = 0, n - 1 DO BEGIN
    stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
    ind = (where(stat.group eq halo[i]))[0]
    vmass[i] = stat[ind].m_tot
    smass[i] = stat[ind].m_star
 ENDFOR
sortind = sort(vmass)
vmass = vmass[sortind]
smass = smass[sortind]
dirs  = dirs[sortind]
files = files[sortind]
halo = halo[sortind]
colors_all = (alog10(vmass) - 9.5)/2.5*256
mean_outflowtime = fltarr(n)
median_outflowtime = fltarr(n)
mode_outflowtime = fltarr(n)
stdv_outflowtime = fltarr(n)

mintime = 0.1
dtime = 0.1
nbins = 50
maxtime = 14;dtime*nbins

outflowtime_all = [0]
outflowmass_all = [0]
outtime_all = [0]
intime_all = [0]

formatplot,outplot = outplot
IF keyword_set(outplot) THEN device,filename = outplot + '_time' + outfile + 'cum.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2  ELSE window, 2, xsize = xsize, ysize = ysize

FOR i = 0, n - 1 DO BEGIN
   readcol,dirs[i] + '/' + readfile + '.grp' + halo[i] + '.txt',outflowtime,outflowmass
   outflowtime[where(outflowtime LT mintime)] = mintime
   readcol,dirs[i] + '/inouttime.grp' + halo[i] + '.txt',outtime,intime
   outflowtime_all = [outflowtime_all,outflowtime]
   outflowmass_all = [outflowmass_all,outflowmass]
   outtime_all = [outtime_all,outtime]
   intime_all = [intime_all,intime]
   IF i EQ 0 THEN histogramp,alog10(outflowtime),weight = outflowmass,min = alog10(mintime), max = alog10(maxtime),nbins = nbins,/normalize,xtitle = 'Log Time Out Of Disk [Gyr]',/ylog,yrange = [0.005,1],/cum,ytitle = textoidl('n < time_{recycling}') ;,yrange = [1e-3,0.45],/ylog;,xrange = [0,5];,ytitle = 'Cumulative Histogram'
    histogramp,alog10(outflowtime),weight = outflowmass,min = alog10(mintime), max = alog10(maxtime),nbins = nbins,color = colors_all[i],/normalize,/overplot,/cum
;    stop
    yarray = weighted_histogram(outflowtime,weight = outflowmass,min = 0.05, max = dtime*nbins,nbins = nbins,locations = xarray)
    temp = max(yarray,maxind)
    mode_outflowtime[i] = xarray[maxind]
    mean_outflowtime[i] = mean(outflowtime)
    median_outflowtime[i] = median(outflowtime)
    stdv_outflowtime[i] = stdev(outflowtime)

    outflowtime2 = fix(outflowtime*1000)/1000.0
    outflowtime2 = outflowtime2[uniq(outflowtime2,sort(outflowtime2))]
;    stop
 ENDFOR
cgColorbar,range=[9.5,12],/vertical,position=[0.90, 0.3, 0.92, 0.9],title = 'Log Virial Mass',divisions = 5,tcharsize = 1.5
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN device,filename = outplot + '_time' + outfile + 'cum2.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2  ELSE window, 2, xsize = xsize, ysize = ysize
FOR i = 0, n - 1 DO BEGIN
   IF keyword_set(expelled) THEN outflowm_all = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_iord.fits',0) ELSE outflowm_all = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reeject.fits',0)
   readcol,dirs[i] + '/' + readfile + '.grp' + halo[i] + '.txt',outflowtime,outflowmass
    outflowtime[where(outflowtime LT mintime)] = mintime
    IF i EQ 0 THEN histogramp,alog10(outflowtime),weight = outflowmass,min = alog10(mintime), max = alog10(maxtime),nbins = nbins,normalize = total(outflowm_all),xtitle = 'Log Time Out Of Disk [Gyr]',/ylog,yrange = [0.005,1],/cum,ytitle = textoidl('n < time_{recycling}');,yrange = [1e-3,0.45],/ylog;,xrange = [0,5];,ytitle = 'Cumulative Histogram'
    histogramp,alog10(outflowtime),weight = outflowmass,min = alog10(mintime), max = alog10(maxtime),nbins = nbins,color = colors_all[i],normalize = total(outflowm_all),/overplot,/cum
;    stop
    yarray = weighted_histogram(outflowtime,weight = outflowmass,min = 0.05, max = dtime*nbins,nbins = nbins,locations = xarray)
    temp = max(yarray,maxind)
    mode_outflowtime[i] = xarray[maxind]
    mean_outflowtime[i] = mean(outflowtime)
    median_outflowtime[i] = median(outflowtime)
    stdv_outflowtime[i] = stdev(outflowtime)

    outflowtime2 = fix(outflowtime*1000)/1000.0
    outflowtime2 = outflowtime2[uniq(outflowtime2,sort(outflowtime2))]
;    stop
 ENDFOR
cgColorbar,range=[9.5,12],position=[0.30, 0.32, 0.9, 0.34],title = 'Log Virial Mass',divisions = 5,/top,tcharsize = 1.5
IF keyword_set(outplot) THEN device, /close ELSE stop

dtime = 0.1069;*2
nbins = 130;/2
IF keyword_set(outplot) THEN  device,filename = outplot + '_time' + outfile + 'log.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
FOR i = 0, n - 1 DO BEGIN
    readcol,dirs[i] + '/' + readfile + '.grp' + halo[i] + '.txt',outflowtime,outflowmass
    IF i EQ 0 THEN histogramp,outflowtime,weight = outflowmass,min = mintime, max = dtime*nbins,nbins = nbins,/normalize,xtitle = 'Time in Outflow [Gyr]',ytitle = '1/M dM/dt',yrange = [1e-3,0.4],/ylog,/xlog,xrange = [mintime,30],xstyle = 1;,yrange = [0,1],/cum,;,ytitle = 'Cumulative Histogram'
    histogramp,               outflowtime,weight = outflowmass,min = mintime, max = dtime*nbins,nbins = nbins,/normalize,color = colors_all[i],/overplot;,/cum
;    readcol,dirs[i] + '/expelltime.grp' + halo[i] + '.txt',exptime,expmass    
;    histogramp,exptime,weight = expmass,min = 0.05, max = dtime*nbins,nbins = nbins,color = colors[i],/normalize,/overplot,/cum,linestyle = 2
 ENDFOR
cgColorbar,range=[9.5,12],/vertical,position=[0.90, 0.3, 0.92, 0.9],title = 'Log Virial Mass',divisions = 5
IF keyword_set(outplot) THEN device, /close ELSE stop

dtime = 0.1069;*2
nbins = 130;/2
IF keyword_set(outplot) THEN  device,filename = outplot + '_time' + outfile + '.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
FOR i = 0, n - 1 DO BEGIN
    readcol,dirs[i] + '/' + readfile + '.grp' + halo[i] + '.txt',outflowtime,outflowmass
    IF i EQ 0 THEN histogramp,outflowtime,weight = outflowmass,min = mintime, max = dtime*nbins,nbins = nbins,/normalize,xtitle = 'Time in Outflow [Gyr]',ytitle = '1/M dM/dt',xrange = [0,4],yrange = [1e-3,0.4],/ylog;,yrange = [0,1],/cum,;,ytitle = 'Cumulative Histogram'
    histogramp,               outflowtime,weight = outflowmass,min = mintime, max = dtime*nbins,nbins = nbins,/normalize,color = colors_all[i],/overplot;,/cum
;    readcol,dirs[i] + '/expelltime.grp' + halo[i] + '.txt',exptime,expmass    
;    histogramp,exptime,weight = expmass,min = 0.05, max = dtime*nbins,nbins = nbins,color = colors[i],/normalize,/overplot,/cum,linestyle = 2
 ENDFOR
cgColorbar,range=[9.5,12],position=[0.2, 0.9, 0.9, 0.92],title = 'Log Virial Mass',divisions = 5,TCHARSIZE = 1.3
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_' + outfile + '_timesmean.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,vmass,mean_outflowtime,xrange = xmrange,yrange = [0,5.5],xtitle = 'Virial Mass [M' + sunsymbol() + ']',ytitle = 'Average Time in Outflow [Gyr]',symsize = symsize,/nodata,/xlog
halomass = 10^(findgen(200)/200*7 + 7)
timeaccrOD = 1.58489d16*halomass^(-0.6)/1d9 ;Oppenheimer and Dave 2008, in Gyr
timeaccrH = 1.8e10*1e10/halomass/1e9 ;Henriques 2013
timeaccrM = 0.52e9*(1 + 0)^(-0.32)*(halomass/1e12)^(-0.45)/1d9
readcol,'~/Datafiles/Oppenheimer10.txt',mhlow,mhhi,trec,nwind,frac
loadct,0
oplot,10^((mhlow+mhhi)/2),trec/1e9,color = 200,linestyle = 3
oplot,halomass,timeaccrH,color = 200,linestyle = 2
oplot,halomass,timeaccrM,color = 200,linestyle = 0
legend,['Oppenheimer & Dave, 2010','Henriques et al., 2013','Mitra et al. 2014'],linestyle = [3,2,0],color = [200,200,200],box = 0,/right
loadct,39
oploterror,vmass,mean_outflowtime,stdv_outflowtime,psym = 3,color = fgcolor,errcolor = fgcolor, thick = 4
oplot,vmass,mean_outflowtime,psym = symcat(symbols[1]),color = colors[1],symsize = symsize, thick = 4
oplot,vmass,mean_outflowtime,psym = symcat(sym_outline(symbols[1])),symsize = symsize, thick = 4
print,'Mean/Median mean time [Gyr]: ',mean(mean_outflowtime),median(mean_outflowtime)
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_' + outfile + '_timesmode.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,vmass,mode_outflowtime,xrange = xmrange,yrange = [0,4],xtitle = 'Virial Mass [M' + sunsymbol() + ']',ytitle = 'Average Time in Outflow [Gyr]',symsize = symsize,/nodata,/xlog
oploterror,vmass,mode_outflowtime,stdv_outflowtime,psym = 3,color = fgcolor,errcolor = fgcolor, thick = 4
oplot,vmass,mode_outflowtime,psym = symcat(symbols[1]),color = colors[1],symsize = symsize, thick = 4
oplot,vmass,mode_outflowtime,psym = symcat(sym_outline(symbols[1])),symsize = symsize, thick = 4
print,'Mean/Median mode time [Gyr]: ',mean(mode_outflowtime),median(mode_outflowtime)
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_' + outfile + '_timesmed.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,vmass,median_outflowtime,xrange = xmrange,yrange = [0,4],xtitle = 'Virial Mass [M' + sunsymbol() + ']',ytitle = 'Median Time in Outflow [Gyr]',symsize = symsize,/nodata,/xlog
halomass = 10^(findgen(200)/200*7 + 7)
;timeaccrOD = 1.58489d16*halomass^(-0.6)/1d9 ;Oppenheimer and Dave 2008, in Gyr
timeaccrH = 1.8e10*1e10/halomass/1e9 ;Henriques 2013
timeaccrM = 0.52e9*(1 + 0)^(-0.32)*(halomass/1e12)^(-0.45)/1d9
readcol,'~/Datafiles/Oppenheimer10.txt',mhlow,mhhi,trec,nwind,frac
loadct,0
oplot,10^((mhlow+mhhi)/2),trec/1e9,color = 100,linestyle = 3,thick = 2
oplot,halomass,timeaccrH,color = 100,linestyle = 2,thick = 2
oplot,halomass,timeaccrM,color = 100,linestyle = 0,thick = 2
legend,['Oppenheimer & Dave, 2010','Henriques et al., 2013','Mitra et al. 2014'],linestyle = [3,2,0],color = [100,100,100],box = 0,/right,thick = [2,2,2]
loadct,39
oploterror,vmass,median_outflowtime,stdv_outflowtime,psym = 3,color = fgcolor,errcolor = fgcolor, thick = 4
oplot,vmass,median_outflowtime,psym = symcat(symbols[1]),color = colors[1],symsize = symsize, thick = 4
oplot,vmass,median_outflowtime,psym = symcat(sym_outline(symbols[1])),symsize = symsize, thick = 4
print,'Mean/Median median time [Gyr]: ',mean(median_outflowtime),median(median_outflowtime)
IF keyword_set(outplot) THEN device, /close ELSE stop

IF 1 THEN BEGIN
;------------------------------ Time in Disk ------------------------------------
mean_disktime = fltarr(n)
median_disktime = fltarr(n)
mode_disktime = fltarr(n)
stdv_disktime = fltarr(n)
IF keyword_set(outplot) THEN  device,filename = outplot + '_timedisk.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
FOR i = 0, n - 1 DO BEGIN
    readcol,dirs[i] + '/disktime.grp' + halo[i] + '.txt',disktime,/silent
    IF i EQ 0 THEN histogramp,disktime,min = 0.05, max = dtime*nbins,nbins = nbins,/normalize,yrange = [1e-4,0.35],xtitle = 'Time in Disk [Gyr]',/ylog,xrange = [0,12];,yrange = [1e-4,1],/ylog
    histogramp,disktime,min = 0.05, max = dtime*nbins,nbins = nbins,color = colors_all[i],/normalize,/overplot
    mean_disktime[i] = mean(disktime)
    median_disktime[i] = median(disktime)
    stdv_disktime[i] = stdev(disktime)
ENDFOR
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_timediskcum.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
FOR i = 0, n - 1 DO BEGIN
    readcol,dirs[i] + '/disktime.grp' + halo[i] + '.txt',disktime,/silent
    IF i EQ 0 THEN histogramp,disktime,min = 0.05, max = dtime*nbins,nbins = nbins,/normalize,yrange = [0,1],xtitle = 'Time in Disk [Gyr]',/cum;,yrange = [1e-4,1],/ylog
    histogramp,disktime,min = 0.05, max = dtime*nbins,nbins = nbins,color = colors_all[i],/normalize,/overplot,/cum
ENDFOR
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_disktimesmean.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,vmass,mean_disktime,xrange = xmrange,yrange = [0,6],xtitle = 'Virial Mass [M' + sunsymbol() + ']',ytitle = 'Average Time in Disk [Gyr]',symsize = symsize,/nodata,/xlog
oploterror,vmass,mean_disktime,stdv_disktime,psym = 3,color = fgcolor,errcolor = fgcolor, thick = 4
oplot,vmass,mean_disktime,psym = symcat(symbols[0]),color = 254,symsize = symsize, thick = 4
IF keyword_set(outplot) THEN device, /close ELSE stop
print,'Mean/Median mean time [Gyr]: ',mean(mean_disktime),median(mean_disktime)

IF keyword_set(outplot) THEN  device,filename = outplot + '_disktimesmed.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,vmass,median_disktime,xrange = xmrange,yrange = [0,5],xtitle = 'Virial Mass [M' + sunsymbol() + ']',ytitle = 'Median Time in Disk [Gyr]',symsize = symsize,/nodata,/xlog
oploterror,vmass,median_disktime,stdv_disktime,psym = 3,color = fgcolor,errcolor = fgcolor, thick = 4
oplot,vmass,median_disktime,psym = symcat(symbols[0]),color = 254,symsize = symsize, thick = 4
IF keyword_set(outplot) THEN device, /close ELSE stop
print,'Mean/Median mean time [Gyr]: ',mean(median_disktime),median(median_disktime)

;------------------------------------- Time expelled -----------------------------------------
mean_exptime = fltarr(n)
median_exptime = fltarr(n)
mode_exptime = fltarr(n)
stdv_exptime = fltarr(n)
IF keyword_set(outplot) THEN  device,filename = outplot + '_timeexp.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
FOR i = 0, n - 1 DO BEGIN
    readcol,dirs[i] + '/expelltime.grp' + halo[i] + '.txt',exptime,expmass,/silent
    IF i EQ 0 THEN histogramp,exptime,weight = expmass,min = 0.05, max = dtime*nbins,nbins = nbins,/normalize,yrange = [1e-4,0.35],xtitle = 'Time in Exp [Gyr]',/ylog;,yrange = [1e-4,1],/ylog
    histogramp,exptime,weight = expmass,min = 0.05, max = dtime*nbins,nbins = nbins,color = colors_all[i],/normalize,/overplot
    mean_exptime[i] = mean(exptime)
    median_exptime[i] = median(exptime)
    stdv_exptime[i] = stdev(exptime)
ENDFOR
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_timeexpcum.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
FOR i = 0, n - 1 DO BEGIN
    readcol,dirs[i] + '/expelltime.grp' + halo[i] + '.txt',exptime,expmass,/silent
    IF i EQ 0 THEN histogramp,exptime,weight = expmass,min = 0.05, max = dtime*nbins,nbins = nbins,/normalize,yrange = [0,1],xtitle = 'Time in Expelled [Gyr]',/cum;,yrange = [1e-4,1],/ylog
    histogramp,exptime,weight = expmass,min = 0.05, max = dtime*nbins,nbins = nbins,color = colors_all[i],/normalize,/overplot,/cum
ENDFOR
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_exptimesmean.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,vmass,mean_exptime,xrange = xmrange,yrange = [0,9],xtitle = 'Virial Mass [M' + sunsymbol() + ']',ytitle = 'Average Time Expelled [Gyr]',symsize = symsize,/nodata,/xlog
oploterror,vmass,mean_exptime,stdv_exptime,psym = 3,color = fgcolor,errcolor = fgcolor, thick = 4
oplot,vmass,mean_exptime,psym = symcat(symbols[0]),color = colors[0],symsize = symsize, thick = 4
oplot,vmass,mean_exptime,psym = symcat(sym_outline(symbols[0])),symsize = symsize, thick = 4
IF keyword_set(outplot) THEN device, /close ELSE stop
print,'Mean/Median mean time [Gyr]: ',mean(mean_exptime),median(mean_exptime)

IF keyword_set(outplot) THEN  device,filename = outplot + '_exptimesmed.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,vmass,median_exptime,xrange = xmrange,yrange = [0,9.5],xtitle = 'Virial Mass [M' + sunsymbol() + ']',ytitle = 'Median Time Expelled [Gyr]',symsize = symsize,/nodata,/xlog
oploterror,vmass,median_exptime,stdv_exptime,psym = 3,color = fgcolor,errcolor = fgcolor, thick = 4
oplot,vmass,median_exptime,psym = symcat(symbols[0]),color = colors[0],symsize = symsize, thick = 4
oplot,vmass,median_exptime,psym = symcat(sym_outline(symbols[0])),symsize = symsize, thick = 4
IF keyword_set(outplot) THEN device, /close ELSE stop
print,'Mean/Median mean time [Gyr]: ',mean(median_exptime),median(median_exptime)
ENDIF
END
