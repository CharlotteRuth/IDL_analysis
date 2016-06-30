;filename = 'Disk_Collapse_1e6.00100.scalez0.1.H2.00020'
;filename = 'Disk_Collapse_1e6.00100.scalez1.H2.00020'
;pfile = '../Disk_Collapse_1e6.param'
;phaseD,filename,pfile = pfile,xrange = [1e-2,1e4], yrange =[10,1e4],/histograms,outplot = filename
PRO phaseD,filename,pfile = pfile,outplot = outplot,xrange = xrange, yrange = yrange,histograms = histograms
formatplot,outplot = outplot
IF NOT keyword_set(xrange) THEN xrange = [1e-10,1e4]
IF NOT keyword_set(yrange) THEN yrange = [1,1e8]


IF NOT keyword_set(pfile) THEN pfile = (findfile('*.param'))[0]
IF NOT keyword_set(pfile) THEN print,'No parameter file found'
units = tipsyunits(pfile)
rtipsy,filename,h,g,d,s
multiplot,[1,1]
IF (KEYWORD_SET(outplot)) THEN device,filename=outplot + '_phaseD.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2  ELSE window,0
plot,g.dens*units.rhounit,g.tempg,/ylog,/xlog,psym = 3,xrange = xrange, yrange = yrange,xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]'
IF (KEYWORD_SET(outplot)) THEN device,/close 

IF keyword_set(histograms) THEN BEGIN
IF (KEYWORD_SET(outplot)) THEN device,filename=outplot + '_phaseHist.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2  ELSE window,1
   multiplot,[2,1],xgap = 0.02,mTitle = filename
   histogramp,alog10(g.dens*units.rhounit),min = alog10(xrange[0]), max = alog10(xrange[1]), xtitle = 'log(Density [amu/cc])',nbins = 100
   multiplot
   histogramp,alog10(g.tempg),min = alog10(yrange[0]), max = alog10(yrange[1]), xtitle = 'log(Temperature [K])',nbins = 100
   multiplot,/reset
IF (KEYWORD_SET(outplot)) THEN device,/close  ELSE stop
ENDIF
END
