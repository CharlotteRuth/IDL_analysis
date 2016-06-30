PRO bulgeSFH,files,munit,lunit,tunit,key,outplot = outplot
loadct,0
yrange1 = [0,2.2]
yout1 = 1.75
;ytickname1 = ['0','5e7','1e8','1.5e8','2e8','2.5e8']
yrange2 = [0,2.2]
yout2 = 1.75
;ytickname2 = ['0','5e7','1e8','1.5e8','2e8','2.5e8']
IF keyword_set(outplot) THEN device,filename=outplot + '_bulgeSFH.eps', /color,xsize = 15, ysize = 12,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 600,ysize = 600
n = n_elements(files); - 1
multiplot,[1,n],mytitle = 'Star Formation Rate [M' + sunsymbol() + textoidl('yr^{-1}]'),mytitoffset = 2.5,mytitsize = 1.5
IF keyword_set(outplot) THEN thick = 4 ELSE thick = 1
FOR i = 0, n - 1 DO BEGIN
    IF i EQ 0 OR i EQ 1 THEN BEGIN
        yrange = yrange1
        yout = yout1
;        ytickname = ytickname1
    ENDIF ELSE BEGIN
        yrange = yrange2
        yout = yout2
;        ytickname = ytickname2
    ENDELSE
    rtipsy,files[i],h,g,d,s
    readarr,files[i] + '.decomp',h,decomp,part = 'star',/ascii
    bulge = where(decomp EQ 5 OR decomp EQ 3)
    disk = where(decomp EQ 1)
    s.mass = s.mass*munit[i]
    s.tform = s.tform*tunit[i]/1e9
    IF i EQ n - 1 THEN xtitle = 'Time [Gyr]' ELSE xtitle = ''
    tbinsize = 14e9/100.0
    histogramp,s.tform,weight = s.mass/tbinsize,min = 1e-6,max = 14,nbins = 100,xtitle = xtitle,yrange = yrange,thick = thick,ytickname = ytickname
;    histogramp,s.tform,weight = s.mass/tbinsize,min = 1e-6,max = 14,nbins = 100,xtitle = xtitle,yrange = [0,1],thick = thick,ytickname = ytickname,/cum,/normalize
;    oplot,[0,14],[0.5,0.5]
    ycum = weighted_histogram(s.tform,weight = s.mass/tbinsize,min = 1e-6,max = 14,nbins = 100,locations = xcum,/cum)
    ycum = ycum/max(ycum)
    asst = spline(ycum[uniq(ycum)],xcum[uniq(xcum)],[0.5])
    print,'Half mass assembled: ',asst,z_from_time(asst*1e9)
    IF i EQ 0 THEN mergertime = [1.4,3.56,5] ELSE mergertime = [1.65,2.56,3.2,4.27,5.34,5.82]
    FOR j = 0, n_elements(mergertime) - 1 DO oplot,[mergertime[j],mergertime[j]],[0,3],color  = 150,thick = 6
    histogramp,s.tform,weight = s.mass/tbinsize,min = 1e-6,max = 14,nbins = 100,xtitle = xtitle,yrange = yrange,thick = thick,ytickname = ytickname
    histogramp,s[bulge].tform,weight=s[bulge].mass/tbinsize,min = 1e-6,max = 14,nbins = 100,linestyle = 2,/overplot,thick = thick
    stop
    IF i EQ 0 THEN legend,['All Stars','Bulge'],linestyle = [0,2],/right,/top,box = 0
    xyouts,0.5,yout,key[i]
    multiplot
    print,total(s[bulge].mass)/total(s.mass)
    print,total(s[bulge].mass)/(total(s[disk].mass)+total(s[bulge].mass))
ENDFOR
multiplot,/reset
IF keyword_set(outplot) THEN device,/close ELSE stop
END
