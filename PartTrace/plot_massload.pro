PRO plot_massload,dirs,files,haloid,outplot = outplot
formatplot,outplot = outplot,thick = formatthick
;loadct,39
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 18
    thick = 2
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 400
    ysize = 400
    thick = 1
ENDELSE
device,decomposed = 0


mvir_all = [0]
ml_eject_all = [0]
ml_lost_all = [0] 
num = intarr(n_elements(dirs))
colors = [30,120,60,60,254,100,230,160,80,210]

FOR i = 0, n_elements(dirs)-1 DO BEGIN
   readcol,dirs[i] + '/grp' + haloid[i] + '.massloading.txt',time,mvir,ml_eject,ml_lost,format = '(F,F,F,F)',/silent
   mvir_all = [mvir_all,mvir]
   ml_eject_all = [ml_eject_all,ml_eject]
   ml_lost_all = [ml_lost_all,ml_lost]
   num[i] = n_elements(ml_lost)
ENDFOR

mvir_all = mvir_all[1:n_elements(mvir_all) - 1]
ml_eject_all = ml_eject_all[1:n_elements(ml_eject_all) - 1]
ml_lost_all = ml_lost_all[1:n_elements(ml_lost_all) - 1]

ml_eject_cut = ml_eject_all[where(mvir_all GT 0 AND ml_eject_all GT 0)]
mvir_cut = mvir_all[where(mvir_all GT 0 AND ml_eject_all GT 0)]

fits_eject = robust_linefit( alog10(mvir_cut), alog10(ml_eject_cut), ml_eject_fit, sigma_eject )
print,fits_eject

;yrange = [min([ml_lost_all[where(ml_lost_all NE 0)],ml_eject_all[where(ml_eject_all NE 0)]]),max([ml_lost_all,ml_eject_all])]
;yrange = [min(ml_eject_all[where(ml_eject_all NE 0)]),max(ml_eject_all)]
yrange = [0.02,200]
xrange = [1e8,1e12]

IF keyword_set(outplot) THEN  device,filename = outplot + '_ml_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,mvir_all,ml_eject_all,xtitle = textoidl('M_{vir} [M') + sunsymbol() + ']',ytitle = textoidl('\eta'),psym = 4,/ylog,/xlog,yrange = yrange,xrange = xrange,/nodata
;oplot,mvir_all,ml_eject_all,psym = symcat(14),color = 254,symsize = 1.5
;oplot,mvir_all,ml_lost_all,psym = symcat(14),color = 60,symsize = 1.5
FOR i = 0, n_elements(dirs) - 1 DO $
   oplot,mvir_all[total(num[0:i]) - num[i]:total(num[0:i]) - 1],ml_eject_all[total(num[0:i]) - num[i]:total(num[0:i]) - 1],psym = symcat(14),color = colors[i],symsize = 1.5
oplot,xrange,10^fits_eject[0]*xrange^(fits_eject[1])
IF keyword_set(outplot) THEN  device,/close ELSE stop
END
