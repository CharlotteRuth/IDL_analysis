PRO tully_fisher_pizagno,vfinal,imagsim,color = color,symbols = symbols,symsizes = symsizes,thicks = thicks,ctables = ctables,obscolor = obscolor,obspsym = obspsym,obssymsize = obssymsize,outplot = outplot,keys = keys

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Code/Datafiles/' ELSE prefix = '/astro/users/christec/code/Datafiles/'

formatplot,outplot = outplot
IF keyword_set(outplot) THEN BEGIN
    l_charsize = 1.0
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    l_charsize = 0.75
    fgcolor = 255
    bgcolor = 0
ENDELSE

nfiles = n_elements(vfinal)
IF keyword_set(color) THEN BEGIN
    IF NOT keyword_set(ctables) THEN ctables = fltarr(nfiles) + 39
    IF color[0] EQ 1 THEN  colors = (findgen(nfiles) + 1)*240/nfiles else colors = color
;    IF NOT keyword_set(thicks) THEN thicks = fltarr(nfiles) + 6;2
    IF NOT keyword_set(symbols) THEN symbols = fltarr(nfiles) + 4
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(nfiles) + 2
    IF NOT keyword_set(obsct) THEN obsct = 39
    IF NOT keyword_set(obscolor) THEN obscolor = fgcolor
    IF NOT keyword_set(obspsym) THEN obspsym = 2
    IF NOT keyword_set(obssymsize) THEN obssymsize = 1
ENDIF ELSE BEGIN
    IF NOT keyword_set(ctables) THEN ctables = fltarr(nfiles)
    colors = fltarr(nfiles) + fgcolor
;    colors = (findgen(nfiles) + 1)*10.0 + 5.0;  fltarr(nfiles) + 5
;    thicks = (findgen(nfiles) + 1)*6/nfiles - 1
;    IF NOT keyword_set(THICKS) THEN thicks = fltarr(nfiles) + 6;2
    IF NOT keyword_set(symbols) THEN symbols = fltarr(nfiles) + 2;(findgen(nfiles)+2)*2
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(nfiles) + 2
    IF NOT keyword_set(obsct) THEN obsct = 0    
    IF NOT keyword_set(obscolor) THEN obscolor = 150
    IF NOT keyword_set(obspsym) THEN obspsym = 2
    IF NOT keyword_set(obssymsize) THEN obssymsize = 1
ENDELSE
IF NOT keyword_set(thicks) THEN thicks = fltarr(nfiles) + 1
IF NOT keyword_set(obsthick) THEN obsthick = 1

readcol,prefix + 'TullyFisher/PizagnoData.txt',name,imag,err_imag,vcirc,err_vcirc,format='(A19 4F)'

IF keyword_set(outplot) THEN device,filename=outplot + '_ptf.eps', /color,xsize = 15, ysize = 18 ELSE window,0,xsize = 400,ysize = 600
plot,vcirc,imag,/xlog,xtitle='Rotational Velocity [km/s]',ytitle = 'SDSS i Magnitude',/nodata,yrange = [-18,-24],xrange = [50,400]
oploterror,vcirc,imag,err_vcirc,err_imag,color = obscolor,errcolor = obscolor,psym = symcat(obspsym),symsize = obssymsize,thick = obsthick
FOR i = 0, nfiles - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[vfinal[i],vfinal[i]],[imagsim[i],imagsim[i]] - 0.37,color = colors[i],psym = symcat(symbols[i]),symsize = symsizes[i],thick = thicks[i]
ENDFOR
IF n_elements(keys) EQ 1 THEN legend,['Pizagno et al., 2007',keys],color = [obscolor,color[0]],psym = [obspsym,symbols[0]],symsize = [obssymsize,symsizes[0]],ctables = [obsct,ctables[0]],thick = [obsthick,thicks[0]],/top,/left,box = 0,charsize = l_charsize ELSE $
  IF n_elements(keys) GT 1 THEN legend,['Pizagno et al., 2007',keys],color = [obscolor,color],psym = [obspsym,symbols],symsize = [obssymsize,symsizes],ctables = [obsct,ctables],thick = [obsthick,thicks],/top,/left,box = 0,charsize = l_charsize
IF keyword_set(outplot) THEN device,/close ELSE stop
END
