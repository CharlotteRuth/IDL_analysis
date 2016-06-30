PRO plot_bulges,filename,n,key = key, color = color, outplot = outplot,symbols = symbols, symsizes = symsizes,ctables = ctables,formatthick = formatthick

readcol,filename,name,M_B,SB_B,log_r_B,n_S,M_D,SB_D,log_r_D,BT,format = '(A,F,F,F,F,F,F,F,F)'

n2 = n_elements(M_B)/n
M_B     = reform(M_B,    n,n2)
log_r_B = reform(log_r_B,n,n2)

formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF keyword_set(color) THEN BEGIN
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n) + 39
    IF color[0] EQ 1 THEN  colors = (findgen(n) + 1)*240/n ELSE colors = color
    IF NOT keyword_set(symbols) THEN symbols = fltarr(n) + 4
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
ENDIF ELSE BEGIN  
    ctables = fltarr(n)
    IF keyword_set(outplot)  THEN colors = fgcolor + (findgen(n) + 1)*10.0 + 5.0 ELSE colors = fgcolor - ((findgen(n) + 1)*10.0 + 5.0)
    IF NOT keyword_set(symbols) THEN symbols = (findgen(n)+2)*2
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
ENDELSE
symthick = 5

IF keyword_set(outplot) THEN device,filename=outplot+'_bulge.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 712,ysize = 712
plot,log_r_B,M_B,psym = 1,xtitle = textoidl('log(r_e) [log(pc)]'),ytitle = textoidl('M_H(Bulge)'),xrange = [1,5],yrange = [-18,-27],/ystyle,/xstyle,/nodata
FOR i = 0, n-1 DO BEGIN
   loadct,ctables[i]
   oplot,log_r_B[i,*],M_B[i,*],psym = symcat(symbols[i],thick = symthick),symsize = symsizes[i],color = colors[i]
ENDFOR
IF keyword_set(key) THEN legend, ['Observed Bulges','Observed Ellipticals',key], psym = [15,15,symbols], color = [50,50,colors], thick = [1,1,fltarr(n) + symthick],ctables = [39,39,ctables], /bottom, /right, box = 0
IF keyword_set(outplot) THEN device,/close
END
