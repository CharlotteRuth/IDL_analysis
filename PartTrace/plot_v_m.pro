PRO plot_v_m,medv_ej,medv_exp,stdv_ej,stdv_exp,dirs,files,halo = halo,colors = colors,outplot = outplot,normalize = normalize,symbols = symbols,formatthick = formatthick,absolute = absolute,vfinal = vfinal
IF NOT keyword_set(finalstep) THEN finalstep = '00512'
n = n_elements(dirs)
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
;    IF keyword_set(stellarmass) THEN outplot = outplot + '_sm' ELSE outplot = outplot + '_vm'
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF keyword_set(stellarmass) THEN BEGIN
   IF NOT keyword_set(absolute) THEN xtitle = 'Stellar Mass [M' + sunsymbol() + ']'
   xrange = [1e7,1e11]
ENDIF ELSE BEGIN
   IF NOT keyword_set(absolute) THEN xtitle = 'Virial Mass [M' + sunsymbol() + ']'
   xrange = [3e9,1e12]
ENDELSE
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symbols) THEN symbols = [15,14]
    symbols_boarder = [6,4]
    IF NOT keyword_set(symsize) THEN symsize = 2
    colore = [254,50]
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2)  
    IF NOT keyword_set(symbols) THEN symbols = [15,14]
    symbols_boarder = [6,4]
    IF NOT keyword_set(symsize) THEN symsize = 2
    colore = [254,50]
ENDELSE

vmass = fltarr(n)
smass = fltarr(n)
vc = fltarr(n)
print,'amiga: '
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group eq halo[i]))[0]
   vmass[i] = stat[ind].m_tot
   smass[i] = stat[ind].m_star
   vc[i] = stat[ind].vc
   vc[i] = sqrt(6.67e-11*stat[ind].m_tot*1.9891e30/(3.08567758e19*stat[ind].rvir))/1000.
   print,stat[ind].vc
ENDFOR
print,'grav: ',vc
IF keyword_set(vfinal) THEN BEGIN
;   plot,vc,vfinal/vc,psym = 4,xtitle = 'V_circ from virial mass and radius',ytitle = 'v_btf/v_circ',yrange = [1,2]
   vc = vfinal
ENDIF

;IF keyword_set(stellarmass) THEN outplotext = '_sm' ELSE outplotext = '_vm'
IF keyword_set(absolute) THEN BEGIN
;    outplotext = outplotext + '_velm_abs.eps'
    ytitle = 'Velocity [km/s]'
    yrange = [0,200];255]
ENDIF ELSE BEGIN
;    outplotext = outplotext + '_velm_scale.eps'
    ytitle = textoidl('V/V_{Circ}')
    yrange = [0,7]
ENDELSE

IF keyword_set(stellarmass) THEN xaxis = smass ELSE xaxis = vmass
plot,xaxis,medv_ej,/xlog,xrange = xrange,xtitle = xtitle,yrange = yrange,ytitle = ytitle,psym = symcat(symbols[0]),symsize = symsize,/nodata
IF keyword_set(absolute) THEN oplot,xaxis,vc,psym = symcat(16),symsize = symsize,color = fgcolor ELSE BEGIN
   oplot,xrange,[1,1],thick = 4; thicks[0]
   oplot,xrange,[sqrt(2),sqrt(2)],thick = 4,linestyle = 2
   oploterror,xaxis,medv_exp,stdv_exp,psym = 3,color = fgcolor,errcolor = colore[1],thick = 4
   oploterror,xaxis,medv_ej,stdv_ej,psym = 3,color = fgcolor,errcolor = colore[0],thick = 4
ENDELSE
oplot,xaxis,medv_exp,psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
oplot,xaxis,medv_ej,psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
IF keyword_set(symbols_boarder) THEN oplot,xaxis,medv_exp,psym = symcat(symbols_boarder[1]),symsize = symsize,color = fgcolor,thick = 1
IF keyword_set(symbols_boarder) THEN oplot,xaxis,medv_ej,psym = symcat(symbols_boarder[0]),symsize = symsize,color = fgcolor,thick = 1
;legend,['Gas Ejected','Star Formation'],psym = [symbols[1],46],color = [colore[1],fgcolor],/bottom,/right,box =0
;IF keyword_set(absolute) THEN legend,['Ejected Gas','Circular Velocity, z = 0'],psym = [symbols[0],16],color = [colore[0],fgcolor],/top,/left,box = 0
IF keyword_set(absolute) THEN legend,['Gas Ejected','Gas Expelled',textoidl('V_{circ}, z = 0')],psym = [symbols,16],color = [colore,fgcolor],/top,/left,box = 0 $
ELSE legend,['Mass ejected','Mass expelled from halo'],psym = [symbols],color = [colore],/top,/right,box = 0
END
