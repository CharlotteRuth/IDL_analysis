PRO plot_z_m,avez,dirs,files,red_avez = red_avez, halo = halo,colors = colors,outplot = outplot,normalize = normalize,symbols = symbols,formatthick = formatthick,absolute = absolute,legend = legend,stellarmass = stellarmass,unscaledz = unscaledz
zsolar  =  0.0130215
IF NOT keyword_set(finalstep) THEN finalstep = '00512'
n = n_elements(dirs)
;formatplot,outplot = outplot,thick = formatthick

IF keyword_set(colors) THEN BEGIN
    distinct_colors,n_colors = 12 ;loadct,39
    white = 13 ;white = 255
    black = 0 ;black = 0
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symbols) THEN symbols = [15]
    symbols_boarder = [6]
    IF NOT keyword_set(symsize) THEN symsize = 2
    colore = [black] ;[254]
    z_colors = [7,2,11,0];[30,80,120,254]
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2)  
    IF NOT keyword_set(symbols) THEN symbols = [15]
    symbols_boarder = [6]
    IF NOT keyword_set(symsize) THEN symsize = 2
;    colore = [254]
    z_colors = [80,120,254,200]
ENDELSE
IF keyword_set(outplot) THEN BEGIN
    fgcolor = black
    bgcolor = white
;    xsize = 18
;    ysize = 12
;    IF keyword_set(stellarmass) THEN outplot = outplot + '_sm' ELSE outplot = outplot + '_vm'
ENDIF ELSE BEGIN
    fgcolor = white
    bgcolor = black
;    xsize = 800
;    ysize = 500
ENDELSE

vmass = fltarr(n)
smass = fltarr(n)
zgas = fltarr(n)
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group eq halo[i]))[0]
   vmass[i] = stat[ind].m_tot
   smass[i] = stat[ind].m_star
   readcol,dirs[i] + '/grp' + halo[i] + '.metals.txt',zmetals,ox,fe,coldg,zmetals_H,ox_H,fe_H,mHI_m,mH2_m,Hgas
   zgas[i] = zmetals[n_elements(z) - 1]/coldg[n_elements(z) - 1]/zsolar
ENDFOR

IF keyword_set(absolute) THEN BEGIN
;    outplotext = outplotext + '_zm_abs.eps'
    ytitle = textoidl('Log(Z_{eject}/Z')+sunsymbol()+')'
    yrange = [-1.8,0]
    showtext = 0
ENDIF ELSE BEGIN
;    outplotext = outplotext + '_zm_scale.eps'
    ytitle = textoidl('Log(Z_{eject}/Z_{ISM})')
    yrange = [-0.2,0.6]
;    yrange = [-0.2,3]
;    yrange = [0.1,0.45]
;    yrange = [-0.05,0.25] range for median
    showtext = 1
 ENDELSE
z_bins_legend = ['2.0','1.0','0.5','0.0']

IF keyword_set(stellarmass) THEN BEGIN
   IF NOT keyword_set(absolute) THEN xtitle = 'Stellar Mass [M' +sunsymbol() + ']'
   xaxis = smass 
   xrange = [500000,5e10]
ENDIF ELSE BEGIN
   IF NOT keyword_set(absolute) THEN xtitle = 'Virial Mass [M' +sunsymbol() + ']'
   xaxis = vmass
   xrange = [1e9,1e12]
ENDELSE

;Read in information to calculate the mass loading at each redshift
nz = 4
reejectmassr  = fltarr(n,nz)
reejectmassmetr  = fltarr(n,nz)
sfmassr = fltarr(n,nz)
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
    readcol,dirs[i] + '/grp' + halo[i] + '.ejectz_quant.txt',z_bins_temp,relostmass_temp,relostmassmet_temp,relostmassr_temp,relostmassmetr_temp,reejectmass_temp,reejectmassmet_temp,reejectmassr_temp,reejectmassmetr_temp,reexpellmass_temp,reexpellmassmet_temp,reexpellmassr_temp,reexpellmassmetr_temp,diskgmass_lost_temp,diskgmass_lostmet_temp,diskgmass_lostr_temp,diskgmass_lostmetr_temp,diskgmass_temp,diskgmassmet_temp,diskgmassr_temp,diskgmassmetr_temp,diskgmass_expell_temp,diskgmass_expellmet_temp,diskgmass_expellr_temp,diskgmass_expellmetr_temp,halogmass_temp,halogmassmet_temp,sfmassr_temp
    reejectmassr[i,*] = reejectmassr_temp
    reejectmassmetr[i,*] = reejectmassmetr_temp
    sfmassr[i,*] = sfmassr_temp
ENDFOR

;IF keyword_set(outplot) THEN  device,filename = outplot + outplotext,/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 0, xsize = xsize, ysize = ysize
plot,xaxis,alog10(avez.meanz),/xlog,xrange = xrange,xtitle = xtitle,yrange = yrange,ytitle = ytitle,psym = symcat(symbols[0]),symsize = symsize,/nodata;,xshowtext = showtext
;IF keyword_set(absolute) THEN oplot,xaxis,alog10(zgas),psym = symcat(16),symsize = symsize,color = fgcolor



IF keyword_set(red_avez) THEN BEGIN
   FOR iz = 0, (size(red_avez))[2] - 1 DO BEGIN
      IF keyword_set(absolute) THEN BEGIN
         IF keyword_set(stellarmass) THEN BEGIN
            oplot,red_avez[*,iz].mstar,alog10(red_avez[*,iz].zeject),psym = symcat(symbols[0]),color = z_colors[iz],symsize = 1.5
            oplot,red_avez[*,iz].mstar,alog10(red_avez[*,iz].zeject),psym = symcat(sym_outline(symbols[0])),color = fgcolor,symsize = 1.5
         ENDIF ELSE BEGIN
            oplot,red_avez[*,iz].mvir,alog10(red_avez[*,iz].zeject),psym = symcat(symbols[0]),color = z_colors[iz],symsize = 1.5
            oplot,red_avez[*,iz].mvir,alog10(red_avez[*,iz].zeject),psym = symcat(sym_outline(symbols[0])),color = fgcolor,symsize = 1.5
            
            massloading = reejectmassr[*,iz]/sfmassr[*,iz]
;            oplot,red_avez[*,iz].mvir,alog10(red_avez[*,iz].zISM + 0.01/massloading/zsolar),psym = symcat(sym_outline(symbols[0]),thick = 4),color = z_colors[iz],symsize = 1.5
            stop
         ENDELSE
      ENDIF ELSE BEGIN
         IF keyword_set(stellarmass) THEN BEGIN
            oplot,red_avez[*,iz].mstar,alog10(red_avez[*,iz].zeject/red_avez[*,iz].zISM),psym = symcat(symbols[0]),color = z_colors[iz],symsize = 1.5
            oplot,red_avez[*,iz].mstar,alog10(red_avez[*,iz].zeject/red_avez[*,iz].zISM),psym = symcat(sym_outline(symbols[0])),color = fgcolor,symsize = 1.5
         ENDIF ELSE BEGIN
            oplot,red_avez[*,iz].mvir,alog10(red_avez[*,iz].zeject/red_avez[*,iz].zISM),psym = symcat(symbols[0]),color = z_colors[iz],symsize = 1.5
            oplot,red_avez[*,iz].mvir,alog10(red_avez[*,iz].zeject/red_avez[*,iz].zISM),psym = symcat(sym_outline(symbols[0])),color = fgcolor,symsize = 1.5

            massloading = reejectmassr[*,iz]/sfmassr[*,iz]
;            oplot,red_avez[*,iz].mvir,alog10((red_avez[*,iz].zISM*zsolar + 0.01/massloading)/(red_avez[*,iz].zISM*zsolar)),psym = symcat(sym_outline(symbols[0]),thick = 4),color = z_colors[iz],symsize = 1.5
            stop
         ENDELSE
      ENDELSE
   ENDFOR
   If keyword_set(legend) THEN legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = symbols[0],color  = z_colors,/right,/bottom,box = 0
ENDIF ELSE BEGIN
   oplot,xaxis,alog10(avez.meanz),psym = symcat(symbols[0]),symsize = symsize,color = fgcolor ;colore[0]
;IF keyword_set(symbols_boarder) THEN oplot,xaxis,alog10(avez.medianz),psym = symcat(symbols_boarder[0]),symsize = symsize,color = fgcolor,thick = 1
ENDELSE
IF NOT keyword_set(absolute) THEN oplot,xrange,[0,0],thick = thicks[0]
;legend,['Gas Ejected','Star Formation'],psym = [symbols[1],46],color = [colore[1],fgcolor],/bottom,/right,box =0
;IF keyword_set(outplot) THEN device, /close ELSE stop
stop
END
