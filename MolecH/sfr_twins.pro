PRO sfr_twins,files,pfiles,outplot = outplot,linestyle = linestyle, colors = colors,thick = thick,scale = scale,keys = keys,vline = vline,formatthick = formatthick
formatplot,outplot = outplot,thick = formatthick

highz_time = 1.93996e+09/1e9 ;dtime = 5.002289e-02 ;72, z = 3.43597, a = 0.22542965
n = N_ELEMENTS(files)
IF KEYWORD_SET(outplot) THEN BEGIN 
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF KEYWORD_SET(colors) THEN BEGIN
    loadct,39
    if colors[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(thick) THEN  $
      IF NOT keyword_set(formatthick) THEN thick = fltarr(n) + 2 $
      ELSE thick = fltarr(n) + 4
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) 
    IF NOT KEYWORD_SET(psym) THEN psym =-1*( fltarr(n) + 4)
    IF NOT KEYWORD_SET(symsize) THEN symsize = fltarr(n) + 2
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT KEYWORD_SET(thick) THEN $
      IF NOT keyword_set(formatthick) THEN  thick = fltarr(n) + 2 $
      ELSE thick = fltarr(n) + 4
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = findgen(n)*2   
    IF NOT KEYWORD_SET(psym) THEN psym = -1*(fltarr(n) + 4)
    IF NOT KEYWORD_SET(symsize) THEN symsize = fltarr(n) + 2
ENDELSE
IF NOT KEYWORD_SET(scale) THEN scale = fltarr(n) + 1.0

ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
z = REVERSE(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
ticktime_in_z = findgen(7)*2*1e9 ;redshift major plot
tickred_in_z = spline(t,z,ticktime_in_z )

tickred_in_t = REVERSE(findgen(5))
ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e9

if (KEYWORD_SET(outplot)) then begin
    device,filename = outplot + '_sfr.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    window,0,xsize = 712,ysize = 392
endelse
;------------------  SFR ---------------------------------
FOR i = 0, n - 1 DO BEGIN
    rtipsy,files[i],h,g,d,s
    units = tipsyunits(pfiles[i])
    if i eq 0 then BEGIN
        sfr,s,massunit = units.massunit*scale[i], timeunit = units.timeunit,linestyle = linestyle[i],thick = thick[i],/redshift,xmargin = [17,17],ymargin = [6,6],binsize = 1e8,yrange = [0,2.5],ystyle = 9,sfh = sfh0,mint = 0.1, maxt = 1.37346e+10
        axis,ystyle = 1,yaxis = 1,ytitle = '8 X SFR [M'+sunsymbol()+' yr!u-1!n]'
        h0 = h
    ENDIF ELSE h1 = h
    if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
    sfr,s,massunit = units.massunit*scale[i], timeunit = units.timeunit,linestyle = linestyle[i],thick = thick[i],color = colors[i],/overplot,binsize = 1e8,sfh = sfh1,mint = 0.1, maxt = 1.37346e+10
ENDFOR
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,/right;,psym = -1.0*psym
if (KEYWORD_SET(outplot)) then device,/close

;readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
;redshift = spline(t,z,time*1e9)
;sfe = fltarr(n_elements(time))
;deltat = 100*1e6
;time0 = time*1e9 - deltat
;FOR i = 0, n_elements(time) - 1 DO sfe[i] = total(s[where(s.tform*units.timeunit GE time0[i] AND s.tform*units.timeunit LT time[i]*1e9)].mass)*units.massunit/deltat
;stop
END
