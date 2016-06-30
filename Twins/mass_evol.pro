
pro mass_evol,dir,files,pfiles,scale = scale,vline = vline,match = match,linestyle = linestyle,thicks = thicks,colors = colors,psym = psym,symsize = symsize,keys = keys,outplot = outplot,formatthick = formatthick
formatplot,outplot = outplot,thick = formatthick
;window,xoffset = 2,yoffset = 2
;colors = [50,245]
;colors = [30,245]

highz_time = 1.93996e+09/1e9 ;dtime = 5.002289e-02 ;72, z = 3.43597, a = 0.22542965
;highz_time = 2.6905329E+00; 100
;104 2.797757 ;e9
ZSOLAR  =  0.0130215
filebase = 'masstrack.dat'
n = N_ELEMENTS(dir)
IF KEYWORD_SET(outplot) THEN BEGIN 
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF KEYWORD_SET(colors) THEN BEGIN
    loadct,39
    IF colors[0] EQ 1 THEN  colors = (findgen(n) + 1)*240/n ELSE colors = colors
    IF NOT KEYWORD_SET(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2    ;1
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) 
    IF NOT KEYWORD_SET(psym) THEN psym =-1*( fltarr(n) + 4)
    IF NOT KEYWORD_SET(symsize) THEN symsize = fltarr(n) + 2
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT KEYWORD_SET(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2;1 ;findgen(n)*2
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

;!p.multi = [0,1,4]
if (KEYWORD_SET(outplot)) then device,filename = outplot + '_mass.eps',/color,bits_per_pixel= 8,/times,ysize = 36,xsize = 18,xoffset =  2,yoffset =  2 else window,0,xsize = 500,ysize = 800
!X.margin = [17,17]
!Y.margin = [6,6]
multiplot,[1,4]

lthick = 4

;------------------ Mvir
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mtot,ytitle = textoidl('M_{vir}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,4e11],ystyle = 9;,xtitle = 'z'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{vir}') + '[M' + sunsymbol() + ']'
        ENDIF
        oplot,time,mtot*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mtot,ytitle = textoidl('M_{vir}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,4e11],ystyle = 9;,xtitle = 'Time [Gyr]'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{vir}') + '[M' + sunsymbol() + ']'
        ENDIF
        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,mtot*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
        if i eq 1 AND keyword_set(match) then oplot,[0,14],[mtot[N_ELEMENTS(mtot) - 1],mtot[N_ELEMENTS(mtot) - 1]],color = colors[i],thick = lthick
    ENDELSE
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,psym = psym,box = 0
ENDFOR
multiplot
;------------------ Mgas
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mgas,ytitle = textoidl('M_{gas}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,3.5e10 - 1e8],ystyle = 9;,xtitle = 'z'
;            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{gas}') + '[M' + sunsymbol() + ']'
        ENDIF
        oplot,time,mgas*scale[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mgas,ytitle = textoidl('M_{gas}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,3.5e10 - 1e8],ystyle = 9;,xtitle = 'Time [Gyr]'
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{gas}') + '[M' + sunsymbol() + ']'
;            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,mgas*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
        if i eq 1  AND keyword_set(match) then oplot,[0,14],[mgas[N_ELEMENTS(mgas) - 1],mgas[N_ELEMENTS(mgas) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
multiplot
;------------------ M cold gas
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,1.4*mcoldg,ytitle = textoidl('M_{cold}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,4.5e9 - 1],ystyle = 9;,xtitle = 'z'
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{cold}') + '[M' + sunsymbol() + ']'
;            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,1.4*mcoldg*scale[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,1.4*mcoldg,ytitle = textoidl('M_{cold}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,4.5e9 - 1],ystyle = 9;,xtitle = 'Time [Gyr]'
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{cold}') + '[M' + sunsymbol() + ']'
;            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,1.4*mcoldg*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
        if i eq 1  AND keyword_set(match) then oplot,[0,14],1.4*[mcoldg[N_ELEMENTS(mcoldg) - 1],mcoldg[N_ELEMENTS(mcoldg) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
multiplot
;------------------ Mstar
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mstar,ytitle = textoidl('M_{star}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,8e9 - 1e7],xtitle = 'z',ystyle = 9
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
;            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,mstar*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mstar,ytitle = textoidl('M_{star}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,8e9 - 1e7],xtitle = 'Time [Gyr]',ystyle = 9
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
;            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,mstar*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
        if i eq 1  AND keyword_set(match) then oplot,[0,14],[mstar[N_ELEMENTS(mstar) - 1],mstar[N_ELEMENTS(mstar) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
multiplot,/reset

;!X.margin = [10,3]
;!Y.margin = [4,2]
;!X.margin = [13,13]
;!Y.margin = [6,6]
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_mvir.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE window,1,xsize = 712,ysize = 392
;------------------ Mvir
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mtot,ytitle = textoidl('M_{vir}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,4e11],ystyle = 9,xmargin = [17,17],ymargin = [6,6];,xtitle = 'z'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{vir}') + '[M' + sunsymbol() + ']'
        ENDIF
        oplot,time,mtot*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mtot,ytitle = textoidl('M_{vir}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,4e11],ystyle = 9,xmargin = [17,17],ymargin = [6,6],xtitle = 'Time [Gyr]'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{vir}') + '[M' + sunsymbol() + ']'
        ENDIF
 ;       if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,mtot*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
 ;       if i eq 1  AND keyword_set(match) then oplot,[0,14],[mtot[N_ELEMENTS(mtot) - 1],mtot[N_ELEMENTS(mtot) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_mgas.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE window,1,xsize = 712,ysize = 392
;------------------ Mgas
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mgas,ytitle = textoidl('M_{gas}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,3.5e10 - 1e8],ystyle = 9,xmargin = [17,17],ymargin = [6,6];,xtitle = 'z'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{gas}') + '[M' + sunsymbol() + ']'
        ENDIF
        oplot,time,mgas*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mgas,ytitle = textoidl('M_{gas}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,3.5e10 - 1e8],ystyle = 9,xmargin = [17,17],ymargin = [6,6],xtitle = 'Time [Gyr]'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{gas}') + '[M' + sunsymbol() + ']'
        ENDIF
 ;       if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,mgas*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
 ;       if i eq 1  AND keyword_set(match) then oplot,[0,14],[mtot[N_ELEMENTS(mtot) - 1],mtot[N_ELEMENTS(mtot) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_mbar.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE window,1,xsize = 712,ysize = 392
;------------------ Mbaryonic
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mgas + mstar,ytitle = textoidl('M_{baryon}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,4e10 - 1e8],ystyle = 9,xmargin = [17,17],ymargin = [6,6];,xtitle = 'z'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{baryon}') + '[M' + sunsymbol() + ']'
        ENDIF
        oplot,time,(mgas + mstar)*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mgas + mstar,ytitle = textoidl('M_{baryon}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,4e10 - 1e8],ystyle = 9,xmargin = [17,17],ymargin = [6,6],xtitle = 'Time [Gyr]'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{baryon}') + '[M' + sunsymbol() + ']'
        ENDIF
 ;       if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,(mgas + mstar)*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
 ;       if i eq 1  AND keyword_set(match) then oplot,[0,14],[mtot[N_ELEMENTS(mtot) - 1],mtot[N_ELEMENTS(mtot) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_cgas.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE window,1,xsize = 712,ysize = 392
;------------------ M Cold Gas
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,1.4*mcoldg,ytitle = textoidl('M_{cold}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,4.5e9 - 1],ystyle = 9,xmargin = [17,17],ymargin = [6,6];,xtitle = 'z'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{cold}') + '[M' + sunsymbol() + ']'
        ENDIF
        oplot,time,1.4*mcoldg*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,1.4*mcoldg,ytitle = textoidl('M_{cold}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,4.5e9 - 1],ystyle = 9,xmargin = [17,17],ymargin = [6,6],xtitle = 'Time [Gyr]'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{cold}') + '[M' + sunsymbol() + ']'
        ENDIF
 ;       if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,1.4*mcoldg*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
 ;       if i eq 1  AND keyword_set(match) then oplot,[0,14],[mtot[N_ELEMENTS(mtot) - 1],mtot[N_ELEMENTS(mtot) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_smass.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
ENDIF ELSE window,1,xsize = 712,ysize = 392
;------------------ Mstar
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mstar,ytitle = textoidl('M_{star}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,8e9 - 1e7],xtitle = 'z',ystyle = 9,xmargin = [17,17],ymargin = [6,6]
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,mstar*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mstar,ytitle = textoidl('M_{star}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,8e9 - 1e7],xtitle = 'Time [Gyr]',ystyle = 9,xmargin = [17,17],ymargin = [6,6]
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
;        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,mstar*scale[i],linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[mstar[N_ELEMENTS(mstar) - 1],mstar[N_ELEMENTS(mstar) - 1]],color = colors[i],thick = lthick
    ENDELSE
 ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,psym = -1.0*psym,box = 0
;!X.margin = [10,3]
;!Y.margin = [4,2]
;!X.margin = [13,13]
;!Y.margin = [6,6]

IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_star_barmass.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE window,1,xsize = 712,ysize = 392
;------------------ Mstar/MB--------------------------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mstar/(mstar + 1.4*mcoldg),ytitle = textoidl('M_{star}/(M_{star} + M_{cold})'),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata,/ylog,yrange = [0.3,1]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,mstar/(mstar + 1.4*mcoldg)
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mstar/(mstar + 1.4*mcoldg),ytitle = textoidl('M_{star}/(M_{star} + M_{cold})'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],/nodata,/ylog,yrange = [0.3,1]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,mstar/(mstar + 1.4*mcoldg),linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[mstar[N_ELEMENTS(mstar) - 1]/(mstar[N_ELEMENTS(mstar) - 1] + 1.4*mcoldg[N_ELEMENTS(mstar) - 1]),mstar[N_ELEMENTS(mstar) - 1]/(mstar[N_ELEMENTS(mstar) - 1] + 1.4*mcoldg[N_ELEMENTS(mstar) - 1])],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,/right,/bottom,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_H2_HI.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE window,2,xsize = 712,ysize = 392
;------------------  H2/HI---------------------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mH2/(mcoldg - mH2),ytitle = textoidl('H_2/HI'),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata,/ylog,yrange =[1e-5,0.03]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,mH2/(mcoldg - mH2)
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mH2/(mcoldg - mH2),ytitle = textoidl('H_2/HI'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],/nodata,/ylog,yrange =[1e-5,0.1]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,mH2/(mcoldg - mH2),linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[mH2[N_ELEMENTS(mstar) - 1]/(mcoldg[N_ELEMENTS(mstar) - 1] - mH2[N_ELEMENTS(mstar) - 1]),mH2[N_ELEMENTS(mstar) - 1]/(mcoldg[N_ELEMENTS(mstar) - 1] - mH2[N_ELEMENTS(mstar) - 1])],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,/right,/bottom,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_H2_HI_inner.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
ENDIF ELSE window,3,xsize = 712,ysize = 392

;------------------  H2/HI inner ---------------------------------

IF 1 THEN BEGIN
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,H2frac,ytitle = textoidl('H_2/HI'),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata,/ylog,yrange =[1e-5,0.1]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,H2frac
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,H2frac,ytitle = textoidl('H_2/HI'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],/nodata,/ylog,yrange =[1e-5,0.1]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,H2frac,linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[mH2[N_ELEMENTS(mstar) - 1]/(mcoldg[N_ELEMENTS(mstar) - 1] - mH2[N_ELEMENTS(mstar) - 1]),mH2[N_ELEMENTS(mstar) - 1]/(mcoldg[N_ELEMENTS(mstar) - 1] - mH2[N_ELEMENTS(mstar) - 1])],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,/right,/bottom,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_metal.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
ENDIF ELSE window,4,xsize = 712,ysize = 392
ENDIF

;------------------  z---------------------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,metal/zsolar,ytitle = textoidl('Z/Z') + sunsymbol(),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata,/ylog,yrange =[1e-4,0.01]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,metal
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,metal/zsolar,ytitle = textoidl('Z/Z') + sunsymbol(),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],/nodata,/ylog,yrange =[3e-3,0.3]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
 ;       if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,metal/zsolar,linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[metal[N_ELEMENTS(metal) - 1],metal[N_ELEMENTS(metal) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,/right,/bottom,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_mox.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE window,5,xsize = 712,ysize = 392

;------------------  ox---------------------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mox,ytitle = textoidl('12 + log(O/H)'),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata,yrange =[1e-4,0.01]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,metal
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mox,ytitle = textoidl('12 + log(O/H)'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],yrange = [6.5,9],/nodata
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
 ;       if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,mox,linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[metal[N_ELEMENTS(metal) - 1],metal[N_ELEMENTS(metal) - 1]],color = colors[i],thick = lthick
    ENDELSE
 ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,/right,/bottom,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_pressure.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE window,6,xsize = 712,ysize = 392

;----------------- Pressure -----------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,pressure,ytitle = textoidl('P/k_B [K cm^{-3}]'),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata,/ylog,yrange = [1e3,1e7]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,pressure
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,pressure,ytitle = textoidl('P/k_B [K cm^{-3}]'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],/nodata,/ylog,yrange = [2e3,3e6]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
  ;      if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,pressure,linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[metal[N_ELEMENTS(metal) - 1],metal[N_ELEMENTS(metal) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,/right,/top,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_sfe.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE window,7,xsize = 712,ysize = 392
;----------------- SFE -----------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,SFR/mcoldg,ytitle = textoidl('Star Formation Efficiency [yr^{-1}]'),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata,/ylog;,yrange = [1e2,1e6]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,pressure
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,SFR/mcoldg,ytitle = textoidl('Star Formation Efficiency [yr^{-1}]'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],/nodata,/ylog,yrange = [3e-12,1e-8]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
  ;      if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,SFR/mcoldg,linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        oplot,time,P_sfe(Pressure),linestyle = 1,thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,/right,/top,psym = -1.0*psym,box = 0

IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_j.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
ENDIF ELSE  window,8,xsize = 712,ysize = 392
;----------------- J -----------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    readcol,'j.dat',j
    spawn,'ls '+dir[i]+'../h*.param',pfilelist
    units = tipsyunits(pfilelist[0])
    g_per_msol = 1.98892d33
    cm_per_kpc = 3.08568025d21
    j = j*units.lengthunit*units.massunit*units.vunit/cm_per_kpc

    IF scale[i] NE 1 THEN scaleL = scale[i]*4 ELSE scaleL = scale[i]
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,j/j[n_elements(j) - 1],ytitle = textoidl('Average Angular Momentum'),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata;,/ylog;,yrange = [1e2,1e6]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,pressure
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,j*scaleL,ytitle = textoidl('<J_z>'),xtitle = 'Time [Gyr]',xstyle = 9,ystyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],/nodata,yrange = [0,0.014];,/ylog
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('32 X <J_z>')
        ENDIF
  ;      if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,j*scaleL,linestyle = linestyle[i],thick = thicks[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[metal[N_ELEMENTS(metal) - 1],metal[N_ELEMENTS(metal) - 1]],color = colors[i],thick = lthick
     ENDELSE
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,/left,/top,psym = -1.0*psym,box = 0
IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot + '_Psfe.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 12,xoffset =  2,yoffset =  2 
ENDIF ELSE window,9,xsize = 712,ysize = 712
;----------------- P/SFE -----------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    if i eq 0 THEN BEGIN
       plot,Pressure,SFR/mcoldg,ytitle = textoidl('Star Formation Efficiency [yr^{-1}]'),xtitle = textoidl('P/k_B [K cm^{-3}]'),xstyle = 1,xmargin = [12,6],ymargin = [12,6],xrange = [1e3,1e6],/nodata,/ylog,/xlog,yrange = [1e-11,1e-8]
    ENDIF
    oplot,Pressure,SFR/mcoldg,linestyle = linestyle[i],thick = thicks[i],psym = symcat(psym[i]),color = colors[i]
 ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,psym = psym,/right,/bottom,box = 0

if (KEYWORD_SET(outplot)) then begin
    device,/close
ENDIF ELSE stop
end
