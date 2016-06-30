

pro mass_evol,dir,files,pfiles,outplot = outplot,linestyle = linestyle, color = color, psym = psym, symsize = symsize, thick = thick,scale = scale,keys = keys,vline = vline
formatplot,outplot = outplot
;window,xoffset = 2,yoffset = 2
color = [50,245]
color = [30,245]

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
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 2;1
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) 
    IF NOT KEYWORD_SET(psym) THEN psym =-1*( fltarr(n) + 4)
    IF NOT KEYWORD_SET(symsize) THEN symsize = fltarr(n) + 2
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 2;1 ;findgen(n)*2
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
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mtot,ytitle = textoidl('M_{vir}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,4e11],ystyle = 9;,xtitle = 'z'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{vir}') + '[M' + sunsymbol() + ']'
        ENDIF
        oplot,time,mtot*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mtot,ytitle = textoidl('M_{vir}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,4e11],ystyle = 9;,xtitle = 'Time [Gyr]'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{vir}') + '[M' + sunsymbol() + ']'
        ENDIF
        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,mtot*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
        if i eq 1 then oplot,[0,14],[mtot[N_ELEMENTS(mtot) - 1],mtot[N_ELEMENTS(mtot) - 1]],color = colors[i],thick = lthick
    ENDELSE
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = [0,0];,psym = -1.0*psym
ENDFOR
multiplot
;------------------ Mgas
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
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
        oplot,time,mgas*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
        if i eq 1 then oplot,[0,14],[mgas[N_ELEMENTS(mgas) - 1],mgas[N_ELEMENTS(mgas) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
multiplot
;------------------ M cold gas
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
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
        oplot,time,1.4*mcoldg*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
        if i eq 1 then oplot,[0,14],1.4*[mcoldg[N_ELEMENTS(mcoldg) - 1],mcoldg[N_ELEMENTS(mcoldg) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
multiplot
;------------------ Mstar
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
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
        oplot,time,mstar*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
        if i eq 1 then oplot,[0,14],[mstar[N_ELEMENTS(mstar) - 1],mstar[N_ELEMENTS(mstar) - 1]],color = colors[i],thick = lthick
    ENDELSE
ENDFOR
multiplot,/reset

;!X.margin = [10,3]
;!Y.margin = [4,2]
;!X.margin = [13,13]
;!Y.margin = [6,6]
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_star_mvir.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else window,1,xsize = 712,ysize = 392
;------------------ Mvir
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mtot,ytitle = textoidl('M_{vir}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,4e11],ystyle = 9,xmargin = [17,17],ymargin = [6,6];,xtitle = 'z'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{vir}') + '[M' + sunsymbol() + ']'
        ENDIF
        oplot,time,mtot*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mtot,ytitle = textoidl('M_{vir}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,4e11],ystyle = 9,xmargin = [17,17],ymargin = [6,6],xtitle = 'Time [Gyr]'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{vir}') + '[M' + sunsymbol() + ']'
        ENDIF
 ;       if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,mtot*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
 ;       if i eq 1 then oplot,[0,14],[mtot[N_ELEMENTS(mtot) - 1],mtot[N_ELEMENTS(mtot) - 1]],color = colors[i],thick = lthick
    ENDELSE
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = [0,0];,psym = -1.0*psym
ENDFOR

if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_smass.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
endif else window,1,xsize = 712,ysize = 392
;------------------ Mstar
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mstar,ytitle = textoidl('M_{star}') + '[M' + sunsymbol() + ']',xstyle =9,xrange = [4,0],/nodata,yrange =[0.1,8e9 - 1e7],xtitle = 'z',ystyle = 9,xmargin = [17,17],ymargin = [6,6]
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,mstar*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mstar,ytitle = textoidl('M_{star}') + '[M' + sunsymbol() + ']',xstyle = 9,xrange = [0,ageUniverse/1e9],/nodata,yrange =[0.1,8e9 - 1e7],xtitle = 'Time [Gyr]',ystyle = 9,xmargin = [17,17],ymargin = [6,6]
            axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
;        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
        oplot,time,mstar*scale[i],linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[mstar[N_ELEMENTS(mstar) - 1],mstar[N_ELEMENTS(mstar) - 1]],color = colors[i],thick = lthick
    ENDELSE
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = [0,0];,psym = -1.0*psym
ENDFOR
;!X.margin = [10,3]
;!Y.margin = [4,2]
;!X.margin = [13,13]
;!Y.margin = [6,6]

if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_star_barmass.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else window,1,xsize = 712,ysize = 392
;------------------ Mstar/MB--------------------------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
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
        oplot,time,mstar/(mstar + 1.4*mcoldg),linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[mstar[N_ELEMENTS(mstar) - 1]/(mstar[N_ELEMENTS(mstar) - 1] + 1.4*mcoldg[N_ELEMENTS(mstar) - 1]),mstar[N_ELEMENTS(mstar) - 1]/(mstar[N_ELEMENTS(mstar) - 1] + 1.4*mcoldg[N_ELEMENTS(mstar) - 1])],color = color[i],thick = lthick
    ENDELSE
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = [0,0],/right,/bottom;,psym = -1.0*psym
ENDFOR

if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_H2_HI.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else window,2,xsize = 712,ysize = 392
;------------------  H2/HI---------------------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,mH2/(mcoldg - mH2),ytitle = textoidl('H_2/HI'),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata,/ylog,yrange =[1e-5,0.1]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,mH2/(mcoldg - mH2)
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,mH2/(mcoldg - mH2),ytitle = textoidl('H_2/HI'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],/nodata,/ylog,yrange =[1e-5,0.1]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
        if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,mH2/(mcoldg - mH2),linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[mH2[N_ELEMENTS(mstar) - 1]/(mcoldg[N_ELEMENTS(mstar) - 1] - mH2[N_ELEMENTS(mstar) - 1]),mH2[N_ELEMENTS(mstar) - 1]/(mcoldg[N_ELEMENTS(mstar) - 1] - mH2[N_ELEMENTS(mstar) - 1])],color = color[i],thick = lthick
    ENDELSE
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = [0,0],/right,/bottom;,psym = -1.0*psym
ENDFOR

if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_H2_HI_inner.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
endif else window,3,xsize = 712,ysize = 392

;------------------  H2/HI inner ---------------------------------

IF 1 THEN BEGIN
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
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
        oplot,time,H2frac,linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[mH2[N_ELEMENTS(mstar) - 1]/(mcoldg[N_ELEMENTS(mstar) - 1] - mH2[N_ELEMENTS(mstar) - 1]),mH2[N_ELEMENTS(mstar) - 1]/(mcoldg[N_ELEMENTS(mstar) - 1] - mH2[N_ELEMENTS(mstar) - 1])],color = color[i],thick = lthick
    ENDELSE
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = [0,0],/right,/bottom;,psym = -1.0*psym
ENDFOR

if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_metal.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
endif else window,4,xsize = 712,ysize = 392
ENDIF

;------------------  z---------------------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
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
        oplot,time,metal/zsolar,linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[metal[N_ELEMENTS(metal) - 1],metal[N_ELEMENTS(metal) - 1]],color = color[i],thick = lthick
    ENDELSE
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = [0,0],/right,/bottom;,psym = -1.0*psym
ENDFOR
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_pressure.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else window,5,xsize = 712,ysize = 392

;----------------- Pressure -----------------------
FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    redshift = spline(t,z,time*1e9)
    
    IF 0 THEN BEGIN
        if i eq 0 THEN BEGIN
            plot,redshift,pressure,ytitle = textoidl('P/k_B [K cm^{-3}]'),xtitle = 'z',xstyle =9,xmargin = [17,17],ymargin = [6,6],xrange = [4,0],/nodata,/ylog,yrange = [1e2,1e6]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'time_from_z_axes',xtitle = 'Time [Gyr]',xtickv = tickred_in_z,xticks = N_ELEMENTS(tickred_in_z) - 1
        ENDIF
        oplot,time,pressure
    ENDIF ELSE BEGIN
        if i eq 0 THEN BEGIN
            plot,time,pressure,ytitle = textoidl('P/k_B [K cm^{-3}]'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],/nodata,/ylog,yrange = [1e2,1e6]
            axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
        ENDIF
  ;      if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[1e-6,1],thick = lthick
        oplot,time,pressure,linestyle = linestyle[i],thick = thick[i],psym = -1.0*symcat(psym[i]),color = colors[i]
;        if i eq 1 then oplot,[0,14],[metal[N_ELEMENTS(metal) - 1],metal[N_ELEMENTS(metal) - 1]],color = color[i],thick = lthick
    ENDELSE
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = [0,0],/right,/bottom;,psym = -1.0*psym
ENDFOR

if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_sfr.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    stop
    window,6,xsize = 712,ysize = 392
endelse

;------------------  SFR ---------------------------------
FOR i = 0, n - 1 DO BEGIN
    rtipsy,files[i],h,g,d,s
    units = tipsyunits(pfiles[i])
    if i eq 0 then BEGIN
        sfr,s,massunit = units.massunit*scale[i], timeunit = units.timeunit,linestyle = linestyle[i],thick = 2,/redshift,xmargin = [20,20],ymargin = [6,6],binsize = 1e8,yrange = [0,2.5],ystyle = 9,sfh = sfh0,mint = 0.1, maxt = 1.37346e+10
        axis,ystyle = 1,yaxis = 1,ytitle = '8 X SFR [M'+sunsymbol()+' yr!u-1!n]'
        h0 = h
    ENDIF ELSE h1 = h
    if i eq 0 AND KEYWORD_SET(vline) then oplot,[highz_time,highz_time],[0,1e12],thick = lthick
    sfr,s,massunit = units.massunit*scale[i], timeunit = units.timeunit,linestyle = linestyle[i],thick = 2,color = colors[i],/overplot,binsize = 1e8,sfh = sfh1,mint = 0.1, maxt = 1.37346e+10
ENDFOR
    IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = [0,0],/right;,psym = -1.0*psym
if (KEYWORD_SET(outplot)) then device,/close

;------------------- SFR ratio ----------------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_sfrratio.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    window,7,xsize = 712,ysize = 392
endelse
nh = N_ELEMENTS(sfh1)
timeaxes = (findgen(nh)*1e8 + 1e8*0.5)/1e9
cumarr=fltarr(nh,nh)
ii=lindgen(nh*nh)
cumarr(where(ii mod nh le ii/nh)) = 1
sfh0cum=reform(sfh0 # cumarr)
sfh1cum=reform(sfh1 # cumarr)
ii=0
cumarr=0
plot,timeaxes,sfh1/sfh0,/ylog,ytitle = textoidl('SFR_{dwarf}/SFR_{spiral}'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],yrange = [1e-2,1e1]
axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1

;------------------ Cum SF ------------------------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_sfrratio.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    window,8,xsize = 712,ysize = 392
endelse
plot,timeaxes,sfh0cum/sfh0cum[nh - 1],ytitle = textoidl('M_*/M_{*,z = 0}'),xtitle = 'Time [Gyr]',xstyle = 9,xmargin = [17,17],ymargin = [6,6],xrange = [0,ageUniverse/1e9],yrange = [0,1],/nodata
axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1
oplot,timeaxes,sfh0cum/sfh0cum[nh - 1],linestyle = linestyle[0],thick = 2,color = colors[0]
oplot,timeaxes,sfh1cum/sfh1cum[nh - 1],linestyle = linestyle[1],thick = 2,color = colors[1]
legend,keys,color = colors,linestyle = [0,0],/right,/bottom

;------------------ FB ---------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_fb.eps',/color,bits_per_pixel= 8,/times,ysize = 16,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
;    stop
    window,9,xsize = 712,ysize = 600
endelse
binsize = 5e8
nbins = ROUND(14.0/(binsize/1e9))

FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    rtipsy,file,h,g,d,s,/justhead
    readarr,file_coolon,h,coolon,part = 'gas',/ascii
    readarr,file_iord,h,iord,part = 'gas',/ascii
;------------------ Outflow -------------------------
    outflowz_all = mrdfits(dir[i] + '../grp1.rvir.outflow_z.fits')
    outflowm_all = mrdfits(dir[i] + '../grp1.rvir.mass_at_outflow.fits')
    outflowiord_all = mrdfits(dir[i] + '../grp1.rvir.outflow_iord.fits')
    
    outflowz = outflowz_all[where(outflowz_all ne 99.0)]
    outflowm = outflowm_all[where(outflowz_all ne 99.0)]
    outflowiord = outflowiord_all[where(outflowz_all ne 99.0)]
    match,iord,outflowiord,ind_iord,ind_iordout
    coolonout = where(coolon[ind_iord] ne 0)

    outflow_z = outflowz[ind_iordout[coolonout]]
    outflow_iord = outflowiord[ind_iordout[coolonout]]
    outflow_mass = outflowm[ind_iordout[coolonout]]

    outflowt = (ageUniverse - wmap3_lookback(outflowz))/1e9

    fbase = strmid(files[i],0,strlen(files[i]) - 7)
    rtipsy,fbase,h,g,d,s,/justhead
    readarr,fbase + '.iord',h,iord,/ascii,part = 'gas'
    readarr,fbase + '.amiga.grp',h,grp,/ascii,part = 'gas'

    match,iord,outflowiord,iordind,ofiordind
    iord0 = iord[iordind];Particles once accreted that still exist
    grp0 = grp[iordind]
    outflowiord0 = outflowiord[ofiordind]
    outflowm0 =    outflowm[ofiordind]
    outflowt0 =    outflowt[ofiordind]
    outflowm00 =   outflowm0[where(grp0 ne 1)]
    outflowt00 =   outflowt0[where(grp0 ne 1)]
    
    outflowmbin = weighted_histogram(outflowt00,weight = outflowm00,min = 0.000001, max = 14.0,nbins = nbins,locations = outflowtbin)
    cumarr=fltarr(nbins,nbins)
    ii=lindgen(nbins*nbins)
    cumarr(where(ii mod nbins le ii/nbins)) = 1
    cumoutflow=reform(outflowmbin # cumarr)
    ii=0
    cumarr=0

; ------------------- Inflow -----------------------------
    inflowz = mrdfits(dir[i] + '../grp1.accrz.fits')
    inflowm = mrdfits(dir[i] + '../grp1.mass_at_accr.fits')
    inflowt = (ageUniverse - wmap3_lookback(inflowz))/1e9     
    inflowmbin = weighted_histogram(inflowt,weight = inflowm,min = 0.000001, max = 14.0,nbins = nbins,locations = outflowtbin)
    cumarr=fltarr(nbins,nbins)
    ii=lindgen(nbins*nbins)
    cumarr(where(ii mod nbins le ii/nbins)) = 1
    cuminflow=reform(inflowmbin # cumarr)
    ii=0
    cumarr=0

  
;-------------------- Bin Outflow and Inflows according to the same
;                     steps in the data file
    a_outflowm = fltarr(N_ELEMENTS(mtot))
    a_outflowm_cum = fltarr(N_ELEMENTS(mtot))
    a_inflowm = fltarr(N_ELEMENTS(mtot))
    a_inflowm_cum = fltarr(N_ELEMENTS(mtot))
    a_time = [0,time]

;IDL> for i=1,10 do begin $
;IDL> & print,i $
;IDL> & endfor

    FOR j = 0, N_ELEMENTS(mtot) - 1 DO BEGIN 
        t0 = a_time[j] 
        t1 = a_time[j + 1]  
        IF (where(outflowt00 GT t0 and outflowt00 LE t1))[0] NE -1 THEN a_outflowm[j] = total(outflowm00[where(outflowt00 GT t0 and outflowt00 LE t1)]) 
        IF (where(outflowt00 LE t1))[0] NE -1 THEN a_outflowm_cum[j] = total(outflowm00[where(outflowt00 LE t1)]) 
        IF (where(inflowt GT t0 and inflowt LE t1))[0] NE -1 THEN a_inflowm[j] = total(inflowm[where(inflowt GT t0 and inflowt LE t1)]) 
        IF (where(inflowt LE t1))[0] NE -1 THEN a_inflowm_cum[j] = total(inflowm[where(inflowt LE t1)]) 
    ENDFOR
    IF i eq 0 THEN a_array0 = [[time],[a_inflowm],[a_inflowm_cum],[a_outflowm],[a_outflowm_cum]] ELSE a_array1 = [[time],[a_inflowm],[a_inflowm_cum],[a_outflowm],[a_outflowm_cum]]

    IF i EQ 0 THEN BEGIN
        plot,outflowtbin,cuminflow*scale[i],ytitle = 'Gas  Mass [M' + sunsymbol() + ']',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = [-2.5e10,5e10],/nodata,ystyle = 9,xmargin = [17,17],ymargin = [6,6]
        axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
        oplot,[0,14],[0,0],thick = 2
    ENDIF
    oplot,outflowtbin,cuminflow*scale[i],linestyle = 2,thick = 2,color = colors[i] ;gas mass in
    oplot,outflowtbin,-1.0*cumoutflow*scale[i],linestyle = 3,thick = 2,color = colors[i];gas mass out
    oplot,outflowtbin,(cuminflow - cumoutflow)*scale[i],linestyle = 0,thick = 2,color = colors[i]
;    oplot,time,mgas*scale[i],linestyle = 2,thick = 2,color = colors[i]
;    oplot,time,mstar*scale[i],linestyle = 2,thick = 2,color = colors[i]
;    oplot,time,mstar*scale[i] + mgas*scale[i],linestyle = 0,thick = 2,color = colors[i]
;    oplot,time,cuminflow*scale[i] - (mstar*scale[i] + mgas*scale[i]),linestyle = 1,thick = 2,color = colors[i]
    IF i EQ 0 THEN outflowratio0 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]] ELSE outflowratio1 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]]
ENDFOR
legend,keys,color = colors,linestyle = [0,0],/bottom
legend,['Gas Accreted','Accreted Gas Still in Halo','Gas Lost'],linestyle = [2,0,3]

;------------------ Outflow Ratio ---------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_fbratio.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    stop
    window,10,xsize = 712,ysize = 392
endelse

FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    IF i EQ 0 THEN BEGIN
        plot,time,a_array0[*,3]/mgas,ytitle = 'Gas Mass Loss/Gas Mass',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = [0,0.8],/nodata,xmargin = [17,17],ymargin = [6,6]
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
        oplot,time,a_array0[*,3]/mgas,thick = 2,color = colors[i]
;        oplot,time,a_array0[*,3]/mcoldg,thick = 2,color = colors[i]
    ENDIF ELSE oplot,time,a_array1[*,3]/mgas,thick = 2,color = colors[i]
;oplot,time,a_array1[*,3]/mcoldg,thick = 2,color = colors[i]
ENDFOR
legend,keys,color = colors,linestyle = [0,0],/top,/right

; Mass of gas lost versus mass of gas accreted so far
;FOR i = 0, n - 1 DO BEGIN
;    IF i EQ 0 THEN BEGIN
;        plot,outflowratio0[*,0],outflowratio0[*,1],ytitle = 'Gas Mass Loss/Gas Mass Accreted',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = [0,0.8],/nodata,xmargin = [17,17],ymargin = [6,6]
;        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
;        oplot,outflowratio0[*,0],outflowratio0[*,1],thick = 2,color = colors[i]
;    ENDIF ELSE oplot,outflowratio1[*,0],outflowratio1[*,1],thick = 2,color = colors[i]
;ENDFOR
;legend,keys,color = colors,linestyle = [0,0],/bottom,/right

;------------------ FB efficiency ---------------------
if (KEYWORD_SET(outplot)) then begin
    device,/close
    device,filename = outplot + '_fbeff.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 
endif else begin
    stop
    window,12,xsize = 712,ysize = 392
endelse

FOR i = 0, n - 1 DO BEGIN
    cd,dir[i]
    readcol,filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    dtime = (time - [0,time[0:N_ELEMENTS(time) - 1]])*1e9
    IF i EQ 0 THEN BEGIN
        plot,time,a_array0[*,3]/dtime/SFR,ytitle = 'Gas Mass Lost Per Year/SFR',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = [0,100],/nodata,xmargin = [17,17],ymargin = [6,6]
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = N_ELEMENTS(ticktime_in_t) - 1        
        oplot,time,a_array0[*,3]/dtime/SFR,thick = 2,color = colors[i],psym = 10
    ENDIF ELSE oplot,time,a_array1[*,3]/dtime/SFR,thick = 2,color = colors[i],psym = 10
ENDFOR
legend,keys,color = colors,linestyle = [0,0],/top,/right
if (KEYWORD_SET(outplot)) then device,/close ELSE stop
end

pro mass_evol_master

dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']

dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']

dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/',$
;         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']
files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512.halo.1',$
;         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
scale = [1,8]
pfiles = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param',$
;         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
maxdistance_photo = 6

;dir   = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/']
;files = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
;scale = [1]
;pfiles = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']


;densRadius,[],[50000.,25000],[1.84793e16,2.310e15],maxdistance = maxdistance_photo ;,outplot = outplot
;mass_evol,dir,files,pfiles,scale = scale,linestyle = [0,0],psym = [14,14],/color,outplot = '~/plots/twins/h603_h516_evol.talk'
mass_evol,dir,files,pfiles,scale = scale,linestyle = [0,0],psym = [14,14],/color,keys = ['Spiral'+textoidl(' Galaxy'),'Dwarf Galaxy'],outplot = '~/plots/twins/h603_h516_evol_flash'
end
