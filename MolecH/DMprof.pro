

pro DMprof, files, dMsolUnit, dKpcUnit, outfile = outfile, keys = keys, color = color,thicks = thicks, linestyle = linestyle,yrange = yrange,maxdistance = maxdistance,vfinal = vfinal,type = type,multiframe = multiframe, label = label,ctables = ctables,formatthick = formatthick
; type: gas, dark, star,all

;*********************************** Units ****************************************
hubble = 73.0
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
sec_per_year = 31556926
grav = 6.67e-8

n = N_ELEMENTS(files)
nbins = 100.0
IF NOT KEYWORD_SET(maxdistance) THEN maxdistance = 8.0
bind = maxdistance/nbins
radius = (findgen(nbins)+1)*bind 
tmass = fltarr(nbins)   
tmassd= tmass
tmassg = tmass
tmasss = tmass
formatplot,outplot = outfile,thick = formatthick
IF KEYWORD_SET(outfile) THEN     device,filename=outfile,/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 ELSE window,0
IF KEYWORD_SET(outfile) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
   IF KEYWORD_SET(multiframe) THEN BEGIN
        xsize = 10*n
        ysize = 12
        mxTitSize = 1.5
        mxTitOffset = 2
    ENDIF ELSE BEGIN
        xsize = 18
        ysize = 18
    ENDELSE
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    IF KEYWORD_SET(multiframe) THEN BEGIN
        xsize = 400*n
        ysize = 475
        mxTitSize = 1.5
        mxTitOffset = 1
    ENDIF ELSE BEGIN
        xsize = 400
        ysize = 400
    ENDELSE
ENDELSE
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if NOT keyword_set(ctables) then ctables = [39,39,39]
    if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) ;REVERSE(findgen(n)*2)
ENDIF ELSE BEGIN
    loadct,0    
    if NOT keyword_set(ctables) then ctables = [0,0,0]
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT KEYWORD_SET(thicks) THEN thicks = (findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = REVERSE(findgen(n)*2)   
ENDELSE
;IF not keyword_set(yrange) then yrange = [0,300]

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
lb_linestyle = intarr(n)

FOR ifile = 0, n - 1 DO BEGIN
    rtipsy,files[ifile],h,g,d,s
    loadct,ctables[ifile]
    a = h.time
    IF (N_ELEMENTS(dKpcUnit) GT 1) THEN kpc_per_syslength = dKpcUnit[ifile]*a ELSE kpc_per_syslength = (dKpcUnit)[0]*a
    IF (N_ELEMENTS(dMsolUnit) GT 1) THEN msol_per_sysmass = dMsolUnit[ifile] ELSE msol_per_sysmass = dMsolUnit[0]
    timeunit = SQRT((cm_per_kpc*kpc_per_syslength)^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year

    distances = SQRT(d.x^2.0 + d.y^2.0 + d.z^2.0)*kpc_per_syslength
    mass = d.mass*msol_per_sysmass

    FOR i = 1, nbins-1 DO BEGIN
        ind = where(distances LE radius[i] AND distances GE radius[i] - bind)
        if (ind[0] ne -1)then tmass[i] = TOTAL(mass[ind])else tmass[i] = 0
    ENDFOR
    rho = tmass/radius^3/1000.0^3

    IF KEYWORD_SET(multiframe) THEN BEGIN
        multiplot,[n,1],/square,mxtitle = 'Radius [kpc]',mxTitSize = mxTitSize, mxTitOffset = mxTitOffset
        plot,radius - bind/2.0,rho,thick = thicks[iFile],linestyle = linestyle[iFile],xtitle = 'Radius [kpc]',ytitle = textoidl('\rho_{DM} [M')+sunsymbol()+ textoidl('pc^{-3}]'),xrange = [0,maxdistance],/ylog
    ENDIF ELSE BEGIN
        IF ifile eq 0 THEN BEGIN
            plot,radius-bind/2.0,rho,thick = thicks[iFile],linestyle = linestyle[iFile],xtitle = 'Radius [kpc]',ytitle = textoidl(' \rho_{DM} [M')+sunsymbol()+ textoidl('pc^{-3}]'),xrange = [0.1,maxdistance],/nodata,title = label,/ylog,/xlog,yrange = [0.001,10.0],ytickname = [textoidl('10^{-3}'),textoidl('10^{-2}'),textoidl('10^{-1}'),textoidl('1'),textoidl('10')]
            eps = MIN(d.eps)*kpc_per_syslength
            oplot,[2.0*eps,2.0*eps],[0.001,100],linestyle = 1
        ENDIF
        oplot,radius-bind/2.0,rho,thick = thicks[iFile],linestyle = linestyle[iFile],color = colors[iFile]

;        IF KEYWORD_SET(keys) THEN BEGIN
;            l_keys = lb_keys
;            l_keys[iFile] = keys[iFile]
;            l_thicks = lb_thicks
;            l_thicks[iFile] = thicks[iFile]
;            l_color = lb_color
;            l_color[iFile] = colors[iFile]
;            l_linestyle = lb_linestyle
;            l_linestyle[iFile] = l_linestyle[iFile]
;            legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/right,box = 0
;        ENDIF
    ENDELSE
        IF KEYWORD_SET(type) THEN $
        IF KEYWORD_SET(multiframe)AND i ne n - 1 THEN multiplot
    ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,ctables = ctables, /right,/top,box = 0
IF KEYWORD_SET(outfile) THEN device,/close else stop
IF KEYWORD_SET(multiframe) THEN multiplot,/reset
END
