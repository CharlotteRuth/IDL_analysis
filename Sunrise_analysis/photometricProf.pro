pro photometricProf, files, filternums = filternums, cameras = cameras, bands = bands, outplot = outplot, keys = keys, colors = colors, thicks = thicks, linestyle= linestyle, maxdistance = maxdistance,yrange = yrange,xrange = xrange, recenter = recenter, label = label,ctables = ctables,formatthick = formatthick
;*********************************** Units ****************************************
n = N_ELEMENTS(files)
If NOT KEYWORD_SET(maxdistance) THEN maxdistance = 8.0
IF not keyword_set(filternums) THEN filternums = intarr(n) + 13    ;index number of filters in fv 
IF not keyword_set(cameras) THEN cameras = intarr(n) + 14     ;index number of desired camera angle in fv
IF NOT keyword_set(bands) THEN bands = intarr(n) + 4     ;filter number starting from zero ;i_SDSS
IF NOT keyword_set(recenter) THEN recenter = intarr(n)

IF N_ELEMENTS(filternums) eq 1 THEN filternums = intarr(n) + filternums
IF N_ELEMENTS(cameras) eq 1 THEN cameras = intarr(n) + cameras
IF N_ELEMENTS(bands) eq 1 THEN bands = intarr(n) + bands
IF N_ELEMENTS(recenter) eq 1 THEN recenter = intar(n) + recenter

formatplot,outplot = outplot,thick = formatthick
IF KEYWORD_SET(outplot) THEN     device,filename=outplot,/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 ELSE window,0
IF KEYWORD_SET(outplot) THEN     fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(outplot) THEN     bgcolor = 255 ELSE bgcolor = 0
IF KEYWORD_SET(colors) THEN BEGIN
    loadct,39
    if colors[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = colors
    if NOT keyword_set(ctables) then ctables = [39,39,39]
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT KEYWORD_SET(symbols) THEN symbols = (fltarr(n)+4)
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    if NOT keyword_set(ctables) then ctables = [0,0,0]
    IF NOT KEYWORD_SET(thicks) THEN thicks = (findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = REVERSE(findgen(n)*2)   
    IF NOT KEYWORD_SET(symbols) THEN symbols = (findgen(n)+2)*2
ENDELSE

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
lb_linestyle = intarr(n)

FOR ifile = 0, n - 1 DO BEGIN
    filternum = filternums[ifile]
    camera = cameras[ifile]
    band = bands[ifile]

;    rtipsy,files[ifile],h,g,d,s,/justhead
;    a = h.time

    temp = mrdfits(files[ifile],2,header)
    fovstr = strsplit(header[13],' ',/extract)
    fov = FLOAT(fovstr[3])

;------------------------------------------------------------------------------
;Get zero point for AB magnitude system
    filters= mrdfits(files[ifile],filternum)
    Lunit = filters[band].ewidth_lambda/filters[band].ewidth_nu ; internal units
;to W/m^2 for Sunrise3
;Lunit = filters[band].l_lambda_to_l_nu; internal units to W/m^2
    Fo=3.631e-23          ;zero point of AB magnitude system, in W/m^2
    units = 4.25e10             ; sr>>> arcsec^2
    sbfactor = Lunit/units/Fo
;ABzero = -2.5*alog10(sbfactor)

;deal with the Sunrise'd image first
    image = mrdfits(files[ifile],camera)
;-----------------------------------------------------------------
;SB profile of the image, using circular annuli
    gal = image(*,*, band)
    makex, gal,x,y              ;set up coordinate grid
    npix = N_ELEMENTS(gal[0,*])
    pixKpc = 1.0*fov/(1.0*npix)
    dummy = max(gal,nn)
    IF (recenter[ifile]) THEN BEGIN             ;center at brightest pixel
        max=array_indices(gal,nn)-npix/2
        x=x-max[0]
        y=y-max[1]
        print,'Center: ',max[0],max[1]
    ENDIF 
    radius = sqrt((x)^2+(y)^2)  ; ciruclar annuli
    binsize = 2.0
;    hist = hist1d(radius,gal,obin=xvals,binsize=binsize,min = 0, max = ceil(max(radius))) 
;    numh = hist1d(radius,binsize=2)
    hist = weighted_histogram(radius,weight = gal,locations=xvals,binsize=binsize,min = 0, max = ceil(max(radius))) 
    numh = weighted_histogram(radius,binsize=binsize)

    SBprofile = hist/numh
    rad_kpc = xvals*pixKpc
    err = SBprofile/sqrt(numh)
    profile = -2.5*alog10(sbfactor*SBprofile)
    loadct,ctables[ifile]
    IF NOT KEYWORD_SET(yrange) THEN yrange = [max(profile[where(rad_kpc lt maxdistance)])+0.8,min(profile[where(rad_kpc lt maxdistance)])-0.8]
    IF ifile eq 0 THEN plot,rad_kpc,profile,xrange=[0,maxdistance],/xstyle,xtit='Radius [kpc]',ytit=textoidl('M/arcsec^2'),yrange=yrange,thick = thicks[iFile],linestyle = linestyle[iFile],/nodata, title = label
    oplot,rad_kpc,profile,thick = thicks[iFile],linestyle = linestyle[iFile],color = colors[iFile]
;    IF KEYWORD_SET(keys) THEN BEGIN
;        l_keys = lb_keys
;        l_keys[iFile] = keys[iFile]
;        l_thicks = lb_thicks
;        l_thicks[iFile] = thicks[iFile]
;        l_color = lb_color
;        l_color[iFile] = colors[iFile]
;        l_linestyle = lb_linestyle
;        l_linestyle[iFile] = l_linestyle[iFile]
;        legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/right,box = 0
;    ENDIF
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thicks,ctables = ctables,/right,/top,box =0
IF KEYWORD_SET(outplot) THEN device,/close else stop
END
