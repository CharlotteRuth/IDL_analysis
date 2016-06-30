;filt = 13
;cam = 14

PRO color_prof, files, filt = filt, cam = cam,outplot = outplot,title = title,key = key,color = color,linestyle = linestyle,thick = thick,dir = dir, r25 = r25,maxdistance = maxdistance

n = N_ELEMENTS(files)
!p.multi = 0
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    black  = 0
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    black = 255
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  color = (findgen(n) + 1)*240/n else colors = color
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 1
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 1 ;findgen(n)*2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = findgen(n)*2   
ENDELSE

if (KEYWORD_SET(outplot)) then begin
    set_plot,'ps'
    device,filename = outplot + '_colorprof.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
endif else begin
    set_plot,'x'
    window,1,xsize = 712,ysize = 318
endelse

;***UPDATE THESE FOR EACH GALAXY***
IF not keyword_set(filt) THEN filt = fltarr(n) + 12    ;index number of filters in fv 
IF not keyword_set(cam) THEN cam = fltarr(n) + 14     ;index number of desired camera angle in fv

;region to fit in (kpc)
rmin =  0
IF KEYWORD_SET(maxdistance) THEN  rmax = maxdistance ELSE rmax = 4 ;8.0

Bband = 21
Vband = 22

FOR i = 0, n-1 DO BEGIN
    if keyword_set(dir) THEN cd,dir[i]
;***END OF USER INPUT***
    temp = mrdfits(files[i],2,header)
    fovstr = strsplit(header[13],' ',/extract)
    fov = FLOAT(fovstr[3])
    
;------------------------------------------------------------------------------
;Get zero point for AB magnitude system
    filters= mrdfits(files[i],filt[i])
    Lunit_B = filters[Bband].ewidth_lambda/filters[Bband].ewidth_nu ; internal units
    Lunit_V = filters[Vband].ewidth_lambda/filters[Vband].ewidth_nu ; internal units
;to W/m^2 for Sunrise3
;Lunit = filters[band].l_lambda_to_l_nu; internal units to W/m^2
    Fo=3.631e-23          ;zero point of AB magnitude system, in W/m^2
    units = 4.25e10             ; sr>>> arcsec^2
    sbfactor_B = Lunit_B/units/Fo
    sbfactor_V = Lunit_V/units/Fo
    
;deal with the Sunrise'd image first
    image = mrdfits(files[i],cam[i])
;-----------------------------------------------------------------
;SB profile of the image, using circular annuli
    gal_B = image(*,*, Bband)
    gal_V = image(*,*, Vband)
    makex, gal_B, x,y           ;set up coordinate grid
    npix = N_ELEMENTS(gal_B[0,*])
    pixKpc = 1.0*fov/(1.0*npix)
    dummy = max(gal_B,nn)
    
    
;for circular annuli
    radius = sqrt((x)^2+(y)^2)  ; ciruclar annuli
    hist_B = weighted_histogram(radius,weight = gal_B,locations=xvals,binsize=2,min = 0, max = ceil(max(radius))) 
    hist_V = weighted_histogram(radius,weight = gal_V,locations=xvals,binsize=2,min = 0, max = ceil(max(radius)))
    numh = weighted_histogram(radius,binsize=2)
    SBprofile_B = hist_B/numh
    SBprofile_V = hist_V/numh
    rad_kpc = xvals*pixKpc
    err_B = SBprofile_B/sqrt(numh)
    err_V = SBprofile_V/sqrt(numh)
    profile_B = -2.5*alog10(sbfactor_B*SBprofile_B)
    profile_V = -2.5*alog10(sbfactor_V*SBprofile_V)
    
    IF i eq 0 THEN plot,rad_kpc,profile_B - profile_V,xtit='Radius [kpc]',ytitle = 'B - V',xrange = [rmin,rmax],yrange = [0,0.8],title = title,linestyle = linestyle[i],thick = thick[i]
    oplot,rad_kpc,profile_B - profile_V,linestyle = linestyle[i],thick = thick[i],color = colors[i]
    if keyword_set(r25) THEN BEGIN
        r25_val = opticalRadii(B_num = Bband)
        oplot,[r25_val,r25_val],[0,0.8],linestyle = linestyle[i],thick = thick[i],color = colors[i]
    ENDIF
ENDFOR
if KEYWORD_SET(key) THEN legend,key,color = colors,linestyle = linestyle,thick = thick,/right,/bottom
if (KEYWORD_SET(outplot)) then device,/close ELSE stop
END


PRO color_prof_master
prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
dir =   prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.1',$
         'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.1'] 
files = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.1/broadband.fits',$
         'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.1/broadband.fits'] 
key = ['DnoH2','DH2']
outplot = '~/plots/h516.cosmo25cmb.paper.colorprof.eps'
filt = [13,13]
color_prof, files,key = key,filt = filt,dir = dir,/r25,outplot = outplot,linestyle = [2,0]
END
