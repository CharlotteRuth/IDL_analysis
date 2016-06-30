;example of fitting an edge-on galaxy using the markwardt package and
;        the code in this directory

;file ='/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512/h603.cosmo50cmb.2304g2bwK.00512.1/broadband.fits'
;file = '/astro/net/scratch1/abrooks/FABIO/h516.cosmo25cmb.2304g2bwK/h516.cosmo25cmb.2304g2bwK.00512/h516.cosmo25cmb.2304g2bwK.00512.1/broadband.fits'


;dir ='/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/steps/h516.cosmo25cmb.1536g2MBWK.00512.dir/h516.cosmo25cmb.1536g2MBWK.00512.1/'
;file = 'broadband.fits'
;camera = 14 ;Extincted, face on
;band = 22 ;V   ; 5 ;SDSS i
;filter = 13
;fitname = 'exp'

;dir ='/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK.00512.1/'
;file = 'broadband.fits'
;camera = 14 ;Extincted, face on
;band = 22 ;V   ; 5 ;SDSS i
;filter = 13
;profilefit, file, filter = filter, camera = camera, band = band,/write_feedme,fitname = fitname,mask = 1

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/h516.cosmo25cmb.2304g14HBWK.00512.1'
;file = 'broadband.fits'
;cam = 16 ;Extincted, face on
;band = 21 ; B
;filt = 13

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.1'
;file = 'broadband.fits'
;camera = 14 ;Extincted, face on
;band = 21 ; B
;band = 22 ; V
;filter = 13
;yrange = [32,20]
;fitname = 'exp'

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00324.dir/h603.cosmo50cmb.3072g14HBWK.00324.1'
;file = 'broadband.fits'
;camera = 14 ;Extincted, face on
;band = 22 ; V
;band = 21; B
;band = 4; SDSS r, Ferah claims that R band is standard but not
;included on the standard set of sunrise filters
;filter = 13
;yrange = yrange

;24, 7

;dir = '/astro/store/nbody3/christensen/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.1'
;file = 'broadband.fits'
;camera = 14 ;FO with dust
;band = 7; 2MASS H
;filter = 13

;dir = '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.00512.dir/h986.cosmo50cmb.3072gs1MbwK.00512.1'
;file = 'broadband.fits'
;camera = 14 ;FO with dust
;band = 7; 2MASS H
;filter = 13


;dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/steps/h239.cosmo50cmb.3072g14HMbwK.00512.dir/h239.cosmo50cmb.3072g14HMbwK.00512.1'
;file = 'broadband.fits'
;camera = 14 ;FO with dust
;band = 7; 2MASS H
;filter = 13

;fitname = 'exp','sersic','double'

pro profilefit, file, filter = filter, camera = camera, band = band, fitname = fitname, mask = mask,yrange = yrange,rmax = rmax,write_feedme = write_feedme,outplot = outplot,prefix = prefix,recenter = recenter;, fov = fov, npix = npix,

;***UPDATE THESE FOR EACH GALAXY***

;IF not keyword_set(fov) THEN fov = 20.     ; field of view in kpc, h603
;IF not keyword_set(npix) THEN npix = 800.   ;number of pixels on each side of the image
IF not keyword_set(filter) THEN filter = 13    ;index number of filters in fv 
IF not keyword_set(camera) THEN camera = 16     ;index number of desired camera angle in fv
IF NOT keyword_set(band) THEN band = 5     ;filter number starting from zero ;i_SDSS
IF NOT keyword_set(fitname) THEN fitname = 'double'
IF NOT keyword_set(mask) THEN mask = 0
IF NOT keyword_set(prefix) THEN prefix = './'

;Initial guesses for fitting parameters of the form:
;init=[Disk Surface Brightness, Disk Scale Length, Bulge Surface Brightness, Bulge Scale Length, Sersic Index]

;***END OF USER INPUT***

temp = mrdfits(file,2,header)
fovstr = strsplit(header[13],' ',/extract)
fov = float(fovstr[3]) ;FOV in kpc
diststr = strsplit(header[12],' ',/extract)
dist = float(diststr[3]) ;Camera distance in kpc

;------------------------------------------------------------------------------
;Get zero point for AB magnitude system
filters= mrdfits(file,filter)
Lunit = filters[band].ewidth_lambda/filters[band].ewidth_nu ;for conversion between W/m^2/m and W/m^2/Hz
Fo=3.631e-23 ;zero point of AB magnitude system for every filter, in W/m^2/Hz
angleunits = 4.25e10 ; sr to arcsec^2
distconv = (1000.*dist/10.)^2 ; To calculate what flux would be at 10 parsecs to get absolute magnitude
sbfactor = Lunit/Fo/angleunits;*distconv

;deal with the Sunrise'd image first
image = mrdfits(file,camera)
;-----------------------------------------------------------------
;SB profile of the image, using circular annuli
gal = image(*,*, band) ;W/m/m^2/s
makex, gal,x,y  ;set up coordinate grid
npix = n_elements(gal[0,*])
kpcpix = 1.0*fov/(1.0*npix);kpc/pix
arcsecpix = atan(fov/dist)*360/(2*!PI)*3600/npix ;arcsec/pix
dummy = max(gal,nn)
IF 0 THEN BEGIN ;center at brightest pixel
    max=array_indices(gal,nn)-npix/2
    x=x-max[0]
    y=y-max[1]
ENDIF ELSE BEGIN
;    x=x-npix/2
;    y=y-npix/2
ENDELSE

;for fitting along a slit
;radius=abs(y)
;slit=where(abs(x) lt 10)
;hist = hist1d(radius[slit],gal[slit],obin=xvals,binsize=10) 
;numh = hist1d(radius[slit],binsize=10)

;for circular annuli
IF keyword_set(recenter) THEN BEGIN
   xhist = weighted_histogram(x,weight = gal,locations = xloc)
   yhist = weighted_histogram(y,weight = gal,locations = yloc)
   temp = max(xhist,xcenind)
   temp = max(yhist,ycenind)
   xcen = xloc[xcenind]
   ycen = yloc[ycenind]
   x = x - xcen
   y = y - ycen
ENDIF ELSE BEGIN
   xcen = 0
   ycen = 0
ENDELSE

radius = sqrt((x)^2+(y)^2) ; ciruclar annuli
binsize = 2.0
medflux = fltarr(ceil(max(radius))/binsize)
meanflux = fltarr(ceil(max(radius))/binsize)
rad_kpcM = indgen(n_elements(medflux))*binsize + binsize/2.0
nbins = n_elements(medflux)
for i=0L,nbins-2 do begin
   rminM = rad_kpcM[i]-(binsize/2.)
   rmaxM = rad_kpcM[i]+(binsize/2.)
   ind = where(radius gt rminM and radius le rmaxM, nind)
;  print, "Number of elements at radius", rad_kpcM[i],n_elements(ind)
;  print, "Median of elements at radius",rad_kpcM[i],median(gal[ind])
   medflux[i] = median(gal[ind])
   meanflux[i] = mean(gal[ind])
;  print,"rad_kpc value is", rad_kpc[i]
ENDfor
rad_kpcM = rad_kpcM * kpcpix
profile_med = -2.5*alog10(sbfactor*medflux)
profile_mean = -2.5*alog10(sbfactor*meanflux)

hist = hist1d(radius,gal,obin=xvals,binsize=binsize,min = 0, max = ceil(max(radius)))
;err = hist
;FOR i = 0L, nbins - 1 DO BEGIN
;    ind = where(radius ge xvals[i] - binsize/2 AND radius lt xvals[i] + binsize/2)
;    err[i] = stdev(gal[ind])
;ENDFOR 

numh = hist1d(radius,binsize=binsize)
SBprofile = hist/(numh)
rad_kpc = xvals*kpcpix
err = SBprofile/sqrt(numh)
profile = -2.5*alog10(sbfactor*SBprofile) ;M/arcsec^2

SBprofile=medflux
profile = profile_med

;region to fit in (kpc)
rmin = mask
IF keyword_set(rmax) THEN $
  IF rmax LT 0 THEN rmax = npix*kpcpix/2
IF NOT keyword_set(rmax) THEN BEGIN
    rmax = opticalRadii(extno = camera,B_num = band,/verbose)
ENDIF
print,'R_max: ',rmax

;Region to plot (kpc)
xmin = mask
xmax = rmax ;8.0

region = where(rad_kpc gt rmin and rad_kpc lt rmax)
IF (KEYWORD_SET(mask)) THEN BEGIN
    xbadpix = where(radius*kpcpix lt rmin) MOD npix
    ybadpix = where(radius*kpcpix lt rmin)/npix
;    xbadpix = where(radius*kpcpix gt rmin and radius*kpcpix lt rmax) MOD npix
;    ybadpix = where(radius*kpcpix gt rmin and radius*kpcpix lt rmax)/npix
    IF (xbadpix[0] ne -1) THEN writecol,'mask.asc',xbadpix,ybadpix
ENDIF

;----------------------- make the fit 
zeropoint = -2.5*alog10(sbfactor*distconv)
IF fitname EQ 'double' THEN begin
;p=[Disk csb, R_s, bulge csb, R_e, n)
;init = [0.39507740, 5.6588991, 5, 0.5, 2.529] ;h603
    init = [3.3842650, 1.0872235, 35.089366, 0.18417694, 2.5 ];1.5] ;h516
    parinfo=replicate({fixed:0},5) 
    parinfo.fixed=[0,0,0,0,1]
    fits = mpfitfun('exp_sersic', rad_kpc[region], SBprofile[region], err[region],init, parinfo=parinfo,/quiet)

    diskfit = abs(fits[0])*exp(-rad_kpc/fits[1])
    diskplot = -2.5*alog10(sbfactor*diskfit)
    k = 1.9992* fits[4] - 0.3271
    bulgefit = abs(fits[2])*exp(-k*((rad_kpc/fits[3])^(1/fits[4])-1))
    bulgeplot = -2.5*alog10(sbfactor*bulgefit)
    fit = bulgefit+diskfit
    fitplot = -2.5*alog10(sbfactor*fit)
    csb = strnsignIF(-2.5*alog10(sbfactor*fits[0]),3)
    scalelength = strnsignIF(fits[1],3)
    u_e= strnsignIF(-2.5*alog10(sbfactor*fits[2]),3)
    r_e = strnsignIF(fits[3],3)
    n= strnsignIF(fits[4],2)
ENDIF
IF fitname EQ 'exp' THEN begin
;p is of form, [disk csb, disk Rl]
;init = [0.39507740, 0.5] ;h603
    init = [3.3842650, 0.18417694] ;h516
    parinfo=replicate({fixed:0},2) 
    parinfo.fixed=[0,0]
    fits = mpfitfun('expon', rad_kpc[region], SBprofile[region], err[region], init, parinfo=parinfo,/quiet)

    diskfit = abs(fits[0])*exp(-rad_kpc/fits[1])
    diskplot = -2.5*alog10(sbfactor*diskfit)
    fit = diskfit
    fitplot = -2.5*alog10(sbfactor*fit)
    csb = strnsignIF(-2.5*alog10(sbfactor*fits[0]),3)
    scalelength = strnsignIF(fits[1],3)
ENDIF
IF fitname EQ 'sersic' THEN begin
;p is of form, [bulge csb, bulge Re, bulge n]
;init = [ 5, 5.6588991, 2.529] ;h603
    init = [35.089366, 1.0872235, 1.5] ;h516
    parinfo=replicate({fixed:0},3) 
    parinfo.fixed=[0,0,1]
    fits = mpfitfun('sersic', rad_kpc[region], SBprofile[region], err[region], init, parinfo=parinfo,/quiet)

    k = 1.9992* fits[2] - 0.3271
    bulgefit = abs(fits[0])*exp(-k*((rad_kpc/fits[1])^(1/fits[2])-1))
    bulgeplot = -2.5*alog10(sbfactor*bulgefit)
    fit = bulgefit
    fitplot = -2.5*alog10(sbfactor*fit)
    u_e= strnsignIF(-2.5*alog10(sbfactor*fits[0]),3)
    r_e = strnsignIF(fits[1],3)
    n= strnsignIF(fits[2],2)
ENDIF
IF fitname EQ 'double' THEN BEGIN
    k = 1.9992* fits[4] - 0.3271
    K_const = 2*!pi*(fits[3]*arcsecpix*npix/fov)^2*exp(k)*k^(-2*fits[4])*Gamma(2*fits[4])
    uo=fits[0]*sbfactor
ENDIF

;------------------ Plot fits ---------------
;region = where(rad_kpc gt rmin and rad_kpc lt rmax)
IF NOT KEYWORD_SET(yrange) THEN BEGIN
    ymin = max(profile[where(rad_kpc LE xmax AND finite(profile))])+0.8
    ymax = min(profile[where(rad_kpc LE xmax AND finite(profile))])-0.8
    
;    ymin = max(profile[region])+0.8
;    ymax = min(profile[region])-0.8
ENDIF ELSE BEGIN
    ymin = yrange[0]
    ymax = yrange[1]
ENDELSE
loadct,39
formatplot,outplot=outplot
IF keyword_set(outplot) THEN BEGIN
    set_plot,'ps'
    device,/portrait,filename='./profilefit_'+strtrim(filters[band].filter,2)+'.eps', /color,bits_per_pixel=8
ENDIF
plot,rad_kpc,profile,xrange=[0,xmax],/xstyle,xtit='Radius (kpc)',ytit=textoidl('M/arcsec^2 [AB]'),yrange=[ymin,ymax],/ystyle
IF keyword_set(mask) THEN BEGIN
    oplot,[mask,mask],[ymin,ymax],linestyle = 1
    oplot,rad_kpc[region],profile[region],thick = 6,color = 254
ENDIF
oplot,rad_kpc,fitplot,linestyle = 3,color = 50
oplot,[rmax,rmax],[ymin,ymax],linestyle = 1
IF (fitname EQ 'sersic' OR  fitname EQ 'double') THEN oplot,rad_kpc,bulgeplot,linestyle=1,color=254,thick = 6
IF (fitname EQ 'exp'    OR  fitname EQ 'double') THEN oplot,rad_kpc,diskplot,linestyle=2,color=150,thick = 6
IF (fitname EQ 'double') THEN BEGIN
    legend, ['Sersic', 'Expdisk', 'Sersic+Expdisk'], linestyle=[1,2,3],color=[250,150,50],/right ;,position=[18,10]
;    legend,[textoidl('R_S') + ':           '+strtrim(fits[3],2),'B/D:         '+strtrim(10^(-(zeropoint - 2.5*alog10(fits[2]*K_const)/2.5))/10^(-(zeropoint - 2.5*alog10(fits[0]*2*!pi*(npix/fov*fits[1])^2)/2.5)),2),'Central SB: '+strtrim(-2.5*alog10(uo),2)],/bottom,/left
ENDIF
IF keyword_set(outplot) THEN device,/close

;Mean Surface Brightness
bulge_ind = where(radius*kpcpix LE fits[3])
meanSB = -2.5*alog10(mean(sbfactor*gal[bulge_ind]))

print,'Sunrise Absolute Magnitude:    ',filters[band].ab_mag0,filters[band].ab_mag0-1.39
print,'Calculated Absolute Magnitude: ',-2.5*alog10(total(gal*sbfactor*arcsecpix^2)*distconv),-2.5*alog10(total(gal*sbfactor*arcsecpix^2)*distconv)-1.39 ;M
IF (fitname EQ 'double') THEN BEGIN
    mb = zeropoint - 2.5*alog10(fits[2]*K_const)
    md = zeropoint - 2.5*alog10(fits[0]*2*!pi*(arcsecpix*npix/fov*fits[1])^2)
    print,'Total Mag:    ',strtrim(-2.5*alog10(10^(-mb/2.5) + 10^(-md/2.5)),2),', ',strtrim(-2.5*alog10(10^(-(mb - 1.39)/2.5) + 10^(-(md - 1.39)/2.5)),2)
    print,'Bulge Mag:    ',strtrim(mb,2),', ',strtrim(mb-1.39,2)
    print,'Mean SB:      ',strtrim(meanSB,2),', ',strtrim(meanSB-1.39,2)
;strtrim(min(bulgeplot),2),', ',strtrim(min(bulgeplot)-1.39,2)
    print,'log(R_e):     ',strtrim(alog10(fits[3]*1000),2)
    print,'Sersic Index: ',n
    print,'Disk Mag:     ',strtrim(md,2),', ',strtrim(md - 1.39,2)
    print,'log(R_h):     ',strtrim(alog10(fits[1]*1000),2)
    print,'B/D:          ',strtrim(10^(-1*mb/2.5)/10^(-1*md/2.5),2)
    print,''
    flux=10^(-1*profile/2.5)
    intmag = -2.5*alog10(total(2*!PI*rad_kpc*(binsize*kpcpix)*arcsecpix^2/kpcpix^2*flux)*distconv)
    print,'Integrated Profile Mag:    ',strtrim(intmag,2),', ',strtrim(intmag - 1.39,2)
    flux=10^(-1*fitplot/2.5)
    fitmag = -2.5*alog10(total(2*!PI*rad_kpc*(binsize*kpcpix)*arcsecpix^2/kpcpix^2*flux)*distconv)
    print,'Integrated Profile Fit Mag:',strtrim(fitmag,2),', ',strtrim(fitmag - 1.39,2)
    flux=10^(-1*bulgeplot/2.5)
    bulgemag = -2.5*alog10(total(2*!PI*rad_kpc*(binsize*kpcpix)*arcsecpix^2/kpcpix^2*flux)*distconv)
    print,'Integrated Bulge Fit Mag:  ',strtrim(bulgemag,2),', ',strtrim(bulgemag - 1.39,2)  
    flux=10^(-1*diskplot/2.5)
    diskmag = -2.5*alog10(total(2*!PI*rad_kpc*(binsize*kpcpix)*arcsecpix^2/kpcpix^2*flux)*distconv)
    print,'Integrated Disk Fit Mag:   ',strtrim(diskmag,2),', ',strtrim(diskmag - 1.39,2)
ENDIF

IF(keyword_set(write_feedme)) THEN BEGIN
    filename = 'gal.'+strtrim(filters[band].filter,2)+'.fits';  ELSE filename = 'gal.fits'
    mwrfits,gal*arcsecpix^2,prefix + filename,/create
    writefeedme,fits,fov,npix,arcsecpix,rmax,-2.5*alog10(Lunit/Fo/angleunits*distconv),filtname = strtrim(filters[band].filter,2),mask = mask, fitname = fitname, diskmag = diskmag,prefix = prefix,xcen = xcen,ycen = ycen
ENDIF
END


pro writefeedme,fits,fov,npix,arcsecpix,rmax,zeropoint, filtname = filtname, fitname = fitname, mask = mask, bulgemag = bulgemag, diskmag = diskmag,prefix = prefix,xcen = xcen, ycen = ycen
IF NOT keyword_set(xcen) THEN xcen = 0
IF NOT keyword_set(ycen) THEN ycen = 0

kpcpix = 1.0*fov/(1.0*npix);kpc/pix
IF KEYWORD_SET (filtname) THEN BEGIN
    filename = prefix + 'gal.'+filtname+'.fits' 
    openw,1,prefix + 'galfit_'+filtname+'.feedme'
ENDIF ELSE BEGIN 
    filename = prefix + 'gal.fits'
    openw,1,prefix + 'galfit.feedme'
ENDELSE
printf,1,'#  Chi^2/nu = 0.010,  Chi^2 = 76.281,  Ndof = 7910'
printf,1,''
printf,1,'================================================================================'
printf,1,'# IMAGE and GALFIT CONTROL PARAMETERS'
IF KEYWORD_SET (filtname) THEN printf,1,'A) '+filename+'            # Input data image (FITS file)' $
ELSE printf,1,'A) '+filename+'            # Input data image (FITS file)' 
printf,1,'B) imgblock.'+filtname+'.fits       # Output data image block'
printf,1,'C) none                # Sigma image name (made from data IF blank or "none") '
printf,1,'D) psf.fits   #        # Input PSF image and (optional) dIFfusion kernel'
printf,1,'E) 1                   # PSF oversampling factor relative to data '
IF KEYWORD_SET (mask) THEN printf,1,'F) mask.asc            # Bad pixel mask (FITS image or ASCII coord list)' $
ELSE printf,1,'F) mask.asc            # Bad pixel mask (FITS image or ASCII coord list)'
printf,1,'G) none                # File with parameter constraints (ASCII file)'
printf,1,'H) '$
       ,strtrim(fix(npix/2.-npix*rmax/fov + xcen),2),' '$
       ,strtrim(fix(npix/2.+npix*rmax/fov + xcen),2),' '$
       ,strtrim(FIX(npix/2.-npix*rmax/fov + ycen),2),' '$
       ,strtrim(FIX(npix/2.+npix*rmax/fov + ycen),2)$
       ,'   # Image region to fit (xmin xmax ymin ymax)'
printf,1,'I) ',strtrim(FIX(npix*rmax/fov*2),2),'   ',strtrim(FIX(npix*rmax/fov*2),2)$
       ,'         # Size of the convolution box (x y)'
printf,1,'J) ',strtrim(zeropoint,2)$
       ,'              # Magnitude photometric zeropoint '
printf,1,'K) '$
       ,strtrim(npix*rmax/fov*2*arcsecpix,2),'  '$
       ,strtrim(npix*rmax/fov*2*arcsecpix,2)$
       ,'        # Plate scale in arcseconds (dx dy) '
printf,1,'O) regular             # Display type (regular, curses, both)'
printf,1,'P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps'
printf,1,''
printf,1,'# INITIAL FITTING PARAMETERS'
printf,1,'#'
printf,1,'#   For object type, the allowed functions are:' 
printf,1,'#       nuker, sersic, expdisk, devauc, king, psf, gaussian, and moffat.'
printf,1,''
printf,1,'# Objtype:      Fit?         Parameters '
printf,1,''
IF (fitname EQ 'double') THEN BEGIN
    k = 1.9992* fits[4] - 0.3271
    K_const = 2*!pi*(fits[3]*arcsecpix*npix/fov)^2*exp(k)*k^(-2*fits[4])*Gamma(2*fits[4])
    printf,1,'# Object number: 1'
    printf,1,' 0)     sersic         #    Object type'
    printf,1,' 1) ',strtrim(fix(npix/2. + xcen),2),' ',strtrim(fix(npix/2. + ycen),2)$
           ,' 0 0  #   position x, y'
    IF keyword_set(bulgemag) THEN $
      printf,1,' 3) ',strtrim(bulgemag,2) ,'    1       #  total magnitude ' ELSE $     
        printf,1,' 3) ',strtrim(zeropoint - 2.5*alog10(fits[2]*K_const),2) ,'    1       #  total magnitude '
    printf,1,' 4) ',strtrim(npix/fov*fits[3],2)$
           ,'    1       #      R_e in pixels'
    printf,1,' 5) ',strtrim(fits[4],2)$
           ,'    1       #  exponent (de Vaucouleurs = 4) '
    printf,1,' 6) 0.0000     0       #     ----- '
    printf,1,' 7) 0.0000     0       #     ----- '
    printf,1,' 8) 1.0000     0       #     -----  '
    printf,1,' 9) 1.0000     0       #  axis ratio (b/a)  '
    printf,1,' 10) 90.0000    0       #  position angle (PA) '
    printf,1,' Z) 0                  #  Output option (0 = residual, 1 = Dont subtract)' 
    printf,1,''
    printf,1,'# Object number: 2'
    printf,1,' 0)    expdisk         #    Object type'
    printf,1,' 1) ',strtrim(fix(npix/2. + xcen),2),' ',strtrim(fix(npix/2. + ycen),2),' 0 0  #   position x, y'
    IF keyword_set(diskmag) THEN $
      printf,1,' 3) ',strtrim(diskmag,2),'    1       #  total magnitude ' ELSE $
        printf,1,' 3) ',strtrim(zeropoint - 2.5*alog10(fits[0]*2*!pi*(arcsecpix*npix/fov*fits[1])^2),2),'    1       #  total magnitude '
    printf,1,' 4) ',strtrim(npix/fov*fits[1],2),'    1       #      Rs in pixels '
    printf,1,' 5) 0.0000     0       #     ----- '
    printf,1,' 6) 0.0000     0       #     ----- '
    printf,1,' 7) 0.0000     0       #     ----- '
    printf,1,' 8) 0.0000     0       #     ----- '
    printf,1,' 9) 1.0000     0       #  axis ratio (b/a)'  
    printf,1,' 10) 90.0000    0       #  position angle (PA) '
;    printf,1,'10) 0     1            #  diskiness(-)/boxiness(+)'
    printf,1,' Z) 0                  #  Output option (0 = residual, 1 = Dont subtract)'
ENDIF
IF (fitname EQ 'sersic') THEN BEGIN
    k = 1.9992* fits[2] - 0.3271
    K_const = 2*!pi*(fits[1]*npix/fov)^2*exp(k)*k^(-2*fits[2])*Gamma(2*fits[2])
    printf,1,'# Object number: 1'
    printf,1,' 0)     sersic         #    Object type'
    printf,1,' 1) ',strtrim(fix(npix/2. + xcen),2),' ',strtrim(fix(npix/2. + ycen),2),' 1 1  #   position x, y'
    printf,1,' 3) ',strtrim(zeropoint - 2.5*alog10(fits[0]*K_const),2),'    1       #  total magnitude '
    printf,1,' 4) ',strtrim(npix/fov*fits[1],2),'     1       #      R_e in pixels'
    printf,1,' 5) ',strtrim(fits[2],2),'    1       #  exponent (de Vaucouleurs = 4) '
    printf,1,' 6) 0.0000     0       #     ----- '
    printf,1,' 7) 0.0000     0       #     ----- '
    printf,1,' 8) 1.0000     0       #     -----  '
    printf,1,' 9) 1.0000     0       #  axis ratio (b/a)  '
    printf,1,'10) 90.0000    0       #  position angle (PA) '
    printf,1,' Z) 0                  #  Output option (0 = residual, 1 = Dont subtract)' 
    printf,1,''
ENDIF
IF   (fitname EQ 'exp') THEN BEGIN
    printf,1,''
    printf,1,'# Object number: 1'
    printf,1,' 0)    expdisk         #    Object type'
    printf,1,' 1) ',strtrim(fix(npix/2. + xcen),2),' ',strtrim(fix(npix/2. + ycen),2),' 0 0  #   position x, y'
    printf,1,' 3) ',strtrim(zeropoint - 2.5*alog10(fits[0]*2*!pi*(npix/fov*fits[1])^2),2),'    1       #  total magnitude '
    printf,1,' 4) ',strtrim(npix/fov*fits[1],2),'    1       #      Rs in pixels'
    printf,1,' 5) 0.0000     0       #     ----- '
    printf,1,' 6) 0.0000     0       #     ----- '
    printf,1,' 7) 0.0000     0       #     ----- '
    printf,1,' 8) 1.0000     0       #     ----- '
    printf,1,' 9) 1.0000     0       #  axis ratio (b/a)'  
    printf,1,' 10) 90.0000    0       #  position angle (PA) '
    printf,1,' Z) 0                  #   Skip this model in output image?  (yes=1, no=0)'
ENDIF
close,1

;print,'fits[0]: ',fits[0],', Area: ',2*!pi*(npix/fov*fits[1])^2,', Zeropoint: ',zeropoint,' sbfactor: ',sbfactor,', ftot: ',zeropoint - 2.5*alog10(fits[0]*2*!pi*(npix/fov*fits[1])^2)
;print,zeropoint - 2.5*alog10(fits[2]*K_const)

END

