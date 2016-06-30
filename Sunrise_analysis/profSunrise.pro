;example of fitting an edge-on galaxy using the markwardt package and
;        the code in this directory

;file ='/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512/h603.cosmo50cmb.2304g2bwK.00512.1/broadband.fits'
;file = '/astro/net/scratch1/abrooks/FABIO/h516.cosmo25cmb.2304g2bwK/h516.cosmo25cmb.2304g2bwK.00512/h516.cosmo25cmb.2304g2bwK.00512.1/broadband.fits'


pro profSunrise_master,outplot = outplot

!p.multi = 0
loadct,39
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps' 
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    device,filename = outplot,/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 
ENDIF ELSE BEGIN
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    window,0
ENDELSE

cam = 14
band = 21;6 - 1
fov =  24 ;20.
filt = 13
npix = 480 ;800.
dirs = ['/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512.1/', $
        '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.1', $
        '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir/h516.cosmo25cmb.1536g6MbwK.00512.1']
file = 'broadband.fits'
colors = [240,140,80]
SPAWN,'pwd',home

cd,dirs[0]
profSunrise, file, npix = npix, fov = fov, filt = filt, cam = cam, band = band

FOR i = 0, N_ELEMENTS(dirs) - 1 DO BEGIN
    cd,dirs[i]
    profSunrise, file, npix = npix, fov = fov, filt = filt, cam = cam, band = band, color = colors[i],/overplot
ENDFOR

IF KEYWORD_SET(outplot) THEN device,/close
cd,home

end


pro profSunrise,file,fov = fov, filt = filt, npix = npix, cam = cam, band = band, color = color, overplot = overplot

;***UPDATE THESE FOR EACH GALAXY***

IF not keyword_set(fov) THEN fov = 20.     ; field of view in kpc, h603
IF not keyword_set(npix) THEN npix = 800.   ;number of pixels on each side of the image
IF not keyword_set(filt) THEN filt = 12    ;index number of filters in fv 
IF not keyword_set(cam) THEN cam = 16     ;index number of desired camera angle in fv
IF NOT keyword_set(band) THEN band = 5     ;filter number starting from zero ;i_SDSS
IF NOT keyword_set(fitname) THEN fitname = 'double'
IF NOT keyword_set(mask) THEN mask = 0

;Region to fit in (kpc)
rmin = mask
rmax = 8.0

;Region to plot (kpc)
xmin=mask;0.0
xmax=10.0

;------------------------------------------------------------------------------
;Get zero point for AB magnitude system
filters= mrdfits(file,filt)
Lunit = filters[band].ewidth_lambda/filters[band].ewidth_nu ; internal units
;to W/m^2 for Sunrise3
;Lunit = filters[band].l_lambda_to_l_nu; internal units to W/m^2
Fo=3.631e-23 ;zero point of AB magnitude system, in W/m^2
units = 4.35e10 ; sr>>> arcsec^2
sbfactor = Lunit/units/Fo
;ABzero = -2.5*alog10(sbfactor)
;deal with the Sunrise'd image first
image = mrdfits(file,cam)

;-----------------------------------------------------------------
;SB profile of the image, using circular annuli
gal = image(*,*, band)
makex, gal,x,y  ;set up coordinate grid
pixKpc = 1.0*fov/(1.0*npix)

dummy = max(gal,nn)
max=array_indices(gal,nn)-npix/2
x=x-max[0]
y=y-max[1]

;for circular annuli
radius = sqrt((x)^2+(y)^2)      ; ciruclar annuli
hist = hist1d(radius,gal,obin=xvals,binsize=2) 
numh = hist1d(radius,binsize=2)

SBprofile = hist/numh
rad_kpc = xvals*pixKpc
err = SBprofile/sqrt(numh)
profile = -2.5*alog10(sbfactor*SBprofile)

loadct,39
region = where(rad_kpc gt rmin and rad_kpc lt rmax)

region = where(rad_kpc gt xmin and rad_kpc lt xmax)
ymin = max(profile[region])+0.8
ymax = min(profile[region])-0.8


IF KEYWORD_SET(overplot) THEN oplot,rad_kpc,profile,color = color, thick = 3 $ 
ELSE plot,rad_kpc,profile,xrange=[0,xmax],/xstyle,charsize=1.25,xtit='Radius [kpc]',ytit=textoidl('SB [M_B/arcsec^2]'),yrange=[ymin,ymax],/ystyle,color = color,thick = 3;,title='Scale Length='+strtrim(npix/fov*fits[1],2)+' kpc'

end


