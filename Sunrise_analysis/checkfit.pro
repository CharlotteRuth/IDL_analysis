;edited for the galfit bootcamp 7/9/10

;plot the profile of the galaxy, with the galfit fits overplotted
;input the galfit fits 'galfit.0x' 
;and fov,npix of your sunrise image, as setin mcrx.stub
;pro checkfit,image,galfit,fov,filter,rad=rad

;camera = 14
;band = 5
;filter = 13
;galfit = 'galfit.04'
;broadbandfile='broadband.fits'

;'/astro/store/nbody3/christensen/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.1'
;camera = 14
;band = 7
;filter = 13
;galfit = 'galfit.05'
;broadbandfile='broadband.fits'

;dir = '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/h986.cosmo50cmb.3072gs1MbwK.00512.dir/h986.cosmo50cmb.3072gs1MbwK.00512.1'
;broadbandfile = 'broadband.fits'
;camera = 14 ;FO with dust
;band = 7; 2MASS H
;filter = 13
;galfit = 'galfit.01'

;checkfit,file,camera,filter,band,'galfit.01',fitname ='exp',yrange = yrange
pro checkfit,broadbandfile,camera,filter,band,galfit,rad = rad,save=save,fitname = fitname,yrange = yrange,recenter = recenter,parameters = parameters

;broadbandfile='broadband.fits'
;camera=14
;filter=13
;galfit='galfit.03'
;fov=50
;npix=500
;rad=14

if (n_params() eq 0) then begin
 print, 'Usage: checkfit, broadband file, index of image, index of filter table, band number, galfit file, field of view(kpc), image size in pixels, radius of fit'
 stop
endif

IF NOT keyword_set(fitname) THEN fitname = 'double'

temp = mrdfits(broadbandfile,2,header)
fovstr = strsplit(header[13],' ',/extract)
fov = float(fovstr[3])
diststr = strsplit(header[12],' ',/extract)
dist = float(diststr[3]) ;Camera distance in kpc

;------------------------------------------------------------------------------
;Get zero point for AB magnitude system
filters= mrdfits(broadbandfile,filter) 
Lunit = filters[band].ewidth_Lambda/filters[band].ewidth_nu; converts internal units to W/m^2
Fo=3.631e-23  ;zero point of AB magnitude system, in W/m^2
angleunits = 4.25e10 ; sr>>> arcsec^2  
distconv = (1000.*dist/10)^2 
sbfactor = Lunit/Fo/angleunits
ABzero = -2.5*alog10(sbfactor*distconv)

;deal with the Sunrise'd image first
image = mrdfits(broadbandfile,camera)  ;14 face on with dust, 19 without dust
;-----------------------------------------------------------------
;SB profile of the image, using circular annuli
gal =image(*,*,band);*distconv
makex, gal,x,y  ;set up coordinate grid
npix = n_elements(gal[0,*])
kpcpix = 1.0*fov/(1.0*npix);kpc/pix
arcsecpix = atan(fov/dist)*360/(2*!PI)*3600/npix ;arcsec/pix

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

radius = sqrt((x)^2+(y)^2) ; total face on
binsize= 2.0 ;two pixels
medflux = fltarr(ceil(max(radius))/binsize)
meanflux = fltarr(ceil(max(radius))/binsize)
rad_kpc = indgen(n_elements(medflux))*binsize + binsize/2.0
nbins = n_elements(medflux)
for i=0L,nbins-2 do begin
   rmin = rad_kpc[i]-(binsize/2.)
   rmax = rad_kpc[i]+(binsize/2.)
   ind = where(radius gt rmin and radius le rmax, nind)
;  print, "Number of elements at radius", rad_kpc[i],n_elements(ind)
;  print, "Median of elements at radius",rad_kpc[i],median(gal[ind])
   medflux[i] = median(gal[ind])
   meanflux[i] = MEAN(gal[ind])
;  print,"rad_kpc value is", rad_kpc[i]
endfor
rad_kpc = rad_kpc * kpcpix
profile = -2.5*alog10(sbfactor*medflux)
profile_mean = -2.5*alog10(sbfactor*meanflux)
;--------------------------------------------------------------------------
zero=where(gal eq 0.)
;The zeros will mess up your plot range, so set them to a minimum value
gal[zero]=replicate(0.001*min(gal[where(gal gt 0)]),n_elements(zero))

;The following is the old checkfit
;It uses the mean flux, but galfit seems to use the median
; hence, the above changes
;hist = hist1d(radius,gal,obin=xvals,binsize=1.0) 
;numh = hist1d(radius,binsize=1.0)
;SBprofile = hist/numh
;rad_kpc = xvals*kpcpix
;profile = -2.5*alog10(sbfactor*SBprofile)

;Now the galfit results
readcol,galfit,dummy,fits,/silent,format='A,F'
readcol,galfit,dummy,com1,com2,/silent,format='A,F,F'
;cy=(com2[3]+com2[12]-npix)/2  ;y positions of bulge and disk

;****************************************************
; galfit fits 

IF fitname eq 'double' then begin
    ftot=10^((-fits[7]+ABzero)/2.5) ;total magnitude of Object 1 (bulge)
    k = 1.9992*fits[9] - 0.3271 ;Sersic index n
;Rb=!PI*(fits[15]+2)/(4*BETA(1/(fits[15]+2),1+1/(fits[15]+2))) ;This param no longer exists in v3
    Rb = 1.0             ;disky/boxy fit, no longer exist in GALFIT v3
    ue=ftot*Rb/(2*!PI*fits[8]^2*fits[9]*exp(k)*k^(-2*fits[9])*GAMMA(2*fits[9])*fits[13])/arcsecpix^2
    ue=ue*sbfactor 
    fits[8]=kpcpix*fits[8];bulge scalelength [kpc]
    bulgefit = abs(ue)*exp(-k*((rad_kpc/fits[8])^(1/fits[9]) - 1))
    bulgeplot = -2.5*alog10(bulgefit)

    ftot=10^((-fits[17]+ABzero)/2.5);total disk flux
;Rb=!PI*(fits[25]+2)/(4*BETA(1/(fits[25]+2),1+1/(fits[25]+2)))  ;obsolete
    uo=ftot*Rb/(2*!PI*fits[18]^2*fits[23])/arcsecpix^2;fits[23] = axis ratio
    uo=uo*sbfactor 
    fits[18]=kpcpix*fits[18];disk scalelength [kpc]
    diskfit = abs(uo)*exp(-(rad_kpc/fits[18]))
    diskplot = -2.5*alog10(diskfit)

    fit=bulgefit+diskfit
    fitplot = -2.5*alog10(fit)
ENDIF 
IF fitname eq 'exp' then begin
    ftot=10^((-fits[7]+ABzero)/2.5)
;Rb=!PI*(fits[25]+2)/(4*BETA(1/(fits[25]+2),1+1/(fits[25]+2)))
;;obsolete
    Rb = 1.0             ;disky/boxy fit, no longer exist in GALFIT v3
    uo=ftot*Rb/(2*!PI*fits[8]^2*fits[13])
    uo=uo*sbfactor 
    fits[8]=kpcpix*fits[8]

    diskfit = abs(uo)*exp(-(rad_kpc/fits[8]))
    diskplot = -2.5*alog10(diskfit)
    fit=diskfit
    fitplot = -2.5*alog10(fit)
ENDIF

;Mean Surface Brightness
n = fits[9]
r_0 = fits[8];(2.17*n - 0.355)^n*fits[8]
bulge_ind = where(radius*kpcpix LE r_0)
meanSB = -2.5*alog10(mean(sbfactor*gal[bulge_ind]))

;************************************************
;plotting

IF keyword_set(save) THEN BEGIN
    formatplot,/outplot
    set_plot,'ps'
    device,filename='./checkfit_'+STRTRIM(filters[band].filter,2)+'.ps', /color   
ENDIF ELSE BEGIN
    formatplot
    window,1
ENDELSE
xmin=0.
if not keyword_set(rad) then xmax=8. else xmax = rad
region=where(rad_kpc gt xmin and rad_kpc lt xmax)
IF NOT KEYWORD_SET(yrange) THEN BEGIN
    temp =  profile[region]
    ymin = max(temp[where(FINITE(temp))])+0.8
    ymax = min(temp[where(FINITE(temp))])-0.8
ENDIF ELSE BEGIN
    ymin = yrange[0]
    ymax = yrange[1]
ENDELSE

pixmin=xmin*(npix/fov)
pixmax=xmax*(npix/fov)

loadct,39
plot,rad_kpc,profile,xrange=[xmin,xmax],/xstyle,thick=4,xtit='Radius (kpc)',yrange=[ymin,ymax],/ystyle,psym=0,ytit=textoidl('M/arcsec^2'),/nodata
oplot,rad_kpc,profile,color = 100,thick=4
;AXIS, XAXIS=1, XRANGE=[pixmin,pixmax], xtitle='Pixels', xticks=1, xtickinterval=50, ticklen=.02, /xstyle,charsize=1.3
;AXIS, XAXIS=0,xrange=[xmin,xmax],/xstyle,charsize=1.3
;AXIS, YAXIS=0,yrange=[ymin,ymax],/ystyle,charsize=1.3,ytitle='Magnitude/sq arcsec'; pos=[.2,.2,.9,.85]
;yrange=[22.49,16.01]
;oplot,rad_kpc,profile_mean,color = 60
IF (fitname EQ 'double' OR fitname EQ 'sersic') THEN oplot,rad_kpc,bulgeplot,linestyle=1,thick = 4,color = 60;color=250
IF (fitname EQ 'double' OR fitname EQ 'exp') THEN oplot,rad_kpc,diskplot,linestyle=2,color = 60,thick = 4;,color=150
oplot,rad_kpc,fitplot,linestyle=3,thick = 4,color = 60;,color=50
oplot,[fits[8],fits[8]],[25,14],linestyle = 1

IF (fitname EQ 'double') THEN BEGIN
;    legend, ['Sersic', 'Expdisk', 'Sersic+Expdisk'], linestyle=[1,2,3],color=[250,150,50], /top, /right
;    legend,[textoidl('R_S') + ':           '+STRTRIM(fits[18],2),'B/D:         '+STRTRIM(10^(-fits[7]/2.5)/10^(-fits[17]/2.5),2),'Central SB: '+STRTRIM(-2.5*alog10(uo),2)],/bottom,/left
ENDIF

IF (fitname eq 'exp') THEN legend,[textoidl('R_S') + ':           '+STRTRIM(fits[8],2)],linestyle=[2],color=[150]
IF (fitname eq 'exp') THEN print,'Disk scalelength =',fits[8],' kpc'

IF (fitname eq 'double') THEN BEGIN
   print,'Total Mag.:    ',strtrim(-2.5*alog10(10^(-0.4*fits[7]) + 10^(-0.4*fits[17])),2),', ',strtrim(-2.5*alog10(10^(-0.4*fits[7]) + 10^(-0.4*fits[17])) - 1.39,2)
   print,'Bulge Mag.:    ',strtrim(fits[7],2),', ',strtrim(fits[7] - 1.39,2)
   print,'Mean SB:       ',strtrim(meanSB,2),', ',strtrim(meanSB-1.39,2)
   print,'log(R_e) [pc]: ',strtrim(alog10(fits[8]*1000),2)
   print,'Sersic Index:  ',strtrim(fits[9],2)
   print,'Disk Mag:      ',strtrim(fits[17],2),', ',strtrim(fits[17] - 1.39,2) 
   print,'log(R_h) [pc]: ',strtrim(alog10(fits[18]*1000),2)
   print,'B/D:           ',strtrim(10^(-fits[7]/2.5)/10^(-fits[17]/2.5),2)
   print,'B/T:           ',strtrim(10^(-fits[7]/2.5)/(10^(-fits[17]/2.5) + 10^(-fits[7]/2.5)),2) 
ENDIF

IF keyword_set(save) THEN BEGIN
    device,/close 
    set_plot,'x' 
ENDIF

parameters = {M:-2.5*alog10(10^(-0.4*fits[7]) + 10^(-0.4*fits[17])),$
              M_B:fits[7],$
              r_e:alog10(fits[8]*1000),$
              r_0:alog10(r_0*1000),$
              mu:float(meanSB),$
              sersic:fits[9],$
              M_D:fits[17],$
              r_h:alog10(fits[18]*1000),$
              B_D:10^(-fits[7]/2.5)/10^(-fits[17]/2.5) $
             }

end
