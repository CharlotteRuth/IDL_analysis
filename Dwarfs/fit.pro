1;example of fitting an edge-on galaxy using the markwardt package and
;        the code in this directory
pro fit_it,file,band
;bands: $
;0 = GALEX FUV,    1 = GALEX NUV,    2 = SDSS u band,   3 = SDSS g band$ 
;4 = SDSS r band,  5 = SDSS i band,  6 = SDSS z band,   7 = HST ACS F435 $
;8 = HST ACS F606, 9 = HST ACS F775, 10 = HST ACS F850, 11 = IRAC1 SIRTF $  
;12 = IRAC2 SIRTF, 13 = MIPS24 SIRTF,14 = MIPS70 SIRTF  15 = MIPS160 SIRTF

image =mrdfits(file,13); face on wo dust. 13 is with dust
filter =mrdfits(file,12)
Lunit =filter[band].L_lambda_to_L_n ; internal units to W/m^2
gal = image(*,*, band)
;gal=gal+1.01*abs(min(gal)) ; use if there are negative gal values
makex, gal,x,y  ;set up coordinate grid
units = 4.35e10  ;sr>>>arcsec^2
Fo = 3.631e-23; AB zero point in W/m^2/Hz
gal=gal 
sb_norm = -2.5*alog10(Lunit)+2.5*alog10(Fo)+2.5*alog10(units)    
stop

radius = sqrt(x^2+y^2) ; total face on
hist = hist1d(radius,gal,obin=xvals,binsize=5) 
numh = hist1d(radius,binsize=5)

SBprofile = hist/numh
rad_kpc = pixKpc*xvals
err = SBprofile/sqrt(numh)
rmin = 4.
rmax = 7.9

loadct,39
region = where(rad_kpc gt rmin and rad_kpc lt rmax)

; plot the profile 
plot,rad_kpc,-2.5*alog10(SBprofile)+sb_norm
; now make the fit:
; make initial guesses
init = [.2,6.,1.5]              ;for sersic profile 


; make the fit 

fit = mpfitfun('sersic', rad_kpc[region], SBprofile[region], err[region], init)

loadct,39

; plot the fit
fit = init
k = 1.9992* fit[2] - 0.3271


bulgefit = abs(fit[0])*exp(-k*((rad_kpc/fit[1])^(1/fit[2])-1))
fitplot = bulgefit   
oplot,rad_kpc,alog10(fitplot),color = 58 ; plot the total fit

print,"bulge csb, Re, n=" ; write out the final fit values to screen
print, fit

stop
end

