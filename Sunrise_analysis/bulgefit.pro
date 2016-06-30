;example of fitting an edge-on galaxy using the markwardt package and
;        the code in this directory
pro bulgefit,file,cam,band

;***UPDATE THESE FOR EACH GALAXY***

fov = 90     ; field of view in kpc
npix = 900   ;number of pixels on each side of the image
;cam = 13     ;index number of desired camera angle in fv
;band = 5     ;filter number starting from zero
filt = 12    ;index number of filters in fv 



;Region to fit in (kpc)

rmin = .1  
rmax = 3

;Region to plot (kpc)

xmin=0.
xmax=15.



;Initial guesses for fitting parameters of the form:
;init=[U_e, R_e, Sersic Index]

init = [19.65,3.,2]


;***END OF USER INPUT***

image =mrdfits(file,cam); face on wo dust.
gal = image(*,*, band)
makex, gal,x,y  ;set up coordinate grid
pixKpc = 1.0*fov/(1.0*npix)
filters= mrdfits(file,filt)
;print,filters,/st
;help,filters,/st
;stop
;Lunit = filters[band].ewidth_lambda/filters[band].ewidth_nu ; internal units
;to W/m^2 for Sunrise3
Lunit = filters[band].l_lambda_to_l_nu; internal units to W/m^2
Fo=3.631e-23
units = 4.35e10 ; sr>>> arcsec^2
sbfactor = Lunit/units/Fo

dummy = max(gal,nn)
max=array_indices(gal,nn)-npix/2
x=x-max[0]
y=y-max[1]

;for fitting along a slit
;radius=abs(y)
;slit=where(abs(x) lt 10)
;hist = hist1d(radius[slit],gal[slit],obin=xvals,binsize=10) 
;numh = hist1d(radius[slit],binsize=10)

;for circular annuli
 radius = sqrt((x)^2+(y)^2) ; ciruclar annuli
 hist = hist1d(radius,gal,obin=xvals,binsize=10) 
 numh = hist1d(radius,binsize=10)

SBprofile = hist/numh
rad_kpc = xvals*pixKpc
err = SBprofile/sqrt(numh)

loadct,39
region = where(rad_kpc gt rmin and rad_kpc lt rmax)


; make initial guesses

parinfo=replicate({fixed:0},3) 
parinfo.fixed=[0,0,0]

; make the fit 

fits = mpfitfun('sersic', rad_kpc[region], SBprofile[region], err[region],init,parinfo=parinfo)

loadct,39

;diskfit = abs(fits[0])*exp(-rad_kpc/fits[1])
;diskplot = -2.5*alog10(sbfactor*diskfit)
k = 1.9992* fits[2] - 0.3271
bulgefit = abs(fits[0])*exp(-k*((rad_kpc/fits[1])^(1/fits[2])-1))
bulgeplot = -2.5*alog10(sbfactor*bulgefit)
fit=bulgefit
profile = -2.5*alog10(sbfactor*SBprofile)
fitplot = -2.5*alog10(sbfactor*fit)
;csb = strnsignif(-2.5*alog10(sbfactor*fits[0]),3)
;scalelength = strnsignif(fits[1],3)
u_e= -2.5*alog10(sbfactor*fits[0])
r_e = fits[1]
n= fits[2]


;set_plot,'ps'
;device,/portrait,filename='./profilefit.ps', /color


region=where(rad_kpc gt xmin and rad_kpc lt xmax)
ymin=max(profile[region])+0.8
ymax=min(profile[region])-0.8
loadct,39

plot,rad_kpc,profile,xrange=[xmin,xmax],/xstyle,thick=10,charsize=1.5,xtit='Radius (kpc)',ytit=textoidl('M/arcsec^2'),yrange=[ymin,ymax],/ystyle,title='Sersic Index='+string(fits[2])+'kpc'

oplot,rad_kpc,bulgeplot,linestyle=1,color=250



;device,/close
;set_plot,'x


print, fits

stop
end

