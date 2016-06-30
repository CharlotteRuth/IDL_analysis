;plot the profile of the galaxy, with the galfit fits overplotted
;input the galfit fits 'galfit.0x' 
;and fov,npix of your sunrise image, as setin mcrx.stub
pro checkfit_JEM,broadbandfile,galfit,fov,filter,rad=rad

if (n_params() eq 0) then begin
 print, 'Usage: checkfit, broadbandfile, galfit, fov, filter, /rad'
 stop
endif

if filter eq 'B' then band=21
if filter eq 'i' then band=5 
;Get zero point for AB magnitude system
filters= mrdfits(broadbandfile,13) ; the 13 needs adjusting depending on no.cameras
Lunit = filters[band].ewidth_Lambda/filters[band].ewidth_nu; converts internal units to W/m^2
Fo=3.631e-23  ;zero point of AB magnitude system, in W/m^2
units = 4.25e10 ; sr>>> arcsec^2  
sbfactor = Lunit/units/Fo
ABzero = -2.5*alog10(sbfactor)
    
;deal with the Sunrise'd image first
image = mrdfits(broadbandfile,14)  ;14 face on with dust, 19 without dust
;SB profile of the image, using circular annuli
gal =image(*,*,band)
makex, gal,x,y  ;set up coordinate grid
npix = 500.
pixKpc = 1.0*fov/(1.0*npix)
radius = sqrt((x)^2+(y)^2) ; total face on

binsize=1.0
totflux = fltarr(ceil(max(radius)))
avgflux = fltarr(ceil(max(radius)))
rad_kpc = indgen(n_elements(totflux))*binsize
nbins = n_elements(totflux)
for i=0,nbins-1 do begin
   rmin = rad_kpc[i]-(binsize/2.)
   rmax = rad_kpc[i]+(binsize/2.)
   ind = where(radius ge rmin and radius le rmax, nind)
 ;  print, "Number of elements at radius", rad_kpc[i],n_elements(ind)
 ;  totflux[i] = total(gal[ind])
 ;  avgflux[i] = total(gal[ind])/(!pi*(rmax^2.-rmin^2.))
endfor
;profile = -2.5*alog10(sbfactor*avgflux)  

zero=where(gal eq 0.)
;The zeros will mess up your plot range, so set them to a minimum value
gal[zero]=replicate(0.001*min(gal[where(gal gt 0)]),n_elements(zero))

hist = hist1d(radius,gal,Obin=xvals,binsize=1.0) 
numh = hist1d(radius,binsize=1.0)
SBprofile = hist/numh
rad_kpc = xvals*pixKpc
profile = -2.5*alog10(sbfactor*SBprofile)

;Now the galfit results
readcol,galfit,dummy,fits,/silent,format='A,F'
readcol,galfit,dummy,com1,com2,/silent,format='A,F,F'
;zeropoint=fits[3]
;The following are not used below.
;cx=(com1[3]+com1[12]-npix)/2  ;x positions of bulge and disk
;cy=(com2[3]+com2[12]-npix)/2  ;y positions of bulge and disk

;****************************************************
; galfit fits 

pi=3.14159265
ftot=10^((-fits[7]+ABzero)/2.5)  ;total magnitude of Object 1 (bulge)
k = 1.9992*fits[9] - 0.3271      ;Sersic index n
;Rb=pi*(fits[15]+2)/(4*BETA(1/(fits[15]+2),1+1/(fits[15]+2))) ;This param no longer exists in v3
Rb = 1.0 ;disky/boxy fit, no longer exist in GALFIT v3
ue=ftot*Rb/(2*pi*fits[8]^2*fits[9]*exp(k)*k^(-2*fits[9])*GAMMA(2*fits[9])*fits[13])
ue=ue*sbfactor 
fits[8]=pixKpc*fits[8]
bulgefit = abs(ue)*exp(-k*((rad_kpc/fits[8])^(1/fits[9])-1))
bulgeplot = -2.5*alog10(bulgefit)

ftot=10^((-fits[17]+ABzero)/2.5)
;Rb=pi*(fits[25]+2)/(4*BETA(1/(fits[25]+2),1+1/(fits[25]+2)))  ;obsolete
uo=ftot*Rb/(2*pi*fits[18]^2*fits[23])
uo=uo*sbfactor 
fits[18]=pixKpc*fits[18]
diskfit = abs(uo)*exp(-(rad_kpc/fits[18]))
diskplot = -2.5*alog10(diskfit)

fit=bulgefit+diskfit
fitplot = -2.5*alog10(fit)
plot,rad_kpc,(fitplot-profile)
;************************************************
;plotting

;set_plot,'ps'
;device,portrait=1,filename='galfit_i.eps',/color,/encapsulate

xmin=0.
if not keyword_set(rad) then xmax=15. else xmax = rad
region=where(rad_kpc gt xmin and rad_kpc lt xmax)
ymin=max(profile[region])+0.8
;ymin=28.
ymax=min(profile[region])-0.8
pixmin=xmin*(npix/fov)
pixmax=xmax*(npix/fov)


loadct,39
plot,rad_kpc,profile,xrange=[xmin,xmax],/xstyle,thick=2,charsize=1.3,xtit='Radius (kpc)',yrange=[ymin,ymax],/ystyle,psym=0, ticklen=0, pos=[.2,.2,.9,.85]
AXIS, XAXIS=1, XRANGE=[pixmin,pixmax], xtitle='Pixels', xticks=1, xtickinterval=50, ticklen=.02, /xstyle,charsize=1.3
AXIS, XAXIS=0,xrange=[xmin,xmax],/xstyle,charsize=1.3
AXIS, YAXIS=0,yrange=[ymin,ymax],/ystyle,charsize=1.3,ytitle='Magnitude/sq arcsec'
;yrange=[22.49,16.01]

oplot,rad_kpc,diskplot,linestyle=1,color=150, thick=5
oplot,rad_kpc,bulgeplot,linestyle=2,color=250, thick=5
oplot,rad_kpc,fitplot,linestyle=3,color=100, thick=5

legend, ['Sersic', 'Expdisk', 'Sersic+Expdisk'], linestyle=[2,1,3],color=[250,150,100], charsize=1.5, /top, /right
;color=[250,150,100]



print,'scalelength =',fits[18],' kpc'
print,'B/D =',10^(-fits[7]/2.5)/10^(-fits[17]/2.5)
print, 'disk central SB = ',-2.5*alog10(uo)
;device,/close
;set_plot,'x'

stop
;endplot
end

