;plot the profile of the galaxy, with the galfit fits overplotted
;input the galfit fits 'galfit.0x' 
;and fov,npix of your sunrise image, as setin mcrx.stub
pro slices,broadbandfile,fov,filter,rad=rad

if (n_params() eq 0) then begin
 print, 'Usage: slices, broadbandfile, fov, filter, /rad'
 stop
endif

if filter eq 'FUV' then band=0
if filter eq 'NUV' then band=1
if filter eq 'u' then band=2 
if filter eq 'g' then band=3 
if filter eq 'r' then band=4 
if filter eq 'i' then band=5 
if filter eq 'z' then band=6 
if filter eq 'H_2MASS' then band=7
if filter eq 'J_2MASS' then band=8
if filter eq 'K_2MASS' then band=9
if filter eq 'IRAC1' then band=10
if filter eq 'IRAC2' then band=11
if filter eq 'IRAC3' then band=12
if filter eq 'IRAC4' then band=13
if filter eq 'MIPS160' then band=14
if filter eq 'MIPS24' then band=15
if filter eq 'MIPS70' then band=16
if filter eq 'F435' then band=17
if filter eq 'F606' then band=18
if filter eq 'F775' then band=19
if filter eq 'F850' then band=20
if filter eq 'B' then band=21
if filter eq 'V' then band=22
if filter eq 'I' then band=23
if filter eq 'H' then band=24
if filter eq 'J' then band=25
if filter eq 'K' then band=26

;Get zero point for AB magnitude system
filters= mrdfits(broadbandfile,13) ; the 13 needs adjusting depending on no.cameras
Lunit = filters[band].ewidth_Lambda/filters[band].ewidth_nu; converts internal units to W/m^2
Fo=3.631e-23  ;zero point of AB magnitude system, in W/m^2
units = 4.25e10 ; sr>>> arcsec^2  
sbfactor = Lunit/units/Fo
ABzero = -2.5*alog10(sbfactor)
    
;deal with the Sunrise'd image first
image = mrdfits(broadbandfile,14)  ;14 face on with dust, 19 without dust
gal =image(*,*,band)

;SB profile of the image, using simple radial distances

;; Before anything else, we need to identify the center of the image
;; This is the roughest of rough kludges, need to refine to exclude
;; CRs

;;; start by reading in image in i-band
i_image = image(*,*,5)

;;; Find brightest pixels, which in i-band should correspond to the
;;; center of the bulge, we hope.

center= max(i_image,brightest)

;;; identify which pixels this corresponds to...
;;; and we have our galaxy center!
width=n_elements(i_image(*,0))
height=n_elements(i_image(0,*))
x_center = (brightest mod 500) 
y_center = (((brightest - (brightest mod 500)) / 500))

;; make an array of size maximum radius from galaxy center in which
;; to store unnormalized flux values, and another to store # of
;; pixels in each of those bins for normalization

;;; First need to find what maximum radius from galaxy center is
if ((x_center-width)^2 gt (x_center)^2) then (maxx = (x_center-width)) else (maxx = x_center)
if ((y_center-height)^2 gt (y_center)^2) then (maxy = (y_center-height)) else (maxy = y_center) 
max_rad = floor(sqrt((maxx)^2 + (maxy)^2))
;;; Now make arrays of size max radius from galaxy center
sum_flux = findgen(max_rad)*0.0
num_pix = findgen(max_rad)*0.0
norm_flux = findgen(max_rad)*0.0
pixKpc = 1.0*fov/(1.0*width)
rad_kpc = indgen(max_rad)* pixKpc ;this holds radius in physical units


;The zeros will mess up your plot range, so set them to a minimum value
zero=where(gal eq 0.)
gal[zero]=replicate(0.001*min(gal[where(gal gt 0)]),n_elements(zero))

;; Start making SB profile
for i=0,width-1 do begin
    for j=0,height-1 do begin

        radius = floor(sqrt((x_center-i)^2 +(y_center-j)^2 )) ;radius from galaxy center
        sum_flux[radius] += gal[i,j]
        num_pix[radius] += 1
 
    endfor
endfor

;; Normalize to # pixels in bin

for i=0,max_rad-1 do begin
    norm_flux[i]=(sum_flux[i]/num_pix[i])
endfor
    
profile = -2.5*alog10(sbfactor*norm_flux)  

;
;;
;;; Time for slice 1: from galcenter increment in y until end
;;
;
slice1rad = height-y_center
slice1flux = findgen(max_rad)*0.0
slice1num_pix = findgen(max_rad)*0.0
slice1norm_flux = findgen(max_rad)*0.0
;slice1rad_kpc = indgen(slice1rad)* pixKpc ;this holds radius in physical units
i=0;
radius=0.;
for i=y_center,height-1 do begin

    radius = floor(sqrt((i-y_center)^2 )) ;radius from galaxy center
    slice1flux[radius] += gal[x_center,i]
    slice1num_pix[radius] += 1
 
endfor

;; Normalize to # pixels in bin

for i=0,slice1rad-1 do begin
    slice1norm_flux[i]=(slice1flux[i]/slice1num_pix[i])
endfor
    
profileslice1 = -2.5*alog10(sbfactor*slice1norm_flux)  
;plot,slice1rad_kpc,slice1norm_flux,xtit='Radius from GalCenter',ytit='Flux'
;PAUSE
;print, 'Moving on'
;
;;
;;; Time for slice 2: from galcenter increment in x until end
;;
;
slice2rad = width-x_center
slice2flux = findgen(max_rad)*0.0
slice2num_pix = findgen(max_rad)*0.0
slice2norm_flux = findgen(max_rad)*0.0
;slice2rad_kpc = indgen(slice2rad)* pixKpc ;this holds radius in physical units
i=0;
radius=0.;
for i=x_center,width-1 do begin

    radius = floor(sqrt((x_center-i)^2)) ;radius from galaxy center
    slice2flux[radius] += gal[i,y_center]
    slice2num_pix[radius] += 1
 
endfor

;; Normalize to # pixels in bin

for i=0,slice2rad-1 do begin
    slice2norm_flux[i]=(slice2flux[i]/slice2num_pix[i])
endfor
    
profileslice2 = -2.5*alog10(sbfactor*slice2norm_flux)  
;plot,slice2rad_kpc,slice2norm_flux,xtit='Radius from GalCenter',ytit='Flux'
;PAUSE
;print, 'Moving on'


 
;
;;
;;; Time for slice 3: from galcenter increment in x and y until end
;;
;
slice3rad = floor(sqrt((width-x_center)^2 + (height-y_center)^2))
slice3flux = findgen(max_rad)*0.0
slice3num_pix = findgen(max_rad)*0.0
slice3norm_flux = findgen(max_rad)*0.0
;slice3rad_kpc = indgen(slice3rad)* pixKpc ;this holds radius in physical units
i=0;
radius=0.;
for i=x_center,width-1 do begin
    for j=y_center, height-1 do begin

    radius = floor(sqrt((x_center-i)^2 +(y_center-j)^2 )) ;radius from galaxy center
    slice3flux[radius] += gal[i,j]
    slice3num_pix[radius] += 1
 
    endfor
endfor

;; Normalize to # pixels in bin

for i=0,slice3rad-1 do begin
    slice3norm_flux[i]=(slice3flux[i]/slice3num_pix[i])
endfor
    
profileslice3 = -2.5*alog10(sbfactor*slice3norm_flux)  
;plot,slice3rad_kpc,slice3norm_flux,xtit='Radius from GalCenter',ytit='Flux'
;PAUSE
;print, 'Moving on'


;************************************************
;plotting

set_plot,'ps'
device,portrait=1,filename='profile_compare.eps',/color,/encapsulate

xmin=0.
xmax=0.
if not keyword_set(rad) then xmax=15. else xmax = rad
region=where(rad_kpc gt xmin and rad_kpc lt xmax)
ymin=max(profile[region])+0.8
;ymin=28.
ymax=min(profile[region])-0.8
pixmin=xmin*(width/fov)
pixmax=xmax*(width/fov)

;window,0,retain=2
loadct,39
plot,rad_kpc,profile,xrange=[xmin,xmax],/xstyle,thick=2,charsize=1.3,xtit='Radius (kpc)',yrange=[ymin,ymax],/ystyle,psym=0, ticklen=0, pos=[.2,.2,.9,.85]
AXIS, XAXIS=1, XRANGE=[pixmin,pixmax], xtitle='Pixels', xticks=1, xtickinterval=50, ticklen=.02, /xstyle,charsize=1.3
AXIS, XAXIS=0,xrange=[xmin,xmax],/xstyle,charsize=1.3
AXIS, YAXIS=0,yrange=[ymin,ymax],/ystyle,charsize=1.3,ytitle='Magnitude/sq arcsec'
;yrange=[22.49,16.01]
oplot,rad_kpc,profileslice1,color=250, linestyle=1,thick=5
oplot,rad_kpc,profileslice2,color=150, linestyle=2,thick=5
oplot,rad_kpc,profileslice3,color=100, linestyle=3,thick=5

;oplot,rad_kpc,diskplot,linestyle=1,color=150, thick=5
;oplot,rad_kpc,bulgeplot,linestyle=2,color=250, thick=5
;oplot,rad_kpc,fitplot,linestyle=3,color=100, thick=5

legend, ['Slice down y', 'Slice across x', 'Slice 45 degrees up'], linestyle=[1,2,3],color=[250,150,100], charsize=1.5, /top, /right
;color=[250,150,100]



;print,'scalelength =',fits[18],' kpc'
;;print,'B/D =',10^(-fits[7]/2.5)/10^(-fits[17]/2.5)
;print, 'disk central SB = ',-2.5*alog10(uo)
device,/close
set_plot,'x'

stop
endplot
return
end
