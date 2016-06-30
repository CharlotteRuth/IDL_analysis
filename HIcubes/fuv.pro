;broadband = '/astro/net/nbody1/abrooks/h799.cosmo25cmb.3072g3bwdK/sunrise/broadband.h799.3072.fits'
pro fuv, broadband

image=mrdfits(broadband, 15)
gal=image(*,*,0) ;picks out FUV image

tvim, gal <0.1  ;displays the image
cam=15 ; face on without dust needs adjusting 
band = 0 ; chose your band (look in the  filters file)
fov=10.  ; adjust this t omatch mcrx input
npixels=500 ; adjust this 
filt=11
;;;;now do various apertures to creat curve of growth;;;;
flux=fltarr(191)
x=250.
y=250.
pixKpc = fov/npixels  ; fov/npixels needs to be adjusted acordingly
filters= mrdfits(broadband,filt) ; the 12 needs adjusting depending on no.cameras
Lunit = filters[band].ewidth_Lambda/filters[band].ewidth_nu; internal units to W/m^2
Lunit=Lunit*1d7 ;convert to ergs/s/m^2
;stop
for i=10,199 do begin

   DIST_ELLIPSE, mask, 500, x, y, 1, 0 ;Create a elliptical image mask
   good = where( mask lt i )   ;Within aperture i
   flux[i-10]=total( gal[good] );/N_ELEMENTS(gal[good])  ;Total pixel values within aperture i
;radii=findgen(200)
;aper, gal, x, y,flux, eflux, 1, [10,20,30,40,50,60,70,80,100],-1,[-32767,80000],/exact,/flux, setsky=0

;   stop
endfor
;steps = findgen(199 - 10) + 10
;contour,gal,levels = steps

;plot, flux, psym=4, xtit='aperture in pixels', ytit='flux'w
window,1
flux=abs(flux*Lunit)
plot, flux, psym=4, xtit='aperture in pixels', ytit='flux'
Fo=3.631e-23
flux=flux/Fo
;mag=(alog10(flux)/(-0.4))+48.6
;plot, mag, psym=4, xtit='aperture in pixels', ytit='mag', yrange=[22,25]
;stop
;;;now compute derivatives;;;;
dr=10.
;dmdr=fltarr(n_elements(flux)-1)
;for i=2, n_elements(flux)+1 do begin
    ;dm=flux[i-1]-flux[i-2]
    ;if dm gt 0 then dmdr[i]=dm/dr
    ;dmdr[i]=dm/dr
;endfor
stop
dm=fltarr(13)
dmdr=fltarr(13)
flux1=fltarr(13)
for i=1, 13 do begin

    dm[i-1]=flux[i*10]-flux[i*10-10]
    flux1[i-1]=flux[i*10]

    dmdr[i-1]=dm[i-1]/dr

endfor
window, 2
plot, dmdr,flux, psym=4


a=linfit(dmdr, flux1)
print, a[0]

print,(alog10(a[0])/(-0.4))+48.6

stop

;;fit by an exponential?
window, 3
r=findgen(200)+10.
light=exp(-(r-10.)/37.5622)


plot, r, exp(-(r-10.)/37.5622)+2.58E11, yrange=[2.5799e11, 2.5801e11]

stop
end
