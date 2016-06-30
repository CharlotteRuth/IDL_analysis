pro initial_funtimes, fitsfile, rad=rad

if (n_params() eq 0) then begin
 print, 'Usage: checkfit, fitsfile without extension, /rad'
 stop
endif

;Lunit = filters[band].ewidth_Lambda/filters[band].ewidth_nu; converts internal units to W/m^2
;Fo=3.631e-23  ;zero point of AB magnitude system, in W/m^2
;units = 4.25e10 ; sr>>> arcsec^2  
;sbfactor = Lunit/units/Fo
;ABzero = -2.5*alog10(sbfactor)
    
;; First, clip those annoying edges of fits file
gal = readfits(strcompress(fitsfile+'.fits'),hdr)
;gal =image[100:900,100:900]

;; try and maybe get rid of annoying sky background?
wg = where (gal ne -7777.0)
box = gal[wg]
sky_bkg = median(box)
print, 'sky background is  ',sky_bkg

n_ele=n_elements(gal)

;for i=0,n_ele-1 do begin
;    if (gal[i] ne -7777) then begin 
;        value = gal[i] - sky_bkg
;        if (value lt 0) then (value = 0)
; ;       print, 'value is',value
;        gal[i] = value
;    endif
;endfor


;writefits,strcompress(fitsfile+'_gal.fits'),gal,hdr
;SB profile of the image, using simple radial distances

;; Before anything else, we need to identify the center of the image

;;; Find brightest pixels, which in i-band should correspond to the
;;; center of the bulge, we hope.
;token = 1
;while (token = 1) do begin

    brightest= max(gal,center)

;;; identify which pixels this corresponds to...
;;; and we have our galaxy center!
    width=n_elements(gal(*,0))
    height=n_elements(gal(0,*))
    x_center = (center mod width) 
    y_center = (((center - (center mod width)) / width))


;; To remove CRs: take brightest point, and average all pixels in a 5
;; pixel radius (that may need adjusting, btw).  If it's more than 3-6 
;; sigma (see what works...) away from mean, reject it.
    sweep=2L
    index=0
    centroid=findgen((2*sweep)^2)*0.0
    for n=(x_center - sweep),(x_center+sweep),1 do begin
        for m=(y_center - sweep),(y_center+sweep),1 do begin
;; excluding center so it won't skew std_dev
            if((n ne x_center) and (m ne y_center)) then begin 
                centroid[index]=gal[n,m]
                index++
            endif
        endfor
    endfor
    ctr_mean=median(centroid)
    ctr_stddev=stddev(centroid)
if (center ge ((6.0*ctr_stddev)+ctr_mean)) then print,'You have a possible CR at centroid'

;endwhile

;; now, see whether brightest point fits critera...

;    if(center le ((6*ctr_stddev)+ctr_mean)) then begin
 ;       print, 'centroid located!'
  ;      print, 'x = ',x_center
   ;     print, 'y = ',y_center
;        token = 0
   ; endif else begin
   ;     print, 'estimated centroid was CR' 
   ;     print, 'flagging pixel'
   ;     gal[x_center,y_center] = -7777.0
   ;     print,'xcenter of',x_center,' ycenter of',y_center,' now set to',gal[x_center,y_center]
   ;     print, 'starting over'
   ; endelse


;;; Now create a pixel mask for the benefit of GALFIT
wb = where(gal eq -7777)
n_ele = n_elements(wb)
maskarray = indgen(n_ele, n_ele)
openw,1,"mask.asc"
for i=0,n_ele-1 do begin
    mask_x = (wb[i] mod width) 
    mask_y = (((wb[i] - (wb[i] mod width)) / width))
    printf,1,mask_x,',',mask_y
    mask_x=0
    mask_y=0
endfor

close,1

print,'X center is probably',x_center
print,'Y center is probably',y_center

;STOP

RETURN
END
