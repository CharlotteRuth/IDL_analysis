;+
; NAME:
;   make_jpeg
;- *** Note that the camera angle needs to be set!!!!!! ***** 
pro make_jpeg,bands=bands,scales=scales,nonlinearity=nonlinearity,outfile=outfile,cam=cam
; default scales = [4.9,5.7,7.8] as set in nw_scale_rgb.pro
; default nonlinearity = 3 as set in nw_arcsinh_fit.pro
if not keyword_set(bands) then bands = [5,4,3] 
if not keyword_set(cam) then cam = 14
;default bands SDSS irg, ie rgb=irg, 5,4,3 in our standard filter list

resizefactor= 1.
;cam=14; face on 
; read data
im= mrdfits('broadband.fits',cam)
tmp= size(im,/dimensions)
nx= tmp[0]
ny= tmp[1]
RGBim= fltarr(nx,ny,3)
RGBim[*,*,0]= im(*,*,bands[0])
RGBim[*,*,1]= im(*,*,bands[1])
RGBim[*,*,2]= im(*,*,bands[2])

; rebin
RGBim= rebin(RGBim,floor(nx*resizefactor),(ny*resizefactor),3)

; scale and set colors
RGBim = nw_scale_rgb(RGBim,scales=scales)
RGBim = nw_arcsinh_fit(RGBim,nonlinearity=nonlinearity)
RGBim = nw_fit_to_box(RGBim,origin=origin)
RGBim = nw_float_to_byte(RGBim)

; write
WRITE_JPEG,outfile,RGBim,TRUE=3,QUALITY=100
return
end
