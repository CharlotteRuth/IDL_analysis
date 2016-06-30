;+
; NAME:
;   make_jpeg
;- *** Note that the camera angle needs to be set!!!!!! ***** 
;13,14,15,16
;bands = 6,5,4

;cd,'/astro/net/scratch2/christensen/Sunrise/MW1lr_v2'
;bands = [5,4,3]
;cam = 14
;outfile = 'edgeon_MW1lr_v2.jpeg'

;cam = 13
;outfile = 'faceon_MW1lr_v2.jpeg'
;make_jpeg,bands = bands, cam=cam, outfile = outfile

;cd,'/astro/net/scratch2/cbrook/analysis/sunrise/h258/h258.cosmo50cmb.1536g2bwK.00264.1';pixels=800
;cd,'/astro/net/scratch1/fabio/h258/h258.cosmo50cmb.1536g2bwK.00264.1.v3'
;bands = [5,4,3]
;bands = [4,3,2]
;resizefactor = .25
;cam = 11
;outfile = 'faceon_h258_v3.jpeg'
;make_jpeg,bands = [5,4,3],resizefactor = 0.25,cam = 11,outfile = 'faceon_h258_v3.jpeg'


;cd,'/astro/net/scratch1/fabio/h258/h258.cosmo50cmb.1536g2bwK.00264.1
;make_jpeg,bands = [4,3,2],resizefactor = 2,cam = 13,outfile ='faceon_h258_v2.jpeg'

pro make_jpeg_h516
;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512.1'
;make_jpeg,cam = 15, scales = [1.5,4.2,1.2]*0.5,outfile = 'angleon.jpeg'

;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir/h516.cosmo25cmb.1536g6MbwK.00512.1'
;make_jpeg,cam = 15, scales = [1.5,4.2,1.2]*0.5,outfile = 'angleon.jpeg'

;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.1'
;make_jpeg,cam = 15, scales = [1.5,4.2,1.2]*0.5,outfile = 'angleon.jpeg'

;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464.dir/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464.1'
;make_jpeg,cam = 15, scales = [1.5,4.2,1.2]*0.5,outfile = 'angleon.jpeg'

;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK.00512.1'
;make_jpeg,cam = 15, scales = [1.5,4.2,1.2]*0.5,outfile = 'h516.cosmo25cmb.1536g14HBWK.00512.angleon.jpeg'

;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00272.dir/h516.cosmo25cmb.2304g14HBWK.00272.1'
;make_jpeg,cam = 15, scales = [1.5,4.2,1.2]*0.5,outfile = 'h516.cosmo25cmb.2304g14HBWK.00512.angleon.jpeg'

;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.2304g/h603.cosmo50cmb.2304g14HBWK/steps/h603.cosmo50cmb.2304g14HBWK.00272.dir/h603.cosmo50cmb.2304g14HBWK.00272.1'
;make_jpeg,cam = 15, scales = [1.5,4.2,1.2]*0.5,outfile = 'h603.cosmo25cmb.2304g14HBWK.00272.angleon.jpeg'

cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir/h516.cosmo25cmb.1536g1MBWK.00512.1'
make_jpeg,cam = 15, scales = [1.5,4.2,1.2]*0.5,outfile = 'h516.cosmo25cmb.1536g1MBWK.00512.angleon.jpeg'

cd,'~/code/'
end

pro make_jpeg_iso
cd,'/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e6g2/MW_disk.00002.sunrise'
bands = [5,4,3]
scales = [4.9,5.7,7.8]
cam = 14
outfile = ' ~/Scratch2/MolecH/results/MW_disk.00002.faceon.jpeg'
make_jpeg,cam = cam, scales = scales, outfile = outfile
cam = 16
outfile = ' ~/Scratch2/MolecH/results/MW_disk.00002.edgeon.jpeg'
make_jpeg,cam = cam, scales = scales, outfile = outfile
cd,'/astro/users/christensen/code/Sunrise_analysis'
end

pro make_jpeg_master

make_jpeg,cam = 14, scales = [1.5,4.2,1.2]*0.5,outfile = 'faceon.jpeg'
make_jpeg,cam = 15, scales = [1.5,4.2,1.2]*0.5,outfile = 'angleon.jpeg'
make_jpeg,cam = 16, scales = [1.5,4.2,1.2]*0.5,outfile = 'edgeon.jpeg'


end

;[17,18,19,20]
;[20,19,17]

;bands = [5,4,3]
;scales = [4.9,5.7,7,8]
; scales = [1.5,4.2,1.2]*0.5
;cam = 14;12
;outfile = '~/h277.jpeg'

pro make_jpeg,bands=bands,scales=scales,nonlinearity=nonlinearity,cam = cam,outfile=outfile,resizefactor=resizefactor,rotateAngle = rotateAngle,range = range,center = center
; default scales = [4.9,5.7,7.8] as set in nw_scale_rgb.pro
; default nonlinearity = 3 as set in nw_arcsinh_fit.pro
if not keyword_set(bands) then bands = [5,4,3] 
;default bands SDSS irg, ie rgb=irg, 5,4,3 in our standard filter list
if not keyword_set(cam) THEN cam =10; face on 
if not keyword_set(outfile) then outfile = 'image'+STRTRIM(cam,2)+'.jpg'
if not keyword_set(resizefactor) THEN resizefactor= 1.0
if not keyword_set(rotateAngle) THEN rotateAngle = 0
if not keyword_set(center) THEN center = [0,0]
; read data
;im= mrdfits('/home/christensen/Storage1/UW/MolecH/Cosmo/h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/steps/h277.cosmo50cmb.3072g14HMbwK.00120.00017.dir/h277.cosmo50cmb.3072g14HMbwK.00120.00017.1/broadband.fits',cam,sr_head)
im=mrdfits('broadband.fits',cam,sr_head)
dxy = sr_head[where(strcmp('CD1_1',sr_head,5))]
dxy = double((strsplit(dxy,' ',/extract))[2])

sr_size = (size(im))[1] ;500 ;480.0
sr_range = dxy*sr_size ;headerHI.NAXIS1*headerHI.CDELT1
sr_axis = indgen(sr_size)
IF KEYWORD_SET(range) THEN BEGIN
    xmin = sr_size/2 + center[0] - range/2.0/dxy
    xmax = sr_size/2 + center[0] + range/2.0/dxy
    ymin = sr_size/2 + center[1] - range/2.0/dxy
    ymax = sr_size/2 + center[1] + range/2.0/dxy
ENDIF ELSE BEGIN
    xmin = 0
    xmax = sr_size
    ymin = 0
    ymax = sr_size
ENDELSE
print,dxy
;tmp= size(im,/dimensions)
;nx= tmp[0]
;ny= tmp[1]

indx = where(sr_axis ge xmin AND sr_axis lt xmax)
indy = where(sr_axis ge ymin AND sr_axis lt ymax)
nx = N_ELEMENTS(indx)
ny = N_ELEMENTS(indy)

RGBim= fltarr(nx,ny,3)
RGBim[*,*,0]= ROTATE(im(indx,indy,bands[0]),rotateAngle)
RGBim[*,*,1]= ROTATE(im(indx,indy,bands[1]),rotateAngle)
RGBim[*,*,2]= ROTATE(im(indx,indy,bands[2]),rotateAngle)
; rebin
RGBim= rebin(RGBim,floor(nx*resizefactor),floor(ny*resizefactor),3)
;RGBim_test = RGBim
for i = 0,floor(ny*resizefactor)-1 Do BEGIN
    sub = RGBim[*,i,*]
    transRGB = transpose(sub,[2,0,1]) 
    color_conversion = [[2.0,0,       0],   $                 ;This row is RGB colors values for the first filter
                        [2.4, 1.6,     0],   $              ;This row is RGB colors values for the second filter
                        [0,       0, 5.0]]                    ;This row is RGB colors values for the third filter

    transRGB = MATRIX_MULTIPLY(color_conversion,transRGB)
    sub = transpose(transRGB)
    RGBim[*,i,*] = sub
ENDFOR

IF NOT keyword_set(scales) THEN scales=[0.4,1,0.4]
print,scales
; scale and set colors
RGBim = nw_scale_rgb(RGBim,scales=scales)
RGBim = nw_arcsinh_fit(RGBim,nonlinearity=nonlinearity)
RGBim = nw_fit_to_box(RGBim,origin=origin)
RGBim = nw_float_to_byte(RGBim)

; write
WRITE_JPEG,outfile,RGBim,TRUE=3,QUALITY=100

return
end
