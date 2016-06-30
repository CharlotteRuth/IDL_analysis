pro getimage_master

cam = 14
band = 6 - 1
filt = 13
dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512.1/'
dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.1'
dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir/h516.cosmo25cmb.1536g6MbwK.00512.1'
file = 'broadband.fits'


cam = 14
band = 8 - 1
filt = 13
dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512.1'
file = 'broadband.fits'

spawn,'pwd',home
cd,dir
getimage,file,   filt = filt, cam = cam, band = band, mask = 1
stop
profilefit, file, filt = filt, cam = cam, band = band, mask = 1, /write_feedme, fitname = 'double'
galfitname = 'galfit.04'
stop
checkfit,file, cam, filt, band, galfitname, rad = rad
cd,home
end

;cd,'/home/christensen/Storage2/UW/MolecH/Cosmo/h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/steps/h285.cosmo50cmb.3072g14HMbwK.00512.dir/h285.cosmo50cmb.3072g14HMbwK.00512.1'
;H band, face-on with dust
;getimage,'broadband.fits',filt = 13, cam = 14, band = 7


pro getimage, file, filt = filt, cam = cam, band = band, mask = mask
IF not keyword_set(filt) THEN filt = 12    ;index number of filters in fv 
IF not keyword_set(cam) THEN cam = 13 ; this camera is face on with dust (adjust number to get this) 
IF not keyword_set(band) THEN band= 5 ; Filter number, starting at 0

temp = mrdfits(file,2,header)
fovstr = strsplit(header[13],' ',/extract)
fov = FLOAT(fovstr[3])

image = mrdfits(file,cam,header)  ; 
filters=mrdfits(file,13)
Lunit = filters[band].ewidth_lambda/filters[band].ewidth_nu   ;m/Hz
Fo=3.631e-23  ;zero point of AB magnitude system, in W/m^2/Hz
;units = 4.25e10 ; sr>>> arcsec^2  
;sbfactor = Lunit/units/Fo
sbfactor = Lunit/Fo   ;units of W/m/m^2
img = image(*,*, band)
img=sbfactor*img
outfile = 'gal.'+STRTRIM(filters[band].filter,2)+'.fits'
mwrfits,img,outfile,/create

gal = image(*,*, band)
npix = n_elements(gal[0,*])
makex, gal,x,y  ;set up coordinate grid

pixKpc = 1.0*fov/(1.0*npix)
dummy = max(gal,nn)
max=array_indices(gal,nn)-npix/2
x=x-max[0]
y=y-max[1]
radius = sqrt((x)^2+(y)^2)*pixKpc ; ciruclar annuli

img = img*0
img[where(radius le 1)] = 1
outfile = 'gal.mask.fits'
mwrfits,img,outfile,/create
end

