pro extractImage, file, cam, band, outfile=outfile

;read, 'Please enter camera (e.g., 14 for face on, 16 for edge on): ', cam
;read, 'Please enter band (e.g., 5 for i band): ', band

;file = 'broadband.fits'
;cam = 14 ; this camera is face on with dust (adjust number to get this)
;band= 5 ; Filter number, starting at 0 (currently SDSS i band)
;24  H

filter_num = 13

image = mrdfits(file,cam)  ;
filters=mrdfits(file,filter_num,header)
;fstr = strsplit(header[filter_num],' ',/extract)
fstr=(strsplit(filters[band].filter,'.',/extract))[0]
Lunit = filters[band].ewidth_lambda/filters[band].ewidth_nu   ;m/Hz
Fo=3.631e-23  ;zero point of AB magnitude system, in W/m^2/Hz
sbfactor = Lunit/Fo   ;units of W/m/m^2

;Units of Sunrised image are W/m/m^2/sr.  Multiplying by sb (below)
;puts you on the zero point of the AB mag system, rather than
;arbitrary values

if NOT keyword_set(outfile) then outfile1 = 'gal.fits' else outfile1 = outfile + '_' + fstr + '.fits' 
img = image(*,*, band)
img=sbfactor*img
mwrfits,img,outfile1

end

