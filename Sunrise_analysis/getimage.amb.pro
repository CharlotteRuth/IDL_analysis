pro getimageamb, file, outfile=outfile, v4=v4

read, 'Please enter camera (e.g., 14 for face on, 16 for edge on, add +1 for v4): ', cam
read, 'Please enter band (e.g., 5 for i band): ', band
;file = The name of the broadband image file
;read, 'Please enter broadband file name: ', file 

;file = 'broadband.fits'
;cam = 14 ; this camera is face on with dust (adjust number to get this) 
;band= 5 ; Filter number, starting at 0 (currently SDSS i band)

image = mrdfits(file,cam)  ; 
;filters=mrdfits(file,11) ;h799.3072 run by Chris, fewer camera angles
filters=mrdfits(file,13)
if keyword_set(v4) then filters=mrdfits(file,14)
Lunit = filters[band].ewidth_lambda/filters[band].ewidth_nu   ;m/Hz
Fo=3.631e-23  ;zero point of AB magnitude system, in W/m^2/Hz
;units = 4.25e10 ; sr>>> arcsec^2  
;sbfactor = Lunit/units/Fo
sbfactor = Lunit/Fo   ;units of W/m/m^2

;Units of Sunrised image are W/m/m^2/sr.  Multiplying by sb (below)
;puts you on the zero point of the AB mag system, rather than 
;arbitrary values
if NOT keyword_set(outfile) then outfile = 'gal.fits'
img = image(*,*, band)
img=sbfactor*img
mwrfits,img,outfile

end

