pro getimage_JEM, file, outfile=outfile

;; This script reads in a broadband file and creates a FITS file 
;; for every bandpass in the Sunrise filter list
;;
;; Enter a general name for outfile (like 'h603.z0.') and code will
;; append the filter name and extension
;;
;; It is currently not set up to iterate through every camera angle,
;; but this can be arranged easily
print, 'Check broadband file for camera extenstion numbers'
read, 'Please enter camera (e.g., 14 for face on, 16 for edge on): ', cam
;read, 'Please enter band (e.g., 5 for i band): ', band
;file = ' '
;read, 'Please enter broadband file name: ', file 

;file = 'broadband.fits'
;cam = 14 ; this camera is face on with dust (adjust number to get this) 
;band= 5 ; Filter number, starting at 0 (currently SDSS i band)
filterlist=['FUV_GALEX','NUV_GALEX','u','g','r','i','z','H_2MASS','J_2MASS','Ks_2MASS','IRAC1', 'IRAC2','IRAC3','IRAC4','MIPS160','MIPS24','MIPS70','F435','F606','F775','F850','B','V','I','H','J','K']

image = mrdfits(file,cam)       ; 
filters=mrdfits(file,13)   ;This obviously changes run to run
;i=band
for i = 0, 26 do begin

;internal units to W/m^2 for Sunrise3
    Lunit = filters[i].ewidth_lambda/filters[i].ewidth_nu 
    Fo=3.631e-23          ;zero point of AB magnitude system, in W/m^2
;units = 4.25e10 ; sr>>> arcsec^2  
;sbfactor = Lunit/units/Fo
    sbfactor = Lunit/Fo
	outfile1=outfile
    outfile1=outfile+filterlist[i]+'.fits'
    outfile1=strcompress(outfile1)
    print, 'Outfile is ',outfile1
    if NOT keyword_set(outfile) then outfile1 = 'gal.fits'
    img = image(*,*, i)
    img=sbfactor*img
    mwrfits,img,outfile1

endfor
    
return
end

