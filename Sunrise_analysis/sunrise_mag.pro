;Returns the magnitude from Sunrise

FUNCTION sunrise_mag,filename,band = band,filtername = filtername,halo_str = halo_str

IF NOT keyword_set(halo_str) THEN halo_str = '1'
;****************** Sunrise Cube **************************\
sfilename = filename+'.' + halo_str + '/broadband.fits'
temp = mrdfits(sfilename,2,header)
fovstr = strsplit(header[13],' ',/extract)
fov = float(fovstr[3]) ;FOV in kpc
diststr = strsplit(header[12],' ',/extract)
distance = float(diststr[3]) ;Camera distance in kpc

IF NOT keyword_set(band) THEN band = 13
IF NOT keyword_set(filtername) THEN filtername = "i_SDSS.res"
filters = mrdfits(sfilename,band,filt_head)
filtind = where(strcmp(filters.filter,filtername,6))
RETURN,filters[filtind].ab_mag0
END
