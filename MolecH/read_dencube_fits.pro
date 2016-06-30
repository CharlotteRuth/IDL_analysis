FUNCTION read_dencube_fits, infile, header, kpcunit=kpcunit, noscale = noscale
if not keyword_set(infile) then begin
    print,"Syntax:"
    print,""
    print,"read_cube_fits(infile, [header, mrdh])"
    print,""
    return,0.
endif

array = mrdfits(infile, 0, mrdh, /silent)
finite = where(FINITE(array),complement = nonfinite)
IF (nonfinite[0] ne -1) then array[nonfinite] = 0
;array = double(10^(double(array)*1.0))
header = {cubedenheader, $
          naxis: 0, $
          naxis1: 0, $
          naxis2: 0, $
          crval1: 0.0, $
          cdelt1: 0.0, $
          crpix1: 0.0, $
          crval2: 0.0, $
          cdelt2: 0.0, $
          crpix2: 0.0, $
          bscale: 0.0, $
          bzero: 0.0, $
          blank: 0.0, $
          kpcunit: 0.0, $
          munit: 0.0 $
         }
header.naxis = sxpar(mrdh, 'NAXIS')
header.naxis1 = sxpar(mrdh, 'NAXIS1')
header.naxis2 = sxpar(mrdh, 'NAXIS2')

header.crval1 = sxpar(mrdh, 'CRVAL1')
header.cdelt1 = sxpar(mrdh, 'CDELT1')
header.crpix1 = sxpar(mrdh, 'CRPIX1')

header.crval2 = sxpar(mrdh, 'CRVAL2')
header.cdelt2 = sxpar(mrdh, 'CDELT2')
header.crpix2 = sxpar(mrdh, 'CRPIX2')

header.bscale = sxpar(mrdh, 'BSCALE')
header.bzero = sxpar(mrdh, 'BZERO')
header.blank = sxpar(mrdh, 'BLANK')

IF NOT KEYWORD_SET(noscale) THEN BEGIN
; scale the array
    array = array * header.bscale + header.bzero
    array = 10.0^(float(array))

; set blanks to 0.
    vals = where(array EQ 10^(float(header.blank * header.bscale + header.bzero)))
    IF (vals[0] NE -1) THEN array[vals] = 0.0
ENDIF

IF keyword_set(kpcunit) THEN BEGIN
    header.crval1 = header.crval1 * kpcunit
    header.crval2 = header.crval2 * kpcunit
    header.cdelt1 = header.cdelt1 * kpcunit
    header.cdelt2 = header.cdelt2 * kpcunit
ENDIF

RETURN,array
END
