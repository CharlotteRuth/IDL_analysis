; PRO SMOOTH_CUBE
;
; Currently only uses a circular beam.
;
; Inputs:
;
;    infile: The input fits file to be convolved.
;
;    outfile: The output, convolved fits file
;
;    pixres: The pixel resolution in kpc
;
;    telres: The telescope FWHM in arcsec.
;
;    distance: The distance to the galaxy. in Mpc
;
;    size: The size of the smoothing kernel to use.
;
; 
; 
function smooth_cube, cube, header, telres1=telres1, telres2 = telres2, size1=size1, size2 = size2, $
                      kpcres1=kpcres1, kpcres2 = kpcres2, distance=distance, outfile = outfile, normalize = normalize
;Normalize: A flag you can turn on to normalize the smoothing gaussian
;such that it has an area of 1 (flux conserving).
;However, in radio astronomy this isn't the case.
;If you don't turn it on, you'll get radio astronomy-esq units
;(Msol/beam as opposed to Msol)
;Originally (before 7/26/11) the output was normalized from this code.
;However, this presents problems for comparing to the noise/sensetivity 

pixres = header.cdelt1

; process input..
if not keyword_set(kpcres1) then begin
    if not (keyword_set(telres1) and keyword_set(distance)) then begin
        print,"KPCRES (or TELRES and DISTANCE) must be set."
        return,cube
    endif
    kpcres1 = telres1 * distance * 1000. / 206265.
;we need to convert kpcres to the standard deviation of the beam (PSF)
;from the full width at half max (FWHM) of the beam 
    kpcres1 = kpcres1/(2.0*SQRT(2.0*alog(2.0)))
endif

;For elliptical gaussians, another axes is set with either kpcres2 or telres2
if keyword_set(telres2) then kpcres2 = telres2*distance*1000. / 206265.
if keyword_set(kpcres2) then kpcres2 = kpcres2/(2.0*SQRT(2.0*alog(2.0)))
;If neither kpcres2 nor telres2 is set, the gaussian is circular with
;the resolution in the second direction the same as in the first
if not keyword_set(kpcres2) then kpcres2 = kpcres1 

if not keyword_set(size1) then size1 = max([2.*ceil(double(kpcres1) / pixres),5])
if not keyword_set(size2) then size2 = max([2.*ceil(double(kpcres2) / pixres),5])

; make the smoothing kernel - just a axisymmetric gaussian for now.
z = shift(dist(2*size1 + 1,2*size2 + 1), size1, size2)
z = z * pixres
;z = exp( -(2. * z / kpcres)^2) Fixed by Adrienne, 7/26/11
z = exp( -(z^2 / 2. / kpcres1/kpcres2))
if keyword_set(normalize) then z = z / total(z)

dim = size(z, /dimension)

zp = reform(z, dim[0], dim[1], 1)

convolved = convol(cube, zp, /edge_truncate)

if keyword_set(outfile) then begin
;    header1 = {smoothheader, $
;                   filename: outfile, $
;                   naxis: 0, $
;                   naxis1: 0, $
;                   naxis2: 0, $
;                   naxis3: 0, $
;                   crval1: 0.0, $
;                   cdelt1: 0.0, $
;                   crpix1: 0.0, $
;                   crval2: 0.0, $
;                   cdelt2: 0.0, $
;                   crpix2: 0.0, $
;                   crval3: 0.0, $
;                   cdelt3: 0.0, $
;                   crpix3: 0.0, $
;                   bscale: 0.0, $
;                   bzero: 0.0, $
;                   blank: 0.0, $
;                   kpcunit: 0.0, $
;                   munit: 0.0, $
;                   kpcres: 0.0 $
;                  }
;        struct_assign,header,header1
;        header1.kpcres = kpcres
       ; stop
        mwrheader = makeheader(header)
        mwrfits,convolved,outfile,mwrheader
endif
return,convolved

end
