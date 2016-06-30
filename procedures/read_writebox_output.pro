;******************************************************************************
;This procedure returns a structure containing the properties of the
;gas particles described in the output of the tipsy procedure
;writebox.
;******************************************************************************
function read_writebox_output, filename

;Open tipsy ascii file.
openr, 2, filename

;Read in the first few lines, which tells you the number of total, gas,
;and star particles to follow.
ntotal=LONG(0) & ngas=LONG(0) & nstar=LONG(0)
time=DOUBLE(0.0)
readf, 2, ntotal, ngas, nstar
readf, 2, dim 
readf, 2, time

ndark = ntotal - ngas - nstar

;Read in the rest of the file.
numlines = (ntotal * 9) + (ngas * 3) + (nstar * 2)
data = dblarr(numlines)
readf, 2, data
close, 2
free_lun, 2

;There are various quantities for each gas, star, and dark particle in
;the tipsy ascii files: mass, x, y, z, vx, vy, vz, glen, rho, T, slen,
;Z.  Declare an array structures for all particles, including 'type'
;to distinguish between 'star', 'dark' and 'gas'
catalog = replicate({mass: 0.d0, x: 0.d0, y: 0.d0, z: 0.d0, vx: 0.d0, vy: 0.d0, vz: 0.d0, rho: 0.d0, temp: 0.d0, metallicity: 0.d0, timeform: 0.d0, type: 'gas'}, ntotal)

;define the types (in the same order as tipsy), but have to beware of
;timesteps before stars have formed
IF (nstar EQ 0) THEN BEGIN
    catalog.type=[REPLICATE('gas',ngas),REPLICATE('dark',ndark)] 
ENDIF ELSE BEGIN
    catalog.type=[REPLICATE('gas',ngas),REPLICATE('dark',ndark),REPLICATE('star',nstar)] 
ENDELSE

;Fill up the array with values defined for each particle type
catalog.mass = data[0 : ntotal - 1]
catalog.x = data[ntotal : 2 * ntotal - 1]
catalog.y = data[2* ntotal : 3 * ntotal - 1]
catalog.z = data[3* ntotal : 4 * ntotal - 1]
catalog.vx = data[4 * ntotal : 5 * ntotal - 1]
catalog.vy = data[5 * ntotal: 6 * ntotal - 1]
catalog.vz = data[6 * ntotal: 7 * ntotal - 1]

;fill in gas only parameters
catalog[0:ngas-1].rho = data[7 * ntotal + ndark + nstar : 8 * ntotal - 1]
catalog[0:ngas-1].temp = data[8 * ntotal : 8 * ntotal + ngas - 1]
catalog[0:ngas-1].metallicity = data[8 * ntotal + 2 * ngas : 8 * ntotal + 3 * ngas - 1]

;fill in star only parameters
IF (nstar GT 0) THEN BEGIN
    catalog[ngas+ndark:ntotal-1].metallicity = data[8 * ntotal + 3 * ngas : 8 * ntotal + 3 * ngas + nstar - 1]
    catalog[ngas+ndark:ntotal-1].timeform = data[8 * ntotal + 3 * ngas + nstar: 8 * ntotal + 3 * ngas + 2 * nstar - 1]
ENDIF

;NOTE: catalog still doesn't include glen(?) and slen(?) cuz I don't
;know what they are...

return, catalog

end
