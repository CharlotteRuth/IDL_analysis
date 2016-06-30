;******************************************************************************
;This function reads in the first three lines of writebox tipsy output
;and returns the numbers of each particle and the expansion factor.
;******************************************************************************
FUNCTION get_dim, filename

;Open tipsy ascii file.
openr, 2, filename

;Read in the first few lines, which tells you the number of total, gas,
;and star particles to follow.
ntotal=LONG(0) & ngas=LONG(0) & nstar=LONG(0)
readf, 2, ntotal, ngas, nstar
readf, 2, dim 
readf, 2, time
close, 2
free_lun, 2
ndark = ntotal - ngas - nstar

;create structure with .ntotal, .ndark, .nstar, .ngas and .time (as a double)
;     note that time is actually 'a' the expansion factor!!!
dim={ntotal: ntotal, ndark: ndark, nstar: nstar, ngas: ngas, time: DOUBLE(time)}
RETURN, dim
END
