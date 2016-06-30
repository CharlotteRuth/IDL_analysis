pro psfconvolve,file

psf = mrdfits('psf/result00_psf.fits')
data = mrdfits(file + '.fits')
psf = rebin(psf,(size(psf))[2]*6,(size(psf))[2]*6)
start = (size(data))[2]/2 - (size(psf))[2]/2 - 1
last = (size(psf))[2] - 1 + start
data = data[start:last,start:last]

window,0,xsize = 800,ysize = 800
ind0 = where(data lt 0)
IF NOT ind0[0] eq -1 THEN data[ind0] = 0
contour,alog10(data),nlevels = 245,/fill,min = -8  ;-12

smoothed = convolve(data,psf)

window,1,xsize = 800,ysize = 800
ind0 = where(smoothed lt 0)
IF NOT ind0[0] eq -1 THEN data[ind0] = 0
contour,alog10(smoothed),nlevels = 245,/fill,min = -8  ;-12

mwrfits,smoothed,file + '.psf.fits'
END
