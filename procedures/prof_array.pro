
;Code to make a profile weighted according to some value
FUNCTION prof_array,r,mass,weight = weight,nbins = nbins, rmax = rmax, rmin = rmin

IF NOT keyword_set(nbins) THEN nbins = 100
IF NOT keyword_set(rmin) THEN rmin = 0
IF NOT keyword_set(rmax) THEN rmax = max(r)
IF NOT keyword_set(weight) THEN weight = fltarr(n_elements(mass)) + 1

dr = (float(rmax) - float(rmin))/nbins
rinner = findgen(nbins)*dr + rmin
rbins = rinner + dr/2.
router = rinner + dr

binarea = !pi*(router^2 - rinner^2)

profile = {rbins: rbins, rinner: rinner, router: router, binarea: binarea, inclosed: fltarr(nbins), mass: fltarr(nbins), sd: fltarr(nbins), mean: fltarr(nbins), stdev: fltarr(nbins)}

profile.inclosed = weighted_histogram(r,weight = weight*mass,max = rmax,min = rmin, nbins = nbins,/cum)
profile.mass = weighted_histogram(r,weight = weight*mass,max = rmax,min = rmin, nbins = nbins)
profile.sd = weighted_histogram(r,weight = weight*mass,max = rmax,min = rmin, nbins = nbins)/binarea

profile.mean = profile.mass/weighted_histogram(r,weight = mass,max = rmax,min = rmin, nbins = nbins)
IF (where(profile.mass EQ 0))[0] NE -1 THEN profile.mean[where(profile.mass EQ 0)] = 0

RETURN,profile
END
