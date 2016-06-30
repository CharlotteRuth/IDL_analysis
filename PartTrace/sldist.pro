PRO sldist,dir,haloid,base,r,sf,nbins = nbins
pfile = dir + '/' + base + '.param'
units = tipsyunits(pfile)
halodat = mrdfits(dir + '/grp' + haloid + '.alignment.fits',1)
halodat.xc = (halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0)*1000.0 - units.lengthunit/2.0 
halodat.yc = (halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0)*1000.0 - units.lengthunit/2.0 
halodat.zc = (halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0)*1000.0 - units.lengthunit/2.0 
filename = dir + '/starlog.' + haloid + '.fits'

sl = mrdfits(filename,1)

IF (abs(median(sl.x)) LT 10 AND abs(median(sl.y)) LT 10 AND abs(median(sl.z)) LT 10 AND stdev(sl.x) LT 50 AND stdev(sl.y) LT 50 AND stdev(sl.z) LT 50) THEN BEGIN 
    xcfit = fltarr(n_elements(sl))
    ycfit = fltarr(n_elements(sl))
    zcfit = fltarr(n_elements(sl))
    uniq_a = fltarr(n_elements(sl)) + 1
ENDIF ELSE BEGIN
    uniq_time = sl[uniq(sl.timeform,sort(sl.timeform))].timeform
    IF max(uniq_time) LT 5 THEN uniq_time = uniq_time*units.timeunit/1e9
    zarr = reverse(10^(findgen(1000)/500) - 1.005)
    tarr = (wmap3_lookback(10000) - wmap3_lookback(zarr))/1e9
    uniq_z = spline(tarr,zarr,uniq_time)
    uniq_a = 1/(1 + uniq_z)
    rvirfit_uniq = spline(halodat.time,halodat.rvir,uniq_time)
    xcfit_uniq = spline(halodat.time,halodat.xc,uniq_time)
    ycfit_uniq = spline(halodat.time,halodat.yc,uniq_time)
    zcfit_uniq = spline(halodat.time,halodat.zc,uniq_time)
    match2,sl.timeform,uniq_time,indsl,induniq
    xcfit = xcfit_uniq[indsl]
    ycfit = ycfit_uniq[indsl]
    zcfit = zcfit_uniq[indsl]
    uniq_a = uniq_a[indsl]
ENDELSE

radius = sqrt((sl.x - xcfit)^2 + (sl.y-  ycfit)^2)*uniq_a
sf = weighted_histogram(radius,weight = sl.massform,locations = r,max = 50,nbins = 500)
plot,r,sf,psym = 10,title = base + '.'  + haloid ,xrange = [0,5];,yrange = [0,1.3e-8]
;stop
END
