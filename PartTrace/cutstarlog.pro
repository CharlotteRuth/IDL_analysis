PRO cutstarlog,filebase,finalid = finalid,molecularH = molecularH,verbose = verbose
IF NOT keyword_set(finalid) THEN finalid = '1'

;spawn,'ls *.starlog',file_starlog 
;spawn,'ls h*.param',pfile
;file_starlog = filebase + '.starlog'
file_starlog = 'starlog.' + finalid + '.fits'
pfile = filebase + '.param'
units = tipsyunits(pfile)
;IF strlen(file_starlog) GT 1 THEN starlog = rstarlog(file_starlog,molecularH = molecularH)
starlog = mrdfits(file_starlog,1)
halodat = mrdfits('grp' + finalid + '.alignment.fits',1)
halodat.xc = (halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0)*1000.0 - units.lengthunit/2.0 
halodat.yc = (halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0)*1000.0 - units.lengthunit/2.0 
halodat.zc = (halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0)*1000.0 - units.lengthunit/2.0 
;center = [[halodat.time],[halodat.z],[halodat.xc],[halodat.yc],[halodat.zc]]
;center[*,2] = center[*,2]*1000.0 - units.lengthunit/2.0 ;"center" is in Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
;center[*,3] = center[*,3]*1000.0 - units.lengthunit/2.0
;center[*,4] = center[*,4]*1000.0 - units.lengthunit/2.0
;a = [[halodat.xa],[halodat.ya],[halodat.za]]

;starlogcut = starlog[0]
starlog.timeform = starlog.timeform*units.timeunit/1e9
starlog.x = starlog.x;*units.lengthunit
starlog.y = starlog.y;*units.lengthunit
starlog.z = starlog.z;*units.lengthunit
uniq_time = starlog[uniq(starlog.timeform,sort(starlog.timeform))].timeform
rvirfit_uniq = spline(halodat.time,halodat.rvir,uniq_time)
xcfit_uniq = spline(halodat.time,halodat.xc,uniq_time)
ycfit_uniq = spline(halodat.time,halodat.yc,uniq_time)
zcfit_uniq = spline(halodat.time,halodat.zc,uniq_time)
match2,starlog.timeform,uniq_time,indsl,induniq
rvirfit = rvirfit_uniq[indsl]
xcfit = xcfit_uniq[indsl]
ycfit = ycfit_uniq[indsl]
zcfit = zcfit_uniq[indsl]
starlogind =intarr(n_elements(starlog)) 

inside = sqrt((starlog.x)^2 + (starlog.y)^2 + (starlog.z)^2) - rvirfit
;help,where(sqrt((starlog.x - xcfit)^2 + (starlog.y - ycfit)^2 + (starlog.z - zcfit)^2) LE rvirfit)
;stop
;FOR i = 0L, n_elements(starlog) -1 DO $
;  IF (sqrt((starlog[i].x - xcfit[i])^2 + (starlog[i].y - ycfit[i])^2 + (starlog[i].z - zcfit[i])^2) LE rvirfit[i]) THEN starlogind[i] = 1
starlogcut = starlog[where(inside LE 0)]
inhalo = where(inside LE 0)
;stop
;histogramp,sqrt((starlogcut.x - xcfit[inhalo])^2 + (starlogcut.y - ycfit[inhalo])^2 + (starlogcut.z - zcfit[inhalo])^2),max = 150
IF keyword_set(verbose) THEN histogramp,sqrt((starlogcut.x)^2 + (starlogcut.y)^2 + (starlogcut.z)^2),max = 150
mwrfits,starlogcut,'starlog.cut.' + finalid + '.fits',/create
END
