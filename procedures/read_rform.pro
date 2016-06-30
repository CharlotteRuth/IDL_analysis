FUNCTION READ_RFORM, tipsyfile
;this function takes the name of the tipsyfile (NOT the rform file)
;and reads in the .rform file and makes sense of it.  

rtipsy,tipsyfile,h,/justhead
ns = h.nstar
ntot = h.n
ngd= h.ngas + h.ndark

filename=tipsyfile+'.rform'
openr,lun,filename,/get_lun
readf,lun,numlines
rgd = fltarr(ngd*3)
readf,lun,rgd
rform = fltarr(ns*3)
readf,lun,rform
indrx = 3*lindgen(ns)
indry = 3*lindgen(ns)+1
indrz = 3*lindgen(ns)+2

r = replicate({x:0.,y:0.,z:0.},ns)

r.x = rform[indrx]
r.y = rform[indry]
r.z = rform[indrz]
return,r
END

