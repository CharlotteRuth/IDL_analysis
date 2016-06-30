FUNCTION READ_OLD_RFORM, tipsyfile
;this function takes the name of the tipsyfile (NOT the rform file)
;and reads in the .rform file and makes sense of it.  

rtipsy,tipsyfile,h,/justhead
ns = h.nstar
ntot = h.n
ngd= h.ngas + h.ndark

filename=tipsyfile+'.rform'
rdfloat,filename,rform,skip=1
nparticle=N_ELEMENTS(rform)/3

r = replicate({x:0.,y:0.,z:0.},ns)
r.x = -rform[ngd:ntot-1]
r.y = rform[ntot+ngd:2*ntot-1]
r.z = rform[2*ntot+ngd:3*ntot-1]
return,r
END

