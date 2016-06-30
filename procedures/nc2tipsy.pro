pro nc2tipsy,file,outdir=outdir,nogas=nogas,nodark=nodark,nostar=nostar

 h = { time:double(0.0), n:0L, ndim:0L, ngas:0L, ndark:0L, nstar:0L }
 if (keyword_set(nogas) EQ 0) then begin
  print,"Reading gas"
  rncarray,file+'/gas/mass',hgm,gm
  g=replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,dens:1.,tempg:1.,h : 1. , zmetal : 1., phi : 1.},hgm.n)
  h.ngas=hgm.n
  h.time=hgm.time
  h.ndim=hgm.ndim
  g.mass=gm
  rncvector,file+'/gas/pos',hgpos,gpos
  g.x=transpose(gpos[0,*])
  g.y=transpose(gpos[1,*])
  g.z=transpose(gpos[2,*])
  h.ndim=hgpos.ndim
  rncvector,file+'/gas/vel',hgvel,gvel
  g.vx=transpose(gvel[0,*])
  g.vy=transpose(gvel[1,*])
  g.vz=transpose(gvel[2,*])
  rncarray,file+'/gas/den',hgd,gd
  g.dens=gd
  rncarray,file+'/gas/temperature',hgt,gast
  g.tempg=gast
  rncarray,file+'/gas/smoothlength',hgt,gh
  g.h=gh
  rncarray,file+'/gas/OxMassFrac',hgox,gox
  rncarray,file+'/gas/FeMassFrac',hgfe,gfe
  g.zmetal=gfe+gox
  rncarray,file+'/gas/pot',hgphi,gphi
  g.phi=gphi
 endif else begin
  print,"Skipping gas"
  g=""
  h.ngas=0
 endelse

 if (keyword_set(nodark) EQ 0) then begin
  print,"Reading dark"
  rncarray,file+'/dark/mass',hdm,dm
  d=replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,eps : 1. , phi : 1.},hdm.n)
  d.mass=dm
  h.ndark=hdm.n
  h.time=hdm.time
  rncvector,file+'/dark/pos',hdpos,dpos
  d.x=transpose(dpos[0,*])
  d.y=transpose(dpos[1,*])
  d.z=transpose(dpos[2,*])
  h.ndim=hdpos.ndim
  rncvector,file+'/dark/vel',hdvel,dvel
  d.vx=transpose(dvel[0,*])
  d.vy=transpose(dvel[1,*])
  d.vz=transpose(dvel[2,*])
  rncarray,file+'/dark/smoothlength',hdt,dh
  d.eps=dh
  rncarray,file+'/dark/pot',hdphi,dphi
  d.phi=dphi
 endif else begin 
  print,"Skipping dark"
  d=""
  h.ndark=0
 endelse

 if (keyword_set(nostar) EQ 0) then begin
  print,"Reading stars"
  rncarray,file+'/star/mass',hsm,sm
  s=replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,metals:1.,tform:1.,eps : 1. , phi : 1.},hsm.n)
  s.mass=sm
  h.nstar=hsm.n
  h.time=hsm.time
  h.ndim=hsm.ndim
  rncvector,file+'/star/pos',hspos,spos
  ; WHAT A PAIN!
  s.x=transpose(spos[0,*])
  s.y=transpose(spos[1,*])
  s.z=transpose(spos[2,*])
  h.ndim=hspos.ndim
  rncvector,file+'/star/vel',hsvel,svel
  s.vx=transpose(svel[0,*])
  s.vy=transpose(svel[1,*])
  s.vz=transpose(svel[2,*])
  rncarray,file+'/star/OxMassFrac',hsox,sox
  rncarray,file+'/star/FeMassFrac',hsfe,sfe
  s.metals=sfe+sox
  rncarray,file+'/star/timeform',hstf,stf
  s.tform=stf
  rncarray,file+'/star/smoothlength',hst,sh
  s.eps=sh
  rncarray,file+'/star/pot',hsphi,sphi
  s.phi=sphi
 endif else begin
  print,"Skipping stars"
  s=""
  h.nstar=0
 endelse

h.n=h.ngas+h.ndark+h.nstar
print,"Writing standard tipsy file"
print,"header=",h
if (keyword_set(outdir) eq 0) then wtipsy,file+".std",h,g,d,s,/standard $
else wtipsy,outdir+"/"+file+".std",h,g,d,s,/standard
end
