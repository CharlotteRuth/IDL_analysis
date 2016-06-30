;
;
;
; calculate angular momentum for the input particles
;
;


pro angmom, p, jvec, lvec, jtot, ltot

rvec = fltarr(n_elements(p),3)
vvec = fltarr(n_elements(p),3)

rvec[*,0] = p.x
rvec[*,1] = p.y
rvec[*,2] = p.z

vvec[*,0] = p.vx
vvec[*,1] = p.vy
vvec[*,2] = p.vz

jvec = crossp_multi(rvec, vvec)

rvec[*,0] = p.x*p.mass
rvec[*,1] = p.y*p.mass
rvec[*,2] = p.z*p.mass

lvec = crossp_multi(rvec, vvec)

jtot = sqrt(jvec[*,0]^2. + jvec[*,1]^2. + jvec[*,2]^2.)

ltot = sqrt(lvec[*,0]^2. + lvec[*,1]^2. + lvec[*,2]^2.)

end
