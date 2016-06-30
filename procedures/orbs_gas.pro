pro orbs_gas,inFile1,inFile2,outFile,p,r
;Puts 2 gals onto parabolic orbits
;provides initial conditions for a merger run
;transcription of orbits.f from Fabio
if (N_PARAMS() eq 0) then begin
    print, "orbs.pro Reads 2 tipsy files and puts them on parabolic orbit: "
    print,"Usage: "
    print,"       orbs,file1,file2,outfile,p,r"
    print,"file1 1st InputTipsy file"
    print,"file2 2nd InputTipsy File"
    print,"outfile output file name"
    print,"p, pericenter distance of merger"
    print,"r, position on orbit: typically r=twice halo radius-->galaxy halos just touching"
    return
endif



;;Read in 1st gal
rtipsy,inFile1,h1,g1,d1,s1
rtipsy,inFile2,h2,g2,d2,s2


;;;Calculate orbital parameters
;Get Mass of each galaxy
m1=total(d1.mass)+total(s1.mass)+total(g1.mass)
m2=total(d2.mass)+total(s2.mass)+total(g2.mass)
v=sqrt(2.*(m1+m2)/r)
f=acos(((2.*p)/r)-1.)
phi=(!pi/2.)+acos(p/r)
theta=phi-f
;;Calculate positions and velocities
xp=r*cos(f)
yp=-r*sin(f)
zp=0
xdotp=v*cos(theta)
ydotp=v*sin(theta)
zdotp=0
print,"xp  yp  xdotp  ydotp m1 m2"
print,xp,yp,xdotp,ydotp,m1,m2


;;Perturb 1st galaxy
d1.x=d1.x+xp
s1.x=s1.x+xp
g1.x=g1.x+xp
d1.y=d1.y+yp
s1.y=s1.y+yp
g1.y=g1.y+yp
d1.z=d1.z+zp
s1.z=s1.z+zp
g1.z=g1.z+zp

d1.vx=d1.vx+xdotp
s1.vx=s1.vx+xdotp
g1.vx=g1.vx+xdotp
d1.vy=d1.vy+ydotp
s1.vy=s1.vy+ydotp
g1.vy=g1.vy+ydotp
d1.vz=d1.vz+zdotp
s1.vz=s1.vz+zdotp
g1.vz=g1.vz+zdotp

;;;Get CM of entire system
;Mass of entire system
totMass=m1+m2
;print,totMass
;Get CM for x,y,z,vx,vy,vz
xcm=(total(d1.x*d1.mass)+total(s1.x*s1.mass)+total(g1.x*g1.mass)+total(d2.x*d2.mass)+total(s2.x*s2.mass)+total(g2.x*g2.mass))/totMass
ycm=(total(d1.y*d1.mass)+total(s1.y*s1.mass)+total(g1.y*g1.mass)+total(d2.y*d2.mass)+total(s2.y*s2.mass)+total(g2.y*g2.mass))/totMass
zcm=(total(d1.z*d1.mass)+total(s1.z*s1.mass)+total(g1.z*g1.mass)+total(d2.z*d2.mass)+total(s2.z*s2.mass)+total(g2.z*g2.mass))/totMass
vxcm=(total(d1.vx*d1.mass)+total(s1.vx*s1.mass)+total(g1.vx*g1.mass)+total(d2.vx*d2.mass)+total(s2.vx*s2.mass)+total(g2.vx*g2.mass))/totMass
vycm=(total(d1.vy*d1.mass)+total(s1.vy*s1.mass)+total(g1.vy*g1.mass)+total(d2.vy*d2.mass)+total(s2.vy*s2.mass)+total(g2.vy*g2.mass))/totMass
vzcm=(total(d1.vz*d1.mass)+total(s1.vz*s1.mass)+total(g1.vz*g1.mass)+total(d2.vz*d2.mass)+total(s2.vz*s2.mass)+total(g2.vz*g2.mass))/totMass
;;Subtract cm values from all positions and velocities
d1.x=d1.x-xcm
d2.x=d2.x-xcm
s1.x=s1.x-xcm
s2.x=s2.x-xcm
g1.x=g1.x-xcm
g2.x=g2.x-xcm

d1.y=d1.y-ycm
d2.y=d2.y-ycm
s1.y=s1.y-ycm
s2.y=s2.y-ycm
g1.y=g1.y-ycm
g2.y=g2.y-ycm

d1.z=d1.z-zcm
d2.z=d2.z-zcm
s1.z=s1.z-zcm
s2.z=s2.z-zcm
g1.z=g1.z-zcm
g2.z=g2.z-zcm

d1.vx=d1.vx-vxcm
d2.vx=d2.vx-vxcm
s1.vx=s1.vx-vxcm
s2.vx=s2.vx-vxcm
g1.vx=g1.vx-vxcm
g2.vx=g2.vx-vxcm

d1.vy=d1.vy-vycm
d2.vy=d2.vy-vycm
s1.vy=s1.vy-vycm
s2.vy=s2.vy-vycm
g1.vy=g1.vy-vycm
g2.vy=g2.vy-vycm

d1.vz=d1.vz-vzcm
d2.vz=d2.vz-vzcm
s1.vz=s1.vz-vzcm
s2.vz=s2.vz-vzcm
g1.vz=g1.vz-vzcm
g2.vz=g2.vz-vzcm

;;doublecheck CM shift OK
xcm=(total(d1.x*d1.mass)+total(s1.x*s1.mass)+total(g1.x*g1.mass)+total(d2.x*d2.mass)+total(s2.x*s2.mass)+total(g2.x*g2.mass))/totMass
ycm=(total(d1.y*d1.mass)+total(s1.y*s1.mass)+total(g1.y*g1.mass)+total(d2.y*d2.mass)+total(s2.y*s2.mass)+total(g2.y*g2.mass))/totMass
zcm=(total(d1.z*d1.mass)+total(s1.z*s1.mass)+total(g1.z*g1.mass)+total(d2.z*d2.mass)+total(s2.z*s2.mass)+total(g2.z*g2.mass))/totMass
vxcm=(total(d1.vx*d1.mass)+total(s1.vx*s1.mass)+total(g1.vx*g1.mass)+total(d2.vx*d2.mass)+total(s2.vx*s2.mass)+total(g2.vx*g2.mass))/totMass
vycm=(total(d1.vy*d1.mass)+total(s1.vy*s1.mass)+total(g1.vy*g1.mass)+total(d2.vy*d2.mass)+total(s2.vy*s2.mass)+total(g2.vy*g2.mass))/totMass
vzcm=(total(d1.vz*d1.mass)+total(s1.vz*s1.mass)+total(g1.vz*g1.mass)+total(d2.vz*d2.mass)+total(s2.vz*s2.mass)+total(g2.vz*g2.mass))/totMass

print,"xcm ycm zcm vxcm vycm vzcm for shifted system"
print,xcm,ycm,zcm,vxcm,vycm,vzcm

;;;Combine new shifted coordinates for both gals into a new
;;;structre-->tipsy output
;;Write new header structure
h1.time=0.0
h1.n=h1.n+h2.n
h1.ndim=3.
h1.ngas=h1.ngas+h2.ngas
h1.ndark=h1.ndark+h2.ndark
h1.nstar=h1.nstar+h2.nstar

;;write out new initial conditions w/ wtipsy
h=h1
g=[g1,g2]
d=[d1,d2]
s=[s1,s2]
help,h,/struc
help,g,/struc
help,d,/struct
help,s,/struct
print,"Writing Output"
wtipsy,outFile,h,g,d,s,/STAN

end
