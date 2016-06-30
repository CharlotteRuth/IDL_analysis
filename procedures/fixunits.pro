PRO fixunits,h,g,d,s,units
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
dens_convert = units.massunit * gm_per_msol * amu_per_gm/units.lengthunit^3/cm_per_kpc^3

g.x = g.x*units.lengthunit*h.time
g.y = g.y*units.lengthunit*h.time
g.z = g.z*units.lengthunit*h.time
g.vx = g.vx*units.vunit*h.time
g.vy = g.vy*units.vunit*h.time
g.vz = g.vz*units.vunit*h.time
g.mass = g.mass*units.massunit
g.dens = g.dens*dens_convert/h.time^3

s.x = s.x*units.lengthunit*h.time
s.y = s.y*units.lengthunit*h.time
s.z = s.z*units.lengthunit*h.time
s.vx = s.vx*units.vunit*h.time
s.vy = s.vy*units.vunit*h.time
s.vz = s.vz*units.vunit*h.time
s.mass = s.mass*units.massunit

d.x = d.x*units.lengthunit*h.time
d.y = d.y*units.lengthunit*h.time
d.z = d.z*units.lengthunit*h.time
d.vx = d.vx*units.vunit*h.time
d.vy = d.vy*units.vunit*h.time
d.vz = d.vz*units.vunit*h.time
d.mass = d.mass*units.massunit

END
