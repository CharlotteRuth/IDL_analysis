;prep_halo,'/astro/net/scratch2/christensen/MolecH/MolecCloud/1E5R/1E5R.std','/astro/net/scratch2/christensen/MolecH/MolecCloud/1E5R/1E5R.gas.std'

;prep_halo,'/astro/net/scratch2/christensen/MolecH/MolecCloud/1E4R/1E4R.std','/astro/net/scratch2/christensen/MolecH/MolecCloud/1E4R/1E4R.gas.std'

;prep_halo,'/astro/net/scratch2/christensen/MolecH/MolecCloud/1E3R/1E3R.std','/astro/net/scratch2/christensen/MolecH/MolecCloud/1E3R/1E3R.gas.std'

;prep_halo,'/astro/net/scratch2/christensen/MolecH/MolecCloud/1E2R/1E2R.std','/astro/net/scratch2/christensen/MolecH/MolecCloud/1E2R/1E2R.gas.std'

PRO prep_halo,infile,outfile
msol_per_sysmass = 232510
kpc_per_syslength = 1.0
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)

zmetal = 0.025
total_He = (0.236 + 2.1*zmetal)/4.0
total_H = 1.0 - total_He*4.0 - zmetal

total_mass = 1e6
tempg = 2e4
density = 1e-6 ;in H atoms per c
density_solmass_kpc = density/gm_per_msol/H_per_gm*cm_per_kpc^3/total_H

rtipsy,infile,h,g,d,s
msol_per_sysmass = total_mass/total(g.mass)
radius = MAX(SQRT(g.x*g.x + g.y*g.y + g.z*g.z))
volume = 4.0/3.0*!PI*radius^3
kpc_per_syslength = (total_mass/volume/density_solmass_kpc)^(1.0/3.0)
print,'dMsolUnit: ',msol_per_sysmass
print,'dKpcUnit: ',kpc_per_syslength

dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3

g.tempg = g.tempg * 0.0 + tempg
g.zmetal = g.zmetal * 0.0 + zmetal
g.dens = g.dens * 0.0; + density/dens_convert/total_H
h.n = h.n - h.ndark
h.ndark = 0
wtipsy,outfile,h,g,d,s,/standard
END
