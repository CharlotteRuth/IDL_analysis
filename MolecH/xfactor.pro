;base ='/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e6g2/MW_disk'
;step = '00004'
;msol_per_sysmass = 1.36e17
;kpc_per_syslength = 1e5
;xfactor,base = base,step = step,msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength
pro xfactor,base = base, step = step, msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength

cm_per_kpc = 3.0857d21
gm_per_msol = 1.989d33
amu_per_gm = 6.022142d23
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
lengthunit =  kpc_per_syslength*cm_per_kpc ;system length unit in cm (=1kpc)
massunit   = msol_per_sysmass*gm_per_msol ;system mass unit in gm
dens_convert =  msol_per_sysmass * gm_per_msol/kpc_per_syslength^3/cm_per_kpc^3 ;massunit*amu_per_gm/molec_weight/lengthunit^3 ;Converts to grams/cm^3
filename = base+'.'+step
rtipsy,filename,h,g,d,s
readcol,filename+".H2",h2_prime,/silent
readcol,filename+".CO",co_prime,/silent
readcol,filename+".I_CO",ico_prime,/silent
readcol,filename+".mach_sheer",mach,/silent
h2 = h2_prime[1:h.ngas]
co = co_prime[1:h.ngas]
ico = ico_prime[1:h.ngas]
length = 1.0/mach[1:h.ngas]*kpc_per_syslength*cm_per_kpc
t = g.tempg
rho = g.dens*dens_convert
zmetal = g.zmetal

sortrho = SORT(rho)
sortz = SORT(zmetal)
plot,rho[sortrho]*5.9753790d23,co[sortrho],psym = 3,/xlog,/ylog
stop

end
