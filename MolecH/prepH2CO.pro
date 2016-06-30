pro setup_prepH2CO
;base ='/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e6g2/MW_disk'
base ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g2HBWK_cool/h516.cosmo25cmb.1536g2HBWK'
step = '00512'
;msol_per_sysmass = 1.36e17
;kpc_per_syslength = 1e5
msol_per_sysmass = 2.310e15 
kpc_per_syslength = 25000.
prepH2CO,base = base,step = step,msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength
!P.thick = 1.5
!P.CHARSIZE = 1.5
!X.Charsize = 1.35
!Y.Charsize = 1.35
end

PRO prepH2CO,base = base, step = step, msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength

set_plot,'x'
loadct,39

cm_per_kpc = 3.0857d21
gm_per_msol = 1.989d33
amu_per_gm = 6.022142d23
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)

c = 1.557e11 
lengthunit =  kpc_per_syslength*cm_per_kpc ;system length unit in cm (=1kpc)
massunit   = msol_per_sysmass*gm_per_msol ;system mass unit in gm
dens_convert =  msol_per_sysmass * gm_per_msol/kpc_per_syslength^3/cm_per_kpc^3 

filename = base+'.'+step
rtipsy,filename,h,g,d,s
t = g.tempg
rho = g.dens*dens_convert
zmetal = g.zmetal

readcol,filename+".HI",HI_prime,/silent
HI = HI_prime[1:h.ngas]
readcol,filename+".H2",h2_prime,/silent
h2 = h2_prime[1:h.ngas]
readcol,filename+'.mach_sheer',mach,/silent
length = 1.0/mach[1:h.ngas]*kpc_per_syslength*cm_per_kpc
readcol,filename+".OxMassFrac",Ofrac_prime,/silent
Ofrac = Ofrac_prime[1:h.ngas]
readcol,filename+".FeMassFrac",Fefrac_prime,/silent
Fefrac = Fefrac_prime[1:h.ngas]
nh2 = h2*2.0
xuv = fltarr(N_ELEMENTS(HI)) + 250.0
co = calcCO(rho*amu_per_gm,g.zmetal,length,xuv)

h2_mass = HI_prime
co1 = HI_prime
t_mass = HI_prime
h2_mass[1:h.ngas] = nH2*g.mass*msol_per_sysmass
co1[1:h.ngas] = reform(co[*,0])*g.mass*msol_per_sysmass/dens_convert
t_mass[1:h.ngas] = g.mass*msol_per_sysmass
t_mass[h.ngas + 1:h.ngas + h.ndark] = d.mass*msol_per_sysmass
t_mass[h.ngas + h.ndark + 1:h.ngas + h.ndark + h.nstar] = s.mass*massunit
openw,1,base+"."+step+".H2_Mass"
openw,2,base+"."+step+".CO_Mass"
openw,3,base+"."+step+".Total_Mass"
printf,1,LONG(HI_prime[0]),FORMAT = '(I)'
printf,2,LONG(HI_prime[0]),FORMAT = '(I)'
printf,3,LONG(HI_prime[0]),FORMAT = '(I)'
for i=1LL,LONG64(N_ELEMENTS(HI_prime)) - 1 do  BEGIN
    printf,1,h2_mass[i]
    printf,2,co1[i]
    printf,3,t_mass[i]
ENDFOR
close,1
close,2
close,3

;plot,rho*amu_per_gm ,co[*,1]/1.73d-26,xrange = [10,1e3],yrange = [1e15,1e19],psym = 1,/xlog,/ylog,ytitle = 'CO Column Density (cm^-2)',xtitle = 'Density'
;plot,rho*amu_per_gm ,co[*,0],xrange = [10,1e3],psym = 1,/xlog,/ylog,ytitle = 'CO Density (cm^-3)',yrange = [1e-4,1],xtitle = 'Density'
;plot,zmetal[sortz]/zmetal_MW,co[sortz,1]/1.73d-26,xtitle = 'Relative Metallicity',/ylog,psym = 1,xrange = [1.69,1.71],ytitle = 'Expected Column Density (cm^-2)'
;plot,nH2_Col,co[*,1]/1.73d-26,psym = 1,/ylog,/xlog,xtitle = 'H2 Column Density',ytitle = 'CO Column Density',yrange=[1e14,1e20],xrange=[1e20,1e22]
;plot,rho*amu_per_gm,h2*length*rho*amu_per_gm,psym = 3,/xlog,/ylog,xtitle = 'Density',ytitle = 'Column Density of H2',xrange = [10,1e3],yrange = [1e16,1e23],ystyle = 1
;oplot,rho[sortrho]*amu_per_gm ,co[sortrho,1]/1.73d-26,psym = 1
;plot,rho*amu_per_gm,(nH2_Col)/co[*,1],psym = 1,/ylog,/xlog,xrange=[10,1000],xtitle = 'Density [cm^3]',ytitle = 'N_{H2}/I_{CO}'

END
