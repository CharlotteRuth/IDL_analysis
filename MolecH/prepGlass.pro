;.r /astro/users/christensen/code/Sunrise_analysis/Helios/Sunrise/wtipsy.pro

PRO prepGlass
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
dens_convert =  gm_per_msol * 5.9753790e+23/cm_per_kpc^3
cd,'/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/'
rtipsy,'shocktube_glass_30_2_2.std',h,g,d,s
initT_highden = 1e8
initT_smallden = 0.6e4
name_highden = 'glass_highd' 
name_smallden  = 'glass_smalld'
mid_smallden = 10.0
mid_highden  =-10.0
init_range =  2.0

tmass_highden = 2.75e8; solar masses
absden_highden = 0.8 ;cm^-3 
;kpc = 1.0/init_range ;system length = 1kpc
time = 1.5e6
zmetal = 0;0.025
nbins = 100.0
dt = (Max(g.x) - Min(g.x))/nbins
x = findgen(nbins)/nbins*(Max(g.x) - Min(g.x)) + Min(g.x)
temp = fltarr(nbins)
dens = fltarr(nbins)

FOR i = 0, nbins - 1 DO BEGIN
    ind = where(g.x GE x[i] AND g.x LT (x[i] + dt))
    temp[i] = MEAN(g[ind].tempg)
    dens[i] = MEAN(g[ind].dens)    
ENDFOR
window,0
plot,x,temp,xtitle = 'x',ytitle = 'Temperature (T)',yrange = [0,3],xrange = [-15,15]
window,1
plot,x,dens,xtitle = 'x',ytitle = 'Density',yrange = [0,1.2],xrange = [-15,15]

g_smallden = g[where(g.x LE mid_smallden + init_range/2.0 AND g.x GE mid_smallden - init_range/2.0)]
g_highden = g[where(g.x LE mid_highden + init_range/2.0 AND g.x GE mid_highden - init_range/2.0)]
msol = tmass_highden/TOTAL(g_highden.mass)
;msol = absden_highden/(MEAN(g_highden.dens)*dens_convert) *kpc*kpc*kpc
kpc = (msol*TOTAL(g_highden.mass)/absden_highden*dens_convert)^(1./3.)
dDelta = time/get_timestep(msol,kpc)

g_highden.x = (g_highden.x - mid_highden)/init_range
g_highden.y =  g_highden.y/init_range
g_highden.z =  g_highden.z/init_range
g_highden.tempg = g_highden.tempg*0.0 + initT_highden
;mass_highden = TOTAL(g_highden.mass)
;msol_highden = 1.0/mass_highden
g_highden.zmetal = g_highden.zmetal * 0.0 + zmetal

g_smallden.x = (g_smallden.x - mid_smallden)/init_range
g_smallden.y =  g_smallden.y/init_range
g_smallden.z =  g_smallden.z/init_range
g_smallden.tempg = g_smallden.tempg*0.0 + initT_smallden
;mass_smallden = TOTAL(g_smallden.mass)
;msol_smallden = 1.0/mass_smallden
g_smallden.zmetal = g_smallden.zmetal * 0.0 + zmetal

h_highden = h
h_smallden = h

h_highden.n = N_ELEMENTS(g_highden)
h_highden.ngas = N_ELEMENTS(g_highden)
h_smallden.n = N_ELEMENTS(g_smallden)
h_smallden.ngas = N_ELEMENTS(g_smallden)

wtipsy,name_highden + '.std',h_highden,g_highden,d,s,/standard
wtipsy,name_smallden + '.std',h_smallden,g_smallden,d,s,/standard
writeprofile,name_highden,msol,kpc,dDelta
writeprofile,name_smallden,msol,kpc,dDelta

END


PRO writeprofile,name,msol,kpc,dDelta
basefile = '/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/base.param'
i = FILE_INFO(name + '.param')
IF i.exists THEN FILE_DELETE,name + '.param'
FILE_COPY,basefile,name + '.param'
openu,1,name + '.param',/APPEND

printf,1,'achInFile	= ',name + '.std'
printf,1,'achOutName	= ',name
printf,1,'dMsolUnit     = ',msol
printf,1,'dKpcUnit      = ',kpc
printf,1,'dDelta	= ',dDelta
close,1
print,'achInFile	= ',name + '.std'
print,'achOutName	= ',name
print,'dMsolUnit     = ',msol
print,'dKpcUnit      = ',kpc
print,'dDelta	= ',dDelta
END
