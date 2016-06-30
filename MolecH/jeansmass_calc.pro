function jeansmass_calc, density, temperature, molec_mass = molec_mass, zmetal = zmetal
;density is actual number density of H

fHII = 0.0
fHI = 1.0
fH2 = (1.0 - fHI - fHII)/2.0

fHeI = 1.0
fHeII = 0.0
fHeIII = 1 - fHeII - fHeII

IF NOT KEYWORD_SET(Zmetal) THEN Zmetal = 0.025 ;solar metallicity
k = 1.380658d-16 ;Boltzmann constant
mh = 1.673534d-24 ;hydrogen mass in grams
grav = 6.67259d-8 ;Graviational constant in cgs
cm_per_kpc = 3.0857d21
gm_per_msol = 1.989d33

if (ZMetal le 0.1) then  yHe = (0.236 + 2.1*ZMetal)/4.0 else yHe = (-0.446*(ZMetal - 0.1)/0.9 + 0.446)/4.0;
yH = 1.0 - yHe*4.0 - ZMetal; 
IF NOT KEYWORD_SET(molec_mass) THEN molec_mass = 1.0/(yH*(fHI + fH2 + 2.0*fHII) + yHe*(fHeI + fHeII*2.0 + fHeIII * 3.0)/4.0 + ZMetal/20) ;molecular weight in cloud with 90% H_2 and 10% He
N = density/yH
mdensity = N*mh*(yH + yHe*4.0 + ZMetal*2.1) ;mass density in gm/cm^3

return,(5.0*k*temperature/grav/molec_mass/mh)^(3.0/2.0)*(3.0/4.0/!PI/mdensity)^(1.0/2.0)/gm_per_msol
END


function jeanslength_calc, density, temperature, molec_mass = molec_mass
;density is actual number density of H

fHII = 0.0
fHI = 1.0
fH2 = (1.0 - fHI - fHII)/2.0

fHeI = 1.0
fHeII = 0.0
fHeIII = 1 - fHeII - fHeII

IF NOT KEYWORD_SET(Zmetal) THEN Zmetal = 0.025 ;solar metallicity
k = 1.380658d-16 ;Boltzmann constant
mh = 1.673534d-24 ;hydrogen mass in grams
grav = 6.67259d-8 ;Graviational constant in cgs
cm_per_pc = 3.0857d18
gm_per_msol = 1.989d33

if (ZMetal le 0.1) then  yHe = (0.236 + 2.1*ZMetal)/4.0 else yHe = (-0.446*(ZMetal - 0.1)/0.9 + 0.446)/4.0;
yH = 1.0 - yHe*4.0 - ZMetal; 
IF NOT KEYWORD_SET(molec_mass) THEN molec_mass = 1.0/(yH*(fHI + fH2 + 2.0*fHII) + yHe*(fHeI + fHeII*2.0 + fHeIII * 3.0)/4.0 + ZMetal/20) ;molecular weight in cloud with 90% H_2 and 10% He
N = density/yH

mdensity = N*mh*(yH + yHe*4.0 + ZMetal*2.1) ;mass density in gm/cm^3

return,(15.0*k*temperature/(4.0*!PI*grav*molec_mass*mh*mdensity))^(1.0/2.0)/cm_per_pc
END


PRO resolveJeans
minlogt = 1
maxlogt = 6
minlogrho = -6
maxlogrho = 3
nel = 100
minparticles = 32
minres_mass = [1e2,1e3,1e4,1e5,1e6]*minparticles
minres_force = [2,10,50,250,1250]

rtipsy,"/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_UV_solar_soft/MW_disk.00010",h,g,d,s
msol_per_sysmass = 1.36e17
kpc_per_syslength = 1e5
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3
total_H = fltarr(N_ELEMENTS(g))
total_He = fltarr(N_ELEMENTS(g))
total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
IF(comple[0] ne -1) THEN total_He[comple] = (-0.446*(g2[comple].zmetal - 0.1)/0.9 + 0.446)/4.0 ;
total_H = 1.0 - total_He*4.0 - g.zmetal ; 
rho_particle = alog10(g.dens * dens_convert * total_H)

temperature = findgen(nel)/(nel - 1)*(maxlogt - minlogt) + minlogt
density = findgen(nel)/(nel - 1)*(maxlogrho - minlogrho) + minlogrho
jeans_input_t = fltarr(nel,nel)
jeans_input_rho = fltarr(nel,nel)
FOR i = 0,nel*nel - 1 DO BEGIN
    jeans_input_t[i MOD nel,i/nel] = temperature[i/nel]
    jeans_input_rho[i MOD nel,i/nel] = density[i MOD nel]
ENDFOR
jeansmass = jeansmass_calc(10.0^jeans_input_rho,10.0^jeans_input_t)
jeanslength = jeanslength_calc(10.0^jeans_input_rho,10.0^jeans_input_t)

;window,0
set_plot,'ps'
device,filename = 'jeansmass.ps',/color,bits_per_pixel= 8,/times
contour,alog10(jeansmass),density,temperature,xtitle="Log(nH)",ytitle = 'Log(T)',/FILL,nlevels = 60,title = 'Log(Jeans Mass)',xrange=[minlogrho,maxlogrho],yrange=[minlogt,maxlogt],xstyle = 1,ystyle = 1
contour,alog10(jeansmass),density,temperature,/overplot,LEVELS = alog10(minres_mass),C_LABELS = fltarr(N_ELEMENTS(minres_mass)) + 1,C_ANNOTATION = STRTRIM(minres_mass/minparticles,2) + 'M'+SUNSYMBOL(),C_CHARSIZE = 1.5
oplot,alog10(g.dens*dens_convert*total_H),alog10(g.tempg),psym = 3
device,/close

;window,1
set_plot,'ps'
device,filename = 'jeanslength.ps',/color,bits_per_pixel= 8,/times
contour,alog10(jeanslength),density,temperature,xtitle="Log(nH)",ytitle = 'Log(T)',/FILL,nlevels = 60,title = 'Log(Jeans Length)',xrange=[minlogrho,maxlogrho],yrange=[minlogt,maxlogt],xstyle = 1,ystyle = 1
contour,alog10(jeanslength),density,temperature,/overplot,LEVELS = alog10(minres_force),C_LABELS = fltarr(N_ELEMENTS(minres_force)) + 1,C_ANNOTATION = STRTRIM(minres_force,2) + 'pc',C_CHARSIZE = 1.5
oplot,alog10(g.dens*dens_convert*total_H),alog10(g.tempg),psym = 3

x_data = alog10(g[where(g.tempg gt 100)].dens*dens_convert*total_H)
y_data = alog10(g[where(g.tempg gt 100)].tempg)

result = poly_fit(x_data,y_data,10)
x_fit = findgen(100)/100.0 * 9 - 6
y_fit = result[0] + x_fit*result[1] + x_fit*x_fit*result[2] + x_fit*x_fit*x_fit*result[3] + x_fit*x_fit*x_fit*x_fit*result[4] + x_fit*x_fit*x_fit*x_fit*x_fit*result[5] + x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[6]+ x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[7] + x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[8]+ x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[9]+ x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*x_fit*result[10]
oplot,x_fit,y_fit
device,/close

jeansmass_curve = jeansmass_calc(10.0^x_fit,10.0^y_fit)
set_plot,'x'
window,2
plot,x_fit,alog10(jeansmass_curve/minparticles),xtitle = "Log(nH)",ytitle = "Log(Jeans Mass) [Solar Mass]",yrange = [0, 7],xrange = [-2,3],xticklen = 1
oplot,[-6,3],[6,6],linestyle = 1
oplot,[-6,3],[5,5],linestyle = 1
oplot,[-6,3],[4,4],linestyle = 1
oplot,[-6,3],[3,3],linestyle = 1
oplot,[-6,3],[2,2],linestyle = 1
oplot,[-6,3],[1,1],linestyle = 1

critical_rho = [-1.4,-0.65,0.05,1.1,2.35,2.8]
critical_T = result[0] + critical_rho*result[1] + critical_rho*critical_rho*result[2] + critical_rho*critical_rho*critical_rho*result[3] + critical_rho*critical_rho*critical_rho*critical_rho*result[4] + critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*result[5] + critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*result[6]+ critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*result[7] + critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*result[8]+ critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*result[9]+ critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*critical_rho*result[10]
print,critical_T
jeanslength = jeanslength_calc(10.0^critical_rho,10.0^critical_T)
print,jeanslength

mvir = 1e6 
rvir =(3.*mvir/(4.*!PI*200*278*0.7^2))^(1./3.)*1000 ;turn from kpc into pc
print,rvir
print,jeanslength/rvir
print,rvir/jeanslength
stop
END
