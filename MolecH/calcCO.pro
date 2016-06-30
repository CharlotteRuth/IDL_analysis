FUNCTION calcCO, rho, zmetal_abs, length, xuv
;taken from Bolatto, Jackson & Ingalls
no0 = 1e4
xuv0 = 250.0
alpha0 = 1.25
delta0 = 435.0
eta0 = 3e-4
zmetal_MW = 0.014791 ;log(z) = 1.42 + log(O/H) from Vila-Costas 92, and 8.75 = 12 + log(O/H) Bolatto 99 referenceing something about Orion

lowz = where(ZMetal_abs le 0.1, complement = highz)
yHe = Zmetal_abs
IF (lowz[0] ne -1) THEN yHe[lowz] = (0.236 + 2.1*ZMetal_abs[lowz])/4.0  
IF (highz[0] ne -1) THEN yHe[highz] = (-0.446*(ZMetal_abs[highz] - 0.1)/0.9 + 0.446)/4.0
yH = 1.0  - yHe*4.0 - ZMetal_abs
nH = yH*rho
alpha = alpha0*(xuv/xuv0)^(0.05)*(rho/no0)^(0.02) ;equation 24 from Bolatto
delta = delta0*(xuv/xuv0)^(-0.13)*(rho/no0)^(1.18) ; equataion 25 from Bolatto

h0 = length/delta
no = 2.0/3.0*rho*eta0*h0
zmetal = zmetal_abs/zmetal_MW
dz = delta*zmetal
a = alpha + 1.0

nco = zmetal
base = zmetal
toosmallz = where(zmetal le (alpha + 1)/delta, complement = correct)
IF (toosmallz[0] ne -1) THEN BEGIN
    nco[toosmallz] = 0
    base[toosmallz] = 0
ENDIF
IF (correct[0] ne -1) THEN BEGIN
    nco[correct] = no[correct]*((a[correct]*a[correct]*a[correct] - 6.0*a[correct]*a[correct]*dz[correct] + 3.0*a[correct]*dz[correct]*dz[correct] + 2.0*dz[correct]*dz[correct]*dz[correct])/(2.0*dz[correct]*dz[correct]) +  3.0*a[correct]*alog10(a[correct]/dz[correct]))
    base[correct] = ((a[correct]*a[correct]*a[correct] - 6.0*a[correct]*a[correct]*dz[correct] + 3.0*a[correct]*dz[correct]*dz[correct] + 2.0*dz[correct]*dz[correct]*dz[correct])/(2.0*dz[correct]*dz[correct]) +  3.0*a[correct]*alog10(a[correct]/dz[correct]))
ENDIF
;IF (correct[0] ne -1) THEN nco1 = no*((alpha[correct]*alpha[correct]*alpha + 3.0*alpha*alpha*(1.0 - 2.0*delta)*zmetal[correct] + 3.0*alpha*(1.0 - 4.0*delta + delta*delta)*zmetal[correct]*zmetal[correct] + (1.0 - 6.0*delta + 3.0*delta*delta + 2.0*delta*delta*delta)*zmetal[correct]*zmetal[correct]*zmetal[correct])/ (2.0*delta*delta*zmetal[correct]*zmetal[correct]) + 3.0*(alpha + zmetal[correct])*alog10((alpha + zmetal[correct])/(delta*zmetal[correct])))
;second formula
I_CO = nco*1.73d-26
rhonco = nco/h0
x = [[rhonco],[I_CO]]
;set_plot,'ps
;device,filename = 'Bolatto99_fig2.eps',/color,bits_per_pixel=8,/times
;loadct,39
;plot,zmetal,nco,/ylog,xrange = [0,1.8],yrange = [2e17,2e20],xstyle = 1,ystyle = 1,xtitle = 'Relative Metallicity (Z)', ytitle = 'Expected Column (cm^-2)'
;oplot,zmetal,nco1,linestyle = 2,color = 240
;plot,rho,(length*rho*yH)/(nco*1.73d-26),psym = 1,xtitle = 'Density',ytitle = 'X Factor, Assuming H2/H = 1',/xlog,/ylog,xrange = [10,1000],yrange = [1e28,1e32]
;plot,rho,nco/h0,/ylog,/xlog,xtitle = 'Density',ytitle = 'CO Column Density',psym = 3
;device,/close
;set_plot,'x'

return,x
END

pro setup_calcCO_master
base ='/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e6g2/MW_disk'
step = '00004'
msol_per_sysmass = 1.36e17
kpc_per_syslength = 1e5
calcCO_master,base = base,step = step,msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength
!P.thick = 1.5
!P.CHARSIZE = 1.5
!X.Charsize = 1.35
!Y.Charsize = 1.35

end

PRO calcCO_master,base = base, step = step, msol_per_sysmass = msol_per_sysmass, kpc_per_syslength = kpc_per_syslength

set_plot,'x'
loadct,39

cm_per_kpc = 3.0857d21
gm_per_msol = 1.989d33
amu_per_gm = 6.022142d23
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)

c = 1.557e11 
zmetal_MW = 0.014791 ;log(z) = 1.42 + log(O/H) from Vila-Costas 92, and 8.75 = 12 + log(O/H) Bolatto 99 referenceing something about Orion
lengthunit =  kpc_per_syslength*cm_per_kpc ;system length unit in cm (=1kpc)
massunit   = msol_per_sysmass*gm_per_msol ;system mass unit in gm
dens_convert =  msol_per_sysmass * gm_per_msol/kpc_per_syslength^3/cm_per_kpc^3 

z_test = fltarr(100) + zmetal_MW ;findgen(100)/100 *1.8
n_test = 10^(findgen(100)/100*2.0 + 1.0)  ;fltarr(100)+1d4
h0_test = fltarr(100)+1d17
xuv = fltarr(100) + 250.0 ;10^(findgen(100)/100.0*3.0 + 1.0) ;fltarr(100) + 250.0
;co_test = calcCO(n_test,1.7*z_test*zmetal_MW,h0_test*435.0,xuv)
;plot,z_test,co_test[*,1]/(1.73d-26),/ylog,xrange = [0,1.8],yrange = [1e17,2e20],xstyle = 1,ystyle = 1,xtitle = 'Relative Metallicity (Z)', ytitle = 'Expected Column (cm^-2)'
;oplot,[(alpha + 1)/delta,(alpha + 1)/delta*zmetal_MW],[1e17,2e20],linestyle = 1
;stop

filename = base+'.'+step
rtipsy,filename,h,g,d,s
t = g.tempg
rho = g.dens*dens_convert
zmetal = g.zmetal
;zmetal = zmetal*0 + 0.0249861

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
nh2_Col = h2*length*rho*amu_per_gm
xuv = fltarr(N_ELEMENTS(HI)) + 250.0
co = calcCO(rho*amu_per_gm,g.zmetal,length,xuv)

nh2 = HI_prime
co1 = HI_prime
co2 = HI_prime
nh2[1:h.ngas] = nH2_Col
co1[1:h.ngas] = reform(co[*,0])
co2[1:h.ngas] = reform(co[*,1])
openw,1,base+"."+step+".NH2"
openw,2,base+"."+step+".CO"
openw,3,base+"."+step+".I_CO"
printf,1,LONG(HI_prime[0]),FORMAT = '(I)'
printf,2,LONG(HI_prime[0]),FORMAT = '(I)'
printf,3,LONG(HI_prime[0]),FORMAT = '(I)'
for i=1LL,LONG64(N_ELEMENTS(HI_prime)) - 1 do  BEGIN
    printf,1,nh2[i]
    printf,2,co1[i]
    printf,3,co2[i]
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
plot,rho*amu_per_gm,(nH2_Col)/co[*,1],psym = 1,/ylog,/xlog,xrange=[10,1000],xtitle = 'Density [cm^3]',ytitle = 'N_{H2}/I_{CO}'

stop
END
