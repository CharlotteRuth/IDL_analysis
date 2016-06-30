PRO derivCheck
set_plot,'x'
loadct,39
dKpcUnit = 25000
expand = 0.77756716
yH = 0.759582
cm_in_kpc = 3.08568025e21

spawn,'ls *.param',pfilelist
pfile = pfilelist[0]
units = tipsyunits(pfile)
readcol,'h516.1536g15HBWK.out',time,temp,rho,ZMetal,correL,dLymanWerner,H2,yH2,dydtH2, $
          Phot_H2, Coll_e_H2, Coll_H2_H2,Coll_HI_H2, Coll_HII_H2,shield,dustform,/silent
file = 'stepstemp/h516.cosmo25cmb.1536g15HBWK.00389.dir/h516.cosmo25cmb.1536g15HBWK.00389.halo.1'
rtipsy,file,h,g,d,s
g.dens = g.dens*units.rhounit/h.time/h.time/h.time
readarr,file+'.HI',h,HI,/ascii,part = 'gas'
readarr,file+'.H2',h,H2,/ascii,part = 'gas'
readarr,file+'.lw',h,lw,/ascii,part = 'gas'
readarr,file+'.smoothlength',h,hl,/ascii,part = 'gas'
hl = hl*units.lengthunit*cm_in_kpc
IF (FILE_TEST(file+".correL")) THEN BEGIN
    readarr,file+'.correL',h,correLfile,/ascii,part = 'gas'
    correLfile = correLfile*units.lengthunit*cm_in_kpc
    totalH = HI + 2.0*H2
    length = correLfile
;    print,'here'
ENDIF
length = hl

yHI = yH - 2.0*yH2
h = dKpcUnit * 3.08568025e21 * expand*correL
;rho = rho*units.rhounit
self_shield = s_H2_file(yH2,h)
dust_shield = s_d_file(yHI,yH2,h,zmetal)
shield = dust_shield*self_shield

steps = findgen(N_ELEMENTS(time))
int = findgen(N_ELEMENTS(time))*0
dydtH21 = findgen(N_ELEMENTS(time))*0
steps[0] = 0
int[0] = yH2[0]
;dydtH21[N_ELEMENTS(time)] = 0
FOR i = 1,N_ELEMENTS(time) - 1 DO steps[i] = steps[i - 1] + time[i]
FOR i = 1,N_ELEMENTS(time) - 1 DO int[i] = time[i]*dydtH2[i] + int[i - 1]
FOR i = 0,N_ELEMENTS(time) - 2 DO dydtH21[i] = (yH2[i + 1] - yH2[i])/time[i]
plot,dydtH2[1:N_ELEMENTS(time) - 1]-dydtH21[1:N_ELEMENTS(time) - 1],psym = 2
stop
plot,steps,shield*Phot_H2
stop
plot,g.dens*length,H2/totalH*2,psym = 3,yrange = [1e-6,1],/xlog,/ylog,xrange = [1e18,1e24],xtitle = 'Surface Density [amu/cm^2]',ytitle = 'H_2/H',title = 'Lyman Werner'
oplot,rho*correL,2.0*yH2/yH,color = 240,psym = -4
stop
END

function S_H2_file,yH2,h
omega_H2 = 0.2
x = yH2*h/5d14
return, (1 - omega_H2)/(1 + x)^2 + omega_H2/SQRT(1 + x)*exp(-0.00085*SQRT(1 + x))
end

function S_d_file,yHI,yH2,h,z
ZSOLAR = 0.0130215 ;0.0177 ;0.0130215
sigma_d = 2d-21
return, exp(-1.0*sigma_d*z/ZSOLAR*(yHI*h + 2.0*yH2*h))
;return, exp(-1.0*sigma_d*(yHI*h + 2.0*yH2*h))
end
