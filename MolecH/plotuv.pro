pro plotuv
loadct,39

h = 4.1357e-15 ;units of eV s
c = 3e18 ;units of Angstroms s^-1
ev_per_erg = 6.24150974e11

openr,1,'haardt_madau_galaxy.dat'
openr,2,'haardt_madau_quasar.dat'
galaxy = fltarr(51,432)
quasar = fltarr(51,432)
readf,1,galaxy
readf,2,quasar
close,/all
zHM = [0, .04921, .1, .1547, .2114, .2709	, .3333	, .3988	, .4675	, 0.5396, .6152	, .6945	, .7778	, .8651	, .9567	, 1.053	, 1.154	, 1.259	, 1.370	, 1.4870, 1.609	, 1.737	, 1.871	, 2.013	, 2.16	, 2.316	, 2.479	, 2.649	, 2.829	, 3.017	, 3.214	, 3.421	, 3.638	, 3.866	, 4.105	, 4.356	, 4.619	, 4.895	, 5.184	, 5.4880, 5.807	, 6.141	, 6.492	, 6.859	, 7.246	, 7.650	, 8.075	, 8.521	, 8.989	, 9.4790]
lambda = REFORM(galaxy[0,*])
dlambda = lambda[0:430] - lambda[1:431]
nu =  3e18/lambda
dnu = nu[0:430] - nu[1:431]
lambda = lambda[0:430]
nu = nu[0:430]

galaxy = galaxy[1:50,0:430]
quasar = quasar[1:50,0:430]

lower_ev = 12.24
upper_ev = 13.51

long_lambda = h*c/lower_ev
short_lambda = h*c/upper_ev

lw = where(lambda le long_lambda AND lambda ge short_lambda)

galaxy_lw = (galaxy[*,lw[0]] + galaxy[*,lw[1]])/2.0
quasar_lw = (quasar[*,lw[0]] + quasar[*,lw[1]])/2.0
total_lw = galaxy_lw + quasar_lw

window,0
plot,zHM,galaxy_lw,linestyle = 2
oplot,zHM,quasar_lw,linestyle = 2
oplot,zHM,total_lw

HI_E0 = 13.6
HI_nu0 = HI_E0/h
epsilon = SQRT(nu/HI_nu0 - 1)
sigmaHI = nu
;Ostriker 89, 2.4
sigmaHI[where(nu ge HI_nu0, complement = invalid)] = 6.30e-18*(HI_nu0/nu[where(nu ge HI_nu0)])^(4.0) $
                                                             *exp(4.0 - 4.0*atan(epsilon)/epsilon)/(1 - exp(-2.0*!PI/epsilon)) 
sigmaHI[invalid] = 0
phot_HI = fltarr(N_ELEMENTS(zHM))
FOR i = 0, N_ELEMENTS(zHM) - 1 DO BEGIN
    phot_HI[i] = TOTAL(4.0*!PI*sigmaHI*galaxy[i,*]*ev_per_erg/h/nu*dnu)
ENDFOR

E0 = 13.61 
sigma0 = 9.492e-16
ya = 1.469
P  = 3.188
yw = 2.039
y0 = 0.4434
y1 = 2.136
HeI_E0 = 24.59
HeI_nu0 = HeI_E0/h
epsilon = SQRT(nu/HeI_nu0 - 1)
sigmaHeI = nu
;Ostriker 89, 2.4
x = h*nu[where(nu ge HeI_nu0)]/E0-y0
y = sqrt(x*x+y1*y1)
Fy = ((x-1.0)^2+yw^2)*y^(0.5*P-5.5)*(1.0+sqrt(y/ya))^(-P)
sigmaHeI[where(nu ge HeI_nu0, complement = invalid)] = sigma0*Fy
sigmaHeI[invalid] = 0
phot_HeI = fltarr(N_ELEMENTS(zHM))
FOR i = 0, N_ELEMENTS(zHM) - 1 DO BEGIN
    phot_HeI[i] = TOTAL(4.0*!PI*sigmaHeI*galaxy[i,*]*ev_per_erg/h/nu*dnu)
ENDFOR

HeII_E0 = 54.4
HeII_nu0 = HeII_E0/h
epsilon = SQRT(nu/HeII_nu0 - 1)
sigmaHeII = nu
;Ostriker 89, 2.4
sigmaHeII[where(nu ge HeII_nu0, complement = invalid)] = 6.30e-18/4.0*(HeII_nu0/nu[where(nu ge HeII_nu0)])^(4.0) $
                                                             *exp(4.0 - 4.0*atan(epsilon)/epsilon)/(1 - exp(-2.0*!PI/epsilon)) 
sigmaHeII[invalid] = 0
phot_HeII = fltarr(N_ELEMENTS(zHM))
FOR i = 0, N_ELEMENTS(zHM) - 1 DO BEGIN
    phot_HeII[i] = TOTAL(4.0*!PI*sigmaHeII*galaxy[i,*]*ev_per_erg/h/nu*dnu)
ENDFOR


window,1
plot,nu/1e16,sigmaHI/1e-18,xrange = [0,6],xtitle = 'nu [10^16 Hz]',ytitle = 'a_nu [10^-18 cm^2]' ;Fig 2.2, Ostriker 89
oplot,nu/1e16,sigmaHI/1e-18,color = 50
oplot,nu/1e16,sigmaHeI/1e-18,color = 120
oplot,nu/1e16,sigmaHeII/1e-18,color = 240
legend,['HI','HeI','HeII'],color = [50,120,240],linestyle = 0,/top,/right

;window,2
;plot,z,phot_HI,/ylog,xrange = [0,8.5],xstyle = 1,xtitle = 'Z',ytitle = 'Rate'

;----------------------------------------------------------------------------

openr,1,'uvtable.txt'
data = fltarr(7,90)
readf,1,data
close,1

window,3
set_plot,'x'
plot,data[0,*],data[1,*],yrange = [1e-17,1e-10],xtitle = 'Z',ytitle = 'Rate, Heat',/ylog,xrange = [0,8.5],xstyle = 1
oplot,data[0,*],data[1,*],color = 50
oplot,data[0,*],data[4,*],color = 50,linestyle = 1
oplot,zHM,phot_HI,color = 40,linestyle = 2
oplot,data[0,*],data[2,*],color = 120
oplot,data[0,*],data[5,*],color = 120,linestyle = 1
oplot,zHM,phot_HeI,color = 110,linestyle = 2
oplot,data[0,*],data[3,*],color = 240
oplot,data[0,*],data[6,*],color = 240,linestyle = 1
oplot,zHM,phot_HeII,color = 230,linestyle = 2
oplot,zHM,total_lw*1.1e8*4.0*!PI,color = 190
oplot,zHM,zHM*0 + 6.4e-13,color = 190,linestyle = 1
legend,['HI','HeI','HeII','H2'],color = [50,120,240,190],linestyle = 0,/bottom,/left
H2_diss_z = REVERSE(spline(zHM,total_lw*1.1e8*4.0*!PI,REVERSE(REFORM(data[0,*]))))
oplot,data[0,*],H2_diss_z,color = 200,linestyle = 2

openw,1,'H2_photoDiss_z.txt'
printf,1,TRANSPOSE([[reform(data[0,*])],[H2_diss_z],[H2_diss_z*0 + 6.4e-13]]) 
close,1

fit = POLY_FIT(alog10(zHM[1:48]),alog10(total_lw[1:48]*1.1d8*4.0*!PI),1)
xfit = 10^(fit[0] + $
              fit[1]*alog10(zHM) );+ $
        ;      fit[2]*zHM*zHM + $ 
        ;      fit[3]*zHM*zHM*zHM + $
        ;      fit[4]*zHM*zHM*zHM*zHM + $
        ;      fit[5]*zHM*zHM*zHM*zHM*zHM + $
        ;      fit[6]*zHM*zHM*zHM*zHM*zHM*zHM + $
        ;      fit[7]*zHM*zHM*zHM*zHM*zHM*zHM*zHM
oplot,zHM,xfit,color = 240
stop

set_plot,'ps'
device,filename = 'UVHeatRate.eps',/color,bits_per_pixel= 8,/times
plot,data[0,*],data[1,*],yrange = [1e-17,1e-10],xtitle = 'Z',ytitle = 'Rate, Heat',/ylog
oplot,data[0,*],data[1,*],color = 50
oplot,data[0,*],data[2,*],color = 50,linestyle = 2
oplot,data[0,*],data[3,*],color = 120
oplot,data[0,*],data[4,*],color = 120,linestyle = 2
oplot,data[0,*],data[5,*],color = 240
oplot,data[0,*],data[6,*],color = 240,linestyle = 2
legend,['HI','HeI','HeII'],color = [50,120,240],linestyle = 0,/bottom,/left
device,/close

set_plot,'x'

if (0) then begin
plot,data[0,*],data[1,*]/MAX(data[1,*]),xtitle = 'Z',ytitle = 'Rate, Heat',/ylog
oplot,data[0,*],data[1,*]/MAX(data[1,*]),color = 50
oplot,data[0,*],data[2,*]/MAX(data[2,*]),color = 50,linestyle = 2
oplot,data[0,*],data[3,*]/MAX(data[3,*]),color = 120
oplot,data[0,*],data[4,*]/MAX(data[4,*]),color = 120,linestyle = 2
oplot,data[0,*],data[5,*]/MAX(data[5,*]),color = 240
oplot,data[0,*],data[6,*]/MAX(data[6,*]),color = 240,linestyle = 2
legend,['HI','HeI','HeII'],color = [50,120,240],linestyle = 1
endif

end
