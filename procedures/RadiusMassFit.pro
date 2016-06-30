;RadiusMassFit.pro
;Charlotte Christensen
;10.7.05

;Ths IDL procedure will fit a curve to a tipsy file filename
;and then output the results of the screen.  It uses curvefit(r,
;m/r^2, weights, b, function_name="nfw", /noderivative), where b is an
;array of the initial guesses for c, m200 (in 1e12 Solar Masses), and
;r200 (in kpc), r an array of
;radiuses, m is an array of masses and nfw.pro is a curve fitting
;function written by Marcus.

PRO RadiusMassFit,datafile
b=[6.2,1e6,167] 

readcol,datafile,r,num,den,m,vc,vr,sr,vth,sth,j,jth,jph
weights = 1.0/(m/r^2)
weights = fltarr(n_elements(m))+1.0
fit = curvefit(r, m/r^2, weights, b, function_name="nfw", /noderivative)

;set_plot,'ps'
;device,filename=datafile+'_plot.ps'
plot,r,(m/r^2),psym=1,/xlog,/ylog,xtitle='Radius (kpc)',ytitle='Mass(R) (1e12 Solar Masses)',title='Galactic Mass as a Function of Radius'
oplot,r,fit
;device,/close
openw,1,datafile+'_fit.txt'
printf,1,b
printf,1,fit
close,1

set_plot,'x'
plot,r,vc,xtit="radius (kpc)",ytit="velocity (km/s)",tit="Rotation Curve"
plot,r,(m/r^2),psym=1,/xlog,/ylog,xtitle='Radius (kpc)',ytitle='Mass(R) (1e12 Solar Masses)',title='Galactic Mass as a Function of Radius '+datafile
print,b
print,fit
oplot,r,fit

END
