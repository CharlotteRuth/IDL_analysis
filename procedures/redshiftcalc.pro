function redshiftcalc,time


end

;stupid, rough calc to get z as a function of t made by fitting to the
;lookback time as a function of z plot
pro redshiftcalc_orig,
loadct,39
ageUniverse = wmap3_lookback(10000)
z = findgen(1000)*10.0/1000.0
time = ageUniverse - wmap3_lookback(z)
t2 = time - ageUniverse
plot,time,z,xtitle = 'Time',ytitle = 'z',yrange = [0.01,10];,/ylog
r = poly_fit(t2,z,7)
zfit = r[0] + $
  r[1]*t2 + $
  r[2]*t2*t2 + $
  r[3]*t2*t2*t2 + $
  r[4]*t2*t2*t2*t2 + $
  r[5]*t2*t2*t2*t2*t2 + $
  r[6]*t2*t2*t2*t2*t2*t2 + $
  r[7]*t2*t2*t2*t2*t2*t2*t2 ; + $
;  r[8]*t2*t2*t2*t2*t2*t2*t2*t2  + $
;  r[9]*t2*t2*t2*t2*t2*t2*t2*t2*t2 + $
;  r[10]*t2*t2*t2*t2*t2*t2*t2*t2*t2*t2+ $
;  r[11]*t2*t2*t2*t2*t2*t2*t2*t2*t2*t2*t2
oplot,time,zfit,color = 240
end
