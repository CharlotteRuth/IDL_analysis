PRO coolingcurve
;openr,density.txt  energy.txt  metalcooling.txt  pdv.txt
;temperature.txt
set_plot,'x'
loadct,39
base = '/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_soft/testingcooling1/'
base = '/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_soft/coolingDataFiles/'
base = '/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_soft_nofit/coolingDataFiles/'
base = '/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_soft/coolingDataFiles/'
base = '/astro/net/scratch1/christensen/MolecH/OneParticle/coolingDataFiles/'
nel = 4000
nel = 1680
nel = 5000
nel = 18520
nel = 32129
nel = 500
ymax = 0.2
ymin = -0.2

density = fltarr(nel)
temperature = fltarr(nel)
bremHII = fltarr(nel)
bremHeII = fltarr(nel)
bremHeIII = fltarr(nel)
radHII = fltarr(nel)
radHeII = fltarr(nel)
radHeIII = fltarr(nel)
photheatHI = fltarr(nel)
photheatHIsheild = fltarr(nel)
photheatHeI = fltarr(nel)
photheatHeII = fltarr(nel)
photheatH2 = fltarr(nel)
photheatH2sheild = fltarr(nel)
mcool = fltarr(nel)
mheat = fltarr(nel)
internal = fltarr(nel)
pdv = fltarr(nel)
edot = fltarr(nel)
energy = fltarr(nel)
temperature2 = fltarr(nel)
h2 = fltarr(nel)
hI = fltarr(nel)
hII = fltarr(nel)

colors = findgen(nel)/nel*240

close,/all
openr,1,base+'density.txt'
readf,1,density
close,1
;density = (density)/(MAX(density) - MIN(density))*(ymax - ymin)/2 + ymin/2

openr,1,base+'temperature.txt'
readf,1,temperature
close,1
;temperature = (temperature)/(MAX(temperature) - MIN(temperature))*(ymax - ymin)/2 + ymin/2

openr,1,base+'bremHII.txt'
readf,1,bremHII
close,1

openr,1,base+'bremHeII.txt'
readf,1,bremHeII
close,1

openr,1,base+'bremHeIII.txt'
readf,1,bremHeIII
close,1

openr,1,base+'radHII.txt'
readf,1,radHII
close,1

openr,1,base+'radHeII.txt'
readf,1,radHeII
close,1

openr,1,base+'radHeIII.txt'
readf,1,radHeIII
close,1

openr,1,base+'photheatHI.txt'
readf,1,photheatHI
close,1

openr,1,base+'photheatHIsheild.txt'
readf,1,photheatHIsheild
close,1

openr,1,base+'photheatHeI.txt'
readf,1,photheatHeI
close,1

openr,1,base+'photheatHeII.txt'
readf,1,photheatHeII
close,1

openr,1,base+'photheatH2.txt'
readf,1,photheatH2
close,1

openr,1,base+'photheatH2sheild.txt'
readf,1,photheatH2sheild
close,1

openr,1,base+'metalcooling.txt'
readf,1,mcool
close,1

openr,1,base+'metalheating.txt'
readf,1,mheat
close,1

openr,1,base+'pdv.txt'
readf,1,pdv
close,1

openr,1,base+'internal.txt'
readf,1,internal
close,1

openr,1,base+'edot.txt'
readf,1,edot
close,1

openr,1,base+'energy.txt'
readf,1,energy
close,1

openr,1,base+'temperature2.txt'
readf,1,temperature2
close,1

openr,1,base+'H2.txt'
readf,1,h2
close,1

openr,1,base+'HI.txt'
readf,1,hI
close,1

openr,1,base+'HII.txt'
readf,1,hII
close,1

start = 5e3 ;0
ending = 8e3 ;nel - 1
start = 31407
ending = nel - 1
start = 0


;set_plot,'ps'
window,0
plot,edot,xstyle = 1,xrange = [start,ending],yrange = [-1,1];,xrange = [2000,2100];,xrange = [0,nel];xrange = [2045,2055],yrange = [-13.5,5.0],yrange = [-0.03,0.03]
;oplot,[2049.5,2049.5],[-13.5,50],linestyle = 1
oplot,internal,linestyle = 1,color = 100
;oplot,pdv,linestyle = 1, color = 150
oplot,-1.0*mcool,linestyle = 2, color = 70
oplot,mheat,linestyle = 2, color = 220
oplot,photheatHI,linestyle = 3, color = 200
oplot,photheatHIsheild,linestyle = 3, color = 190
legend,['Edot','internal','mcool','mheat','photheatHI','photheatHIsheild'],color = [0,100,70,220,200,190],linestyle = [0,1,2,2,3,3]

window,1
plot,alog10(temperature), linestyle = 2, xrange = [start,ending], xstyle = 1,yrange = [-4,4.5];, xrange = [0,nel], xrange = [2045,2055]
oplot,alog10(temperature), linestyle = 2, color = 100
oplot,[0,nel],[2,2], color = 100
oplot,alog10(density), linestyle = 2, color = 150
legend,['temperature','density'], color = [100,150], linestyle = [2,2]


plot,temperature, linestyle = 2, xrange = [start,ending], xstyle = 1
oplot,temperature, linestyle = 2, color = 100
stop
ind = where(temperature gt 100)
;window,2
;plot,energy[ind],edot[ind],psym = 2,/xlog,xrange = [7e9,7.5e9],yrange = [-0.0072 ,-0.0068 ];,yrange = [-2.5,0.5]
;oplot,energy[ind],edot[ind],psym = 2,color = 70
;oplot,[7.00541e+09,7.00541e+09],[-0.00703148,-0.00703148],color = 240,psym = 4
;print,where(energy gt 7e9 AND energy lt 7.5e9 AND edot gt -0.0072 AND edot -0.0068)

ind = where(temperature lt 100)

window,2
plot,h2*2,xrange = [start,ending], xstyle = 1,yrange = [0, 0.7],linestyle = 1
oplot,h2*2,color = 50,linestyle = 1
oplot,ind,h2[ind],color = 50,psym = 2
oplot,hI,color = 120,linestyle = 1
oplot,ind,hI[ind],color = 120,psym = 2
oplot,hII,color = 240,linestyle = 1
oplot,ind,hII[ind],color = 240,psym = 2
legend,['H2','HI','HII'],color = [50,120,240],linestyle = [0,0,0]

;device,filename = 'coolingcurve.ps',/color
window,3
;start = 5e3 ;0
;ending = 8e3 ;nel - 1
colors = findgen(ending-start)/(ending-start)*240
t_short = temperature[start:ending]
d_short = density[start:ending]
plot,d_short,t_short,xtitle = 'Density',ytitle = 'Temperature',/xlog,/ylog,psym = 3,xrange = [3,100],yrange=[10,400];xrange = [0.1,300],yrange = [10,1e4]
for i= start, ending - 1 do oplot,[density[i],density[i+1]],[temperature[i],temperature[i+1]],color = colors[i-start],psym = 3
ind = where(t_short lt 90)
for i= 0, N_ELEMENTS(ind) -2  do oplot,[d_short[ind[i]],d_short[ind[i+1]]],[t_short[ind[i]],t_short[ind[i+1]]],color = colors[ind[i]],psym = 2
;;stop
;device,/close

;yaxes = fltarr(nel)
;yaxes[ind] = 1
;plot,yaxes,yrange = [0,1.2],xrange = [5e3,8e3]

;;window,1
;plot,energy,temperature,xtitle = 'Energy',ytitle = 'Temperature',/xlog,/ylog,psym = 3,yrange = [10,1e4]
;for i= 0, nel-2 do oplot,[energy[i],energy[i+1]],[temperature[i],temperature[i+1]],color = colors[i],psym = 3
;for i= 0, N_ELEMENTS(ind) -2  do oplot,[energy[ind[i]],energy[ind[i+1]]],[temperature[ind[i]],temperature[ind[i+1]]],color = colors[ind[i]],psym = 2
;;stop
;device,filename = 'pdv.ps',/color

;;window,1
;plot,pdv,temperature,xtitle = 'PdV',ytitle = 'Temperature',/xlog,/ylog,psym = 3,yrange = [10,1e4]
;for i= 0, nel-2 do oplot,[pdv[i],pdv[i+1]],[temperature[i],temperature[i+1]],color = colors[i],psym = 3
;for i= 0, N_ELEMENTS(ind) -2  do oplot,[pdv[ind[i]],pdv[ind[i+1]]],[temperature[ind[i]],temperature[ind[i+1]]],color = colors[ind[i]],psym = 2
;;stop
;device,/close

;device,filename = 'metalcool.ps',/color
;;window,2
;plot,mcool,temperature,xtitle = 'Metal Cooling',ytitle = 'Temperature',/xlog,/ylog,psym = 3,yrange = [10,1e4]
;for i= 0, nel-2 do oplot,[mcool[i],mcool[i+1]],[temperature[i],temperature[i+1]],color = colors[i],psym = 3
;for i= 0, N_ELEMENTS(ind) -2  do oplot,[mcool[ind[i]],mcool[ind[i+1]]],[temperature[ind[i]],temperature[ind[i+1]]],color = colors[ind[i]],psym = 2
;device,/close

stop
END


;-12.5920 + 5.25788*Tlog -0.949444*nHlog + 1.02849*Tlog*nHlog - 0.647718*Tlog^2 - 0.137914*Tlog^2*nHlog+ 0.000204617*Tlog^2*nHlog^2 - 0.0751646*Tlog*nHlog^2 + 0.247727*nHlog^2;

;-13.3866 + 3.62845*10^(-10)*Tlog + nHlog -3.12419*10^(-11)*Tlog*nHlog - 6.09512*10^(-11)*Tlog^2 + 5.49561*10^(-12)*Tlog^2*nHlog + 7.33096*10^(-12)*Tlog^2*nHlog^2 - 3.80744*10^(-11)*Tlog*nHlog^2 + 2.65849*10^(-8)*nHlog^2
