pro plotkenn,sfrgas

;sims in which the number of particles per mass remains the same
col=[0,254,64]
m=2.362e5
t=1e9
sfr=sfrgas.sfr ;starformation rate
gas=sfrgas.gas ;gas 

loadct,39;load color table

fbind=[0,1,2,4]
plot,gas,sfr,psym=1,xrange=[0.8,max(gas)*2],yrange=[1e-6,max(sfr)*2],xtitle=textoidl(' \Sigma')+"!lgas!n [M!lsolar!n pc!u-2!n]",ytitle=textoidl(' \Sigma')+"!lSFR!n [M!lsolar!n kpc!u-2!n yr!u-1!n]",xstyle=1,/xlog,/ylog,xmargin=[6.4,1],ymargin=[3.7,0.6],/nodata
oplot,gas,sfr,psym=1,color=0
;oplot,gas[1,fbind],sfr[1,fbind],psym=1,color=254
oplot,gas[*,1],sfr[*,1],psym=1,color=254
oplot,gas[*,2],sfr[*,2],psym=1,color=64
xline=findgen(2)+min(gas[where(gas GT 0)])
xline[1]=max(gas)*2
yline=2.5e-4*xline^1.4
print,xline
print,yline
oplot,xline,yline
 legend,['1E6 Mass per Particle','1E7 Mass per Particle','1E8 Mass per Particle'],psym=[1,1,1],color=col,/bottom,/right

END
