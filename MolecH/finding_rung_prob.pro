pro finding_rung_prob
cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6HBWK'
units = tipsyunits('h516.cosmo25cmb.1536g6HBWK.param')
rtipsy,'h516.cosmo25cmb.1536g6HBWK',h,g,d,s
g.dens = g.dens*units.RHOUNIT
readarr,'h516.cosmo25cmb.1536g6HBWK.rung',h,rung,/ascii,part = 'gas'
readarr,'h516.cosmo25cmb.1536g6HBWK.SPHdt',h,sphdt,/ascii,part = 'gas'
sphdt = sphdt*units.timeunit
readarr,'h516.cosmo25cmb.1536g6HBWK.c',h,c,/ascii,part = 'gas'
c = c*units.vunit
readarr,'h516.cosmo25cmb.1536g6HBWK.dt',h,dt,/ascii,part = 'gas'
dt = dt*units.timeunit
readarr,'h516.cosmo25cmb.1536g6HBWK.HI',h,HI,/ascii,part = 'gas'
readarr,'h516.cosmo25cmb.1536g6HBWK.H2',h,H2,/ascii,part = 'gas'
readarr,'h516.cosmo25cmb.1536g6HBWK.eDot',h,eDot,/ascii,part = 'gas'
readarr,'h516.cosmo25cmb.1536g6HBWK.eCool',h,eCool,/ascii,part = 'gas'
readarr,'h516.cosmo25cmb.1536g6HBWK.eHeat',h,eHeat,/ascii,part = 'gas'

ind = where(rung eq MAX(rung))
g_rung = g[where(rung eq MAX(rung))]
sphdt_rung = sphdt[where(rung eq MAX(rung))]
c_rung = c[where(rung eq MAX(rung))]
dt_rung = dt[where(rung eq MAX(rung))]
H2_rung = H2[where(rung eq MAX(rung))]
HI_rung = HI[where(rung eq MAX(rung))]
eDot_rung = eDot[ind]
eHeat_rung = eHeat[ind]
eCool_rung = eCool[ind]

y = histogram(alog10(dt),locations = x)
plot,x,y,psym = 10,/ylog
stop
plot,g.dens,g.tempg,/xlog,/ylog,psym = 3
oplot,g_rung.dens,g_rung.tempg,psym = 2
stop
end
