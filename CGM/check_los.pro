;check_los,'../h603.cosmo50cmb.3072g14HBWK.00512.halo.1','h603.g14HBWK.1_11.5_9.90_0x0.00kpc_y0.00kpc','../h603.cosmo50cmb.3072g14HBWK.param'

PRO check_los,tfilename,cgmfilename,pfile
bind  = binzread('binzfile.'   + cgmfilename)
partd = partzread('partzfile.' + cgmfilename)

units = tipsyunits(pfile)
rtipsy,tfilename,h,g,d,s
vir_rad = max(sqrt(g.x*g.x + g.y*g.y + g.z*g.z))
indbin = where(bind.bcoord GE -1*vir_rad AND bind.bcoord LE vir_rad)
bind = bind[indbin]
readarr,tfilename + '.hsm',h,smooth,part = 'gas',/ascii
dist = sqrt(g.x*g.x + g.y*g.y)
ind = where(dist LE 2*smooth)

empty = where(bind.temp EQ 0)


histogramp,smooth[ind]*units.lengthunit,xtitle = 'Smoothing Length [kpc]',nbins = 100
stop

histogramp,g[ind].z*units.lengthunit,xtitle = 'z coordinate [kpc]',nbins = 100
stop

msunpc2_gcm2 = 0.000208908219
amucm2_msuncp2 = 1.25e20
amucm3_msuncp3 = 40.77
msuncp3_amucm3 = 1/40.77
amu_gram = 6.022e23
cm_pc = 3.08567758e18
plot,bind.bcoord*units.lengthunit,bind.temp,/ylog,xtitle = 'Radius [kpc]',ytitle = 'Temperature [K]',yrange = [1e3,1e7]
stop
plot,bind.bcoord*units.lengthunit,bind.mass,/ylog,xtitle = 'Radius [kpc]',ytitle = 'Mass'
stop
plot,bind.bcoord*units.lengthunit,bind.mass/(bind.bsize^3*units.lengthunit^3*1000.0^3),/ylog,xtitle = 'Radius [kpc]',ytitle = 'Rho [Msol/pc^3]'
stop
plot,bind.bcoord*units.lengthunit,bind.rho*(bind.bsize*units.lengthunit*1000.0*cm_pc)/msunpc2_gcm2,/ylog,xtitle = 'Radius [kpc]',ytitle = 'Surface Density [Msun/pc^2]'
;gm/cm^3*binsize in sysunits*kpc_per_sysunit*pc_per_kpc*cm_per_pc*amu_per_gm
oplot,bind.bcoord*units.lengthunit,bind.mass/(bind.bsize*units.lengthunit*1000.0)^2,color = 100
stop
END
