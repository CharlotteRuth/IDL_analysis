;What are the characteristics of the accreting gas


;filename = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00072.dir/h603.cosmo50cmb.3072g14HBWK.00072.halo.2'
;pfile =    '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param'
;filename2 = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1'
;pfile2    = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param'

PRO accreting,filename,pfile,filename2,pfile2
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23

units = tipsyunits(pfile)
rtipsy,filename,h,g,d,s
readarr,filename + '.eDot',h,edot,part = 'gas',/ascii
;dens_convert = units.massunit * gm_per_msol * amu_per_gm/units.lengthunit^3/cm_per_kpc^3/h.time^3
;g.dens = g.dens*dens_convert
fixunits,h,g,d,s,units
edot_amu = edot/g.dens
angmom, g, jvecg, lvec, jmagg, ltot
;readcol,filename + '.eDot',eall
;eall1 = eall
;eall1[1:h.ngas] = edot_amu
;writecol,filename + '.eDotamu',eall1

units2 = tipsyunits(pfile2)
rtipsy,filename2,h2,g2,d2,s2
readarr,filename2 + '.eDot',h2,edot2,part = 'gas',/ascii
;dens_convert = units2.massunit * gm_per_msol * amu_per_gm/units2.lengthunit^3/cm_per_kpc^3/h.time^3
;g2.dens = g2.dens*dens_convert
fixunits,h2,g2,d2,s2,units2
edot_amu2 = edot2/g2.dens
angmom, g2, jvecg2, lvec2, jmagg2, ltot2

loadct,39

window,0
diskg = where(g.tempg le 2e4,complement = halog)
;plot,-1.0*edot_amu,jvecg[*,2],psym = 3,/xlog,/ylog,yrange = [1e-9,1e-2],xrange = [1e-2,1e5],xtitle = 'Cooling Rate',ytitle = 'Angular Momentum'
;oplot,-1.0*edot_amu[diskg],jvecg[diskg,2],psym = 3,color = 50
;oplot,-1.0*edot_amu[halog],jvecg[halog,2],psym = 3,color = 50

plot,-1.0*edot_amu[halog],jvecg[halog,2],psym = 3,/xlog,/ylog,yrange = [1e-9,1e-2],xrange = [1e-2,1e5]

plot,-1.0*edot,jvecg[*,2],psym = 3,/xlog,/ylog,yrange = [1e-9,1e-2]
oplot,-1.0*edot[halog],jvecg[halog,2],psym = 3,color = 50

plot,g.dens,jvecg[*,2],psym = 3,/xlog,/ylog,yrange = [1e-9,1e-2],xrange = [1e-8,1e3],xstyle = 1
oplot,g.dens,jvecg[halog,2],psym = 3,color = 50

plot,-1.0*edot_amu,jvecg[*,2],psym = 3,/xlog,/ylog,yrange = [1e-2,1e4],xrange = [1e-4,1e3],xstyle = 1,xtitle = 'Cooling Rate',ytitle = 'Angular Momentum'
oplot,-1.0*edot_amu[diskg],jvecg[diskg,2],psym = 3,color = 80
oplot,-1.0*edot_amu[halog],jvecg[halog,2],psym = 3
stop

window,1
diskg2 = where(g2.tempg le 2e4,complement = halog2)
plot,g2.dens,jvecg2[*,2],psym = 3,/xlog,/ylog,yrange = [1e-2,1e4],xrange = [1e-8,1e3],xstyle = 1
oplot,g2.dens,jvecg2[halog2,2],psym = 3,color = 50

plot,-1.0*edot_amu2,jvecg2[*,2],psym = 3,/xlog,/ylog,yrange = [1e-2,1e4],xrange = [1e-4,1e3],xstyle = 1,xtitle = 'Cooling Rate',ytitle = 'Angular Momentum'
oplot,-1.0*edot_amu2[diskg2],jvecg2[diskg2,2],psym = 3,color = 80
oplot,-1.0*edot_amu2[halog2],jvecg2[halog2,2],psym = 3
stop
END
