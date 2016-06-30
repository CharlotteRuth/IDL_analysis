PRO  lawPaper
rtipsy,'h277.cosmo50cmb.3072g14HMbwK.00120.00017.halo.1',h,g,d,s
readarr,'h277.cosmo50cmb.3072g14HMbwK.00120.00017.halo.1.HI',h,HI,part = 'gas',/ascii
readarr,'h277.cosmo50cmb.3072g14HMbwK.00120.00017.halo.1.H2',h,H2,part = 'gas',/ascii
;readarr,'h277.cosmo50cmb.3072g14HMbwK.00120.00017.massform',h,massform,part = 'star',/ascii
units = tipsyunits('../../h277.cosmo50cmb.3072g14HMbwK.param')

r = 0.00015
deltat = 100e6
tcurrent = MAX(s.tform)
indd = where(SQRT(d.x*d.x + d.y*d.y + d.z*d.z) LE r )
indg = where(SQRT(g.x*g.x + g.y*g.y + g.z*g.z) LE r )
indghot = where(SQRT(g.x*g.x + g.y*g.y + g.z*g.z) LE r AND g.dens*units.rhounit GE 0.1 AND g.tempg GE 2e4)
inds = where(SQRT(s.x*s.x + s.y*s.y + s.z*s.z) LE r )
indsnow = where(SQRT(s.x*s.x + s.y*s.y + s.z*s.z) LE r AND s.tform GE tcurrent - deltat/units.timeunit)

z = 1/h.time - 1
Hmass = total(g[indg].mass*HI[indg] + g[indg].mass*H2[indg]*2.0)*1.3*units.massunit
Smass = total(s[inds].mass)*units.massunit
Dmass = total(d[indd].mass)*units.massunit

SFR = N_ELEMENTS(indsnow)*units.istarmass/deltat
SFRden = SFR/!PI/(r*h.time*units.lengthunit)^2

Hden = (SFRden/2.5e-4)^(1/1.4)
HmassSK = Hden*!PI*(r*h.time*units.lengthunit)^2*1e6

arms = selectregion_tipsy(s[indsnow].x*h.time*units.lengthunit,s[indsnow].y*h.time*units.lengthunit,area = area)
SFRbar = N_ELEMENTS(arms)*units.istarmass/deltat/area

;histogramp,g.vz*units.vunit*h.time,nbins = 100,/overplot,color = 100,/normalize,min = -300,max = 300
sigz = histogram(g[indghot].vz*units.vunit*h.time,nbins = 100,min = -300,max = 300,locations = sigz_x)
fit = gaussfit(sigz_x,sigz,coeff)
plot,sigz_x,sigz,psym = 10,xtitle = 'V_z [km/s]'
oplot,sigz_x,fit,color = 100
print,coeff[2]
END
