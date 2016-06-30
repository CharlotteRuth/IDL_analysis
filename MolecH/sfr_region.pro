PRO sfr_region,s,units
deltat = 1e6
tcurrent = MAX(s.tform*unit.timeunit)
snew = s[where(s.tform*units.timeunit ge tcurrent - 1e7)] 

sarmind =  selectregion_tipsy(snew.x*units.lengthunit,snew.y*units.lengthunit,area = area)
sarm = snew[sarmind]
sarm = sarm[where(sarm.tform*units.timeunit ge tcurrent - deltat)]
sfr = N_ELEMENTS(sarm)*units.istarmass/deltat/area
END
