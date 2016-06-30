;Calculates the gas mass of the set of simulations

FUNCTION tully_fisher_obs_gasmass, files, msol_per_sysmass, smass = smass, gmassall = gmassall
gmass = fltarr(N_ELEMENTS(files))
smass = fltarr(N_ELEMENTS(files))
gmassall = fltarr(N_ELEMENTS(files))
FOR i = 0, N_ELEMENTS(files)  - 1 DO BEGIN
    rtipsy,files[i],h,g,d,s
    a = h.time
    readarr,files[i]+".HI",h,HI_frac,part = 'gas',/ascii
    HI = HI_frac*g.mass*msol_per_sysmass[i]
    He = (0.236 + 2.1*g.zmetal)/4.0*g.mass*msol_per_sysmass[i]
    gmass[i] = TOTAL(HI + He)
    print,gmass[i]
    gmass[i] = TOTAL(HI)*1.4;scaled up by primordial HeI (Sanders and McGaugh 02); Geha 06
    smass[i] = TOTAL(s.mass)*msol_per_sysmass[i]
    gmassall[i] = TOTAL(g.mass)
ENDFOR
RETURN, gmass
END
