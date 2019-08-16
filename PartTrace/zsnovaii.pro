FUNCTION zsnovaii,tmin,tmax,s,sfmass,which_imf
;This calculates the metals produced by supernovas that have gone off
;within a given time frame for star particles of a range of ages

;(For a Kroupa 93 IMF) 1e4 Msol should produce 48.6 SNII, ejecting 611 Msol of material, 2.2
;solar masses of iron and 47.5 solar masses of oxygen
z_outflow = {massloss: 0.0, femassloss: 0.0, oxmassloss: 0.0, zmassloss: 0.0}

zmin = 7e-5
zmax = 3e-2
zsol = 0.02
masses = sfmass
dMSNrem = 1.4 ;mass of SNII remnants

dMetals = s.metals
IF (where(dMetals LT zmin))[0] NE -1 THEN dMetals[where(dMetals LT zmin)] = zmin
IF (where(dMetals GT zmax))[0] NE -1 THEN dMetals[where(dMetals GT zmax)] = zmax

agemin = tmin - s.tform
agemax = tmax - s.tform

IF (where(s.tform GT tmin))[0] NE -1 THEN agemin[where(s.tform GT tmin)] = 0 ;For all stars formed after the timebin started, set their minimum age to zero
IF (where(s.tform GT tmax))[0] NE -1 THEN agemax[where(s.tform GT tmax)] = 0 ;For all stars formed after the timebin bin, set their minimum age to zero. This should result in zero mass producing SNII after the integration

mmax = mass_from_lifetime(agemin,dMetals)
mmin = mass_from_lifetime(agemax,dMetals)
;print,''
;print,'Min Mass: ',mmin,', Max Mass: ',mmax
IF keyword_set(debug) THEN BEGIN
    masses = 1e4
    mmax = 100
    mmin = 1
ENDIF

ind_massive = where(mmax GE 8.0 AND mmin LE 40.0) ;The indicies for the stars which will have a type II supernova 

IF (ind_massive[0] EQ -1) THEN BEGIN
    z_outflow = {massloss: 0.0, femassloss: 0.0, oxmassloss: 0.0, zmassloss: 0.0, stellarRemMass: 0.0,stellarRemZ: 0.0}
    RETURN,z_outflow ;If there are no star particles young enough to have produced SNII, return
ENDIF

masses_upper = mmax[ind_massive]
masses_lower = mmin[ind_massive]

;Set the minimum integration bound to the maximum of 8 or the minmum
;mass that would have produced supernovae
IF (where(masses_lower LT 8))[0] NE -1 THEN masses_lower[where(masses_lower LT 8)] = 8 
; Only stars less that 40 solar masses produce SNII
IF (where(masses_upper GT 40))[0] NE -1 THEN masses_upper[where(masses_upper GT 40)] = 40
IF (where(masses_lower GT 40))[0] NE -1 THEN masses_lower[where(masses_lower GT 40)] = 40 

; I haven't tested the integration of the IMF with anything
; other than Kroupa 93 and Kroupa 01

mass_sn = fltarr(n_elements(mmax))
num_sn = fltarr(n_elements(mmax))
IF (which_imf EQ 0) THEN mass_sn[ind_massive] = int_kroupa(masses_upper,masses_lower) ;ELSE mass_sn[ind_massive] = int_ms(masses_upper,masses_lower)
IF (which_imf EQ 0) THEN num_sn[ind_massive] = int_kroupa(masses_upper,masses_lower,/nstar) ;ELSE num_sn[ind_massive] = int_ms(masses_upper,masses_lower,/nstar)
IF (which_imf EQ 1) THEN mass_sn[ind_massive] = int_kroupa01(masses_upper,masses_lower) ;ELSE mass_sn[ind_massive] = int_ms(masses_upper,masses_lower)
IF (which_imf EQ 1) THEN num_sn[ind_massive] = int_kroupa01(masses_upper,masses_lower,/nstar) ;ELSE num_sn[ind_massive] = int_ms(masses_upper,masses_lower,/nstar)
mass_sn = mass_sn*sfmass
num_sn = num_sn*sfmass

;metallicity (Ox and Fe fractions) for SN ejecta as calculated from
;Raiteri, Villata and Navarro, 1996 (eqn 6-8)
ox_frac = fltarr(n_elements(mmax))
fe_frac = fltarr(n_elements(mmax))
IF (which_imf EQ 0) THEN dMassLoss = (0.7682*int_kroupa(masses_upper,masses_lower,/nstar,mod_exp = 1.056)) ;As calculated by Raiteri
IF (which_imf EQ 1) THEN dMassLoss = (0.7682*int_kroupa01(masses_upper,masses_lower,/nstar,mod_exp = 1.056)) ;As calculated by Raiteri
ind_massloss = where(dMassLoss NE 0)
IF dMassLoss[0] NE -1 THEN $
   IF (which_imf EQ 0) THEN ox_frac[ind_massive[ind_massloss]] = 4.586e-4*int_kroupa(masses_upper[ind_massloss],masses_lower[ind_massloss],/nstar,mod_exp = 2.721)/dMassLoss[ind_massloss] $ ;this should be updated to allow for other IMFs
   ELSE IF (which_imf EQ 1) THEN ox_frac[ind_massive[ind_massloss]] = 4.586e-4*int_kroupa01(masses_upper[ind_massloss],masses_lower[ind_massloss],/nstar,mod_exp = 2.721)/dMassLoss[ind_massloss]
IF dMassLoss[0] NE -1 THEN $
   IF (which_imf EQ 0) THEN fe_frac[ind_massive[ind_massloss]] = 2.802e-4*int_kroupa(masses_upper[ind_massloss],masses_lower[ind_massloss],/nstar,mod_exp = 1.864)/dMassLoss[ind_massloss] $
   ELSE IF (which_imf EQ 1) THEN fe_frac[ind_massive[ind_massloss]] = 2.802e-4*int_kroupa01(masses_upper[ind_massloss],masses_lower[ind_massloss],/nstar,mod_exp = 1.864)/dMassLoss[ind_massloss]
z_frac =  1.06*fe_frac + 2.09*ox_frac

IF keyword_set(debug)  THEN dMassLoss = dMassLoss*1e4 ELSE dMassLoss = mass_sn - num_sn*dMSNrem

;print,'Mass of SNII: ',total(mass_sn)
;print,'Number of SNII: ',total(num_sn)
;print,'Mass Lost: ',total(dMassLoss),total(dMassLoss);/1.66500e+09
oxMassLoss = dMassLoss*ox_frac
;print,'Ox Mass Lost: ',total(oxMassLoss),total(oxMassLoss)/total(dMassLoss)
feMassLoss = dMassLoss*fe_frac
;print,'Fe Mass Lost: ',total(feMassLoss),total(feMassLoss)/total(dMassLoss)
zMassLoss = dMassLoss*z_frac
;print,'Metal Mass Lost: ',total(zMassLoss),total(zMassLoss)/total(dMassLoss)
;Mass of star particles now in the form of stellar remnants
stellarRemMass = num_sn*dMSNrem
;Mass of metals now in the form of stellar remnants
stellarRemZ = num_sn*dMSNrem*dMetals

z_outflow = {massloss: float(total(dMassLoss)), femassloss: float(total(feMassLoss)), oxmassloss: float(total(oxMassLoss)), zmassloss: float(total(zMassLoss)),stellarRemMass: float(total(stellarRemMass)),stellarRemZ: float(total(stellarRemZ))}
RETURN,z_outflow
END
