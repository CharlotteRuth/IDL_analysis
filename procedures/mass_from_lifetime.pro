;Charlotte Christensen
;8/19/2016
;Calculates the stellar mass [in solar masses] corresponding to stellar
;lifetime dStarTime [in yr]. 
;This code is adapted from startime.C in ChaNGA,
;which in turn references the Raiteri, Villata, and Navarro (A&A, 315, 105, 1996) fit to Padova group stellar models
;dMetals is defined as the total mass fraction in elements heavier than
;He, is in the range 0.0004-0.05.
;Both dStarLtime and dMetals can be scalars or arrays 

FUNCTION mass_from_lifetime, dStarLtime, dMetals
StarMass = fltarr(n_elements(dStarLtime))
zmin = 7e-5
zmax = 3e-2
dMetals = float(dMetals)

IF (where(dMetals LT zmin))[0] NE -1 THEN dMetals[where(dMetals LT zmin)] = zmin
IF (where(dMetals GT zmax))[0] NE -1 THEN dMetals[where(dMetals GT zmax)] = zmax
 
a0 = 10.13 + 0.07547*alog10(dMetals) - 0.008084*(alog10(dMetals))^2
a1 = -4.424 - 0.7939*alog10(dMetals) - 0.1187*(alog10(dMetals))^2
a2 = 1.262 + 0.3385*alog10(dMetals) + 0.05417*(alog10(dMetals))^2

c = a0
c -= alog10(dStarLtime)
b = a1
a = a2

logStarMass = (-b - sqrt(b*b - 4*a*c))/(2*a)
StarMass = 10.^logStarMass

IF (where(b*b - 4*a*c LT 0.0))[0] NE -1 THEN $
  StarMass(where(b*b - 4*a*c LT 0.0)) = 1000
; time is too small for fitting formula */

IF (where(finite(StarMass) EQ 0))[0] NE -1 THEN $
  StarMass(where(finite(StarMass) EQ 0)) = 1000
; time is too small for fitting formula */

RETURN, StarMass
END
