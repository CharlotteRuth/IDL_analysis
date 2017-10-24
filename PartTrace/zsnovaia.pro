FUNCTION zsnovaia,tmin,tmax,s, sfmass
;This calculates the metals produced by supernovae Ia that have gone off
;within a given time frame for star particles of a range of ages

;Should produce 7.8 SNIa per 1e4 solar mass star particle in the
;lifetime of a one solar mass star for Raitari scaling of Aprime

;1e4 Msol should ejecting 11 Msol of material, 4.92
;solar masses of iron and 1.0 solar masses of oxygen

zmin = 7e-5
zmax = 3e-2
zsol = 0.02
dMSNrem = 1.4 ;mass of SNII remnants

dMetals = s.metals
IF (where(dMetals LT zmin))[0] NE -1 THEN dMetals[where(dMetals LT zmin)] = zmin
IF (where(dMetals GT zmax))[0] NE -1 THEN dMetals[where(dMetals GT zmax)] = zmax
;stop

agemin = tmin - s.tform
agemax = tmax - s.tform
IF (where(s.tform GT tmin))[0] NE -1 THEN agemin[where(s.tform GT tmin)] = 0 ;For all stars formed after the timebin started, set their minimum age to zero
IF (where(s.tform GT tmax))[0] NE -1 THEN agemax[where(s.tform GT tmin)] = 0 ;For all stars formed after the timebin bin, set their minimum age to zero. This should result in zero mass producing SNII after the integration

mmax = mass_from_lifetime(agemin,dMetals) ; The maximum mass of stars that may have ended their lives during this time
mmin = mass_from_lifetime(agemax,dMetals) ; the minimum mass of stars that may have ended their lives during this time
;print,''
;print,'Mmin: ',mmin[0],' Mmax: ',mmax[0];,agemin,agemax

;Debug
;mmin = findgen(700,increment = 0.01,start = 1)
;mmax = masses_lower*0 + 8
;if tmax GT 2.7e7 then stop
ind_lowmass = where(mmin LT 8.0) ;The indicies for the stars which will have a type Ia supernova 

IF (ind_lowmass[0] EQ -1) THEN BEGIN
    z_outflow = {massloss: 0.0, femassloss: 0.0, oxmassloss: 0.0, zmassloss: 0.0}
    RETURN,z_outflow ;If there are no star particles old enough to have produced SNIa, return
ENDIF

masses_upper = mmax[ind_lowmass]
masses_lower = mmin[ind_lowmass]
IF (where(masses_upper GT 3./2.))[0] EQ -1 THEN BEGIN
    z_outflow = {massloss: 0.0, femassloss: 0.0, oxmassloss: 0.0, zmassloss: 0.0}
    RETURN,z_outflow ;If there are no star particles old enough to have produced SNIa, return
ENDIF

; Only stars less that 8 solar masses produce SNIa
IF (where(masses_upper GT 8))[0] NE -1 THEN masses_upper[where(masses_upper GT 8)] = 8
IF (where(masses_upper LT 3./2.))[0] NE -1 THEN masses_upper[where(masses_upper LT 3./2.)] = 3./2. 
IF (where(masses_lower LT 3./2.))[0] NE -1 THEN masses_lower[where(masses_lower LT 3./2.)] = 3./2. 


num_sn = fltarr(n_elements(mmax))
indlow = where(masses_lower LT 3./2., complement = indhigh)

;Which of the following integrals you do depends on what version of
;gasoline/Changa you ran the code in. Compare against supernovaia.c

;normalization constant, set by the total number of binaries
dFracBinSNIa = 0.16; 0.05 #Changes depending on the version. Do a grep
imf1to8PreFactor = 0.3029

;Aprime = 0.04847 ;from Raiteri 1996
Aprime = dFracBinSNIa*24*imf1to8PreFactor/alog(10)^2
Aprime = dFracBinSNIa*imf1to8PreFactor;*.85
IF indlow[0] NE -1 THEN num_sn[ind_lowmass[indlow]] = Aprime/3.7*( $
                           3^(-4.7)*((3./2.+fltarr(n_elements(indlow)))^3 - masses_lower[indlow]^3) $
                           - 2^(-3.7)/0.7*(masses_upper[indlow]^(-0.7) - (3./2.+fltarr(n_elements(indlow)))^(-0.7)) $
                           + (1.42857*masses_upper[indlow]^2 + 13.4454*masses_upper[indlow] + 39.8382)*(masses_upper[indlow] + 8)^(-2.7) - (1.42857*masses_lower[indlow]^2 + 13.4454*masses_lower[indlow] + 39.8382)*(masses_lower[indlow] + 8)^(-2.7))
IF indhigh[0] NE -1 THEN num_sn[ind_lowmass[indhigh]] = Aprime/3.7*( $
                           (-1)*2^(-3.7)/0.7*(masses_upper[indhigh]^(-0.7) - masses_lower[indhigh]^(-0.7)) $
                           + ((1.42857*masses_upper[indhigh]^2 + 13.4454*masses_upper[indhigh] + 39.8382)*(masses_upper[indhigh] + 8)^(-2.7) - (1.42857*masses_lower[indhigh]^2 + 13.4454*masses_lower[indhigh] + 39.8382)*(masses_lower[indhigh] + 8)^(-2.7)))


;from Gasoline
IF 0 THEN BEGIN
Aprime = dFracBinSNIa*24*imf1to8PreFactor
IF indlow[0] NE -1 THEN num_sn[ind_lowmass[indlow]] = Aprime/4.7*( $
                           3^(-5.7)*((3./2.+fltarr(n_elements(indlow)))^3 - masses_lower[indlow]^3) $
                           - 2^(-4.7)/1.7*(masses_upper[indlow]^(-1.7) - (3./2.+fltarr(n_elements(indlow)))^(-1.7)) $
                           + (0.588235*masses_upper[indlow]^2 + 3.48584*masses_upper[indlow] + 7.53695)*(masses_upper[indlow] + 8)^(-3.7) $
                           - (0.588235*masses_lower[indlow]^2 + 3.48584*masses_lower[indlow] + 7.53695)*(masses_lower[indlow] + 8)^(-3.7))

IF indhigh[0] NE -1 THEN num_sn[ind_lowmass[indhigh]] = Aprime/4.7*( $
                           (-1)*2^(-4.7)/1.7*(masses_upper[indhigh]^(-1.7) - masses_lower[indhigh]^(-1.7)) $
                           + (0.588235*masses_upper[indhigh]^2 + 3.48584*masses_upper[indhigh] + 7.53695)*(masses_upper[indhigh] + 8)^(-3.7) $
                           - (0.588235*masses_lower[indhigh]^2 + 3.48584*masses_lower[indhigh] + 7.53695)*(masses_lower[indhigh] + 8)^(-3.7))
ENDIF


;plot,mmin,num_sn*1e4 ;Debug

num_sn = num_sn*sfmass
dMassLoss = dMSNrem*num_sn
ind_lowmass = where(dMassLoss NE 0) ;because sometimes a lower mass can be just below 8 and because of rounding num_sn = 0

;metallicity (Ox and Fe fractions) for SN ejecta as calculated from
;Raiteri, Villata and Navarro, 1996 (eqn 6-8)
ox_frac = fltarr(n_elements(mmax))
fe_frac = fltarr(n_elements(mmax))
ox_frac[ind_lowmass] = num_sn[ind_lowmass]*0.13/dMassLoss[ind_lowmass]
fe_frac[ind_lowmass] = num_sn[ind_lowmass]*0.63/dMassLoss[ind_lowmass]
z_frac =  1.06*fe_frac + 2.09*ox_frac
oxMassLoss = dMassLoss*ox_frac
feMassLoss = dMassLoss*fe_frac
zMassLoss = dMassLoss*z_frac

;print,'Number of SNIa: ',total(num_sn)
;print,'Mass Lost: ',total(dMassLoss),total(dMassLoss);/1.66500e+09
;print,'Ox Mass Lost: ',total(oxMassLoss),total(oxMassLoss)/total(dMassLoss)
;print,'Fe Mass Lost: ',total(feMassLoss),total(feMassLoss)/total(dMassLoss)
;print,'Metal Mass Lost: ',total(zMassLoss),total(zMassLoss)/total(dMassLoss)
z_outflow = {massloss: float(total(dMassLoss)), femassloss: float(total(feMassLoss)), oxmassloss: float(total(oxMassLoss)), zmassloss: float(total(zMassLoss))}
;print,float(total(dMassLoss)),float(total(feMassLoss)),float(total(oxMassLoss)),float(total(zMassLoss))
;print,z_outflow
;plot,mmin[where(dMetals GT 0.0007 AND dMetals LT 0.0009)],num_sn[where(dMetals GT 0.0007 AND dMetals LT 0.0009)],psym = 3,xrange = [0.8,1]
;colors = (alog10(dMetals)+4.2)/2.3*254
;for j=0,n_elements(color)-1 DO oplot,[mmin[j],mmin[j]],[num_sn[j],num_sn[j]],psym = 3,color = color[j]
RETURN,z_outflow
END

