;Very similar to calc_sn but slightly cleaner

 ; .r calc_sn

FUNCTION int_kroupa01,mass1,mass2,nstar = nstar,mod_exp = mod_exp
;Integratees the IMF to give total number of stars between bounding
;masses, mass1 and mass2, where mass1 is the greater mass

IF n_elements(mass1) NE n_elements(mass2) THEN BEGIN
    print,'mass1 and mass2 must have the same number of elements'
    RETURN,0*mass1
ENDIF

bmasses = [0.08,0.5,100] ;bounding masses and mass of break in IMF
exp = [-0.3,-1.3] ;exponents of IMF to produce mass of stars in given range. IMF to produce # will have exp-1
a = [2.0*alog(10.0),alog(10.0)]/14.4981 ;coefficients for IMF
IF keyword_set(nstar) THEN exp = exp-1  ;so that it is no longer in terms of mass but is in number
IF keyword_set(mod_exp) THEN exp = exp + mod_exp ;Enables one to include a function of stellar mass (for instance, oxygen mass released) inside of the integral

IF (where(mass1 GT bmasses[2]))[0] NE -1 THEN mass1[where(mass1 GT bmasses[2])] = bmasses[2]
IF (where(mass2 GT bmasses[2]))[0] NE -1 THEN mass2[where(mass2 GT bmasses[2])] = bmasses[2]
IF (where(mass1 LT bmasses[0]))[0] NE -1 THEN mass1[where(mass1 LT bmasses[0])] = bmasses[0]
IF (where(mass2 LT bmasses[0]))[0] NE -1 THEN mass2[where(mass2 LT bmasses[0])] = bmasses[0]

cumM1 = fltarr(n_elements(mass1)) ;Create an array to hold the total mass between mass1 and the maximum mass

ind1 = where(mass1 GT bmasses[1] AND mass1 LE bmasses[2]) ;Anything that is in the high-mass part of the IMF
ind2 = where(mass1 GE bmasses[0] AND mass1 LE bmasses[1]) ;Anything that is in the low-mass part of the IMF

IF (ind1[0] NE -1) THEN cumM1[ind1] = a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - mass1[ind1]^(exp[1] + 1))
IF (ind2[0] NE -1) THEN cumM1[ind2] = a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - bmasses[1]^(exp[1] + 1)) + $
                                      a[0]/(exp[0] + 1)*(bmasses[1]^(exp[0] + 0) - mass1[ind2]^(exp[0] + 1))

cumM2 = fltarr(n_elements(mass1)) ;Create an array to hold the total mass between mass2 and the maximum mass

ind1 = where(mass2 GT bmasses[1] AND mass2 LE bmasses[2]) ;Anything that is in the high-mass part of the IMF
ind2 = where(mass2 GE bmasses[0] AND mass2 LE bmasses[1]) ;Anything that is in the low-mass part of the IMF

IF (ind1[0] NE -1) THEN cumM2[ind1] = a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - mass2[ind1]^(exp[1] + 1))
IF (ind2[0] NE -1) THEN cumM2[ind2] = a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - bmasses[1]^(exp[1] + 1)) + $
                                      a[0]/(exp[0] + 1)*(bmasses[1]^(exp[0] + 0) - mass2[ind2]^(exp[0] + 1))

RETURN,cumM2 - cumM1
END


FUNCTION int_kroupa,mass1,mass2,nstar = nstar,mod_exp = mod_exp
;Integratees the IMF to give total number of stars between bounding
;masses, mass1 and mass2, where mass1 is the greater mass

IF n_elements(mass1) NE n_elements(mass2) THEN BEGIN
    print,'mass1 and mass2 must have the same number of elements'
    RETURN,0*mass1
ENDIF

bmasses = [0.08,0.5,1.0,100] ;bounding masses and mass of break in IMF
exp = [-0.3,-1.2,-1.7] ;exponents of IMF to produce mass of stars in given range. IMF to produce # will have exp-1
a = [2^(0.9),1.0,1.0]*0.3029 ;coefficients for IMF
IF keyword_set(nstar) THEN exp = exp - 1  ;so that it is no longer in terms of mass but is in number
IF keyword_set(mod_exp) THEN exp = exp + mod_exp ;Enables one to include a function of stellar mass (for instance, oxygen mass released) inside of the integral

IF (where(mass1 GT bmasses[3]))[0] NE -1 THEN mass1[where(mass1 GT bmasses[3])] = bmasses[3]
IF (where(mass2 GT bmasses[3]))[0] NE -1 THEN mass2[where(mass2 GT bmasses[3])] = bmasses[3]
IF (where(mass1 LT bmasses[0]))[0] NE -1 THEN mass1[where(mass1 LT bmasses[0])] = bmasses[0]
IF (where(mass2 LT bmasses[0]))[0] NE -1 THEN mass2[where(mass2 LT bmasses[0])] = bmasses[0]

cumM1 = fltarr(n_elements(mass1)) ;Create an array to hold the total mass between mass1 and the maximum mass

ind0 = where(mass1 GT bmasses[2] AND mass1 LE bmasses[3]) ;Anything that is in the high-mass part of the IMF
ind1 = where(mass1 GT bmasses[1] AND mass1 LE bmasses[2]) ;Anything that is in the med-mass part of the IMF
ind2 = where(mass1 GE bmasses[0] AND mass1 LE bmasses[1]) ;Anything that is in the low-mass part of the IMF
IF (ind0[0] NE -1) THEN cumM1[ind0] = a[2]/(exp[2] + 1)*(bmasses[3]^(exp[2] + 1) - mass1[ind0]^(exp[2] + 1))
IF (ind1[0] NE -1) THEN cumM1[ind1] = a[2]/(exp[2] + 1)*(bmasses[3]^(exp[2] + 1) - bmasses[2]^(exp[2] + 1)) + $
                                      a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - mass1[ind1]^(exp[1] + 1))
IF (ind2[0] NE -1) THEN cumM1[ind2] = a[2]/(exp[2] + 1)*(bmasses[3]^(exp[2] + 1) - bmasses[2]^(exp[2] + 1)) + $
                                      a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - bmasses[1]^(exp[1] + 1)) + $
                                      a[0]/(exp[0] + 1)*(bmasses[1]^(exp[0] + 1) - mass1[ind2]^(exp[0] + 1))

cumM2 = fltarr(n_elements(mass1)) ;Create an array to hold the total mass between mass2 and the maximum mass


ind0 = where(mass2 GT bmasses[2] AND mass2 LE bmasses[3]) ;Anything that is in the high-mass part of the IMF
ind1 = where(mass2 GT bmasses[1] AND mass2 LE bmasses[2]) ;Anything that is in the med-mass part of the IMF
ind2 = where(mass2 GE bmasses[0] AND mass2 LE bmasses[1]) ;Anything that is in the low-mass part of the IMF
IF (ind0[0] NE -1) THEN cumM2[ind0] = a[2]/(exp[2] + 1)*(bmasses[3]^(exp[2] + 1) - mass2[ind0]^(exp[2] + 1))
IF (ind1[0] NE -1) THEN cumM2[ind1] = a[2]/(exp[2] + 1)*(bmasses[3]^(exp[2] + 1) - bmasses[2]^(exp[2] + 1)) + $
                                      a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - mass2[ind1]^(exp[1] + 1))
IF (ind2[0] NE -1) THEN cumM2[ind2] = a[2]/(exp[2] + 1)*(bmasses[3]^(exp[2] + 1) - bmasses[2]^(exp[2] + 1)) + $
                                      a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - bmasses[1]^(exp[1] + 1)) + $
                                      a[0]/(exp[0] + 1)*(bmasses[1]^(exp[0] + 1) - mass2[ind2]^(exp[0] + 1))
RETURN,cumM2 - cumM1
END

FUNCTION zsnovaii,tmin,tmax,s,sfmass,which_imf
;This calculates the metals produced by supernovas that have gone off
;within a given time frame for star particles of a range of ages

;1e4 Msol should produce 48.6 SNII, ejecting 611 Msol of material, 2.2
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
IF (where(s.tform GT tmax))[0] NE -1 THEN agemax[where(s.tform GT tmin)] = 0 ;For all stars formed after the timebin bin, set their minimum age to zero. This should result in zero mass producing SNII after the integration

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
; other than Kroupa 93

mass_sn = fltarr(n_elements(mmax))
num_sn = fltarr(n_elements(mmax))
IF (which_imf EQ 0) THEN mass_sn[ind_massive] = int_kroupa(masses_upper,masses_lower) ;ELSE mass_sn[ind_massive] = int_ms(masses_upper,masses_lower)
IF (which_imf EQ 0) THEN num_sn[ind_massive] = int_kroupa(masses_upper,masses_lower,/nstar) ;ELSE num_sn[ind_massive] = int_ms(masses_upper,masses_lower,/nstar)
mass_sn = mass_sn*sfmass
num_sn = num_sn*sfmass

;metallicity (Ox and Fe fractions) for SN ejecta as calculated from
;Raiteri, Villata and Navarro, 1996 (eqn 6-8)
ox_frac = fltarr(n_elements(mmax))
fe_frac = fltarr(n_elements(mmax))
dMassLoss = (0.7682*int_kroupa(masses_upper,masses_lower,/nstar,mod_exp = 1.056)) ;As calculated by Raiteri
ind_massloss = where(dMassLoss NE 0)
IF dMassLoss[0] NE -1 THEN ox_frac[ind_massive[ind_massloss]] = 4.586e-4*int_kroupa(masses_upper[ind_massloss],masses_lower[ind_massloss],/nstar,mod_exp = 2.721)/dMassLoss[ind_massloss] ;this should be updated to allow for other IMFs
IF dMassLoss[0] NE -1 THEN fe_frac[ind_massive[ind_massloss]] = 2.802e-4*int_kroupa(masses_upper[ind_massloss],masses_lower[ind_massloss],/nstar,mod_exp = 1.864)/dMassLoss[ind_massloss]
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


PRO test_metal_production
dirbase = '/nobackupp8/crchrist/MolecH/' ;PFE
;dirbase = '/home/christensen/Storage2/UW/MolecH/Cosmo/';ymir

steps = ['00024','00048','00072','00096','00120','00144','00168','00192','00216','00240','00264','00288','00312','00336','00360','00384','00408','00432','00456','00480','00504','00512']
steps = ['00024','00096','00168','00240','00312','00360','00480','00512']
steps = ['00512']

;fileroot_short = 'h516.cosmo25cmb.1536g'
;fileroot = 'h516.cosmo25cmb.1536g14HBWK' 
;fileroot_short = 'h516.cosmo25cmb.3072g'
;fileroot = 'h516.cosmo25cmb.3072g14HBWK'
;fileroot_short = 'h603.cosmo50cmb.3072g'
;fileroot = 'h603.cosmo50cmb.3072g14HBWK'
fileroot_short = 'h799.cosmo25cmb.3072g/'
fileroot = 'h799.cosmo25cmb.3072g14HBWK/'

dir = dirbase + fileroot_short + '/' + fileroot
filebase = fileroot + '.'+ steps + '/' + fileroot + '.' + steps
pfile = fileroot + '.param'
;dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/'
;filebase = 'h516.cosmo25cmb.3072g14HBWK.' + steps + '.dir/h516.cosmo25cmb.3072g14HBWK'
;pfile = '../h516.cosmo25cmb.3072g14HBWK.param'

dDelta = 0.000690912 ;4.075253e-04
dDeltaStarForm = 1e6

IF 0 THEN BEGIN
   dir = '/home/christensen/Code/gasoline-radiation/test/onestar/'
   pfile = 'onestar.param'
;steps = ['00020','00031','00032','00033','00034','00035','00036','00037','01000']
;steps = ['00001','00002','00003','00004','00005','00006','00007','00008','00009','00010','00011','00012','00013','00014','00015','00016','00017','00018','00019','00020','00021','00022','00023','00024','00025','00026','00027','00028','00029','00030','00031','00032','00033','00034','00035','00036','00037','00038','00039','00040']
   steps = ['01000','02000','03000','04000','05000','06000','07000','08000','09000','10000']
;steps = ['00100','00200','00300','00400','00500','00600','00700','00800','00900','01000']
;steps = ['00010','00020','00030','00040','00050','00060','00070','00080','00090','00100']
   filebase = 'onestar.' + steps
;filebase = 'onestar.z0.01.' + steps
   dDelta = 7.812500e-04
ENDIF

cd,dir
units = tipsyunits(pfile)
dDelta = dDelta*units.timeunit
IF keyword_set(dDeltaStarForm) THEN dDeltaStarForm = dDelta
totalOxMass = fltarr(n_elements(steps))
totalFeMass = fltarr(n_elements(steps))
totalOxOut = fltarr(n_elements(steps))
totalFeOut = fltarr(n_elements(steps))

FOR i=0,n_elements(steps)- 1 DO BEGIN
    step = steps[i]
    filename = filebase[i]
    rtipsy,filename,h,g,d,s
    readarr,filename + '.FeMassFrac',h,fe_gas,/ascii,part = 'gas',type = 'float'
    readarr,filename + '.FeMassFrac',h,fe_star,/ascii,part = 'star',type = 'float'
    readarr,filename + '.OxMassFrac',h,ox_gas,/ascii,part = 'gas',type = 'float'
    readarr,filename + '.OxMassFrac',h,ox_star,/ascii,part = 'star',type = 'float'
    readarr,filename + '.massform',h,starmass,/ascii,part = 'star',type = 'float'

;sl = rstarlog('../../h516.cosmo25cmb.1536g14HBWK.starlog')
;s.mass = sl.massform ;make sure these line up

    s.tform = s.tform*units.timeunit
    sfmass = starmass*units.massunit
    s.mass = s.mass*units.massunit

    time = max(s.tform)
;time = 0.354154*units.timeunit    
    time = dDelta*float(step)
;stop
    
    which_imf = 0               ;kroupa93
    tmax = time + dDelta
    tmin = time                 ;- dDelta;0 ;time ;time
    tmin = 0
    
    print,''
    print,i,tmin,tmax,time/units.timeunit
    z_outflow_snii = zsnovaii(tmin, tmax, s, sfmass, which_imf)
    z_outflow_snia = zsnovaia(tmin, tmax, s, sfmass)
 
IF pfile eq 'onestar.param' THEN BEGIN
    zstart_ox =  total(sfmass*ox_star)
    zstart_fe =  total(sfmass*fe_star)
ENDIF ELSE BEGIN
    zstart_ox = 0
    zstart_fe = 0
ENDELSE 

    print,''
    print,'Totals'
    print,'Total Ox Mass: ',(total(g.mass*ox_gas)*units.massunit + total(s.mass*ox_star)),', Stars: ',total(s.mass*ox_star),', Gas: ',total(g.mass*ox_gas)*units.massunit
    totalOxMass[i] = (total(g.mass*ox_gas)*units.massunit + total(s.mass*ox_star))
    print,'Total Ox Out: ',z_outflow_snii.oxmassloss + z_outflow_snia.oxmassloss +  zstart_ox,', from SNII: ',z_outflow_snii.oxmassloss,', from SNIa: ',z_outflow_snia.oxmassloss
    totalOxOut[i] = z_outflow_snii.oxmassloss + z_outflow_snia.oxmassloss +  zstart_ox
    print,'Total Fe Mass: ',(total(g.mass*fe_gas)*units.massunit + total(s.mass*fe_star)),', Stars: ',total(s.mass*fe_star),', Gas: ',total(g.mass*fe_gas)*units.massunit
    totalFeMass[i] = (total(g.mass*fe_gas)*units.massunit + total(s.mass*fe_star))
    print,'Total Fe Out: ',z_outflow_snii.femassloss + z_outflow_snia.femassloss + zstart_fe,', from SNII: ',z_outflow_snii.femassloss,', from SNIa: ',z_outflow_snia.femassloss
    totalFeOut[i] = z_outflow_snii.femassloss + z_outflow_snia.femassloss + zstart_fe
ENDFOR
totalzMass =  1.06*totalFeMass + 2.09*totalOxMass
totalzOut = 1.06*totalFeOut + 2.09*totalOxOut

loadct,39
window,0
plot,float(steps)*dDelta,totalzMass
oplot,float(steps)*dDelta,totalzOut,linestyle = 2
oplot,float(steps)*dDelta,totalOxMass,color = 60
oplot,float(steps)*dDelta,totalOxOut,color = 60,linestyle = 2
oplot,float(steps)*dDelta,totalFeMass,color = 254
oplot,float(steps)*dDelta,totalFeOut,color = 254,linestyle = 2

window,1
plot,float(steps)*dDelta,totalzOut/totalzMass
oplot,float(steps)*dDelta,totalOxOut/totalOxMass,color = 60
oplot,float(steps)*dDelta,totalFeOut/totalFeMass,color = 254
stop
END

;print,(total(g.mass*femassfrac)*units.massunit + total(s.mass*femassfrac))
;      280469.

;print,(total(g.mass*oxmassfrac)*units.massunit + total(s.mass*oxmassfrac))
;  2.37371e+06
