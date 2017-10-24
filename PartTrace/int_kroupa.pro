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

