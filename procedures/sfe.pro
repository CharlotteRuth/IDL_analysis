;Calculates the star formation efficiency of an array of gas particles

FUNCTION sfe,g,H2 = H2,HI = HI,cstar = cstar, tempcut = tempcut, denscut = denscut
gm_per_H = 1.673534d-24
grav = 6.6725985e-8 ;am^3/g/s
s_per_yr = 3.15569e7

IF keyword_set(H2) THEN BEGIN
    IF NOT keyword_set(cstar) THEN cstar = 0.1
    IF NOT keyword_set(tempcut) THEN tempcut = 1e3
    IF NOT keyword_set(denscut) THEN denscut = 0.1
ENDIF ELSE BEGIN
    IF NOT keyword_set(cstar) THEN cstar = 0.1
    IF NOT keyword_set(tempcut) THEN tempcut = 1e4
    IF NOT keyword_set(denscut) THEN denscut = 10
ENDELSE
deltat = 1e6*s_per_yr           ;yr
sfeff = fltarr(n_elements(g))
indsf = where(g.tempg LE tempcut AND g.dens GE denscut)
tdyn = 1.0/sqrt(4.0*!PI*grav*g.dens*gm_per_H)
IF keyword_set(H2) THEN $
  IF indsf[0] NE -1 THEN sfeff[indsf] = 1.0 - exp(-1.0*cstar*deltat/tdyn[indsf]*2.0*H2[indsf]/(2.0*H2[indsf] + HI[indsf])) $
ELSE $
  IF indsf[0] NE -1 THEN sfeff[indsf] = 1.0 - exp(-1.0*cstar*deltat/tdyn[indsf])

RETURN,sfeff
END
