FUNCTION IMF_KROUPA01, m

a1 = 0.22038*alog(10.0)*2
b1 = -0.3
m1 = 0.08

a2 = 0.22038*alog(10.0)
b2 = -1.3
m2 = 0.5

m3 = 100

imf_m = fltarr(n_elements(m))
ind1 = where(m GE m1 AND m LT m2)
ind2 = where(m GE m2 AND m LE m3)

IF (ind1[0] NE -1) THEN imf_m[ind1] = a1*m[ind1]^b1
IF (ind2[0] NE -1) THEN imf_m[ind2] = a2*m[ind2]^b2

RETURN,imf_m
END

FUNCTION FRACTION_REM, m

f_rem = fltarr(n_elements(m))

ind1 = where(m LT 1)
ind2 = where(m GE 1 AND m LT 8)
ind3 = where(m GE 8)

IF (ind1[0] NE -1) THEN f_rem[ind1] = 1
IF (ind2[0] NE -1) THEN f_rem[ind2] = 1 - (0.86 - exp(-m[ind2]/1.1))
IF (ind3[0] NE -1) THEN f_rem[ind3] = 1.4/m[ind3]

RETURN,f_rem
END

FUNCTION MaxStarMass, dStarLtime, dMetals

;finds stellar mass in solar masses corresponding to stellar lifetime
;dStarTime in yr 

StarMass = fltarr(n_elements(dStarLtime))

IF (where(dStarLtime LE 0.0))[0] THEN StarMass[where(dStarLtime LE 0.0)] = 1000
  
zmax = 3e-2
zmin = 7e-5

IF (where(dMetals LT zmin))[0] NE -1 THEN dMetals[where(dMetals LT zmin)] = zmin
IF (where(dMetals GT zmax))[0] NE -1 THEN dMetals[where(dMetals GT zmax)] = zmax

a00 = 10.13
a01 = 0.07547
a02 = 0.008084
a10 = -4.424
a11 = -0.7939
a12 = -0.1187
a20 = 1.262
a21 = 0.3385
a22 = 0.05417

logZ = alog10(dMetals)                       
a0 = a00 + a01*logZ + a02*logZ*logZ         
a1 = a10 + a11*logZ + a12*logZ*logZ         
a2 = a20 + a21*logZ + a22*logZ*logZ         

c = a0                          
c = c - alog10(dStarLtime)      
b = a1                          
a = a2                          
; time is too small for fitting formula */
IF (where(b*b - 4*a*c LT 0.0))[0] NE -1 THEN StarMass[where(b*b - 4*a*c LT 0.0)] = 1000
logStarMass = (-b - sqrt(b*b - 4*a*c))/(2*a) 
StarMass[where(dStarLtime GT 0.0 AND (b*b - 4*a*c GE 0.0))] = 10.^logStarMass[where(dStarLtime GT 0.0 AND (b*b - 4*a*c GE 0.0))]
RETURN,StarMass                             
END

FUNCTION INTEGRATE_SSP, age, zmetal
ms_turnoff = MaxStarMass(age, zmetal)

m_min = 0.08
m_max = 100.0

turnoff_arr = 10^(findgen(1000)/1000*3 - 1)
star_frac_arr = fltarr(n_elements(turnoff_arr))

n = 1000
bin = (m_max - m_min)/n
mi = (findgen(n) + 1)*bin


FOR i = 0, n_elements(turnoff_arr) - 1 DO BEGIN
   mf_ms = mi[where(mi LT turnoff_arr[i])]
   mf_rem = mi[where(mi GT turnoff_arr[i])] 
   
   total_m_ms = total(imf_kroupa01(mf_ms)*bin)
   total_m_rem = total(imf_kroupa01(mf_rem)*fraction_rem(mf_rem)*bin)
   star_frac_arr[i] = total_m_ms/(total_m_rem + total_m_ms)
ENDFOR
;stop
;print,'Main Sequence Turn-off Mass [M_sol]: ',m_turnoff

;mf_ms = mi[where(mi LT m_turnoff)]
;mf_rem = mi[where(mi GT m_turnoff)] 

;total_m_ms = total(imf_kroupa01(mf_ms)*bin)
;total_m_rem = total(imf_kroupa01(mf_rem)*fraction_rem(mf_rem)*bin)
;print,'Main Sequence Mass Fraction: ',total_m_ms/total(imf_kroupa01(mi)*bin)

ms_turnoff_uniq = ms_turnoff[uniq(ms_turnoff,sort(ms_turnoff))]

star_frac_uniq = spline(turnoff_arr,star_frac_arr,ms_turnoff_uniq)
IF (where(star_frac_uniq GT 1))[0] NE -1 THEN star_frac_uniq[where(star_frac_uniq GT 1)] = 1
IF (where(finite(star_frac_uniq,/nan)))[0] NE -1 THEN star_frac_uniq[where(finite(star_frac_uniq,/nan))] = 1
match2,ms_turnoff,ms_turnoff_uniq,ind1,ind2
RETURN,star_frac_uniq[ind1]
;total_m_ms/(total_m_rem + total_m_ms)
END
