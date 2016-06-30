PRO krumholtz,z,x,y
;Krumholz et al, 2009
;ON THE ABSENCE OF HIGH METALLICITY-HIGH COLUMN DENSITY DAMPED
;Ly\u03b1 SYSTEMS: MOLECULE FORMATION IN A TWO-PHASE INTERSTELLAR
;MEDIUM
;eq 1

nel = 1000.0
range = 5.0
min = -1.
x = 10^(findgen(nel)/nel*range + min)

e = 2.71828 ;Euler's number
zsol = 0.0177 ;solar metallicity
surf_den_convert = 1.258075e20 ;unit conversion from amu cm^-2 to msol pc^-2

sigmad = 1e-21 ;eq 3

phi = 1.0  ;3.0 ;eq 18
phimol = 9.6 ;eq 18

chi = 2.3*(1 + 3.1*(z/zsol)^0.365)/phi  ;eq 7
psi = chi*(2.5 + chi)/(2.5 + chi*e) ;eq 9 psi = 1.6 in solar neighborhood

tr = x*surf_den_convert*sigmad ; eq 14

a = 0.2 ;eq 25
trd = tr + a*psi ;eq 25
tc = tr*(phimol + 3.0*psi/(4.0*trd)*(1 - phimol)) ;eq 26
tc[where(tc lt 0)] = 0

y = 1 - (3.0*psi/4.0/tc)/(1 + 4*a*psi*phimol/(4*tc + 3.0*(phimol - 1.0)*psi)) ;eq 34
;s = x*(z/zsol)/psi                             ;eq 36 and 37
;y = (1 + (s/11.0)^3*((125.0 + s)/(96.0 + s))^3)^(1.0/3.0)

n20 = x*surf_den_convert/1e20
zprime = z/zsol
s = alog(1.0 + 0.6*chi)/(0.045*n20*zprime)
delta = 0.0712*(0.1*s^(-1) + 0.675)^(-2.8)
y = 1.0 - (1.0 + (3.0/4.0*s/(1.0 + delta))^(-5.0))^(-1.0/5.0)

IF (where(y lt 0))[0] ne -1 THEN y[where(y lt 0)] = 0

END

