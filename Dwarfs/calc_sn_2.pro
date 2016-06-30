PRO sn_quick,filename
;sn_quick,'200pc.00300'
;This function will do an analysic calculation using gasoline functions:
;startime.c and millerscalo.c
;It reads in the collection of star particles and from their intial
;mass and formation time, determins the number of supernovas over time

;maxtime = 3.e9
;tbinsize=1.e6
;ntbins = maxtime/tbinsize
timeunit=1.e9
massunit=2.325e5 ; massunit in terms of solar masses
;tbins = findgen(ntbins)*tbinsize
;tbins_next =  findgen(ntbins)*tbinsize + tbinsize
;snbins = fltarr(ntbins)*tbinsize

s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
base = '/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/'
file = base+filename
;file='o'+m[mct]+'_1.00300'
rtipsy,file,h,g,d,s
nstars = N_ELEMENTS(s)
s.tform=s.tform*timeunit
s.mass=s.mass*massunit
print,'Number of Stars: ',nstars

z = MEAN(s.metals)
a = 1.262 + 0.3385*alog10(Z) - 0.1187*(alog10(Z))^2
b = -4.424 - 0.7939*alog10(Z) - 0.1187*(alog10(Z))^2
c = 10.13 + 0.07547*ALOG10(Z) - 0.008084*(ALOG10(Z))^2
time_lower = findgen(9.)*1.e6
time_lower[0] = 0.1
time_upper = findgen(9.)*1.e6 + 1.e6

time_upper = ALOG10(time_upper)
time_lower = ALOG10(time_lower)

c_upper = c - time_upper
c_lower = c - time_lower

masses_upper = 10.^((-1*SQRT(b*b - 4.*a*c_upper) - b)/(2.0*a))
masses_lower = 10.^((-1*SQRT(b*b - 4.*a*c_lower) - b)/(2.0*a))
masses_lower[0] = 150.0

nstars = int_krupa(masses_lower,masses_upper)

print,nstars*MIN(s.mass)
print,nstars*MEAN(s.mass)
stop


END


FUNCTION krupa,smass
;For stars with masses greater than 1 solar mass (parameters from A&A,
;315, 1996)
a = 0.3029
b = -1.7
nstars = a * smass^b
RETURN,nstars
END

FUNCTION int_krupa,mass1,mass2
;Integratees the IMF to give total number of stars between bounding
;masses
a = 0.3029
b = -1.7
RETURN,a/b*(mass1^b - mass2^b)
END

FUNCTION nsnova,time,z
;This calculates the number of supernovas that have gone off at a
;given amount of time in a star particle of general age
c = 10.13 + 0.07547*ALOG10(Z) - 0.008084*(ALOG10(Z))^2
b = -4.424 - 0.7939*alog10(Z) - 0.1187*(alog10(Z))^2
a = 1.262 + 0.3385*alog10(Z) - 0.1187*(alog10(Z))^2
log_age = ALOG10(time)
c = c - log_age
masses = (-1*SQRT(b*b - 4.*a*c) - b)/(2.0*a)
masses = 10.^masses             ;An array of masses with lifetimes t
ind_massive = where(masses gt 8.0) ;The indicies for hte stars which will have a type II supernova 
time_lower = time-1e6
masses_lower = masses[0,masses[0:N_ELEMENTS(masses)-2]]
imf = krupa(masses)
imf_lower = krupa(masses_lower)
num_sn = fltarr(N_ELEMENTS(masses))
num_sn[ind_massive] = (imf[ind_massive] + imf_lower[ind_massive])/2.0*1e6
print
stop
RETURN,num_sn
END

PRO calc_sn,filename
;calc_sn,'200pc.00300'
;This function will do an analysic calculation using gasoline functions:
;startime.c and millerscalo.c
;It reads in the collection of star particles and from their intial
;mass and formation time, determins the number of supernovas over time

maxtime = 3.e9
tbinsize=1.e6
ntbins = maxtime/tbinsize
timeunit=1.e9
massunit=2.325e5 ; massunit in terms of solar masses
tbins = findgen(ntbins)*tbinsize
tbins_next =  findgen(ntbins)*tbinsize + tbinsize
snbins = fltarr(ntbins)*tbinsize

s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
base = '/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/10M/200pc/'
file = base+filename
;file='o'+m[mct]+'_1.00300'
rtipsy,file,h,g,d,s
nstars = N_ELEMENTS(s)
s.tform=s.tform*timeunit
s.mass=s.mass*massunit
print,'Number of Stars: ',nstars
FOR istar = 0, nstars-1 DO BEGIN
    ind = WHERE(tbins GE s[istar].tform)
    IF (s[istar].metals eq 0.0) THEN z = 0.0004 ELSE z = s[istar].metals
    IF(ind[0] ne -1) THEN snbins[ind] = snbins[ind] + s[istar].mass*(nsnova(tbins_next[ind] - s[istar].tform,z))
    IF(istar MOD 1000 EQ 0) THEN BEGIN 
        print,istar
        stars = istar + 1
        plot,tbins/1e9,snbins,xtitle = 'time',ytitle = 'Number of Supernovas (per million years)',title = stars+'Stars'
    endif
ENDFOR

plot,tbins/1e9,snbins,xtitle = 'time',ytitle = 'Number of Supernovas (per million years)'
stop
END

