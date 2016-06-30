FUNCTION kroupa,smass,nstar = nstar
;For stars with masses greater than 1 solar mass (parameters from A&A,
;315, 1996) -- normalized so that one solar mass of stars are formed
ind1 = WHERE(smass gt 1.0)
ind2 = WHERE(smass gt 0.5 AND smass le 1.0)

a = smass*0 + 0.3029*1.86606 ;for masses less than 0.08
b = smass*0 - 0.3
IF (ind2[0] ne -1) THEN BEGIN
    a[ind2] = 0.3029
    b[ind2] = -1.2
ENDIF
IF (ind1[0] ne -1) THEN BEGIN
    a[ind1] = 0.3029
    b[ind1] = -1.7
ENDIF
IF(KEYWORD_SET(nstar)) THEN b = b-1 ;so that it is no longer in terms of mass but in number of stars
nstars = a * smass^b
RETURN,nstars
END

FUNCTION ms,smass,nstar = nstar
;For stars with masses greater than 1 solar mass (parameters from A&A,
;315, 1996) -- normalized so that one solar mass of stars are formed
;Miller-Scalo IMF
ind1 = WHERE(smass gt 10.0)
ind2 = WHERE(smass gt 1.0 AND smass le 10.0)
c = 118.643 ;It sure looks like MS is scalled to 75 solar masses rather than 1
;c = 1.0

a = smass*0 + 42.0 ;for masses less than 0.08
b = smass*0 - 0.4
IF (ind2[0] ne -1) THEN BEGIN
    a[ind2] = 42.0
    b[ind2] = -1.5
ENDIF
IF (ind1[0] ne -1) THEN BEGIN
    a[ind1] = 240.0
    b[ind1] = -2.3
ENDIF
IF(KEYWORD_SET(nstar)) THEN b = b-1 ;so that it is no longer in terms of mass but in number of stars
nstars = a * smass^b
RETURN,nstars/c
END

FUNCTION int_kroupa,mass1,mass2,nstar = nstar
;Integratees the IMF to give total number of stars between bounding
;masses
ind1 = WHERE(mass1 gt 1.0 AND mass1 le 100.0)
ind2 = WHERE(mass1 gt 0.5 AND mass1 le 1.0)

a = mass1*0 + 0.3029*1.86606 ;for masses less than 0.08
b = mass1*0 - 0.3
IF (ind2[0] ne -1) THEN BEGIN
    a[ind2] = 0.3029
    b[ind2] = -1.2
ENDIF
IF (ind1[0] ne -1) THEN BEGIN
    a[ind1] = 0.3029
    b[ind1] = -1.7
ENDIF
IF(KEYWORD_SET(nstar)) THEN b = b-1  ;so that it is no longer in terms of mass but is in number
;of stars
;stop
RETURN,a/(b+1)*(mass1^(b+1) - mass2^(b+1))
END

FUNCTION int_ms,mass1,mass2,nstar = nstar
ind1 = WHERE(mass1 gt 10.0 AND mass1 le 100.0)
ind2 = WHERE(mass1 gt 1.0 AND mass1 le 10.0)
c = 118.643 ;It sure looks like this is normalized for 75 solar masses
;c = 1.0

a = mass1*0 + 42.0 ;for masses less than 0.08
b = mass1*0 - 0.4
IF (ind2[0] ne -1) THEN BEGIN
    a[ind2] = 42.0
    b[ind2] = -1.5
ENDIF
IF (ind1[0] ne -1) THEN BEGIN
    a[ind1] = 240.0
    b[ind1] = -2.3
ENDIF
IF(KEYWORD_SET(nstar)) THEN b = b-1 ;so that it is no longer in terms of mass but in muber of stars
RETURN,(a/(b+1)*(mass1^(b+1) - mass2^(b+1)))/c
END

FUNCTION nsnova,time,z,which_imf
;This calculates the number of supernovas that have gone off at a
;given amount of time in a star particle of general age
c = 10.13 + 0.07547*ALOG10(Z) - 0.008084*(ALOG10(Z))^2
b = -4.424 - 0.7939*alog10(Z) - 0.1187*(alog10(Z))^2
a = 1.262 + 0.3385*alog10(Z) + 0.05417*(alog10(Z))^2
log_age = ALOG10(time)
c = c - log_age
masses = (-1*SQRT(b*b - 4.*a*c) - b)/(2.0*a)
if((where((b*b - 4.*a*c) lt 0))[0] ne -1) then masses[where((b*b - 4.*a*c) lt 0)] = 2.0
masses = 10.^masses             ;An array of masses with lifetimes t
;stop
ind_massive = where(masses gt 8.0) ;The indicies for the stars which will have a type II supernova 
if (ind_massive[0] ne -1) then ind_massive = [ind_massive, MAX(ind_massive) + 1]else ind_massive = 0
time_lower = time-1e6
IF (N_ELEMENTS(masses) GT 1) THEN masses_lower = [100,masses[0:N_ELEMENTS(masses)-2]] ELSE masses_lower = [100]
IF(which_imf EQ 0) THEN imf = kroupa(masses) ELSE imf = ms(masses);numer of stars in each of the mass bins
IF(which_imf EQ 0) THEN imf_lower = kroupa(masses_lower) ELSE imf_lower = ms(masses_lower)
num_sn = fltarr(N_ELEMENTS(masses))
;print,'First try:    ',(imf[ind_massive] + imf_lower[ind_massive])/2.0*(masses_lower[ind_massive] - masses[ind_massive])
IF(which_imf EQ 0) THEN num_sn[ind_massive] = int_kroupa(masses_lower[ind_massive],masses[ind_massive]) ELSE num_sn[ind_massive] = int_ms(masses_lower[ind_massive],masses[ind_massive])
;print,'Number of SN: ',num_sn[where(num_sn ne 0)]
;stop
RETURN,num_sn
END

function calc_sn,filename,imf
;calc_sn,'200pc.00300'
;This function will do an analysic calculation using gasoline functions:
;startime.c and millerscalo.c
;It reads in the collection of star particles and from their intial
;mass and formation time, determins the number of supernovas over time

;The parameter filename is the name of your tipsy output
;The parameter imf is either 0 (kroupa) or 1 (Miller Scalo)

maxtime = 3.e9
tbinsize=5.e6
ntbins = maxtime/tbinsize
timeunit=1.e9
massunit=2.325e5 ; massunit in terms of solar masses
tbins = findgen(ntbins)*tbinsize ; time bins starting at 0
tbins_next =  findgen(ntbins)*tbinsize + tbinsize ; time bins starting at dtime
snbins = fltarr(ntbins)*tbinsize ;Histogram Number of SN will go in

s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/'
file = base+filename
rtipsy,file,h,g,d,s
nstars = N_ELEMENTS(s)
s.tform=s.tform*timeunit
s.mass=s.mass*massunit
print,'Number of Stars: ',nstars
;Read in data and change units

FOR istar = 0L, nstars-1 DO BEGIN
;for each star particle
    ind = WHERE(tbins GE s[istar].tform)
    ;The bins with time later than formation time
    IF (s[istar].metals eq 0.0) THEN z = 0.0004 ELSE z = s[istar].metals
    ;No zero metallicity stars allowed
    IF(ind[0] ne -1) THEN snbins[ind] = snbins[ind] + s[istar].mass*(nsnova(tbins_next[ind] - s[istar].tform,z,imf))
    ;For those times greater than the
    ;formation time, add in the sn from that particle
    IF(istar MOD 1000 EQ 0) THEN BEGIN 
        stars = istar + 1
        plot,tbins/1e9,snbins,xtitle = 'time',ytitle = 'Mass of Supernovas (per million years)',title = strtrim(stars)+'Stars'
    endif
ENDFOR

plot,tbins/1e9,snbins,xtitle = 'time',ytitle = 'Mass of Supernovas (per million years)'
return,snbins
END
 
PRO comp_IMF_sn
loadct,39
set_plot,'x'
maxtime = 3.e9
tbinsize=5.e6
ntbins = maxtime/tbinsize
timeunit=1.e6
tbins = findgen(ntbins)*tbinsize ; time bins starting at 0

;file1 = '1E5R/10M_k/o10M_1.00300'
;file2 = '1E5R/10M_ms/o10M_1.00300'
file1 = '1E5R/12M_k/o12M_1.00300'
file2 = '1E5R/12M_ms/o12M_1.00300'
sn_k = calc_sn(file1,0)
sn_k_other = calc_sn(file1,1)
sn_ms = calc_sn(file2,1)
sn_ms_other = calc_sn(file2,0)
maxval = MAX([MAX(sn_k),MAX(sn_ms)])

plot,tbins/1e9,sn_ms,xtitle = 'time',ytitle = 'Number of Supernovas (per million years)',title='10^12 Solar Mass Halo',yrange = [0,maxval]
oplot,tbins/1e9,sn_k, color = 50
oplot,tbins/1e9,sn_ms, color = 240
;oplot,tbins/1e9,sn_k_other, color = 50,linestyle = 1
;oplot,tbins/1e9,sn_ms_other, color = 240,linestyle = 1
legend,['Miller-Scalo'+STRTRIM(TOTAL(sn_ms)),'Kroupa'+STRTRIM(TOTAL(sn_k))],linestyle = [0,0],color = [240,50],/right
;legend,['Miller-Scalo'+STRTRIM(TOTAL(sn_ms)),'Kroupa'+STRTRIM(TOTAL(sn_k)),'Miller-Scalo Galaxy, Kroupa Feedback'+STRTRIM(TOTAL(sn_ms_other)),'Kroupa Galaxy, MS Feedback'+STRTRIM(TOTAL(sn_k_other))],linestyle = [0,0,1,1],color = [240,50,240,50],/right

set_plot,'ps'
device,filename = '../results/compIMF_sn_12M.eps',/color,bits_per_pixel=8
plot,tbins/1e9,sn_ms,xtitle = 'time',ytitle = 'Number of Supernovas (per 5 million years)',title='10^12 Solar Mass Halo',yrange = [0,maxval]
oplot,tbins/1e9,sn_k, color = 50
oplot,tbins/1e9,sn_ms, color = 240
;oplot,tbins/1e9,sn_k_other, color = 50,linestyle = 1
;oplot,tbins/1e9,sn_ms_other, color = 240,linestyle = 1
legend,['Miller-Scalo'+STRTRIM(TOTAL(sn_ms)),'Kroupa'+STRTRIM(TOTAL(sn_k))],linestyle = [0,0],color = [240,50],/right
;legend,['Miller-Scalo'+STRTRIM(TOTAL(sn_ms)),'Kroupa'+STRTRIM(TOTAL(sn_k)),'Miller-Scalo Galaxy, Kroupa Feedback'+STRTRIM(TOTAL(sn_ms_other)),'Kroupa Galaxy, MS Feedback'+STRTRIM(TOTAL(sn_k_other))],linestyle = [0,0,1,1],color = [240,50,240,50],/right
device,/close
END

PRO plot_IMF
set_plot,'x'
max_mass = 100
min_mass = 0.08
dmass = 0.02
nmass = (max_mass - min_mass)/dmass + 1
masses = findgen(nmass)*dmass + min_mass
n_stars_k = kroupa(masses,/nstar)
n_stars_ms = ms(masses,/nstar)
m_stars_k = kroupa(masses)
m_stars_ms = ms(masses)
plot,masses,m_stars_ms,linestyle = 0,ytitle = 'Mass of stars formed at a given mass',xtitle = 'Mass',title = 'Comparing IMFs',/ylog,/xlog
oplot,masses,m_stars_k,linestyle = 1
legend,['Miller-Scalo','Kroupa'],linestyle = [0,1],/right

set_plot,'ps'
device,filename = '../results/plotIMF_m.eps'
plot,masses,m_stars_ms,linestyle = 0,ytitle = 'Mass of stars formed at a given mass',xtitle = 'Mass',title = 'Comparing IMFs',/ylog,/xlog
oplot,masses,m_stars_k,linestyle = 1
legend,['Miller-Scalo','Kroupa'],linestyle = [0,1],/right
device,/close
stop

set_plot,'x'
plot,masses,n_stars_ms,linestyle = 0,ytitle = 'Mass of stars formed at a given mass',xtitle = 'Mass',title = 'Compareing IMFs',/ylog,/xlog
oplot,masses,n_stars_k,linestyle = 1
legend,['Miller-Scalo','Kroupa'],linestyle = [0,1],/right

set_plot,'ps'
device,filename = '../results/plotIMF_n.eps'
plot,masses,n_stars_ms,linestyle = 0,ytitle = 'Number of stars formed at a given mass',xtitle = 'Mass',title = 'Compareing IMFs',/ylog,/xlog
oplot,masses,n_stars_k,linestyle = 1
legend,['Miller-Scalo','Kroupa'],linestyle = [0,1],/right
device,/close

END
