; Written by Alyson Brooks, updated Dec 7, 2011
;
; This version takes the input file and finds *any* resolved halos (non-satellites)
; where resolved is defined as N_dark > 3500 particles
;
; If z < 0.1 (h.time gt 0.9), this code finds where 2/3 of the stellar radius 
; is and gets the gas phase O abundance of the gas in that radius
; This is to mimic the SDSS fibers.
;
; If z > 0.1, all cold gas is included

;metals = mzr('h516.cosmo25cmb.3072g14HBWK.00024',/obs)
function mzr, filename, obs = obs, onehalo = onehalo
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
grav = 6.6725985e-8 ;cm^3/g/s
s_per_yr = 3.15569e7
filebase = filename
IF keyword_set(onehalo) THEN filename = filename + '.halo.' + strtrim(onehalo)
rheader, filename, h
a = h.time
;Units
;name = strmid(filename,26,7, /reverse_offset)
;if name eq 'cosmo25' or name eq '.cosmo2' or name eq 'osmo25c' then begin
IF strpos(filename,'cosmo25') NE -1 THEN BEGIN
  boxsize = 18250.
  massunit=2.310e15
  vunit = 630.28
endif
;if name eq 'cosmo50' or name eq 'osmo50c' or name eq 'smo50cm' or name eq '50cmb.2' then begin
IF strpos(filename,'cosmo50') NE -1 THEN BEGIN
  boxsize = 36500.
  massunit=1.84793e16
  vunit = 1261.05
endif
h0=0.73
lengthunit = boxsize/h0
timeunit=SQRT((lengthunit*3.086d21)^3/(6.67d-8*massunit*1.99d33))/(3600.*24.*365.24)
lengthunit=lengthunit*a
rhounit=massunit * gm_per_msol * H_per_gm/lengthunit^3/cm_per_kpc^3

;Filenames
grpfile = filename+'.amiga.grp'
gtpfile = filename+'.amiga.gtp'
fefile = filename+'.FeMassFrac'
oxfile = filename+'.OxMassFrac'
statfile = filename+'.amiga.stat'
HIfile = filename+'.HI'
HeIfile = filename+'.HeI'
HeIIfile = filename+'.HeII'
H2file = filename+'.H2'
doH2 = file_test(H2file)

;Read in files
rtipsy, filename, h, g, d, s
IF NOT keyword_set(onehalo) THEN BEGIN
    rtipsy, gtpfile, h2, g2, d2, halo
    grp = read_lon_array(grpfile)
    grpallstar = grp[h.ngas+h.ndark:h.n-1]
    grp = grp[0:h.ngas-1]
ENDIF ELSE BEGIN
    grp = fltarr(h.ngas) + onehalo
    grpallstar = fltarr(h.nstar) + onehalo
ENDELSE
;read_tipsy_arr,fefile,h,fe,part = 'gas'
;read_tipsy_arr,oxfile,h,ox,part = 'gas'
read_tipsy_arr,fefile,h,fe,type = 'float'
festar = fe[h.ngas+h.ndark:h.n-1]
fe = fe[0:h.ngas-1]
read_tipsy_arr,oxfile,h,ox,type = 'float'
oxstar = ox[h.ngas+h.ndark:h.n-1]
ox = ox[0:h.ngas-1]
read_tipsy_arr,HIfile,h,HI,part = 'gas'
;readarr, fefile,   h, fe, part = 'gas',/ascii
;readarr, oxfile,   h, ox, part = 'gas',/ascii
;readarr, HIfile,   h, HI, part = 'gas',/ascii
;readarr, HeIfile,  h, HeI, part = 'gas',/ascii
;readarr, HeIIfile, h, HeII, part = 'gas',/ascii
IF doH2 THEN BEGIN
    read_tipsy_arr,H2file,h,H2,part = 'gas'
;    readarr, H2file, h, H2, part = 'gas',/ascii
    IF (max(H2) EQ 0) THEN doH2 = 0
ENDIF ELSE H2 = HI*0
IF file_test(HeIfile) THEN read_tipsy_arr,HeIfile,h,HeI,part = 'gas' ELSE HeI = (HI + 2.0*H2)*0.4/4
IF file_test(HeIIfile) THEN read_tipsy_arr,HeIIfile,h,HeII,part = 'gas' ELSE HeII = fltarr(n_elements(HI))
wrongHeI = where(1 - HI - 2.0*H2 - 4.0*(HeI - HeII) - g.zmetal LT 0)
HeI[wrongHeI] = (1 - HI[wrongHeI] - 2.0*H2[wrongHeI] - 4.0*HeII[wrongHeI] - g[wrongHeI].zmetal)/4.0
IF NOT keyword_set(onehalo) THEN stat = read_stat_struc_amiga(statfile)

;Select the cold/hot gas
;For the cold/hot gas, create an array of Ox and Fe metallicities and HI mass

;Weight the gas by the star formation efficiency
IF doH2 THEN BEGIN
    cstar = 0.1
    tempcut = 1e3
    denscut = 0.1
    deltat = 1e6*s_per_yr ;yr
    sfeff = fltarr(h.ngas)
    indsf = where(g.tempg LE tempcut AND g.dens*rhounit GE denscut)
    tdyn = 1.0/sqrt(4.0*!PI*grav*g.dens*rhounit*gm_per_H)
    IF indsf[0] NE -1 THEN sfeff[indsf] = 1.0 - exp(-1.0*cstar*deltat/tdyn[indsf]*2.0*H2[indsf]/(2.0*H2[indsf] + HI[indsf]))
ENDIF ELSE BEGIN
    cstar = 0.1
    tempcut = 1e4
;    denscut = 100.0
    denscut = 10
    deltat = 1e6*s_per_yr ;yr
    sfeff = fltarr(h.ngas)
    indsf = where(g.tempg LE tempcut AND g.dens*rhounit GE denscut)
    tdyn = 1.0/sqrt(4.0*!PI*grav*g.dens*rhounit*gm_per_H)
    IF indsf[0] NE -1 THEN sfeff[indsf] = 1.0 - exp(-1.0*cstar*deltat/tdyn[indsf])
ENDELSE

;All gas Ox
allhe_mass = g.mass*HeI*massunit + g.mass*HeII*massunit
allh2_plus_h_mass = (1-g.zmetal)*g.mass*massunit
allh_mass = allh2_plus_h_mass - (4*allhe_mass)
allnh = (allh_mass/1.6726d-24)*(2.0d33)
alln_ox = (ox*g.mass*massunit)*2.0d33/2.66d-23
all_ox_abund = alln_ox/allnh
allyeff = ox   

;Cold gas Ox
coolgas = where(g.tempg LE tempcut and grp gt 0) ;Cold Gas
grpcoolgas = grp[coolgas]
he_mass = g[coolgas].mass*HeI[coolgas]*massunit + g[coolgas].mass*HeII[coolgas]*massunit
he_plus_h_mass = (1-g[coolgas].zmetal)*g[coolgas].mass*massunit
h_mass = he_plus_h_mass - (4*he_mass)
nh = (h_mass/1.6726d-24)*(2.0d33)
n_ox = (ox[coolgas]*g[coolgas].mass*massunit)*2.0d33/2.66d-23
cool_ox_abund = n_ox/nh
yeff = ox[coolgas]                                  

;Hot gas Ox
hotgas = where(g.tempg GT tempcut and grp gt 0)  ;Hot gas
grphotgas = grp[hotgas]
he_mass = g[hotgas].mass* HeI[hotgas]* massunit + g[hotgas].mass* HeII[hotgas]* massunit
he_plus_h_mass = (1- g[hotgas].zmetal)* g[hotgas].mass*massunit
h_mass = he_plus_h_mass - (4*he_mass)
nh = (h_mass/1.6726d-24)*(2.0d33)
n_ox = (ox[hotgas]*g[hotgas].mass*massunit)*2.0d33/2.66d-23
hot_ox_abund = n_ox/nh

;Atomic gas Ox
;atomicgas = where(grp gt 0)
;grpatomicgas = grp[atomicgas]
;atomicgas_mass = (HI[atomicgas] + 2.0*H2[atomicgas])/max(HI + 2.0*H2)*g[atomicgas].mass ;Weight the gas mass by the atomic fraction
;h_mass = (HI[atomicgas] + 2.0*H2[atomicgas])*g[atomicgas].mass*massunit
;nh = (h_mass/1.6726d-24)*(2.0d33)
;n_ox = (ox[atomicgas]*atomicgas_mass*massunit)*2.0d33/2.66d-23
;atomic_ox_abund = n_ox/nh

;Dense gas
densegas = where(g.dens*rhounit GT 0.1 and grp gt 0)  ;dense gas
grpdensegas = grp[densegas]
g[densegas].mass = g[densegas].mass                        ;Mass of hot gas
he_mass = g[densegas].mass*HeI[densegas]*massunit + g[densegas].mass*HeII[densegas]*massunit
he_plus_h_mass = (1-g[densegas].zmetal)*g[densegas].mass*massunit
h_mass = he_plus_h_mass - (4*he_mass)
nh = (h_mass/1.6726d-24)*(2.0d33)
n_ox = (ox[densegas]*g[densegas].mass*massunit)*2.0d33/2.66d-23
dense_ox_abund = n_ox/nh

;Note original 2007 paper didn't exclude satellites
;resolved = where(stat.npart GE 3500 and stat.sat ne 'yes')
IF NOT keyword_set(onehalo) THEN BEGIN
    resolved = where(stat.npart GE 3500)
    unique = stat[resolved].group
ENDIF ELSE BEGIN
    resolved = [onehalo]
    unique = [onehalo]
ENDELSE

metals = replicate({maxdist:dblarr(1), mvir:dblarr(1), grp:dblarr(1), ntot:dblarr(1), ngas:dblarr(1), nstar:dblarr(1), nhot:lonarr(1), ncold:lonarr(1), inner_gasmass:dblarr(1), cool_inner_gasmass:dblarr(1), coolgasmass:dblarr(1), hotgasmass:dblarr(1), massHI:dblarr(1), ninner:lonarr(1), ox_inner:dblarr(1), fe_inner:dblarr(1), ox_innerAtom:dblarr(1), fe_innerAtom:dblarr(1), ox_cold:dblarr(1), fe_cold:dblarr(1), ox_hot:dblarr(1), fe_hot:dblarr(1), ox_sfr:dblarr(1), fe_sfr:dblarr(1), ox_stellar:dblarr(1), fe_stellar:dblarr(1), ox_stellar_rband:dblarr(1), fe_stellar_rband:dblarr(1), smass:dblarr(1), gmass:dblarr(1), yeff:dblarr(1), yeff_inner:dblarr(1), smassO:dblarr(1), sat:lonarr(1)},n_elements(unique))

;Calculate the gas masses, metallicities etc for each of the satellities
FOR i=0L,n_elements(unique)-1 DO BEGIN
  grpind = unique[i]
  ind = where(grpcoolgas EQ grpind, nind)
  ind_star = where(grpallstar EQ grpind, nind_star)

  IF keyword_set(obs) THEN BEGIN
;      IF 0 THEN BEGIN ;
      IF 0 THEN BEGIN ;file_test(filebase + '.' + strtrim(grpind,2) + '/broadband.fits') THEN BEGIN
          lums = mrdfits(filebase + '.' + strtrim(grpind,2) + '/broadband.fits', 13, header,/silent)
          IF n_elements(lums) LE 1 THEN lums = mrdfits(filebase + '.' + strtrim(grpind,2) + '/broadband.fits', 14, header,/silent)
          Bn = where(strcmp(strtrim(lums.filter,2),'B_Johnson.res') EQ 1)
          Bmag = lums[Bn].AB_mag1
          Kn = where(strcmp(strtrim(lums.filter,2),'K_Johnson.res') EQ 1)
          Kmag = lums[Kn].AB_mag1
          Vn = where(STRCMP(strtrim(lums.filter,2),'V_Johnson.res') eq 1)
          Vmag = lums[Vn].AB_mag1
     ENDIF ;ELSE BEGIN
          mags = mrdfits(filebase + '.amiga_virbk.halos.star99_K_ab.Mv.fits',1,/silent)
          ind2grp = where(mags.id EQ grpind)
          magB = mags[ind2grp].B
          magK = mags[ind2grp].K
          mags = mags[ind2grp]
          mags_bes = transpose([[mags.u],[mags.b],[mags.v],[mags.r],[mags.i]])
          mags_errs = fltarr(5,n_elements(mags.u))+0.02
          mags_ivar=1./mags_errs^2
          dmod = 19.4576        ;distance modulus for z = 0 in this model
          z = fltarr(n_elements(mags.u))
          mgy =(10.D)^(-(0.4D)*(mags_bes + dmod))
          mgy_ivar = mags_ivar/(0.4*alog(10.)*mgy)^2.
;      ENDELSE
  ENDIF
  mags_stars = mrdfits(filebase + '.stars.star99_K_ab.Mv.fits',1,/silent)
;  IF file_test(filebase + '.' + strtrim(grpind,2) + '/broadband.fits') THEN BEGIN
;      print,'K (sunrise/idl): ',Kmag,magK
;      print,'B (sunrise/idl): ',Bmag,magB
;  ENDIF
  ;Fill in stat file values
;IF 0 THEN BEGIN
  IF NOT keyword_set(onehalo) THEN BEGIN
      ind2grp = where(stat.group EQ grpind)
      IF stat[ind2grp].sat EQ 'yes' THEN metals[i].sat = 1
      metals[i].smass = stat[ind2grp].m_star
      metals[i].gmass = stat[ind2grp].m_gas
      metals[i].mvir = stat[ind2grp].m_tot
      metals[i].grp = stat[ind2grp].group
      metals[i].ntot = stat[ind2grp].npart
      metals[i].ngas = stat[ind2grp].ngas 
      metals[i].nstar = stat[ind2grp].nstar
   ENDIF ELSE BEGIN
      metals[i].smass = total(s.mass)*massunit
      metals[i].gmass = total(g.mass)*massunit
      metals[i].mvir = (total(s.mass) + total(g.mass) + total(d.mass))*massunit
      metals[i].grp = onehalo
      metals[i].ntot = n_elements(s) + n_elements(g) + n_elements(d)
      metals[i].ngas = n_elements(g)
      metals[i].nstar = n_elements(s)
   ENDELSE

  IF ind[0] NE -1 THEN BEGIN
      metals[i].ncold = n_elements(ind)
      metals[i].ox_cold = 12.+alog10(mean(cool_ox_abund[ind]))
      metals[i].fe_cold = mean(fe[coolgas(ind)])
      metals[i].coolgasmass = total(g[coolgas(ind)].mass)*massunit
      gasfrac = metals[i].coolgasmass/(metals[i].smass+metals[i].coolgasmass)
      metals[i].yeff = alog10(mean(yeff[ind])/alog(1./gasfrac))
  ENDIF

  indhot = where(grphotgas eq grpind)
  metals[i].nhot = n_elements(indhot)
  IF indhot[0] NE -1 THEN BEGIN
      metals[i].ox_hot = 12.+alog10(mean(hot_ox_abund[indhot]))
      metals[i].fe_hot = mean(fe[hotgas(indhot)])
      metals[i].hotgasmass = total(g[hotgas[indhot]].mass)*massunit
  ENDIF

  metals[i].ox_stellar = total(oxstar[ind_star]*s[ind_star].mass)/total(s[ind_star].mass)
  metals[i].fe_stellar = total(festar[ind_star]*s[ind_star].mass)/total(s[ind_star].mass)
  metals[i].ox_stellar_rband = total(oxstar[ind_star]*10^(mags_stars[ind_star].r/(-2.5)))/total(10^(mags_stars[ind_star].r/(-2.5)))
  metals[i].fe_stellar_rband = total(festar[ind_star]*10^(mags_stars[ind_star].r/(-2.5)))/total(10^(mags_stars[ind_star].r/(-2.5)))

  IF keyword_set(obs) THEN BEGIN
      IF 0 THEN BEGIN
;------- get stellar masses from K and B band luminosities using the
;        output from get_halo_mags
;        Follows Lee 06
          Bmag = magB
          Kmag = magK
          a_prime = -0.851
          b_prime =  0.254
          mk = 3.28
          ml = a_prime + b_prime*(Bmag - Kmag)
          Metals[i].smassO = 10^(ml + 0.4*(mk - Kmag))
      ENDIF ELSE BEGIN
          IF 0 THEN BEGIN ;file_test(filebase + '.' + strtrim(grpind,2) + '/broadband.fits') THEN BEGIN
              L_B_sol = 4.67e25 ;Watts, from B+M, p 53, Table 2.1
              MLa = -0.942
              MLb = 1.69
              Bmag = Bmag + 0.09                ; to Vega
              Vmag = Vmag + 0.03                ; to Vega
              L_B_band = lums[Bn].L_lambda_eff1 ;Watts/m
              B_band_width = lums[Bn].ewidth_lambda ;m
              L_B = L_B_band*B_band_width           ;Watts
              B_V = Bmag - Vmag
              MLfit = 10^(MLa + MLb*B_V)
              Metals[i].smassO = MLfit*L_B/L_B_sol
          ENDIF ELSE BEGIN
              kcorrect, mgy, mgy_ivar, z, kcor, mass = star, mtol = mtol_k, absmag = absmag_k, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']
              Metals[i].smassO = star
;              print,Metals[i].smassO,star,
          ENDELSE
      ENDELSE
  ENDIF ELSE Metals[i].smassO = 0
;  kcorrect, mgy, mgy_ivar, z, kcor, mass = star, mtol = mtol_k, absmag = absmag_k, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']
  print,'Stellar Mass (actual/observed/kcorrect): ',strtrim(Metals[i].smass,2),strtrim(Metals[i].smassO,2),strtrim(star,2)
  IF abs(alog10(Metals[i].smass) - alog10(Metals[i].smassO)) GT 3 AND Metals[i].smass NE 0 THEN stop
  indall = where(grp EQ grpind, nind)
  IF nind EQ 0 THEN CONTINUE
  metals[i].ox_sfr = 12.+alog10(total(all_ox_abund[indall]*sfeff[indall])/total(sfeff[indall]))
  metals[i].fe_sfr = total(fe[indall]*sfeff[indall])/total(sfeff[indall])
  metals[i].massHI = total(HI(where(grp EQ grpind))*g(where(grp EQ grpind)).mass)*massunit
ENDFOR

;if h.time gt 0.9 then begin 
;Create an array that for each gas particle in the simulation contains 
;its distance from its parent halo
;Create an array for each star particle in the simulation that contains 
;its distance from its parent halo

inner_gas = -1
sorted = grp[sort(grp)]   ;sort the gas grp array 
uniqueg = sorted[uniq(sorted)]  ;find the unique grp values for the gas
sorted = grpallstar[sort(grpallstar)] ;sort the start grp array
uniques = sorted[uniq(sorted)]  ; find the unique grp values for the stars
IF keyword_set(onehalo) THEN ibaryons = 1 ELSE ibaryons = binfind(uniqueg, uniques) ;searches uniqueg for  uniques 
wbaryons = uniqueg[ibaryons[where(ibaryons ne -1)]] ;finds all the unique halos with gas and stars
IF NOT keyword_set(onehalo) THEN BEGIN
    ibaryons = binfind(stat[resolved].group, wbaryons) ;find only those unique halos that are resolved
    wbaryons = wbaryons[where(ibaryons ne -1)]         ;removes extraneous halos
ENDIF
innerrad = fltarr(n_elements(wbaryons)) ;an array for each of the halos that can contain the inner radius

for i=0,n_elements(wbaryons)-1 do begin
;Loop through halos
  inds = where(grpallstar eq wbaryons[i], ninds) ;find the stars that belong to that halo 
  IF keyword_set(onehalo) THEN stardist = sqrt(s[inds].x^2 + s[inds].y^2 + s[inds].z^2) $
  ELSE BEGIN
      indh = where(stat.group eq wbaryons[i]) ;find the stat id that corresponds to it
      stardist = sqrt((halo[indh].x-s[inds].x)^2.+(halo[indh].y-s[inds].y)^2.+(halo[indh].z-s[inds].z)^2.) ;find the center of the stars from the center of the halo
  ENDELSE
  k=0.
  test = 0 
;   To math SDSS, Consider the gas that is within the radius that contains 2/3 of the total stellar mass
  while test LE (2./3.)*total(s[inds].mass) do begin
	k=k+1. 
	inner = where(stardist*lengthunit LE k*0.1, ninner)
	if ninner ne 0 then test = total(s[inds(inner)].mass) else continue
  endwhile
  innerrad[i] = k*0.1

;Alternatively, find the radius that contains 95% of the recently
;(last 100 Myr) formed stars
  percent = 0.95
  deltat = 100*1e6
  tcurrent = max(s[inds].tform)*timeunit
  IF (tcurrent)[0] NE -1 THEN BEGIN
     newstarind = where(s[inds].tform*timeunit GE tcurrent - deltat)
     stardistsort = stardist[newstarind[sort(stardist[newstarind])]]
;     print,innerrad[i],stardistsort[round(n_elements(stardistsort)*percent) - 1]*lengthunit
     innerrad[i] = stardistsort[round(n_elements(stardistsort)*percent) - 1]*lengthunit
  ENDIF

  indg = where(grp eq wbaryons[i], nindg) ;Find the gas that belongs to the halo
  IF keyword_set(onehalo) THEN distance = sqrt(g[indg].x^2 + g[indg].y^2 + g[indg].z^2) $
  ELSE distance = sqrt((halo[indh].x-g[indg].x)^2.+(halo[indh].y-g[indg].y)^2.+(halo[indh].z-g[indg].z)^2.) ;Find the distance of that gas from the center of the halo
;  ind2g = where(g[indg].tempg LE tempcut and distance LE innerrad[i]/lengthunit, nind2g)
  ind2g = where(distance LE innerrad[i]/lengthunit, nind2g) 
  if nind2g ne 0 then inner_gas = [inner_gas, indg[ind2g]]
endfor
inner_gas = inner_gas[where(inner_gas GE 0)] ;remove negetive one from start of array

;Cool inner gas
cool_inner_gas = where(g[inner_gas].tempg le tempcut)
grp_cool_inner_gas = grp[inner_gas[cool_inner_gas]]
he_mass = 4.*(g[inner_gas[cool_inner_gas]].mass*HeI[inner_gas[cool_inner_gas]]*massunit + g[inner_gas[cool_inner_gas]].mass*HeII[inner_gas[cool_inner_gas]]*massunit)
he_plus_h_mass = (1-g[inner_gas[cool_inner_gas]].zmetal)*g[inner_gas[cool_inner_gas]].mass*massunit
h_mass = he_plus_h_mass - he_mass
nh = (h_mass/1.6726d-24)*(2.0d33)
n_ox = (ox[inner_gas[cool_inner_gas]]*g[inner_gas[cool_inner_gas]].mass*massunit)*2.0d33/2.66d-23
ox_abund = n_ox/nh 
yeff = ox[inner_gas[cool_inner_gas]]

;Atomic inner gas
;atomicgas_mass = g[inner_gas].mass
atomicgas_mass = (HI[inner_gas] + 2.0*H2[inner_gas])/max(HI + 2.0*H2)*g[inner_gas].mass ;Weight the gas mass by the atomic fraction
grpatomicgas_inner = grp[inner_gas]
h_mass = (HI[inner_gas] + 2.0*H2[inner_gas])*g[inner_gas].mass*massunit
nh = (h_mass/1.6726d-24)*(2.0d33)
n_ox = (ox[inner_gas]*atomicgas_mass*massunit)*2.0d33/2.66d-23
atomic_ox_abund = n_ox/nh

;grpgas_inner = 
;sorted = grpgas_inner(sort(grpgas_inner))
;unique = sorted(uniq(sorted))

;sorted = grp_cool_inner_gas(sort(grp_cool_inner_gas))
;unique = sorted(uniq(sorted))
;for i=0L,n_elements(unique)-1 do begin
;endfor

for i=0L,n_elements(wbaryons)-1 do begin
  grpind = where(metals.grp eq wbaryons[i], ngrpind)
  IF ngrpind EQ 0 THEN CONTINUE

  ind = where(grp_cool_inner_gas eq wbaryons[i], nind)
  IF nind NE 0 THEN BEGIN
      metals[grpind].maxdist = innerrad[i]
      metals[grpind].ox_inner = 12.+alog10(mean(ox_abund[ind]))
      metals[grpind].ninner = n_elements(ind) 
      metals[grpind].fe_inner = mean(fe[inner_gas[cool_inner_gas[ind]]]) 
      metals[grpind].cool_inner_gasmass = total(g[inner_gas[cool_inner_gas[ind]]].mass)*massunit 
      gasfrac = metals[grpind].cool_inner_gasmass/(metals[grpind].smass+metals[grpind].cool_inner_gasmass)
      metals[grpind].yeff_inner = alog10(mean(yeff[ind])/alog(1./gasfrac))
  ENDIF

  ind = where(grpatomicgas_inner eq wbaryons[i], nind)
  IF nind NE 0 THEN BEGIN
      metals[grpind].ox_innerAtom = 12.+alog10(mean(atomic_ox_abund[ind]))
      metals[grpind].fe_innerAtom = total(fe[inner_gas[ind]]*atomicgas_mass[ind]*massunit)/total(atomicgas_mass[ind]*massunit) 
      metals[grpind].inner_gasmass = total(atomicgas_mass[ind])*massunit 
  ENDIF
ENDFOR
return, metals
  
END

PRO mzr_calc,filebase 
spawn,'pwd',home
home = '/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/'
home = '/home/christensen/Storage2/UW/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/'
home = '/astro/store/nbody3/christensen/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/'
home = '/home/christensen/Storage1/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/'
spawn,'ls -d ' + home + filebase + '*.dir',dirnames
n = n_elements(dirnames)

FOR i = 1, n - 1 DO BEGIN
    cd,dirnames[i]
    startpos = strpos(dirnames[i],'/',/reverse_search) 
    IF startpos NE -1 THEN filename = strmid(dirnames[i],startpos + 1,strlen(dirnames[i]) - startpos - 1 - 4) ELSE filename = strmid(dirnames[i],0,strlen(dirnames[i]) - 4)
    print,dirnames[i],', ',filename
    IF file_test(filename + '.amiga.stat') THEN BEGIN
        IF NOT file_test(filename + '.metals.fits') THEN BEGIN
;            IF file_test(filename + '.amiga_vir.halos.star99_K_ab.Mv.fits') THEN metals = mzr(filename,/obs) ELSE 
            metals = mzr(filename) 
            mwrfits,metals,filename + '.metals.fits',/create
        ENDIF
    ENDIF
    cd,home
ENDFOR

IF (keyword_set(outfile)) THEN device,/close ;ELSE stop
END
