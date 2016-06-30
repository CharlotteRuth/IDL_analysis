;
;
;
;  generate profiles for the given set of particles
;
;  The types of profiles generated depends on the kind of particles
;  you input. The procedure returns a structure containing all the
;  profiles. 
;
;
;  For gas, the profiles are:  
;  density (rho), 
;  angular momentum (j),
;  z-component of ang. mom (j_z),
;  mass-weighted density (norm_rho),
;  standard deviation of mass-weighted density (stddev_rho),
;  temperature (temp),
;  mass enclosed (mass_enc)
;  metal surface density (metal_den)
;
;
;  For stars, the profiles are:
;  density (rho),
;  mean age per bin (age),
;  star formation rate (sfr),
;  angular momentum (j),
;  z-component of ang. mom (j_z),
;  mass enclosed (mass_enc)
;
;
;  For dark matter, the profiles are:
;  density (rho),
;  angular momentum (j)
;  z-component of ang. mom (j_z),
;  mass encolsed (mass_enc)
;
;
;  OPTIONAL KEYWORDS:
;
;  nbins - number of bins
;  rmax  - maximum radius
;  rmin  - minimum radius
;  rbins - return radial bin array
;  num_in_bin - return number of particles per bin
;  sph   - spherical bins; default is cylindrical bins
;
;
;  USAGE:
;
;  prof = profile(g, 'gas', nbins = 50, rbins = rbins)
;  plot, rbins, prof.norm_rho, /ylog
;
; ***********************************************************************

function prof, p, part, time, NBINS = nbins, RMAX = rmax, RMIN = rmin, SPH = sph


if keyword_set(sph) then $
  r = sqrt(p.x^2+p.y^2+p.z^2) $
else r = sqrt(p.x^2 + p.y^2)

if keyword_set(rmax) eq 0 then begin
    rmax = max(r)
    print, 'warning no rmax specified, using max(r)'
endif

if keyword_set(rmin) eq 0 then begin
    rmin = 0.0
    print, 'warning no rmin specified, using 0.0'
endif

if keyword_set(nbins) eq 0 then begin 
    nbins = 100.
    print, 'warning no nbins specified, using 100'
endif

; ********************************* 
; arrays holding positions of bins
; ********************************* 

dr = (rmax - rmin)/nbins
rbins = findgen(nbins)*dr + rmin
rcentbins = rbins + dr/2.

; **********************************
; initialize the profile structures
; **********************************

if part eq 'gas' then begin
    prof = {rbins:rbins, rho:fltarr(nbins), j:fltarr(nbins), j_z:fltarr(nbins),$
            norm_rho:fltarr(nbins), stddev_rho:fltarr(nbins), temp:fltarr(nbins), $
            mass_enc:fltarr(nbins), num_in_bin: fltarr(nbins), metal_den:fltarr(nbins) } 
endif 

if part eq 'star' then begin
    prof = {rbins:rbins, rho:fltarr(nbins), j:fltarr(nbins), j_z:fltarr(nbins), $
            mass_enc:fltarr(nbins), age:fltarr(nbins), sfr:fltarr(nbins), $
            num_in_bin:fltarr(nbins)}
endif 

if part eq 'dark' then begin
    prof = {rbins:rbins, rho:fltarr(nbins), j:fltarr(nbins), j_z:fltarr(nbins), $
            mass_enc:fltarr(nbins), num_in_bin:fltarr(nbins)}
endif

; ******************************************************
; need angular momentum regardless of the particle type
; ******************************************************

angmom, p, jvec, lvec, jtot, ltot

for i = 0L, nbins-1 do begin
    ind = where((r ge rbins[i] and r le rbins[i]+dr), count)
    
    prof.num_in_bin[i] = count

; *****************************************************
; generate the profiles that are used by all particles
; *****************************************************

binmass = fltarr(nbins)

if count ne 0 then begin

    binmass[i] = total(p[ind].mass)

    if keyword_set(sph) eq 0 then $  
      binarea = !pi*((rbins[i]+dr)^2.-rbins[i]^2.) $    
    else $
      binarea = 4./3.*!pi*((rbins[i]+dr)^3.-rbins[i]^3.)

; ****
; rho
; ****
    prof.rho[i] = binmass[i]/binarea
 
; ********
; J & J_z
; ********
    prof.j[i] = total(jtot[ind])
    prof.j_z[i] = total(jvec[ind,2])
    
; **************
; mass enclosed
; **************

    if i gt 0 then $
      prof.mass_enc = total(binmass[indgen(i)]) $
    else $
      prof.mass_enc = binmass[0]
        

; *************
; gas profiles
; *************

    if part eq 'gas' then begin
        
; *********
; norm_rho
; *********
        rho_part = p[ind].mass*p[ind].dens
        prof.norm_rho[i] = total(rho_part)/binmass[i]


; ***********
; stddev_rho
; ***********
        if count gt 1 then $
          stats = moment(rho_part, sdev = prof.stddev_rho[i]) $
        else $
          prof.stddev_rho[i] = 0.0

; *****
; temp
; *****
        prof.temp[i] = total(p[ind].mass*p[ind].tempg)/binmass[i]

; *****
; metal surface density
; *****
        prof.metal_den[i] = total(p[ind].mass*p[ind].zmetal)/binarea
    endif


    if part eq 'star' then begin
        
; ****
; age
; ****
        prof.age[i] = total(p[ind].tform*p[ind].mass)/binmass[i]

; ****
; sfr
; ****
        new = where(p[ind].tform gt time - 0.01, nsfr)
        new = where(p[ind].tform gt time - 0.045, nsfr) ;Slightly different definition of SFR, based on that from Saitoh 2008, figure 7
        if (nsfr gt 0) then $
          prof.sfr[i] = total(p[new].mass)/binarea/0.01/1e9 $
        else $
          prof.sfr[i] = 0.0

    endif
endif

endfor
return, prof

end

