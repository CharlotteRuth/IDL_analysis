function isohalomag, file, mass_unit=mass_unit, timeunit=timeunit

; NAME:
;   getisohalomag
; PURPOSE:
;   Get the absolute magnitudes of each star in a simulation output
;   and total for galaxy
; CALLING SEQUENCE:
;  IDL> getisohalomag, mass_unit
;
;   This must be done *in* the directory containing the tipsy
;   files and the SKID output files you want to analyze
;
;   The code only connects  SKID and TIPSY files that have the same
;   prefix (i.e. tipsyfile.grp) or have only a g added
;   (i.e. tipsyfileg.grp)
;
; INPUTS:
;   mass_unit - mass unit of simulation
; OUTPUTS:
;   tipsyfile.halo.Mv.fits
;   tipsyfile.stars.Mv.fits
;   
;   For each tipsyfile with SKID output the code identifies
;   in the working directory
; EXTERNAL ROUTINES:
;     rtipsy
;     
;     You can get both from me
; EXTERNAL DATA:
;     Padova SPS grids.  I have hard coded their path into the code.
;     Not ideal, but live with it for now.
; COMMENTS:
;     I apologize for use of a common block.
;     This has to be changed
;
;     All stars with Z = 0.000 are assigned Z = .0005.  
;     This is hard coded for now.
;
;     I hard coded a time unit of 4e10
;
;     The interpolation I use is stupid, but it will be improved.
; 
;     To allow my stupid interploation:
;     Right now I set 'minimum' and 'maximum' log_ages as [8.5,10.1]
;     Right now I set 'minimum' and 'maximum' metals as [.0004,.1]

;     The "halo" output fits file contains a structure with:
;      r (is 0.00 for starless halo)
;      v (is 0.00 for starless halo)
;      b (is 0.00 for starless halo)
;      i (is 0.00 for starless halo)
;      k (is 0.00 for starless halo)
;      u (is 0.00 for starless halo)
;      rprime (SDSS; these are *not* properly normalized now)
;      iprime (SDSS; these are *not* properly normalized now)
;      Z (average metallicity, lum weighted)
;      age (average, lum weighted)
;
;      To read:  IDL> halo = MRDFITS(halo_file,1)
;      halo.r etc .. contain data
;
;     The "stars" output fits file contains a structure with:
;      r
;      v
;      b
;      i
; REVISION HISTORY:
;   2005-Jun-24  Beth Willman

common halo_stats, group, npart, vc, vc_half, radius, r_vc_max, r_m_half, xcenter, ycenter, zcenter, xVcm, yVcm, zVcm, M_star, m_tot, m_gas, vdisp

IF (keyword_set (timeunit) EQ 0) THEN timeunit=1e9
IF (keyword_set(mass_unit) EQ 0) THEN mass_unit=2.326e5
datadir='/net/grads-1/stinson/trace/'

;print, tipsy_files
outfile = file + '.Mv.fits'
outfile2 = file + '.stars.Mv.fits'

;-----------------;
; Read SPS grids  ;
;-----------------;
Z = [.0004,.004,.008,.02,.05,.1]
n_metallicity = 26
n_age = 18
V_grid = fltarr(26,18)
vr_grid = fltarr(26,18)
ri_grid = fltarr(26,18)
bv_grid = fltarr(26,18)
vk_grid = fltarr(26,18)
openr, 1, datadir+'V_grid'
openr, 2, datadir+'vr_grid'
openr, 3, datadir+'ri_grid'
openr, 4, datadir+'bv_grid'
openr, 5, datadir+'vk_grid'
openr, 6, datadir+'ub_grid'
readf, 1, V_grid
readf, 2, vr_grid
readf, 3, ri_grid
readf, 4, bv_grid
readf, 5, vk_grid
readf, 6, ub_grid
close, /all

age_vector = fltarr(n_age)
metallicity_vector = fltarr(n_metallicity)
age_max = 10.2
age_min = 8.5
age_bin = (age_max-age_min)/(n_age - 1.)
metal_max = .1
metal_min = .0004
metal_bin = (alog10(metal_max)-alog10(metal_min))/(n_metallicity - 1.)

for i = 0, n_age-1 do begin
 age_vector[i] = age_min + i*age_bin
endfor
for i = 0, n_metallicity-1 do begin
 metallicity_vector[i] = 10.^(alog10(metal_min)+i*metal_bin)
endfor

;-------------------;
; loop through each ;
; timestep          ;
;-------------------;
    tipsy_file = file

    temp = FILE_SEARCH(tipsy_file)

    IF (strmatch(temp[0], tipsy_file)) then begin
	print, tipsy_file
        ;xc = 7.168929e-2
        ;yc = -7.45516e-2
        ;zc = -7.94397e-2

        RTIPSY, tipsy_file, h, gas, dark, star

        halo = {r:0., v:0., b:0., i:0., k:0., u:0., rprime:0., iprime:0., Z:0., age:0.}
        star_mags = replicate({r:0., v:0., b:0., i:0., k:0., u:0.},h.nstar)

        ; For runs with no metals, this assigns minimal metals
        ;print, ''
        ;print, tipsy_file
        print, 'metallicity range:', MINMAX(star.metals)
        if MAX(star.metals) eq 0 then star.metals = REFORM(fltarr(n_elements(star))+.0004)

        age_star = (h.time - star.tform) * timeunit
        log_age = alog10(age_star)
        print, 'age range:', MINMAX(log_age)
	print, ''
        log_mass = alog10(mass_unit*star.mass)

        temp = where (star.metals ge 0.1)
        if temp[0] ne -1 then star(temp).metals = 0.1
        temp = where (star.metals le .0004)
        if temp[0] ne -1 then star(temp).metals = 0.0004
        temp = where (log_age ge 10.12)
        if temp[0] ne -1 then log_age(temp) = 10.12
        temp = where (log_age le 8.5)
        if temp[0] ne -1 then log_age(temp) = 8.5

;-----------------------;
; calc phot params for  ;
;   stars               ;
;-----------------------;
        metal_ref = (alog10(star.metals) - alog10(metal_min))/metal_bin
        metal_ref_diff = metal_ref - FIX(metal_ref)
        temp = where (metal_ref_diff lt 0.5)
        if temp[0] ne -1 then metal_ref(temp) = metal_ref(temp)
        temp = where (metal_ref_diff ge 0.5)
        if temp[0] ne -1 then metal_ref(temp) = metal_ref(temp)+1

        age_ref = (log_age - age_min)/age_bin
        age_ref_diff = age_ref - FIX(age_ref)
        age_temp = where (age_ref_diff lt 0.5)
        age_temp2 = where (age_ref_diff ge 0.5)
        if age_temp[0] ne -1 then age_ref(age_temp) = age_ref(where (age_ref_diff lt 0.5))
        if age_temp2[0] ne -1 then age_ref(age_temp2) = age_ref(where (age_ref_diff ge 0.5))+1

        star_mags.v = V_grid[metal_ref,age_ref]
        star_mags.r = V_grid[metal_ref,age_ref]-vr_grid[metal_ref,age_ref]
        star_mags.b = V_grid[metal_ref,age_ref]+bv_grid[metal_ref,age_ref]
        star_mags.i = star_mags.r-ri_grid[metal_ref,age_ref]
        star_mags.k = V_grid[metal_ref,age_ref]-vk_grid[metal_ref,age_ref]
        star_mags.u = ub_grid[metal_ref,age_ref] + star_mags.b

        star_mags.v = star_mags.v - 2.51*log_mass
        star_mags.r = star_mags.r -  2.51*log_mass
        star_mags.b = star_mags.b -  2.51*log_mass
        star_mags.i = star_mags.i -  2.51*log_mass
        star_mags.k = star_mags.k -  2.51*log_mass
        star_mags.u = star_mags.u -  2.51*log_mass

;------------------------;
; Get Magnitudes of Each ;
; Halo                   ;
;------------------------;
        npartstar = h.nstar
        mass = TOTAL(10.^log_mass)
        lumtemp = TOTAL(10.^(-0.4*star_mags.k))

        Z_temp = TOTAL(star.metals)/npartstar
        age_temp = alog10(TOTAL(10.^(log_age)))/npartstar
        halo.Z = TOTAL(star.metals * 10.^(-0.4*star_mags.k))/lumtemp
        halo.age = alog10(TOTAL(10.^log_age * 10.^(-0.4*star_mags.k))/lumtemp)
        ;print, 'Age test:', age_temp, halo[i-1].age
        ;print, 'Z test:', Z_temp, halo[i-1].Z
        halo.r = TOTAL(10.^(-0.4*star_mags.r))
        rlum = TOTAL(10.^(-0.4*(star_mags.v-4.46)))
        halo.r = -2.5*alog10(halo.r)
print,'absolute r:'+string(halo.r)
        halo.v = TOTAL(10.^(-0.4*star_mags.v))
        vlum = TOTAL(10.^(-0.4*(star_mags.v-4.83)))
print,'v luminosity:'+string(vlum)
        halo.v = -2.5*alog10(halo.v)
print,'absolute v:'+string(halo.v)
        halo.b = TOTAL(10.^(-0.4*star_mags.b))
        blum = TOTAL(10.^(-0.4*(star_mags.b-5.45)))
        halo.b = -2.5*alog10(halo.b)
print,'absolute b:'+string(halo.b)
        halo.i = TOTAL(10.^(-0.4*star_mags.i))
        halo.i = -2.5*alog10(halo.i)
print,'absolute i:'+string(halo.i)
        halo.k = TOTAL(10.^(-0.4*star_mags.k))
        halo.k = -2.5*alog10(halo.k)
print,'absolute k:'+string(halo.k)
        halo.u = TOTAL(10.^(-0.4*star_mags.u))
        halo.u = -2.5*alog10(halo.u)
print,'absolute u:'+string(halo.u)

;-------------;
; in sdss mag ;
;-------------;
        temp = where (halo.r - halo.i lt 1.15)
        temp2 = where (halo.r - halo.i ge 1.15)
        r_minus_i = fltarr(n_elements(halo))

        halo.rprime = (halo.v - .81*(halo.v - halo.r) + .13)
        if temp[0] ne -1 then r_minus_i(temp) = (halo.r - halo.i) - 0.21
        if temp2[0] ne -1 then r_minus_i(temp2) = 1.42*(halo.r - halo.i) - 0.69
        halo.iprime = halo.rprime - r_minus_i

;---------------------------;
; calculate halo mags based ;
; on average Z and age and  ;
; compare with the summed   ;
; values                    ;
;---------------------------;
;    b_halo_avg = fltarr(n_elements(r_halo))
;    r_halo_avg = fltarr(n_elements(r_halo))
;    rprime_avg = fltarr(n_elements(r_halo))

;    metal_ref = (alog10(Z_halo(where (Z_halo ne 0.))) - alog10(metal_min))/metal_bin
;    metal_ref_diff = metal_ref - FIX(metal_ref)
;    metal_ref(where (metal_ref_diff lt 0.5)) = metal_ref(where (metal_ref_diff lt 0.5))
;    metal_ref(where (metal_ref_diff ge 0.5)) = metal_ref(where (metal_ref_diff ge 0.5))+1

;    age_ref = (age_halo(where (Z_halo ne 0.)) - age_min)/age_bin
;    age_ref_diff = age_ref - FIX(age_ref)
;    age_temp = where (age_ref_diff lt 0.5)
;    age_temp2 = where (age_ref_diff ge 0.5)
;    if age_temp[0] ne -1 then age_ref(age_temp) = age_ref(where (age_ref_diff lt 0.5))
;    if age_temp2[0] ne -1 then age_ref(age_temp2) = age_ref(where (age_ref_diff ge 0.5))+1

;        Z = Z_halo(temp)
;        age = age_halo(temp)
;        b_halo_avg = (V_grid[metal_ref,age_ref]+bv_grid[metal_ref,age_ref])-2.51*alog10(mass_unit*M_star(temp))
;        r_halo_avg = (V_grid[metal_ref,age_ref]-vr_grid[metal_ref,age_ref])-2.51*alog10(mass_unit*M_star(temp))

;    print, r_halo(temp[0:9])
;    print, r_halo_avg(temp[0:9])

;    print, ''
;    print, b_halo(temp[0:9])
;    print, b_halo_avg(temp[0:9])

;    print, ''
;    print, b_halo(temp[0:9])-r_halo(temp[0:9])
;    print, b_halo_avg(temp[0:9])-r_halo_avg(temp[0:9])

;-----------------------------;
; output magnitudes to a file ;
; one file for each timestep, ;
; one entry for each halo     ;
;-----------------------------;
        MWRFITS, halo, outfile, /create
        MWRFITS, star_mags, outfile2, /create
        ENDIF ELSE print, STRING('No matching tipsy file found for: ' + file)

RETURN,blum
;RETURN,rlum
END

temp = SORT(b_halo)
for i = 0,10 do begin
 print, V_halo(temp[i]), b_halo(temp[i]) - r_halo(temp[i]), alog10(mass(temp[i]))
endfor

;-----------------;
;-----------------;

RETURN
END
