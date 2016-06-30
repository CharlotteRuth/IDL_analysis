pro get_magnitudes

common halo_stats, group, vc, vc_half, radius, r_vc_max, r_m_half, xcenter, ycenter, zcenter, xVcm, yVcm, zVcm, M_star, m_tot

readcol, '/host/toast/parkay/willman/nbody/tipsy/IClight_runs/time_vs_z.log', timearr, zarr

;step_array = [1024,998,973,925,896,852,839,768,698,641,590,505,437,384,371]
;rvir_array = [.02127,.02155,.0215,.0201,.0192,.0176,.0173,.01635,.01675,.01635,.0141,.01255,.0124,.0112,.01075]
;tipsy_prefix = '/net/grads-1/willman/coma/cl2c6.288gSF/TIPSY_files/cl2c6.288gSF.0'
;skid_prefix = '/net/grads-1/willman/coma/cl2c6.288gSF/SKID_files/cl2c6.288gSF.tau125.dark.0'
;outfile_prefix = '/net/grads-1/willman/coma/cl2c6.288gSF/cluster_mags_uncorr.tau125.0'

;step_array = [1024,998,973,925,896,852,839,768,698,641,590,505,437,384,371]
;rvir_array = [.02127,.02155,.0215,.0201,.0192,.0176,.0173,.01635,.01675,.01635,.0141,.01255,.0124,.0112,.01075]
step_array = [1024]
rvir_array = [.02127]
tipsy_prefix = '/net/grads-1/willman/coma/cl2c6.288gSF/TIPSY_files/cl2c6.288gSF.0'
skid_prefix = '/net/grads-1/willman/coma/cl2c6.288gSF/SKID_files/cl2c6.288gSF.tau187.dark.0'
outfile_prefix = '/net/grads-1/willman/coma/cl2c6.288gSF/cluster_mags_uncorr.0'

;step_array = [1024,998,973,925,903,860,839,820,800,768,698,640,590,512,505,437]
;rvir_array = [.0215,.02145,.0214,.0201,.0193,.0179,.0173,.01665,.01635,.0162,.01675,.01635,.0141,.01245,.01255,.0124]
;tipsy_prefix = '/net/grads-1/willman/coma/cl2c6.288gSF.05k/TIPSY_files/cl2c6.288gSF.z3.05k.0'
;skid_prefix = '/net/grads-1/willman/coma/cl2c6.288gSF.05k/SKID_files/cl2c6.288gSF.z3.05k.tau187.dark.0'
;outfile_prefix = '/net/grads-1/willman/coma/cl2c6.288gSF.05k/cluster_mags_uncorr.0'

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
openr, 1, './data/V_grid'
openr, 2, './data/vr_grid'
openr, 3, './data/ri_grid'
openr, 4, './data/bv_grid'
openr, 5, './data/vk_grid'
openr, 6, './data/ub_grid'
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

print, age_vector
print, ''
print, metallicity_vector
print, ''

;print, bv_grid+vr_grid

;age_temp = [10.02,9.8,10.165]
age_temp = [10.03,10.03,9.7,9.6,9.6,10.1,10.1,10.1]
;Z_temp = [.0001,.009,.0179]
Z_temp = [0.009,.018,0.018,.002,.01,.002,.008,.018]
metal_ref = (alog10(Z_temp) - alog10(metal_min))/metal_bin
metal_ref_diff = metal_ref - FIX(metal_ref)
metal_ref(where (metal_ref_diff lt 0.5)) = metal_ref(where (metal_ref_diff lt 0.5))
metal_ref(where (metal_ref_diff ge 0.5)) = metal_ref(where (metal_ref_diff ge 0.5))+1

age_ref = (age_temp - age_min)/age_bin
age_ref_diff = age_ref - FIX(age_ref)
age_temp = where (age_ref_diff lt 0.5)
age_temp2 = where (age_ref_diff ge 0.5)
if age_temp[0] ne -1 then age_ref(age_temp) = age_ref(where (age_ref_diff lt 0.5))
if age_temp2[0] ne -1 then age_ref(age_temp2) = age_ref(where (age_ref_diff ge 0.5))+1

b_halos = (V_grid[metal_ref,age_ref]+bv_grid[metal_ref,age_ref])
r_halos = (V_grid[metal_ref,age_ref]-vr_grid[metal_ref,age_ref])

rprime_halos = V_grid[metal_ref,age_ref] -0.44*bv_grid[metal_ref,age_ref] + 0.12

print, rprime_halos
print, r_halos

print, b_halos - r_halos

RETURN

;-------------------;
; loop through each ;
; timestep          ;
;-------------------;
FOR j = 0, n_elements(step_array)-1 do begin
    if j eq 0 then tipsy_file = STRCOMPRESS( (tipsy_prefix+strn(step_array[j],format='(I4)')) )
    if j gt 0 then tipsy_file = STRCOMPRESS( (tipsy_prefix+'0'+strn(step_array[j],format='(I3)')) )
    if j eq 0 then grp_file = STRCOMPRESS( (skid_prefix+strn(step_array[j],format='(I4)')+'.grp') )
    if j gt 0 then grp_file = STRCOMPRESS( (skid_prefix+'0'+strn(step_array[j],format='(I3)')+'.grp') )
    if j eq 0 then stat_file = STRCOMPRESS( (skid_prefix+strn(step_array[j],format='(I4)')+'.stat') )
    if j gt 0 then stat_file = STRCOMPRESS( (skid_prefix+'0'+strn(step_array[j],format='(I3)')+'.stat') )
    if j eq 0 then outfile = STRCOMPRESS( (outfile_prefix+strn(step_array[j],format='(I4)')) )
    if j gt 0 then outfile = STRCOMPRESS( (outfile_prefix+'0'+strn(step_array[j],format='(I3)')) )
    READ_STAT, file = stat_file
    nhalos = n_elements(M_star)
    RTIPSY, tipsy_file, h, gas, dark, star
    OPENR, 2, grp_file
    grp_data = intarr(h.ndark+h.nstar+h.ngas+1)
    readf, 2, grp_data
    close, 2
    halo_index = grp_data[h.ngas+h.ndark+1:h.ndark+h.nstar+h.ngas]

    z_star_form = INTERPOL(zarr, timearr, star.tform)
    age_star = fltarr(n_elements(z_star_form))
    for i = 0l, n_elements(z_star_form) - 1 do age_star(i) = lookback(z_star_form(i))
    log_age = alog10(age_star)
    log_mass = alog10(13.6e16*star.mass)

    temp = where (star.metals ge 0.1)
    if temp[0] ne -1 then star(temp).metals = 0.1
    temp = where (star.metals le .0004)
    if temp[0] ne -1 then star(temp).metals = 0.0004
    temp = where (log_age ge 10.1)
    if temp[0] ne -1 then log_age(temp) = 10.2
    temp = where (log_age le 8.5)
    if temp[0] ne -1 then log_age(temp) = 8.5

;-----------------------;
; calc phot params for  ;
;   stars               ;
;-----------------------;
    metal_ref = (alog10(star.metals) - alog10(metal_min))/metal_bin
    metal_ref_diff = metal_ref - FIX(metal_ref)
    metal_ref(where (metal_ref_diff lt 0.5)) = metal_ref(where (metal_ref_diff lt 0.5))
    metal_ref(where (metal_ref_diff ge 0.5)) = metal_ref(where (metal_ref_diff ge 0.5))+1

    age_ref = (log_age - age_min)/age_bin
    age_ref_diff = age_ref - FIX(age_ref)
    age_temp = where (age_ref_diff lt 0.5)
    age_temp2 = where (age_ref_diff ge 0.5)
    if age_temp[0] ne -1 then age_ref(age_temp) = age_ref(where (age_ref_diff lt 0.5))
    if age_temp2[0] ne -1 then age_ref(age_temp2) = age_ref(where (age_ref_diff ge 0.5))+1

    v_star = V_grid[metal_ref,age_ref]
    r_star = V_grid[metal_ref,age_ref]-vr_grid[metal_ref,age_ref]
    b_star = V_grid[metal_ref,age_ref]+bv_grid[metal_ref,age_ref]
    i_star = r_star-ri_grid[metal_ref,age_ref]
    k_star = V_grid[metal_ref,age_ref]-vk_grid[metal_ref,age_ref]
    u_star = ub_grid[metal_ref,age_ref] + b_star

    v_star = v_star - 2.51*log_mass
    r_star = r_star -  2.51*log_mass
    b_star = b_star -  2.51*log_mass
    i_star = i_star -  2.51*log_mass
    k_star = k_star -  2.51*log_mass
    u_star = u_star -  2.51*log_mass

;------------------------;
; Get Magnitudes of Each ;
; Halo                   ;
;------------------------;
    b_halo = fltarr(nhalos)
    v_halo = fltarr(nhalos)
    r_halo = fltarr(nhalos)
    i_halo = fltarr(nhalos)
    k_halo = fltarr(nhalos)
    u_halo = fltarr(nhalos)
    Z_halo = fltarr(nhalos)
    age_halo = fltarr(nhalos)
    big_halo = where (M_star eq MAX(M_star))
    for i = 1, nhalos do begin
        halo_member = where (halo_index eq i)
        if n_elements(halo_member) gt 3 then begin
            v_stars = v_star(halo_member)
            r_stars = r_star(halo_member) 
            b_stars = b_star(halo_member) 
            i_stars = i_star(halo_member) 
            k_stars = k_star(halo_member) 
            u_stars = u_star(halo_member) 

;            if big_halo+1 eq i then begin                
;                ops, file = '/net/grads-1/willman/coma/cl2c6.288gSF/tau187_plots/age_vs_z_cd.ps', /landscape 
;                plot, star(halo_member).metals, log_age(halo_member), psym = 3, yrange = [8.5,10.15]
;                cps, /noprint
;            endif

            Z_halo[i-1] = TOTAL(star(halo_member).metals)/n_elements(halo_member)
            age_halo[i-1] = alog10(TOTAL(10.^(log_age(halo_member)))/n_elements(halo_member))
            r_halo[i-1] = TOTAL(10.^(-0.4*r_stars))
            r_halo[i-1] = -2.5*alog10(r_halo[i-1])
            v_halo[i-1] = TOTAL(10.^(-0.4*v_stars))
            v_halo[i-1] = -2.5*alog10(v_halo[i-1])
            b_halo[i-1] = TOTAL(10.^(-0.4*b_stars))
            b_halo[i-1] = -2.5*alog10(b_halo[i-1])
            i_halo[i-1] = TOTAL(10.^(-0.4*i_stars))
            i_halo[i-1] = -2.5*alog10(i_halo[i-1])
            k_halo[i-1] = TOTAL(10.^(-0.4*k_stars))
            k_halo[i-1] = -2.5*alog10(k_halo[i-1])
            u_halo[i-1] = TOTAL(10.^(-0.4*u_stars))
            u_halo[i-1] = -2.5*alog10(u_halo[i-1])
        endif
    endfor

;-------------;
; in sdss mag ;
;-------------;
    rprime = fltarr(n_elements(r_halo))
    r_minus_i = fltarr(n_elements(r_halo))
    iprime = fltarr(n_elements(r_halo))

    temp = where (r_halo - i_halo lt 1.15)
    temp2 = where (r_halo - i_halo ge 1.15)

    rprime = (v_halo - .81*(v_halo - r_halo) + .13)
    if temp[0] ne -1 then r_minus_i(temp) = (r_halo - i_halo) - 0.21
    if temp2[0] ne -1 then r_minus_i(temp2) = 1.42*(r_halo - i_halo) - 0.69
    iprime = rprime - r_minus_i

;---------------------------;
; calculate halo mags based ;
; on average Z and age and  ;
; compare with the summed   ;
; values                    ;
;---------------------------;
    b_halo_avg = fltarr(n_elements(r_halo))
    r_halo_avg = fltarr(n_elements(r_halo))
    rprime_avg = fltarr(n_elements(r_halo))

    metal_ref = (alog10(Z_halo(where (Z_halo ne 0.))) - alog10(metal_min))/metal_bin
    metal_ref_diff = metal_ref - FIX(metal_ref)
    metal_ref(where (metal_ref_diff lt 0.5)) = metal_ref(where (metal_ref_diff lt 0.5))
    metal_ref(where (metal_ref_diff ge 0.5)) = metal_ref(where (metal_ref_diff ge 0.5))+1

    age_ref = (age_halo(where (Z_halo ne 0.)) - age_min)/age_bin
    age_ref_diff = age_ref - FIX(age_ref)
    age_temp = where (age_ref_diff lt 0.5)
    age_temp2 = where (age_ref_diff ge 0.5)
    if age_temp[0] ne -1 then age_ref(age_temp) = age_ref(where (age_ref_diff lt 0.5))
    if age_temp2[0] ne -1 then age_ref(age_temp2) = age_ref(where (age_ref_diff ge 0.5))+1

    r_halo_avg(where (Z_halo ne 0.)) = (V_grid[metal_ref,age_ref]-vr_grid[metal_ref,age_ref]) - 2.51*(alog10(13.6d16*M_star(where (Z_halo ne 0.))))
    b_halo_avg(where (Z_halo ne 0.)) = (V_grid[metal_ref,age_ref]+bv_grid[metal_ref,age_ref]) - 2.51*(alog10(13.6d16*M_star(where (Z_halo ne 0.))))
    rprime_avg(where (Z_halo ne 0.)) = V_grid[metal_ref,age_ref] - 0.81*vr_grid[metal_ref,age_ref] + 0.13 - 2.51*(alog10(13.6d16*M_star(where (Z_halo ne 0.))))

    temp = (SORT(r_halo))

    print, ''
    print, r_halo(temp[0:9])
    print, r_halo_avg(temp[0:9])

    print, ''
    print, rprime(temp[0:9])
    print, rprime_avg(temp[0:9])

    print, ''
    print, b_halo(temp[0:9])
    print, b_halo_avg(temp[0:9])

    print, ''
    print, b_halo(temp[0:9])-r_halo(temp[0:9])
    print, b_halo_avg(temp[0:9])-r_halo_avg(temp[0:9])
RETURN

;-----------------------------;
; output magnitudes to a file ;
; one file for each timestep, ;
; one entry for each halo     ;
;-----------------------------;
    openw, 1, outfile
    for i = 0l, nhalos - 1 do printf, 1, format = '(8(F7.2), F8.5, F6.2)', r_halo[i], v_halo[i], b_halo[i], i_halo[i], k_halo[i], u_halo[i], rprime[i], iprime[i], Z_halo[i], age_halo[i]
    close, 1

ENDFOR

;-----------------;
;-----------------;
temp = where (b_halo lt -5)
metal_ref = (alog10(Z_halo(temp)) - alog10(metal_min))/metal_bin
metal_ref_diff = metal_ref - FIX(metal_ref)
metal_ref(where (metal_ref_diff lt 0.5)) = metal_ref(where (metal_ref_diff lt 0.5))
metal_ref(where (metal_ref_diff ge 0.5)) = metal_ref(where (metal_ref_diff ge 0.5))+1

age_ref = (age_halo(temp) - age_min)/age_bin
age_ref_diff = age_ref - FIX(age_ref)
age_temp = where (age_ref_diff lt 0.5)
age_temp2 = where (age_ref_diff ge 0.5)
if age_temp[0] ne -1 then age_ref(age_temp) = age_ref(where (age_ref_diff lt 0.5))
if age_temp2[0] ne -1 then age_ref(age_temp2) = age_ref(where (age_ref_diff ge 0.5))+1

Z = Z_halo(temp)
age = age_halo(temp)
b_halos = (V_grid[metal_ref,age_ref]+bv_grid[metal_ref,age_ref])-2.51*alog10(13.6d16*M_star(temp))
r_halos = (V_grid[metal_ref,age_ref]-vr_grid[metal_ref,age_ref])-2.51*alog10(13.6d16*M_star(temp))

temp = SORT(b_halos)
for i = 0,10 do begin
 print, b_halos(temp[i]), r_halos(temp[i]), b_halos(temp[i]) - r_halos(temp[i]), age(temp[i]), Z(temp[i])
endfor

RETURN
END
