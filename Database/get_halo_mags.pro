pro get_halo_mags, mass_unit, length_unit, halo_finder = halo_finder, radius_file = radiusfile, iso_set = iso_set, mag_sys = mag_sys, init_mass = init_mass, snaptime = snaptime, multiple = multiple, tipsy_file = tipsy_file, virial_radius = virial_radius

; NAME:
;   get_halo_mags
; PURPOSE:
;   Get the absolute magnitudes of each star in a simulation output,
;   as well as of each halo ID'd by SKID or in an entire simulation
; CALLING SEQUENCE:
;  IDL> get_halo_mags, mass_unit, length_unit
;
;  Examples:
;  IDL>
;  get_halo_mags,3.171e15,28.5714,iso_set=0,init_mass=11.92e-12,snaptime= 0.333091,mag_sys='ab',tipsy_file='galaxy.std'

;  IDL>  get_halo_mags,2.310e15,  25,iso_set=0,init_mass=4.331e-13,mag_sys='ab',/multiple
;  IDL>  get_halo_mags,1.69875e16,50,iso_set=0,init_mass=2.331e-10,snaptime=0.333091,mag_sys='ab',/multiple ;call for multiple amiga files
;  IDL>  get_halo_mags,1.69875e16,50,iso_set=0,init_mass=2.331e-10,snaptime=0.333091,mag_sys='ab',tipsy_file = 'h1201-h1755-X2X2g1bwK.00512.182.std', call for when the magnitude of the entire simulation is desired
;  IDL>  get_halo_mags,3.17111799e15,28.571,iso_set=0,init_mass=7.630752e-10,  snaptime = 0.3330907, mag_sys = 'ab', tipsy_file = '/astro/net/scratch2/christensen/Sunrise/MW1lr_v2/MW1lr_001.std'
;  IDL>  get_halo_mags,1.84854092e16,50,iso_set = 0, init_mass =  2.7048e-11,tipsy_file = 'galaxy.std',mag_sys = 'ab',snaptime = 0.457870
;  IDL>  get_halo_mags,9.9669e16,87.671,iso_set = 2,init_mass = 1.46042e-12,mag_sys = 'ab',/multiple
;
;  IDL>  get_halo_mags,1.84854092e16,50,iso_set = 0, init_mass=4.330647e-13,snaptime=0.333091,mag_sys='ab',/multiple
;
;   If you use the "multiple" option, this must be done *in* the directory containing both the tipsy
;   files and the SKID output files you want to analyze
;
;   The code only connects  SKID and TIPSY files that have the same
;   prefix (i.e. tipsyfile.grp) or have only a g added
;   (i.e. tipsyfileg.grp)
;
; INPUTS:
;   mass_unit - mass unit of simulation in solar masses
;   length_unit - length unit of the simulation in Mpc
;
;   halo_finder - the name of the halo finder used (either 'skid' or 'amiga'), 
;     amiga is default  
;   radius_file is an optional file to define the radius for each halo
;      you want used when determining stellar membership, rather than
;      using the particles selected by amiga
;   iso_set -- 0 for Starburst99 with Kroupa IMF, 2 for Starburst99
;              with MS IMF, 3 for Leo's, 4 for Beth Willman's  --
;      Starburst99 is default
;   mag_sys can be either 'ab' or 'vega', 'vega' is default
;   init_mass is the initial mass of a star particle in simulation units,
;      if not provided, it is set to the mass of the largest particle
;   snaptime is the age of the simulation in simulation units,
;      if not provided, is is set to the formation time of the
;      youngest particle
;   multiple -- This flag must be turned on if you want to find the magnitudes of
;      multiple timesteps or multiple halos in a simulation
;      If it is not turned on, a tipsy filename of the simulation must
;      be provided
;      If you do want to find the magnitude of multiple files or halos, the
;      program must be run from a directory that contains all of them
;   tipsy_name -- If multiple is not turned on and you want the
;      magnitude of the entire simulation, use this option give the
;      name of the simulation
;                 
;
; OUTPUTS:
;   tipsyfile.halos.Mv.fits
;   tipsyfile.stars.Mv.fits
;   
;   For each tipsyfile with SKID output the code identifies
;   in the working directory
; EXTERNAL ROUTINES:
;     rtipsy
;     read_stat_struc_amiga.pro
;     get_timestep     
;     
;     You can get all three from me
; EXTERNAL DATA:
;     Leo Girardi, Starburst99 and Padova SSP grids.  I have hard coded their path into the code.
;     Not ideal, but live with it for now.
;
;     If you are working on a different computer system or would like
;     to use your own copies of the grids, you will need to update
;     grid path on line 70
;
; COMMENTS:
;     All stars with Z = 0.000 are assigned the minimum metallicity of
;     the models  
;
;     The interpolation is done using the IDL procedure TRI_SURF.
;     It's ok, not great

;     The "halos" output fits file contains a structure with:
;      r (is 0.00 for starless halos)
;      v (is 0.00 for starless halos)
;      b (is 0.00 for starless halos)
;      i (is 0.00 for starless halos)
;      k (is 0.00 for starless halos)
;      u (is 0.00 for starless halos)
;      rprime (SDSS; these are *not* properly normalized now)
;      iprime (SDSS; these are *not* properly normalized now)
;      Z (average metallicity, lum weighted)
;      age (average, lum weighted)
;
;      To read:  IDL> halos = MRDFITS(halo_file,1)
;      halos.r etc .. contain data
;
;     The "stars" output fits file contains a structure with:
;      r
;      v
;      b
;      i
; REVISION HISTORY:
;   2005-Jun-24  Beth Willman
;   2008-Jul-23  Charlotte Christensen
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN gridpath = '/home/christensen/Code/Datafiles/Database/' ELSE $
  IF strcmp(hostname,'pfe',3) THEN gridpath = '/home6/crchrist/Datafiles/Database/' ELSE gridpath = '/astro/users/christensen/code/Datafiles/Database/'

close,/all
sec_to_year = 365.*24.*3600.
timestep = get_timestep(mass_unit,length_unit*1000.0)
;NOT in system units, mass in Solar mass, tfom in years, age is log age
IF (keyword_set(iso_set)) THEN isochrone = iso_set ELSE isochrone = 0
CASE isochrone OF
    0:  iso_name = 'star99_K'
    1:  iso_name = 'star99_MS'
    2:  iso_name = 'leo'
    3:  iso_name = 'bw'
ENDCASE
IF (NOT keyword_set(mag_sys)) THEN mag_sys = 'vega'

grp_files = ['']
stat_files = ['']
IF (keyword_set(multiple)) THEN BEGIN
    IF (keyword_set(halo_finder)) THEN halo_finder = halo_finder ELSE halo_finder = 'amiga'
    IF halo_finder eq 'skid' then begin
        grp_files = file_search('*.gd.skid.grp')
        stat_files = STRARR(n_elements(grp_files))
        tipsy_files = STRARR(n_elements(grp_files))
        FOR j = 0, n_elements(grp_files)-1 do begin
            IF (strmatch(grp_files[j], '*gd.skid.grp')) then begin
                tipsy_files[j] = (STREGEX(grp_files[j], '(.*)\.gd.skid.grp', /extr, /sub))[1]
                stat_files[j] = (STREGEX(grp_files[j], '(.*)\.grp', /extr, /sub))[1] + '.stat'
            ENDIF
        ENDFOR
    ENDIF
    IF halo_finder eq 'amiga' then begin
        stat_files = file_search('*.amiga.stat')
        grp_files = STRARR(n_elements(stat_files))
        tipsy_files = STRARR(n_elements(stat_files))
        FOR j = 0, n_elements(stat_files)-1 do begin
            IF (strmatch(stat_files[j], '*amiga.stat')) then begin
                tipsy_files[j] = (STREGEX(stat_files[j], '(.*)\.amiga.stat', /extr, /sub))[1]
                grp_files[j] = (STREGEX(stat_files[j], '(.*)\.stat', /extr, /sub))[1] + '.grp'
            ENDIF
        ENDFOR
    ENDIF
    IF(n_elements(grp_files) gt 27) THEN BEGIN
        grp_files = grp_files[27:n_elements(grp_files)-1]
        stat_files = stat_files[27:n_elements(stat_files)-1]
        tipsy_files = tipsy_files[27:n_elements(tipsy_files)-1]
    ENDIF
ENDIF ELSE tipsy_files = [tipsy_file]

;looping through the different tipsteps of outputs
IF (mag_sys eq 'vega') THEN outfiles_stars = tipsy_files + '.starsbk.'+iso_name+'_vega.Mv.fits' ELSE outfiles_stars = tipsy_files + '.stars.'+iso_name+'_ab.Mv.fits'
IF (keyword_set(halo_finder)) THEN BEGIN
    IF halo_finder eq 'skid' then begin
        IF (keyword_set(RADIUSFILE) OR keyword_set(virial_radius)) THEN outfiles = tipsy_files + '.halos_r200.'+iso_name ELSE outfiles = tipsy_files + '.halos_virbk.'+iso_name
    ENDIF
    IF halo_finder eq 'amiga' then begin
        IF (keyword_set(RADIUSFILE) OR keyword_set(virial_radius)) THEN outfiles = tipsy_files + '.amiga_r200.halos.'+iso_name ELSE outfiles = tipsy_files + '.amiga_virbk.halos.'+iso_name
    ENDIF
ENDIF ELSE outfiles = STRMID(tipsy_file,0,STRLEN(tipsy_file)-4)+'.'+iso_name
IF (mag_sys eq 'vega') THEN outfiles = outfiles + '_vega.Mv.fits' ELSE outfiles = outfiles + '_ab.Mv.fits'

;-----------------;
; Read SPS grids  ;
;-----------------; 
IF (isochrone eq 0) THEN BEGIN
;Starburst99 Isochrones with Kroupa IMF ----
    print,'Starburst 99 Isochrones with Kroupa IMF'
    n_metallicity = 5
    readcol,gridpath+'data_star99_K/ages',age_vector
    age_vector = ALOG10(age_vector)
    age_min = MIN(age_vector)
    age_max = MAX(age_vector)
    n_age = n_elements(age_vector)
    metallicity_vector = [0.0004,0.004,0.008,0.02,0.05]

    openr, 1, gridpath+'data_star99_K/V_grid'
    openr, 2, gridpath+'data_star99_K/vr_grid'
    openr, 3, gridpath+'data_star99_K/ri_grid'
    openr, 4, gridpath+'data_star99_K/bv_grid'
    openr, 5, gridpath+'data_star99_K/vk_grid'
    openr, 6, gridpath+'data_star99_K/ub_grid'
    openr, 7, gridpath+'data_star99_K/mb_grid'
    mb_grid = fltarr(n_metallicity,n_age)
    readf, 7, mb_grid
ENDIF
IF (isochrone eq 1) THEN BEGIN
;Starburst99 Isochrones with MS IMF ----
    printf,'Starburst 99 Isochrones with Miller-Scalo IMF'
    n_metallicity = 5
    readcol,gridpath+'data_star99_MS/ages',age_vector
    age_vector = ALOG10(age_vector)
    age_min = MIN(age_vector)
    age_max = MAX(age_vector)
    n_age = n_elements(age_vector)
    metallicity_vector = [0.0004,0.004,0.008,0.02,0.05]

    openr, 1, gridpath+'data_star99_MS/V_grid'
    openr, 2, gridpath+'data_star99_MS/vr_grid'
    openr, 3, gridpath+'data_star99_MS/ri_grid'
    openr, 4, gridpath+'data_star99_MS/bv_grid'
    openr, 5, gridpath+'data_star99_MS/vk_grid'
    openr, 6, gridpath+'data_star99_MS/ub_grid'
    openr, 7, gridpath+'data_star99_MS/mb_grid'
    mb_grid = fltarr(n_metallicity,n_age)
    readf, 7, mb_grid
ENDIF
IF (isochrone eq 2) THEN BEGIN
;Leo Girardi Isochrones ---
    printf,'Leo Girardis Isochrones'
    n_age = 42
    n_metallicity = 15 
    age_bin = 0.1
    age_min = 6.0
    age_vector = findgen(n_age)*age_bin + age_min
    metallicity_vector = [0.0004,0.0005,0.0007,0.0010,0.0013,0.0017,0.0023,0.0030,0.0040,0.0054,0.0073,0.0098,0.0130,0.0234,0.03]

    V_grid = fltarr(n_age,n_metallicity)
    vr_grid = fltarr(n_age,n_metallicity)
    ri_grid = fltarr(n_age,n_metallicity)
    bv_grid = fltarr(n_age,n_metallicity)
    vk_grid = fltarr(n_age,n_metallicity)
    ub_grid = fltarr(n_age,n_metallicity)
    openr, 1, gridpath+'data_leo/V_grid'
    openr, 2, gridpath+'data_leo/vr_grid'
    openr, 3, gridpath+'data_leo/ri_grid'
    openr, 4, gridpath+'data_leo/bv_grid'
    openr, 5, gridpath+'data_leo/vk_grid'
    openr, 6, gridpath+'data_leo/ub_grid'
ENDIF 
IF (isochrone eq 3) THEN BEGIN
;Beth Willam's Isochrones ---
    printf,'Beth Willams Isochrones'
    n_age = 18
    age_max = 10.2
    age_min = 8.5
    age_bin = (age_max-age_min)/(n_age - 1.)
    age_vector = findgen(n_age)*age_bin + age_min
    n_metallicity = 26          ;bw
    metal_max = .1
    metal_min = .0004
    metal_bin = (alog10(metal_max)-alog10(metal_min))/(n_metallicity - 1.)
    metallicity_vector = 10.^(findgen(n_metallicity)*metal_bin + alog10(metal_min))

    openr, 1, gridpath+'data_bw/V_grid'
    openr, 2, gridpath+'data_bw/vr_grid'
    openr, 3, gridpath+'data_bw/ri_grid'
    openr, 4, gridpath+'data_bw/bv_grid'
    openr, 5, gridpath+'data_bw/vk_grid'
    openr, 6, gridpath+'data_bw/ub_grid'
ENDIF

print,'SSP Range of Metallicities: ',MINMAX(metallicity_vector)
print,'SSP Range of Ages:          ',MINMAX(age_vector)

V_grid =  fltarr(n_metallicity,n_age)
vr_grid = fltarr(n_metallicity,n_age)
ri_grid = fltarr(n_metallicity,n_age)
bv_grid = fltarr(n_metallicity,n_age)
vk_grid = fltarr(n_metallicity,n_age)
ub_grid = fltarr(n_metallicity,n_age)

readf, 1, V_grid
readf, 2, vr_grid
readf, 3, ri_grid
readf, 4, bv_grid
readf, 5, vk_grid
readf, 6, ub_grid
close, /all
 
m_low = MIN(metallicity_vector)
m_hi = MAX(metallicity_vector)
a_low = MIN(age_vector)
a_hi = MAX(age_vector)
n_ind_m = 26
n_ind_a = 830       
m_spacing = (m_hi - m_low)/(n_ind_m-1)
m_vector = findgen(n_ind_m)*m_spacing + m_low
a_spacing = (a_hi - a_low)/(n_ind_a-1)
a_vector = findgen(n_ind_a)*a_spacing + a_low
V_grid2 = TRI_SURF(V_grid, xvalues =  metallicity_vector, yvalues = age_vector,nx = n_ind_m, ny = n_ind_a)
vr_grid2 = TRI_SURF(vr_grid, xvalues =  metallicity_vector, yvalues = age_vector,GS=[m_spacing, a_spacing], BOUNDS=[m_low,a_low,m_hi,a_hi])
bv_grid2 = TRI_SURF(bv_grid, xvalues =  metallicity_vector, yvalues = age_vector,GS=[m_spacing, a_spacing], BOUNDS=[m_low,a_low,m_hi,a_hi])
ri_grid2 = TRI_SURF(ri_grid, xvalues =  metallicity_vector, yvalues = age_vector,GS=[m_spacing, a_spacing], BOUNDS=[m_low,a_low,m_hi,a_hi])
vk_grid2 = TRI_SURF(vk_grid, xvalues =  metallicity_vector, yvalues = age_vector,GS=[m_spacing, a_spacing], BOUNDS=[m_low,a_low,m_hi,a_hi])
ub_grid2 = TRI_SURF(ub_grid, xvalues =  metallicity_vector, yvalues = age_vector,GS=[m_spacing, a_spacing], BOUNDS=[m_low,a_low,m_hi,a_hi])
mb_grid2 = TRI_SURF(mb_grid, xvalues =  metallicity_vector, yvalues = age_vector,GS=[m_spacing, a_spacing], BOUNDS=[m_low,a_low,m_hi,a_hi])

;-------------------;
; loop through each ;
; timestep          ;
;-------------------;

FOR j = 0, n_elements(tipsy_files)-1 do begin
    tipsy_file = tipsy_files[j]
    grp_file = grp_files[j]
    stat_file = stat_files[j]
    outfile = outfiles[j]
    outfile_stars = outfiles_stars[j]

    temp = FILE_SEARCH(tipsy_file)
    IF (strmatch(temp[0], tipsy_file)) THEN BEGIN
;Looping Through Timesteps
        RTIPSY, tipsy_file, h, gas, dark, star
        star = star[where(star.tform gt 0)] ;selects the stars but not the Black Holes
        h.nstar = n_elements(star)
        star_mags = replicate({age:0.0, metallicity:0.0, r:0.0, v:0.0, b:0.0, i:0.0, k:0.0, u:0.0, mb:0.0},h.nstar)
        IF (keyword_set(snaptime)) THEN maxtime = snaptime ELSE maxtime = max(star.tform)
        age_star = (maxtime - star.tform) * timestep
        star_mags.age = age_star
        star_mags.metallicity = star.metals 
        IF((where(age_star le 1))[0] ne -1) THEN age_star[where(age_star le 1)] = 1
        log_age = alog10(age_star)
        print, 'Range of SSP Metallicities:  ', MINMAX(star.metals)
        print, 'Range of SSP Ages:            ', MINMAX(log_age)
	print, ''
        IF (keyword_set(init_mass)) THEN log_mass = alog10(mass_unit*init_mass) ELSE log_mass = alog10(mass_unit*MAX(star.mass)) ;The original mass of the star particle (is about equal to the mass of the largest star particle)
        print,'Age of Simulation:               ',maxtime*timestep
        print,'Log Initial Mass of Star Particle (solar masses): ',log_mass
        print,''

        temp = where (star.metals ge m_hi)
        if temp[0] ne -1 then star(temp).metals = m_hi
        temp = where (star.metals le m_low)
        if temp[0] ne -1 then star(temp).metals = m_low
        temp = where (log_age ge a_hi)
        if temp[0] ne -1 then log_age(temp) = a_hi
        temp = where (log_age le a_low)
        if temp[0] ne -1 then log_age(temp) = a_low

;-----------------------;
; calc phot params for  ;
;   stars               ;
;-----------------------;
        metal_ref = (star.metals - m_low)/m_spacing
        metal_ref_diff = metal_ref - FIX(metal_ref)
        temp = where (metal_ref_diff ge 0.5*m_spacing)
        if temp[0] ne -1 then metal_ref[temp] = metal_ref[temp]+1
        metal_ref = FIX(metal_ref)

        age_ref = (log_age - a_low)/a_spacing
        age_ref_diff = age_ref - FIX(age_ref)
        age_temp = where (age_ref_diff ge 0.5*a_spacing)
        if age_temp[0] ne -1 then age_ref[age_temp] = age_ref[age_temp]+1
        age_ref = FIX(age_ref)

        star_mags.v = V_grid2[metal_ref,age_ref]
        star_mags.r = star_mags.v - vr_grid2[metal_ref,age_ref]
        star_mags.b = bv_grid2[metal_ref,age_ref] + star_mags.v
        star_mags.i = star_mags.r - ri_grid2[metal_ref,age_ref]
        star_mags.k = star_mags.v - vk_grid2[metal_ref,age_ref]
        star_mags.u = ub_grid2[metal_ref,age_ref] + star_mags.b

        IF(isochrone gt 1) THEN BEGIN
            star_mags.v = star_mags.v -  2.5*(log_mass)
            star_mags.r = star_mags.r -  2.5*(log_mass)
            star_mags.b = star_mags.b -  2.5*(log_mass)
            star_mags.i = star_mags.i -  2.5*(log_mass)
            star_mags.k = star_mags.k -  2.5*(log_mass)
            star_mags.u = star_mags.u -  2.5*(log_mass)
        ENDIF ELSE BEGIN  
            star_mags.mb = mb_grid2[metal_ref,age_ref]
;The Starburst99 isochrones are for a 10^6 solar mass star particle
            star_mags.v = star_mags.v -  2.5*(log_mass  - 6.0)
            star_mags.r = star_mags.r -  2.5*(log_mass  - 6.0)
            star_mags.b = star_mags.b -  2.5*(log_mass  - 6.0)
            star_mags.i = star_mags.i -  2.5*(log_mass  - 6.0)
            star_mags.k = star_mags.k -  2.5*(log_mass  - 6.0)
            star_mags.u = star_mags.u -  2.5*(log_mass  - 6.0)
            star_mags.mb = star_mags.mb -  2.5*(log_mass  - 6.0)
        ENDELSE
        IF (mag_sys EQ 'ab') THEN BEGIN
            print,'AB Magnitude System'
            star_mags.r = star_mags.r + 0.055
            star_mags.v = star_mags.v - 0.044
            star_mags.b = star_mags.b - 0.163
            star_mags.i = star_mags.i + 0.309
        ENDIF ELSE print, 'Vega Magnitude System'
;------------------------;
; Get Magnitudes of Each ;
; Halo                   ;
;------------------------;
        IF (keyword_set(multiple)) THEN BEGIN
            if halo_finder eq 'skid' then input_halos = READ_STAT_STRUC_SKID(stat_file)
            if halo_finder eq 'amiga' then input_halos = READ_STAT_STRUC_AMIGA(stat_file)
            openr, 2, grp_file
            grp_data = intarr(h.ndark+h.nstar+h.ngas+1)
            readf, 2, grp_data
            close, 2
            halo_index = grp_data[h.ngas+h.ndark+1:h.ndark+h.nstar+h.ngas]
            IF keyword_set(virial_radius) THEN BEGIN
                galid = input_halos.group
                radii = input_halos.rvir
            ENDIF
         ENDIF ELSE input_halos = replicate({xc:0, yc:0, zc:0,group:0},1)
        nhalos = n_elements(input_halos)
        halos = replicate({id:-1.,mb:0.0 ,u:0.0, b:0.0, v:0.0, r:0.0, i:0.0, k:0.0, rprime:0.0, iprime:0.0, Z:0.0, age:0.0, mass:0.0},nhalos)
        npartstar = fltarr(nhalos)
        mass = fltarr(nhalos)
        IF (keyword_set(RADIUSFILE)) THEN BEGIN
            readcol,RADIUSFILE,galid,radii,format = 'F3,F3'
        ENDIF
        
        FOR i = 0, nhalos - 1 DO BEGIN
            IF (keyword_set(radius_file) OR keyword_set(virial_radius)) THEN BEGIN
;                readcol, read_col
                ind = (where(galid EQ input_halos[i].group))[0]
 ;               print,'*** Galaxy ',input_halos[i].group,' *** ';
                IF (ind NE -1) THEN BEGIN
                    halo_member = where (sqrt((length_unit*star.x-input_halos[i].xc + length_unit/2)^2.0 + (length_unit*star.y - input_halos[i].yc + length_unit/2)^2.0 + (length_unit*star.z - input_halos[i].zc + length_unit/2)^2.0) LE radii[ind]/1000.0) 
;                    print,'Radius [kpc]: ',radii[ind]
;                    print,'20*V200/220 Halo: ',n_elements(halo_member)
;                    amiga_matches = where (halo_index EQ input_halos[i].group)
;                    print,'Amiga Halo:       ',n_elements(amiga_matches)
;                    IF(n_elements(amiga_matches)*2.0 LE n_elements(halo_member) OR (n_elements(halo_member)*2.0 le n_elements(amiga_matches))) THEN BEGIN
;                        print,"Very different number of stars found"
;                        halo_member = where (halo_index EQ input_halos[i].group)
;                    ENDIF
                 ENDIF ELSE BEGIN
                    halo_member = where (halo_index EQ input_halos[i].group)
;                    print,'Amiga Halo:       ',n_elements(halo_member)
                ENDELSE
            ENDIF ELSE BEGIN
                IF (keyword_set(multiple)) THEN halo_member = where (halo_index EQ input_halos[i].group) ELSE halo_member = lindgen(long(h.nstar))
            ENDELSE
            hnew = h
            hnew.nstar = n_elements(halo_member)
            npartstar[i] = n_elements(halo_member)
            IF n_elements(halo_member) GT 3 THEN BEGIN
                halos[i].id = input_halos[i].group
                mass[i] = TOTAL(10.^log_mass(halo_member))

                lumtemp = TOTAL(10.^(-0.4*star_mags(halo_member).k))
                Z_temp = TOTAL(star[halo_member].metals)/n_elements(halo_member)
                age_temp = alog10(TOTAL(10.^(log_age[halo_member]))/n_elements(halo_member))

                halos[i].Z = TOTAL(star[halo_member].metals * 10.^(-0.4*star_mags[halo_member].k))/lumtemp
                halos[i].age = alog10(TOTAL(10.^log_age[halo_member] * 10.^(-0.4*star_mags[halo_member].k))/lumtemp)
                halos[i].mass = mass[i]

                halos[i].r = -2.5*alog10(TOTAL(10.^(-0.4*star_mags[halo_member].r)))
                halos[i].v = -2.5*alog10(TOTAL(10.^(-0.4*star_mags[halo_member].v)))
                halos[i].b = -2.5*alog10(TOTAL(10.^(-0.4*star_mags[halo_member].b)))
                halos[i].i = -2.5*alog10(TOTAL(10.^(-0.4*star_mags[halo_member].i)))
                halos[i].k = -2.5*alog10(TOTAL(10.^(-0.4*star_mags[halo_member].k)))
                halos[i].u = -2.5*alog10(TOTAL(10.^(-0.4*star_mags[halo_member].u)))   
                halos[i].mb = -2.5*alog10(TOTAL(10.^(-0.4*star_mags[halo_member].mb)))  

                IF i LT 4 THEN BEGIN
                    print,'Magnitudes',i
                    print,'Mbol:',halos[i].mb
                    print,'U:',halos[i].u
                    print,'B:',halos[i].b
                    print,'V:',halos[i].v
                    print,'I:',halos[i].i
                    print,'R:',halos[i].r
                    print,'K:',halos[i].k
                ENDIF
            ENDIF
          ENDFOR

;-------------;
; in sdss mag ;
;-------------;
        temp  = where (halos.r - halos.i LT 1.15)
        temp2 = where (halos.r - halos.i GE 1.15)
        r_minus_i = fltarr(n_elements(halos))

        halos.rprime = (halos.v - .81*(halos.v - halos.r) + .13)
        IF temp[0]  NE -1 THEN r_minus_i[temp]  =      (halos[temp].r  - halos[temp].i)  - 0.21
        IF temp2[0] NE -1 THEN r_minus_i[temp2] = 1.42*(halos[temp2].r - halos[temp2].i) - 0.69
        halos.iprime = halos.rprime - r_minus_i

;-----------------------------;
; output magnitudes to a file ;
; one file for each timestep, ;
; one entry for each halo     ;
;-----------------------------;
        MWRFITS, halos, outfile, /create
        MWRFITS, star_mags, outfile_stars, /create
        ENDIF ELSE print, STRING('No matching tipsy file found for: ' + tipsy_file)
    ENDFOR
close,5
RETURN
END
