;.r highres_grp_stat.pro
;.r ~/code/procedures/wtipsy.pro
;filename = '/astro/net/scratch2/fabio/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512.D/h603.cosmo50cmb.2304g2bwK.00512'
;highres_grp_stat,filename,'.z0.000',18250,munit = 1.84793e16,vunit = 4.052167 ,h0 = 0.73
;highres_grp_stat,'/astro/net/scratch2/christensen/LG_dm/LG_dm.00512','.z0.000',63999.83,munit = 9.9669d16,vunit = 2211.1567 ,h0 = 0.73

; This program takes the output from AMIGA and creates .grp, .gtp and 
; .stat files similar to that produced by SKID. 
;
; Author: Alyson Brooks
; Last updated: Dec 6, 2006

pro highres_grp_stat, tipsyfile, redshift, boxsize, munit = munit, vunit = vunit, lunit = lunit, h0 = h0, multiplefile = multiplefile

; INPUTS
;
; tipsyfile is the name of the tipsy file for the run.  If you are not running 
; 	this code in the simulation directory where tipsyfile is located, then 
; 	tipsyfile must include the full directory path
;
; redshift: when AMIGA outputs .AHF_halos files .AHF_particles files, it appends 
; 	the redshift of the timestep to the base filename.  IMPORTANT: this code 
;	assumes that the output file name from AMIGA is simply the tipsy file name 
;	+ addenda from AMIGA.  If this is not the case, this code will not find your 
;	AMIGA results. 
; 	The redshift MUST be entered to three decimal places, just as AMIGA outputs.  
;	It MUST be ascii! Include the 'z', e.g., '.z0.000' is an appropriate input, 
;	with the quotes to designate ascii.  
;
; boxsize: this is the boxsize as found by AMIGA.  This should be the comoving 
;	length unit (boxsize) in kpc/h.   
;	IMPORTANT: in some cases, where an old version of tipsy2snap was used, this 
;	is not simply the length unit of the simulation.  If you are using .gadget 
;	files from before Sept 2007, (HIGHLY improbable, as they should have been 
;	deleted) verify the box length (grep boxsize: amigoutput.log)  
;
; munit: the mass unit of the simulation.  If no input is given here, munit will 
; 	default to 3.171d15 
;
; vunit: the velocity unit of the simulation.  If no input is given here, vunit will 
; 	default to 691.0 km/s
;
; h0: the value of h used in the simulation.  If no input is given, h0 will 
;	default to 0.7
;
; multiple: if you wish to run this code on more than one input, multiplefile is a 
;	file that contains a line for each input.  In this case, all inputs must 
;	be given (i.e., munit, vunit, and h0 are not optional inputs) for each 
;	input filex

if (n_params() eq 0) then begin
  print, "highres_grp_stat.pro -- Creates a .grp, .gtp, and .stat file from AMIGA outputs"
  print, " "
  print, "Usage: "
  print, "highres_grp_stat, tipsyfile, redshift, boxsize, munit = munit, vunit = vunit, h0 = h0, multiplefile = multiplefile"
  print, "Please read highres_grp_stat.pro for more info on the inputs"
  print, " "
  print, " "
endif

IF n_elements(multiplefile) EQ 0 THEN BEGIN
  ; Set mass unit and h0 if not given as inputs
  if n_elements(munit) eq 0 then munit = 3.171d15
  if n_elements(vunit) eq 0 then vunit = boxsize/2.894405;vunit = 691.0
  if n_elements(h0) eq 0 then h0 = 0.7  
ENDIF ELSE BEGIN
  readcol, multiplefile, tipsyfile, redshift, boxsize, munit, vunit, h0, format='a,a,d,d,d,d'
ENDELSE
  lunit = boxsize/h0

; start a big loop over the number of input simulations
FOR j=0,n_elements(tipsyfile)-1 DO BEGIN

rtipsy, tipsyfile[j], h,g,d,s
halosfile = strcompress(tipsyfile[j]+redshift[j]+'.AHF_halos')
particlefile = strcompress(tipsyfile[j]+redshift[j]+'.AHF_particles') 

; The format of the particlefile is a long array.  For each halo, the number 
; of particles in the halo is listed, and that many next lines contain the 
; index of the particles that belong to that halo.  When the first halo is 
; done, then the number of particles in the next halo is given and that many 
; next lines contain the indices to the particles in that halo, and so on.  So, 
; we must grab all the indices without treating npart of each halo as an index. 
readcol, halosfile, npart, nvpart, xc, yc, zc, vxc, vyc, vzc, mvir, rvir, vmax, rmax, sigv, spin, format='l,l,f,f,f,f,f,f,f,f,f,f,f,f'
;delta = nvpart-npart
;real = where(delta lt 100, comp=heavy)
real = where(npart ge nvpart) ; If nvpart GT npart, then the halo has definitely 
; been contaminated by low res (heavy) dm particles
; NOTE that if there are a significant number of baryonic particles (which have 
; low particle masses) in a halo that has been mildly contaminated by heavy dm 
; particles, it is still possible to have nvpart LE npart.  Or, if a low res halo 
; is purely dm, then nvpart = npart.  Thus, heavies must be searched for on a per 
; halo basis.  The .stat file contains a column to identify these contaminated 
; halos
heavy = where(nvpart GT npart) ; Positively contaminated

; **************************************************************************
; First, create the .grp file
; **************************************************************************
if h.n lt 4.2e6 then begin 
 rdfloat, particlefile, index
 index = long(index)
endif else readcol, particlefile, index, format='l'
init=0
startind = lonarr(n_elements(npart))
endind = lonarr(n_elements(npart))
for i=0L, n_elements(npart)-1 do begin
  addind = i+1 
  startind[i] = addind+init
  endind[i] = startind[i]+npart[i]-1 
  init = init+npart[i]
endfor

realstartind = startind(real)
realendind = endind(real)
IF heavy[0] ne -1 then begin
  heavy_startind = startind(heavy)
  heavy_endind = endind(heavy)
ENDIF

; Now we have the starting and ending point of each halo's particles in the 
; particlefile.  Let's put them into an array consecutively rather than 
; separated by npart.
particles = lonarr(n_elements(index))
grp = lonarr(n_elements(index))
for i=0L, n_elements(realstartind)-1 do begin
  particles[realstartind(i):realendind(i)] = index[realstartind(i):realendind(i)]
  grp[realstartind(i):realendind(i)] = i+1
endfor
IF heavy[0] ne -1 then begin
  for i=0L, n_elements(heavy_startind)-1 do begin
    particles[heavy_startind(i):heavy_endind(i)] = index[heavy_startind(i):heavy_endind(i)]
    grp[heavy_startind(i):heavy_endind(i)] = 0
  endfor
ENDIF
nonzero = where(particles ne 0)
particles = particles(nonzero)
grp = grp(nonzero) 

; Now put the particle indices in ascending order, and get the grp number 
; that each particle belongs to.
ntot = long(h.n)
sortind = sort(particles)
sortedgrp = grp(sortind)
totarray = lonarr(ntot) 
for i=0L,n_elements(sortind)-1 do begin 
  ind = particles(sortind[i])
  totarray[ind] = sortedgrp(i) 
endfor
; Totarray now contains the grp that each particle belongs to, and we 
; can write it to a .grp file.

; NOTE: some particles are assigned by AMIGA to more than one halo!  This
; means that some of the ind values that go into totarray get written over
; after an intial assignment of grp number.  This is ok if the halos are
; read in by order of mass, as that guarantees that the reassigned value
; is a substructure grp number.  However, AMIGA currently orders the halos
; by npart rather than my mass.  However, it is likely that the oddball
; halos that are out of mass order are lowres halos that accidentally made
; it into the catalog.  So, for now I will not rearrange the input by mass.
; This may become necessary at a later time, though.


; Write the .grp file
get_lun, lun
filename = strcompress(tipsyfile[j]+'.amiga.grp')
openw, lun, filename
printf, lun, ntot
for i=0L,n_elements(totarray)-1 do printf, lun, long(totarray[i])
close, lun
free_lun, lun

; **************************************************************************
; Now, create the .stat file
; **************************************************************************
stat = replicate({grp:0L, contam:strarr(1), false:strarr(1), n_part:0L, n_gas:0L, n_star:0L, n_dark:0L, m_vir:dblarr(1), r_vir:dblarr(1), gas_mass:dblarr(1), star_mass:dblarr(1), dark_mass:dblarr(1), v_max:dblarr(1), r_max:dblarr(1), sig_v:dblarr(1), x_c:dblarr(1), y_c:dblarr(1), z_c:dblarr(1), vx_c:dblarr(1), vy_c:dblarr(1), vz_c:dblarr(1), amiga_orig_id:0L}, n_elements(real))

sorted = sortedgrp(sort(sortedgrp))
unique = sorted(uniq(sorted))
; Take care of the easy stuff first. 
; By definition, there should not be grp 0 unless there are heavies, which are given grp=0
; See lines 127-128
IF heavy[0] ne -1 then begin
for i=0L,n_elements(real)-1 do stat[i].grp = unique[i+1]
endif else begin
for i=0L,n_elements(real)-1 do stat[i].grp = unique[i]
endelse
stat.n_part = npart(real)
stat.amiga_orig_id = real 
for i=0L,n_elements(real)-1 do stat[i].m_vir = double(mvir(real[i])/h0[j]) 
for i=0L,n_elements(real)-1 do stat[i].r_vir = double(rvir(real[i])/h0[j])
for i=0L,n_elements(real)-1 do stat[i].v_max = double(vmax(real[i]))
for i=0L,n_elements(real)-1 do stat[i].r_max = double(rmax(real[i])/h0[j])
for i=0L,n_elements(real)-1 do stat[i].sig_v = double(sigv(real[i]))
for i=0L,n_elements(real)-1 do stat[i].x_c = double(xc(real[i])/h0[j])
for i=0L,n_elements(real)-1 do stat[i].y_c = double(yc(real[i])/h0[j])
for i=0L,n_elements(real)-1 do stat[i].z_c = double(zc(real[i])/h0[j])
for i=0L,n_elements(real)-1 do stat[i].vx_c = double(vxc(real[i]))
for i=0L,n_elements(real)-1 do stat[i].vy_c = double(vyc(real[i]))
for i=0L,n_elements(real)-1 do stat[i].vz_c = double(vzc(real[i]))
; Now figure out info on gas and stars, if they exist!
IF h.ngas+h.nstar GT 0 THEN BEGIN 
  gmass = double(g.mass[0])
  if h.nstar GT 0 then smass = double(s.mass[0])
  dmass = double(d.mass[0])
  FOR i=0L, n_elements(real)-1 do begin
    particles_considered = index[realstartind(i):realendind(i)]
    gas_ind = where(particles_considered LT h.ngas)
    if gas_ind[0] NE -1 then stat[i].n_gas = n_elements(gas_ind) else stat[i].n_gas = 0
    if gas_ind[0] NE -1 then stat[i].gas_mass = total(gmass[particles_considered(gas_ind)]) else stat[i].gas_mass = 0
; stat[i].gas_mass = 0
  IF h.nstar GT 0 then begin
    star_ind = where(particles_considered GE h.ngas+h.ndark and particles_considered LT h.ngas+h.ndark+h.nstar)
    if star_ind[0] NE -1 then stat[i].n_star = n_elements(star_ind) else stat[i].n_star = 0 
    if star_ind[0] NE -1 then stat[i].star_mass = total(smass[particles_considered(star_ind)-(h.ngas+h.ndark)]) else stat[i].star_mass = 0 
; stat[i].star_mass = 0
  ENDIF else stat[i].star_mass = 0
    dark_ind = where(particles_considered GE h.ngas and particles_considered LT h.ngas+h.ndark)
    if dark_ind[0] NE -1 then begin
	stat[i].n_dark = n_elements(dark_ind);else stat[i].n_dark = 0 
	stat[i].dark_mass = total(dmass[particles_considered(dark_ind)-h.ngas]);else stat[i].dark_mass = 0 
	test = min(dmass[particles_considered(dark_ind)-h.ngas])
	test2 = where(dmass[particles_considered(dark_ind)-h.ngas] GT test)
	if test2[0] NE -1 then stat[i].contam = 'contam' else stat[i].contam = 'clean'
    endif else begin
	stat[i].n_dark = 0
	stat[i].dark_mass = 0
        stat[i].contam = 'false'
    endelse	
  ENDFOR
ENDIF ELSE BEGIN
  dmass = double(d.mass[0])
  FOR i=0L, n_elements(real)-1 do begin
    particles_considered = index[realstartind(i):realendind(i)]
    stat[i].gas_mass = 0
    stat[i].n_gas = 0
    stat[i].star_mass = 0
    stat[i].n_star = 0
    dark_ind = where(particles_considered GE h.ngas and particles_considered LT h.ngas+h.ndark)
;   if dark_ind[0] NE -1 then 
    stat[i].n_dark = n_elements(dark_ind) 
; else stat[i].n_dark = 0 
;   if dark_ind[0] NE -1 then 
    stat[i].dark_mass = total(dmass[particles_considered(dark_ind)-h.ngas]) 
    test = min(dmass[particles_considered(dark_ind)-h.ngas])
    test2 = where(dmass[particles_considered(dark_ind)-h.ngas] GT test)
    if test2[0] NE -1 then stat[i].contam = 'contam' else stat[i].contam = 'clean'
; else stat[i].dark_mass = 0 
  ENDFOR
ENDELSE

; This next is designed to detect if the main halo has been broken down into smaller 
; false halos at the center, a common problem.
test_xc = stat.x_c-stat[0].x_c[0]
test_yc = stat.y_c-stat[0].y_c[0]
test_zc = stat.z_c-stat[0].z_c[0]
test_xc = test_xc[1:n_elements(real)-1]
test_yc = test_yc[1:n_elements(real)-1]
test_zc = test_zc[1:n_elements(real)-1]
distance = ((test_xc^2.+test_yc^2.+test_zc^2.)^0.5)*1000. ;change from Mpc to kpc
pot_false = where(distance lt 0.025*stat[0].r_vir[0], nfalse)
if nfalse ne 0 then begin 
    stat[pot_false+1].false = 'false' ;subtracted off 1 to ignore main halo, now add back
    stat[0].false = 'false?'
endif


; Write the .stat file
get_lun, lun
filename = strcompress(tipsyfile[j]+'.amiga.stat')
openw, lun, filename
 printf, lun, format='(A5,2x,A7,2x,A7,2x,A7,2x,A7,2x,A16,2x,A10,2x,A16,2x,A16,2x,A16,2x,A10,2x,A9,2x,A9,2x,A9,2x,A9,2x,A9,2x,A13,2x,A13,2x,A13,2x,A5,2x,A7,2x,A7)', 'Grp', 'N_tot', 'N_gas', 'N_star', 'N_dark', 'Mvir(M_sol)', 'Rvir(kpc)', 'GasMass(M_sol)', 'StarMass(M_sol)', 'DarkMass(M_sol)', 'V_max', 'R@V_max', 'VelDisp', 'Xc', 'Yc', 'Zc', 'VXc', 'VYc', 'VZc', 'ID_A', 'Contam', 'False?'

for i=0L,n_elements(real)-1 do printf, lun, format='(I5,2x,I7,2x,I7,2x,I7,2x,I7,2x,A16,2x,D10.4,2x,A16,2x,A16,2x,A16,2x,D10.5,2x,D9.5,2x,D9.5,2x,D9.5,2x,D9.5,2x,D9.5,2x,D13.7,2x,D13.7,2x,D13.7,2x,I5,2x,A7,2x,A7)', stat[i].grp, stat[i].n_part, stat[i].n_gas, stat[i].n_star, stat[i].n_dark, stat[i].m_vir, stat[i].r_vir , stat[i].gas_mass*munit[j], stat[i].star_mass*munit[j], stat[i].dark_mass*munit[j], stat[i].v_max, stat[i].r_max, stat[i].sig_v, stat[i].x_c, stat[i].y_c, stat[i].z_c, stat[i].vx_c, stat[i].vy_c, stat[i].vz_c, stat[i].amiga_orig_id, stat[i].contam, stat[i].false

close, lun
free_lun, lun

; **************************************************************************
; Now, create the .gtp file
; **************************************************************************

statfile = strcompress(tipsyfile[j]+'.amiga.stat')
readcol, statfile, grp, npart, ngas, nstar, ndark, mvir, rvir, gmass, smass, dmass, vmax, rmax, sigv, xc, yc, zc, vxc, vyc, vzc, format='l,l,l,l,l,d,d,d,d,d,d,d,d,d,d,d,d,d,d'
  
star = replicate({mass:0.0, x:0.0, y:0.0, z:0.0, vx:0.0, vy:0.0, vz:0.0, metals:0.0, tform:0.0, eps:0.0, phi:0.0},n_elements(grp))
  for i=0L,n_elements(grp)-1 do begin
	star[i].mass = (gmass[i]+smass[i]+dmass[i])/munit[j]
	star[i].x = (xc[i]*1000./(lunit[j]))-0.5
	star[i].y = (yc[i]*1000./(lunit[j]))-0.5
	star[i].z = (zc[i]*1000./(lunit[j]))-0.5
	star[i].vx = vxc[i]/vunit[j]
	star[i].vy = vyc[i]/vunit[j]
	star[i].vz = vzc[i]/vunit[j]
	star[i].eps = rvir[i]/lunit[j]
  endfor
  star[*].metals = 0.0  
  star[*].phi = 0.0  
  star[*].tform = h.time
  header = h
  header.n = n_elements(grp)
  header.nstar = n_elements(grp)
  header.ngas = 0
  header.ndark = 0
  outfile = strcompress(tipsyfile[j]+'.amiga.gtp')
  wtipsy, outfile, header, gas, dark, star, /standard

ENDFOR

end

