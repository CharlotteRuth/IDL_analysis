;******************************************************************************
;This function analyzes the history of gas and star particles by
;taking a list of input particles and outputing a .fits file containing
;their properties at a given timestep.  
;
;REQUIREMENTS:
;
;     This procedure requires a minimum of one array if only gas is 
;     being traced (orig_iord), but three arrays if stars or gas+stars
;     are being traced (orig_iord, ind_stars, igasord_stars)  
;
;     orig_iord: the iord values of the objects you'd like to track, as 
;     given in the *iord files.  If you are tracing both gas and stars, 
;     this array will include the iord values for all particles 
;
;     igasord_stars: if you are tracing star particles, this is the 
;     igasorder value of those stars.  If you are tracing both gas and 
;     stars simultaneously, this array is ONLY for the stars.  If you 
;     are tracing only gas, do not include this input
;
;     ind_stars: If you are tracing stars (or gas+stars), this array 
;     is the index (the order in the list) of the STARS ONLY.  That is,
;     this array has a minimum value of 0 and maximum value of h.nstar.
;     Ignore h.ngas+h.ndark.  This is not the index to the entire 
;     (gas+dm+star) list, but ONLY the star list.  This is true even 
;     if you are tracing gas and stars simultaneously.  If you are 
;     tracing gas only, do not include this input
;
;     The program also deals with AMIGA .grp outputs (.stat could be 
;     added).  The output structure will contain the ID of the halo  
;     (history.haloid) that a star/gas particle belonged to at that 
;     timestep.  (Thus, AMIGA should be run on the timestep of choice 
;     before running this script.)
;
;     The program writes a fits file for each timestep using the 
;     name given in the simoutputlist file plus a unique string 
;     'tonamefiles' to identify a specific trace from others (e.g., gas 
;     vs stars of a MW type galaxy, traced separately).
;
; OPTIONAL KEYWORDS
;     The above assumes that you are tracing stars or gas particles that
;     still exist at the final timestep.  In some cases, you may wish to 
;     follow the history of gas particles that get deleted by the final 
;     step.  In that case, use keyword DEL.  This will write the properties 
;     of particles as in the standard trace, but when a particle no longer
;     exists it is left as zeros.
;
;     A number of our AMIGA outputs have found halos down to a minimum of 
;     16 particles.  However, our mass function may only be resolved down 
;     to 64 particles.  If you'd like to set a minimum number of particles 
;     for AMIGA to use, set MINPART to the minimum number of particles. 
;     If keyword TYPE isn't also set, then MINPART will require that there be 
;     a minimum number of MINPART dm particles.  Other options are 'gas', 'star', 
;     'baryon', and 'tot.'
;
;     This program can trace particles from multiple halos and split
;     them up in the output. If you would like to have it do so,
;     include GRP_ARR equal to an array of the same length as the
;     particles you would like to trace with each value equal to the
;     final grp value of the particle
;******************************************************************************
pro find_history, simoutputlist, tonamefiles, tipsyunitsfile, ORIG_IORD=orig_iord, IND_STARS=ind_stars, IGASORD_STARS=igasord_stars, MINPART=minpart, DEL=del, TYPE=type, GRP_ARR=GRP_ARR
; tonamefiles is a string to append to the output files, so multiple searches can 
; be done on multiple criteria and kept separate.
;
; tipsyunitsfile is tipsy.units.idl: this is the same as tipsy.units, but with the 
; ()'s removed so that idl can read it 
; 
; ind_stars should be an array at starts at h.ngas+h.ndark = 0 (that is, the maximum 
; length of ind_stars is the length of h.nstar for the timesteps they are pulled
; from). 
;
; orig_iord includes the iord values of both gas and star particles, if both are being 
; traced
;
; igasord_stars is the igasorder values for stars being traced 
;
; The number of elements in igasord_star and ind_stars should be equal, but 
; orig_iord will be longer if you are tracing both gas and stars  


readcol, tipsyunitsfile, lengthunit, massunit, velunit, format='d,d,f'
lengthunit = lengthunit[0]
massunit = massunit[0]
velunit = velunit[0]
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
H_per_gm = 5.9753790d23
densityunit = massunit * gm_per_msol * H_per_gm/lengthunit^3/cm_per_kpc^3
root = simoutputlist 

;get number of timesteps and desired particles
ntimesteps=N_ELEMENTS(root)
norig=N_ELEMENTS(orig_iord)

;create structure to hold information of particles over time:
;history structure is subscripted by number of desired particles, many
;elements within it are subscripted by number of timesteps
history=REPLICATE({iord:0l, igasorder:0l, mark:0l, mass:DBLARR(ntimesteps), x:DBLARR(ntimesteps), y:DBLARR(ntimesteps), z:DBLARR(ntimesteps), vx:DBLARR(ntimesteps), vy:DBLARR(ntimesteps), vz:DBLARR(ntimesteps), rho:DBLARR(ntimesteps), temp:DBLARR(ntimesteps), metallicity:DBLARR(ntimesteps), haloid:LONARR(ntimesteps)}, norig)
grp_arr_final = fltarr(norig)

history.iord=orig_iord
orig_iord = 0 ;Conserve memory
IF n_elements(igasord_stars) NE 0 THEN history.igasorder=igasord_stars 
FOR i = 0l, ntimesteps - 1 DO BEGIN
  haloind=read_lon_array(root[i]+'.amiga.grp')
  IF keyword_set(MINPART) THEN BEGIN
      readcol, root[i]+'.amiga.stat', grp, ntot, ngas, nstar, ndark, format='l,l,l,l,l', /silent
      IF keyword_set(TYPE) THEN BEGIN
;For all particles in halos that are smaller than minpart, set their
;grp to zero
         IF type EQ 'tot' THEN mingrp = grp(where(ntot LT minpart, ngrp)) 
         IF type EQ 'gas' THEN mingrp = grp(where(ngas LT minpart, ngrp)) 
         IF type EQ 'star' THEN mingrp = grp(where(ngas LT minpart, ngrp)) 
         IF type EQ 'baryon' THEN mingrp = grp(where(ngas+nstar LT minpart, ngrp)) 
         IF type EQ 'dm' THEN mingrp = grp(where(ndark LT minpart, ngrp)) 
         IF ngrp NE 0 THEN BEGIN
            FOR j=0,n_elements(mingrp)-1 DO BEGIN
               ind = where(haloind EQ mingrp[j])
               haloind[ind] = 0
            ENDFOR
         ENDIF
      ENDIF ELSE BEGIN
         test = where(ndark LT minpart, ngrp)
         IF ngrp NE 0 THEN BEGIN
            mingrp = grp(where(ndark LT minpart, ngrp)) 
            FOR j=0,n_elements(mingrp)-1 DO BEGIN
               ind = where(haloind EQ mingrp[j])
               haloind[ind] = 0
            ENDFOR
         ENDIF
      ENDELSE
  ENDIF

  ;Read in valuable information from the tipsy file
  print,root[i]
  rtipsy, root[i], h,g,d,s
;Conserve memory by getting rid of dark matter particles
;  s = s[where(s.tform GT 0)] Deleting black hole particles causes
;  problems when matching iord, oxmassfrac, femassfrac
  ;Read in the .iord file.
  iord=read_lon_array(root[i]+'.iord')
  
  IF file_test(root[i]+'.OxMassFrac') THEN BEGIN
     read_tipsy_arr,root[i]+'.OxMassFrac',h,Ox,type='float'
     ox_gas = ox[0:h.ngas-1]
     ox_star = ox[h.ngas+h.ndark:h.n-1]
  ENDIF ELSE BEGIN
     ox_gas = fltarr(h.ngas)
     ox_star = fltarr(h.nstar)
  ENDELSE

  IF file_test(root[i]+'.FeMassFrac') THEN BEGIN
     read_tipsy_arr,root[i]+'.FeMassFrac',h,Fe,type='float'
     fe_gas = fe[0:h.ngas-1]
     fe_star = fe[h.ngas+h.ndark:h.n-1]
  ENDIF ELSE BEGIN
     fe_gas = fltarr(h.ngas)
     fe_star = fltarr(h.nstar)
 ENDELSE

  haloind_gas = haloind[0:h.ngas-1]
  haloind_star = haloind[h.ngas+h.ndark:h.n-1]
 
  nstars = n_elements(ind_stars)
  ngas = norig-nstars
  IF ngas NE 0 THEN BEGIN
    orig_iord_gas = history[0:ngas-1].iord
    gas = indgen(ngas, /long)
  ENDIF
  IF nstars NE 0 THEN BEGIN
    orig_iord_star = history[ngas:norig-1].iord
    stars = indgen(nstars, /long)+ngas
  ENDIF
;  stop
  ;For stars,
  IF nstars GT 0 THEN BEGIN
    ; Find which are still stars and which are gas at this step
    gasprog = where(orig_iord_star GT max(iord), nprog, comp=oldstars, ncomplement=nold)
    ind2stars = ind_stars[oldstars]	;Found the stars, now find gas indices
    IF nprog NE 0 THEN BEGIN
       ind2gas = ind_stars[gasprog] ;For indexing into the history array
       igasorder = igasord_stars[gasprog]
       iord_star = iord[h.ngas+h.ndark:h.n-1]
       iord=iord[0:h.ngas]	;shorten for faster reading
;       gasind = FINDEX(iord, igasorder) 
       match2,iord,igasorder,temp,gasind
       IF n_elements(where(gasind NE -1)) NE n_elements(igasorder) THEN stop
    ENDIF
    ;Now fill up the history structure
    ;first for stars
    IF nold NE 0 THEN BEGIN
        history[stars(0:nold-1)].mark[i]=long(h.ngas+h.ndark+ind2stars)+1
;    FOR j=0L,nold-1 DO BEGIN
        history[ngas : ngas + nold - 1].mass[i]=s[ind2stars].mass*massunit
        history[ngas : ngas + nold - 1].x[i]=s[ind2stars].x*lengthunit*h.time
        history[ngas : ngas + nold - 1].y[i]=s[ind2stars].y*lengthunit*h.time
        history[ngas : ngas + nold - 1].z[i]=s[ind2stars].z*lengthunit*h.time
        history[ngas : ngas + nold - 1].vx[i]=s[ind2stars].vx*velunit*h.time
        history[ngas : ngas + nold - 1].vy[i]=s[ind2stars].vy*velunit*h.time
        history[ngas : ngas + nold - 1].vz[i]=s[ind2stars].vz*velunit*h.time
        history[ngas : ngas + nold - 1].rho[i]=0.
        history[ngas : ngas + nold - 1].temp[i]=0.
        history[ngas : ngas + nold - 1].metallicity[i]=2.09*ox_star[ind2stars]+1.06*fe_star[ind2stars] ;s[ind2stars].metals
        history[ngas : ngas + nold - 1].haloid[i]=haloind_star[ind2stars]
        IF keyword_set(grp_arr) THEN grp_arr_final[ngas : ngas + nold - 1] = grp_arr[oldstars] ELSE grp_arr_final[ngas : ngas + nold - 1] = 1
 ;   ENDFOR
;    FOR j=0L,nold-1 DO BEGIN
;       history[stars(j)].mass[i]=s[ind2stars(j)].mass*massunit 
;       history[stars(j)].x[i]=s[ind2stars(j)].x*lengthunit*h.time
;       history[stars(j)].y[i]=s[ind2stars(j)].y*lengthunit*h.time
;       history[stars(j)].z[i]=s[ind2stars(j)].z*lengthunit*h.time
;       history[stars(j)].vx[i]=s[ind2stars(j)].vx*velunit*h.time
;       history[stars(j)].vy[i]=s[ind2stars(j)].vy*velunit*h.time
;       history[stars(j)].vz[i]=s[ind2stars(j)].vz*velunit*h.time
;       history[stars(j)].rho[i]=0.
;       history[stars(j)].temp[i]=0.
;       history[stars(j)].metallicity[i]=s[ind2stars(j)].metals
;       history[stars(j)].haloid[i]=haloind[ind2stars(j)+h.ngas+h.ndark]
;    ENDFOR
 ;   ;now for stars that are still gas at this step
    ENDIF
    IF nprog NE 0 THEN BEGIN
       history[stars(nold:nstars-1)].mark[i]=long(gasind)+1
       history[ngas + nold :ngas + nold + nprog - 1].mass[i]=g[gasind].mass*massunit
       history[ngas + nold :ngas + nold + nprog - 1].x[i]=g[gasind].x*lengthunit*h.time
       history[ngas + nold :ngas + nold + nprog - 1].y[i]=g[gasind].y*lengthunit*h.time
       history[ngas + nold :ngas + nold + nprog - 1].z[i]=g[gasind].z*lengthunit*h.time
       history[ngas + nold :ngas + nold + nprog - 1].vx[i]=g[gasind].vx*velunit*h.time
       history[ngas + nold :ngas + nold + nprog - 1].vy[i]=g[gasind].vy*velunit*h.time
       history[ngas + nold :ngas + nold + nprog - 1].vz[i]=g[gasind].vz*velunit*h.time
       history[ngas + nold :ngas + nold + nprog - 1].rho[i]=g[gasind].dens*densityunit/h.time^3
       history[ngas + nold :ngas + nold + nprog - 1].temp[i]=g[gasind].tempg
       history[ngas + nold :ngas + nold + nprog - 1].metallicity[i]=2.09*ox_gas[gasind]+1.06*fe_gas[gasind] ;g[gasind].zmetal
       history[ngas + nold :ngas + nold + nprog - 1].haloid[i]=haloind_gas[gasind]
       IF keyword_set(grp_arr) THEN grp_arr_final[ngas + nold :ngas + nold + nprog - 1] = grp_arr[gasprog] ELSE grp_arr_final[ngas + nold :ngas + nold + nprog - 1] = 1
;       FOR j=0l,nprog-1 DO BEGIN
;          history[stars(j+nold)].mass[i]=g[gasind(j)].mass*massunit
;          history[stars(j+nold)].x[i]=g[gasind(j)].x*lengthunit*h.time 
;          history[stars(j+nold)].y[i]=g[gasind(j)].y*lengthunit*h.time
;          history[stars(j+nold)].z[i]=g[gasind(j)].z*lengthunit*h.time
;          history[stars(j+nold)].vx[i]=g[gasind(j)].vx*velunit*h.time
;          history[stars(j+nold)].vy[i]=g[gasind(j)].vy*velunit*h.time
;          history[stars(j+nold)].vz[i]=g[gasind(j)].vz*velunit*h.time
;          history[stars(j+nold)].rho[i]=g[gasind(j)].dens*densityunit/h.time^3
;          history[stars(j+nold)].temp[i]=g[gasind(j)].tempg
;          history[stars(j+nold)].metallicity[i]=g[gasind(j)].zmetal
;          history[stars(j+nold)].haloid[i]=haloind[gasind(j)]
;       ENDFOR     
    ENDIF
  ENDIF 

  IF ngas GT 0 THEN BEGIN
    iord=iord[0:h.ngas]	;shorten for faster reading
;    gasind = binfind(iord, orig_iord_gas) 
    match2,iord,orig_iord_gas,temp,gasind
;    stop
    del = where(gasind EQ -1)
    exist = where(gasind NE -1)
    gasind = gasind[where(gasind NE -1)]
    history[exist].mark[i]=long(gasind)+1
    IF keyword_set(del) THEN BEGIN
       history[exist].mass[i]=g[gasind].mass*massunit
       history[exist].x[i]=g[gasind].x*lengthunit*h.time
       history[exist].y[i]=g[gasind].y*lengthunit*h.time
       history[exist].z[i]=g[gasind].z*lengthunit*h.time
       history[exist].vx[i]=g[gasind].vx*velunit*h.time
       history[exist].vy[i]=g[gasind].vy*velunit*h.time
       history[exist].vz[i]=g[gasind].vz*velunit*h.time
       history[exist].rho[i]=g[gasind].dens*densityunit/h.time^3
       history[exist].temp[i]=g[gasind].tempg
       history[exist].metallicity[i]=2.09*ox_gas[gasind]+1.06*fe_gas[gasind];g[gasind].zmetal
       history[exist].haloid[i]=haloind_gas[gasind]
       IF keyword_set(grp_arr) THEN grp_arr_final = grp_arr ELSE grp_arr_final[exist] = 1
;        FOR j=0L,n_elements(gasind)-1 DO BEGIN
;            history[exist(j)].mass[i]=g[gasind(j)].mass*massunit
;            history[exist(j)].x[i]=g[gasind(j)].x*lengthunit*h.time
;            history[exist(j)].y[i]=g[gasind(j)].y*lengthunit*h.time
;            history[exist(j)].z[i]=g[gasind(j)].z*lengthunit*h.time
;            history[exist(j)].vx[i]=g[gasind(j)].vx*velunit*h.time
;            history[exist(j)].vy[i]=g[gasind(j)].vy*velunit*h.time
;            history[exist(j)].vz[i]=g[gasind(j)].vz*velunit*h.time
;            history[exist(j)].rho[i]=g[gasind(j)].dens*densityunit/h.time^3
;            history[exist(j)].temp[i]=g[gasind(j)].tempg
;            history[exist(j)].metallicity[i]=g[gasind(j)].zmetal
;            history[exist(j)].haloid[i]=haloind[gasind(j)]
;        ENDFOR
    ENDIF ELSE BEGIN
       history.mass[i]=g[gasind].mass*massunit
       history.x[i]=g[gasind].x*lengthunit*h.time
       history.y[i]=g[gasind].y*lengthunit*h.time
       history.z[i]=g[gasind].z*lengthunit*h.time
       history.vx[i]=g[gasind].vx*velunit*h.time
       history.vy[i]=g[gasind].vy*velunit*h.time
       history.vz[i]=g[gasind].vz*velunit*h.time
       history.rho[i]=g[gasind].dens*densityunit/h.time^3
       history.temp[i]=g[gasind].tempg
       history.metallicity[i]=2.09*ox_gas[gasind]+1.06*fe_gas[gasind];g[gasind].zmetal
       history.haloid[i]=haloind_gas[gasind]
       IF keyword_set(grp_arr) THEN grp_arr_final = grp_arr[gasind] ELSE grp_arr_final = fltarr(n_elements(history)) + 1
;        FOR j=0L,ngas-1 DO BEGIN
;            history[j].mass[i]=g[gasind(j)].mass*massunit
;            history[j].x[i]=g[gasind(j)].x*lengthunit*h.time
;            history[j].y[i]=g[gasind(j)].y*lengthunit*h.time
;            history[j].z[i]=g[gasind(j)].z*lengthunit*h.time
;            history[j].vx[i]=g[gasind(j)].vx*velunit*h.time
;            history[j].vy[i]=g[gasind(j)].vy*velunit*h.time
;            history[j].vz[i]=g[gasind(j)].vz*velunit*h.time
;            history[j].rho[i]=g[gasind(j)].dens*densityunit/h.time^3
;            history[j].temp[i]=g[gasind(j)].tempg
;            history[j].metallicity[i]=g[gasind(j)].zmetal
;            history[j].haloid[i]=haloind[gasind(j)]       
; ;           IF j EQ 7175807 THEN stop
 ;       ENDFOR
    ENDELSE
ENDIF
  halos = grp_arr_final[uniq(grp_arr_final,sort(grp_arr_final))]
  IF n_elements(halos) NE n_elements(tonamefiles) THEN BEGIN
      print,'Wrong number of files to write to given for the number of halos'
      stop
  ENDIF
;  plot,history.x,history.y,psym = 3
  halocolors = (alog10(grp_arr_final + 2))/max(alog10(grp_arr_final + 2))*254
  FOR iwrite = 0, n_elements(tonamefiles) - 1 DO BEGIN
      outfile = root[i]+'.'+tonamefiles[iwrite]+'.history.fits'
      mwrfits, history[where(grp_arr_final EQ halos[iwrite])], outfile, /create
      prehalos = reform(history[where(grp_arr_final EQ halos[iwrite])].haloid[i])
      print,'Halo ',halos[iwrite],': ',prehalos[uniq(prehalos,sort(prehalos))]
;      oplot,history[where(grp_arr_final EQ halos[iwrite])].x[i],history[where(grp_arr_final EQ halos[iwrite])].y[i],psym = 3,color = (iwrite + 1)*20
;      stop
  ENDFOR
  
ENDFOR
orig_iord = history.iord
END

