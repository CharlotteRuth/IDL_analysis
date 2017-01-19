;Charlotte Christensen
;5/25/12
;This program uses grp1.allgas.entropy.fits and alignment.txt to
;establish when the gas is in the main halo and when it is in the disk
;of the main halo

;run from the main simulation directory

;  dir = '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK'
;  dir = '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK'
;  dir = '/astro/store/nbody3/christensen/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK'
  
;  dir = '/astro/store/nbody3/christensen/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK'

PRO gas_gal_disk_history,dir = dir,finalid = finalid,laststep = laststep,centralhalo = centralhalo,plots = plots,nowrite = nowrite,outplot = outplot,all_eject = all_eject,rvircut_half = rvircut_half,rvircut_fifth = rvircut_fifth,debug = debug,steps_debug = steps_debug,Hcut = Hcut,cooloncut = cooloncut
  formatplot,outplot = outplot
  device, decomposed=0
  loadct,39
  debug_cum = 1

  IF keyword_set(rvircut_fifth) OR keyword_set(rvircut_half) OR keyword_set(all_eject) THEN limitedwrite = 1 ELSE limitedwrite = 0
  CASE 1 OF
     keyword_set(all_eject): BEGIN
        ejectext = 'reeject_all'
        accrext  = 'reaccrdisk_all'
     END
     keyword_set(rvircut_half): BEGIN
        ejectext = 'reeject_rvir0.5'
        accrext  = 'reaccrdisk_rvir0.5'
     END
     keyword_set(rvircut_fifth): BEGIN
        ejectext = 'reeject_rvir0.2'
        accrext  = 'reaccrdisk_rvir0.2'
     END
     ELSE:  BEGIN
        ejectext = 'reeject'
        accrext  = 'reaccrdisk'
     END
  ENDCASE

  IF keyword_set(cooloncut) THEN coolext = 'coolon_' ELSE coolext = ''
;-------------------- Read in simulation data for final step
  IF NOT keyword_set(finalid) THEN finalid = '1' 
  IF NOT keyword_set(laststep) THEN laststep = '512'
  IF NOT keyword_set(dir) THEN dir = '.'
  cd,dir
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.coolontime',file_coolon
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord
;  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.OxMassFrac',file_ox
;  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.FeMassFrac',file_fe
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file
  rtipsy,file[0],h,g,d,s,/justhead
  readarr,file_coolon[0],h,coolon,part = 'gas',/ascii
  readarr,file_iord[0],h,iord,part = 'gas',/ascii,type = 'long' 
;  readarr,file_ox[0],h,oxmassfrac,part = 'gas',/ascii
;  readarr,file_fe[0],h,femassfrac,part = 'gas',/ascii
;  stop
;  readarr,'',
  spawn,'ls ' + dir + '/h*.param',pfile
  units = tipsyunits(pfile[0])

;-------------------- Read in the information on the halos
  halodat = mrdfits('grp' + finalid + '.alignment.fits',1)
  nsteps = n_elements(halodat.file)
  halodat.xc = halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.yc = halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.zc = halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.vxc = halodat.vxc*units.vunit
  halodat.vyc = halodat.vyc*units.vunit
  halodat.vzc = halodat.vzc*units.vunit
  center = [[halodat.time],[halodat.z],[halodat.xc],[halodat.yc],[halodat.zc]]
  center[*,2] = center[*,2]*1000.0 - units.lengthunit/2.0 ;"center" is in kpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
  center[*,3] = center[*,3]*1000.0 - units.lengthunit/2.0
  center[*,4] = center[*,4]*1000.0 - units.lengthunit/2.0
  a = [[halodat.xa],[halodat.ya],[halodat.za]]
  diskmass = fltarr(n_elements(halodat.file))

;-------------------- Read in the information about the gas particles at each timestep
  gpart = mrdfits('grp' + finalid + '.allgas.entropy.fits',1)
  print,'Read grp' + finalid + '.allgas.entropy.fits'
  npart = n_elements(gpart)

;-------------------- Need to exclude 'early' particles to compare to the outflow files
;  early = mrdfits('early.iord.fits',0)
;  all = gpart.iord
;  test = binfind(all,early)
;  early = 0           ;delete early array
;  all[test] = -1      ;mark all particles in the early array
;  exclude = where(all eq -1, comp=keep) ;find all particles in the early array
;  exclude = 0         ;delete exclude array
;  gpart = gpart[keep] ;gpart now only contains the particles that were not accreted early

;-------------------- Make sure the halo information and gpart are for equivalent outputs
  IF NOT nsteps EQ n_elements(gpart[0].mass) THEN BEGIN
     print,'Files are for different time steps'
     stop
  ENDIF

  indiskarr = intarr(npart) ; List of gas particles that are ever in the disk

  IF keyword_set(Hcut) $
  THEN gloc_filename = 'grp' + finalid + '.' + ejectext + '.gloc.fits' $
  ELSE gloc_filename = 'grp' + finalid + '.' + ejectext + '.dcut.gloc.fits'

  IF NOT file_test(gloc_filename) OR keyword_set(steps_debug) THEN BEGIN
     stop
     gloc = intarr(npart,nsteps) ;fill this matrix with numbers representing whether the gas is gone (-1), not in the halo (0), in the halo but not the disk (1), or in the disk (2)
  
;-------------------- Step through outputs to find the location of the particles
     FOR i = 0, nsteps - 1 DO BEGIN
        rtipsy,halodat[i].file,h,g,d,s,/justhead
        readarr,halodat[i].file + '.iord',h,iord,part = 'gas',/ascii,type = 'long'
;        readarr,halodat[i].file + '.OxMassFrac',h,oxmassfrac,part = 'gas',/ascii
;        readarr,halodat[i].file + '.FeMassFrac',h,femassfrac,part = 'gas',/ascii
;        stop
        match,gpart.iord,iord,ind1,ind2
        iord = iord[ind2]
;        oxmassfrac = oxmassfrac[ind2]
;        femassfrac = femassfrac[ind2]
        IF keyword_set(Hcut) THEN BEGIN
           read_tipsy_arr,halodat[i].file + '.HI',h,HI,part = 'gas'
           read_tipsy_arr,halodat[i].file + '.H2',h,H2,part = 'gas'
           HI = HI[ind2]
           H2 = H2[ind2]
        ENDIF
        match2,gpart.iord,iord,ind1,ind2
 ;       oxmassfrac = oxmassfrac[ind1[where(ind1 NE -1)]]
 ;       femassfrac = femassfrac[ind1[where(ind1 NE -1)]]
        IF keyword_set(Hcut) THEN BEGIN
           Hall = fltarr(n_elements(gpart.iord))
           Hall[where(ind1 NE -1)] = (HI[ind1[where(ind1 NE -1)]] + 2*H2[ind1[where(ind1 NE -1)]])/max(HI + 2*H2)
        ENDIF        

        stat = read_stat_struc_amiga(halodat[i].file + '.amiga.stat')
        main = where(halodat[i].haloid EQ stat.group)
        satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir AND stat.ngas gt 0)
;     satellites2 = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
;     IF keyword_set(debug) THEN print,halodat[i].file,satellites,satellites2
        z = center[i,1] 
        scale = 1.0/(1.0 + z)
        
;-------------------- Convert the position information gpart to the frame of the halo
        az = reform(a[i,*])
        az0 = az[0]
        az1 = az[1]
        az2 = az[2]
        ax = [az2/sqrt(az0*az0 + az2*az2),0,-1.0*az0/sqrt(az0*az0 + az2*az2)]
        ay = crossp(az,ax)      ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
        basis = [[ax],[ay],[az]]
        gpos = [[gpart.x[i] - center[i,2]*scale],[gpart.y[i] - center[i,3]*scale],[gpart.z[i] - center[i,4]*scale]] ;*scale
        gpos = transpose(transpose(basis)#transpose(gpos))
        gpart.x[i] = gpos[*,0]
        gpart.y[i] = gpos[*,1]
        gpart.z[i] = gpos[*,2]
        gvel = [[gpart.vx[i] - halodat[i].vxc*scale],[gpart.vy[i] - halodat[i].vyc*scale],[gpart.vz[i] - halodat[i].vzc*scale]]
        gvel = transpose(transpose(basis)#transpose(gvel))
        gpart.vx[i] = gvel[*,0]
        gpart.vy[i] = gvel[*,1]
        gpart.vz[i] = gvel[*,2]
        
;-------------------- Gas particles that have been deleted are set to -1
        gone = where(gpart.grp[i] EQ -1)
        IF gone[0] NE -1 THEN gloc[gone,i] = -1
        
;-------------------- Gas particles that are in the halo/satellites of the halo are set to 2 
;     distance = sqrt(gpart.x[i]*gpart.x[i] + gpart.y[i]*gpart.y[i] + gpart.z[i]*gpart.z[i]) ;should be the same thing to gpart.radius[i]
;     inhalo = where(distance LE halodat[i].rvir)
;     inhalo = where(gpart.grp[i] EQ halodat[i].haloid)
        inhalo = [0]
        IF keyword_set(centralhalo) THEN BEGIN
           test = where(gpart.grp[i] EQ stat[main].group, ntest)
           IF ntest NE 0 THEN inhalo = [inhalo,test]
        ENDIF ELSE BEGIN
           FOR j = 0, n_elements(satellites) - 1 DO BEGIN
              test = where(gpart.grp[i] EQ stat[satellites[j]].group, ntest)
              IF ntest NE 0 THEN inhalo = [inhalo,test]
           ENDFOR
        ENDELSE
        IF n_elements(inhalo) NE 1 THEN BEGIN
           inhalo = inhalo[1:n_elements(inhalo) - 1]
           gloc[inhalo,i] = 2
        ENDIF ELSE BEGIN
           print,'Empty Halo'
;        stop
        ENDELSE
;-------------------- Gas particles that are part of the disk 
        IF keyword_set(Hcut) THEN $
           disk = where(gloc[*,i] EQ 2 AND abs(gpart.z[i]) LE 3.0 AND Hall GE 0.5) ELSE $
              disk = where(gloc[*,i] EQ 2 AND gpart.rho[i] GE 0.1 AND gpart.temp[i] LE 1.2e4 AND abs(gpart.z[i]) LE 3.0) ;1.2    
;(in the galaxy, with density GE 0.1, temperature LE 2e4, and abs(z) LE 10kpc) are set to 3
;     disk = where(gpart.grp[i] EQ halodat[i].haloid AND gpart.rho[i] GE 0.1 AND gpart.temp[i] LE 2e4 AND abs(gpart.z[i]) LE 10.0)
        IF disk[0] NE -1 THEN BEGIN
           gloc[disk,i] = 3
           indiskarr[disk] = 1
        ENDIF
        
 ;------------------- Gas particles that are in the halo but unbound
;                    from the disk are set to 1
        IF disk[0] NE -1 THEN diskmass[i] = total(gpart[disk].mass[i]) + stat[main].m_star ELSE diskmass[i] = stat[main].m_star ;solarmasses
        inhalo = where(gloc[*,i] EQ 2)
        IF inhalo[0] NE -1 THEN BEGIN
           pos = [[gpart[inhalo].x[i]],[gpart[inhalo].y[i]],[gpart[inhalo].z[i]]]
           magpos = sqrt(pos[*,0]*pos[*,0] + pos[*,1]*pos[*,1] + pos[*,2]*pos[*,2])
           pos = pos/[[magpos],[magpos],[magpos]]
           tanvel = gpart[inhalo].vx[i]*pos[*,0] + gpart[inhalo].vy[i]*pos[*,1] + gpart[inhalo].vz[i]*pos[*,2] ;km/s
           ind_unbound = where(tanvel - 1.1*sqrt(2*6.67384e-11*diskmass[i]*2e30/(3.086e19*magpos))/1000. GE 0)
           CASE 1 OF
              keyword_set(all_eject): ind_out = indgen(n_elements(inhalo))
              keyword_set(rvircut_half): ind_out = where(magpos GE 0.5*stat[main].rvir*scale)
              keyword_set(rvircut_fifth): ind_out = where(magpos GE 0.2*stat[main].rvir*scale)
              ELSE: ind_out = ind_unbound
           ENDCASE
           gloc[inhalo[ind_out],i] = 1
           print,'Halo/Halo-out/Disk-unbound (' + strtrim(i,2) + '): ',n_elements(inhalo),n_elements(ind_out),n_elements(ind_unbound)
        ENDIF
        
        IF keyword_set(steps_debug) THEN BEGIN
;     IF 1 THEN BEGIN
           IF n_elements(inhalo) LE 1 THEN BEGIN
              print,halodat[i].file
              print,'No particles in halo'
              stop
           ENDIF ELSE BEGIN
              print,halodat[i].file,minmax(gpart[inhalo].grp[i])
              window,0,xsize = 400,ysize = 400
              plot,gpart[inhalo].x[i],gpart[inhalo].z[i],psym = 3,title = halodat[i].file,xtitle = 'X',ytitle = 'Y' ;,yrange = [-1*stat[main].rvir,stat[main].rvir],xrange = [-1*stat[main].rvir,stat[main].rvir]
              oplot,[0,0],[-1000,1000],linestyle = 2
              oplot,[-1000,1000],[0,0],linestyle = 2
              IF ind_out[0] NE -1 THEN oplot,gpart[inhalo[ind_out]].x[i],gpart[inhalo[ind_out]].z[i],psym = 3,color = 80
              IF disk[0] NE -1 THEN oplot,gpart[disk].x[i],gpart[disk].z[i],psym = 3,color = 245
              if 1 then begin
                 window,3,xsize = 600,ysize = 400
                 plot,gpart[inhalo].rho[i],gpart[inhalo].temp[i],psym = 3,/xlog,/ylog,xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]',xrange = [1e-7,1e4],yrange = [10,1e8],title = halodat[i].file,xstyle = 1
                 IF disk[0] NE -1 THEN oplot,gpart[disk].rho[i],gpart[disk].temp[i],psym = 3,color = 245
                 IF (ind_out)[0] NE -1 THEN oplot,gpart[inhalo[ind_out]].rho[i],gpart[inhalo[ind_out]].temp[i],psym = 3,color = 80
                 window,2,xsize = 600,ysize = 400
                 plot,magpos,tanvel,xtitle = '|Pos|',ytitle = 'Tan Vel',psym = 3
                 IF (ind_out)[0] NE -1 THEN oplot,magpos[ind_out],tanvel[ind_out],psym = 3,color = 80
                 oplot,findgen(200),1.1*sqrt(2*6.67384e-11*diskmass[i]*2e30/(3.086e19*findgen(200)))/1000.0,color = 100
                 
                 window,1,xsize = 600,ysize = 400
                 histogramp,abs(gpart[inhalo].z[i]),nbins = 100
                 IF (ind_out)[0] NE -1 THEN histogramp,abs(gpart[inhalo[ind_out]].z[i]),/overplot,color = 80,nbins = 100        
                 IF i GE 98 THEN stop
              ENDIF
           ENDELSE
        ENDIF
     ENDFOR
     IF NOT keyword_set(nowrite) THEN $
        IF keyword_set(Hcut) THEN mwrfits,gloc,gloc_filename,/create ELSE mwrfits,gloc,gloc_filename,/create
  ENDIF ELSE BEGIN
     gloc = mrdfits(gloc_filename)
     gloc_temp = gloc
     gloc_temp[where(gloc NE 3 AND gloc NE -1)] = 0
     gloc_temp[where(gloc EQ 3 OR  gloc EQ -1)] = 1    
     stepsarr = fltarr(nsteps) + 1
     indiskarr = gloc_temp#stepsarr
     indiskarr[where(indiskarr NE 0)] = 1
  ENDELSE
;====================================================================
  
  coolong = fltarr(npart)
  match,gpart.iord,iord,indg,inda ;Match the tracked particles to the last gas array
  IF keyword_set(cooloncut) THEN BEGIN
     coolong[indg] = coolon[inda]*units.timeunit/1e9 
     indcoolon = where(coolong NE 0, complement = indcooloff) 
  ENDIF ELSE BEGIN
     coolong = coolong + 14
     indcoolon = lindgen(n_elements(gpart))
     indcooloff = [-1]
  ENDELSE
  ycoolong = coolong[indcoolon]
  ncoolpart = n_elements(indcoolon)

;------------------- Re-accretion of gas onto the halo (reaccr_iord.fits)
  gloc_temp = gloc
  gloc_temp[where(gloc EQ  1 OR gloc EQ 2)] = 3 ;Set halo and disk gas to 3
  gloc_temp[where(gloc EQ -1)] = 3 ;Set gas that formed stars to disk
  gloc_diff = gloc_temp[*,0:nsteps - 2] - gloc_temp[*,1:nsteps - 1]
  ind_accr = where(gloc_diff EQ -3) ;Select for instances particles were first out of disk and then in it
  istep = fix(ind_accr/npart)
  ipart = ind_accr - istep*npart
  indstart = where(gloc_temp[*,0] EQ 3)
  mass = gpart[ipart].mass
  arraybase = lindgen(n_elements(ipart))*n_elements(gpart[0].mass)
  IF indstart[0] NE -1 THEN BEGIN
     reaccr_iord = [gpart[indstart].iord,gpart[ipart].iord] 
     reaccr_z = [(indstart*0 + halodat[0].z),halodat[istep + 1].z]
     reaccr_m = [reform((gpart[indstart].mass)[0,indstart]),mass[arraybase + istep + 1]]
  ENDIF ELSE BEGIN
     reaccr_iord = gpart[ipart].iord
     reaccr_z = halodat[istep + 1].z
     reaccr_m = mass[arraybase + istep + 1]
  ENDELSE
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,reaccr_iord,'grp' + finalid + '.' + coolext + 'reaccr_iord.fits',/create ;iord of the reaccreted particle
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,reaccr_z,'grp' + finalid + '.' + coolext + 'reaccr_z.fits',/create ;Redshift at which particle is first in the halo
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,reaccr_m,'grp' + finalid + '.' + coolext + 'mass_at_reaccr.fits',/create ;Mass of particle when first in the halo
  IF keyword_set(debug) THEN BEGIN
     window,1
     haloaccr_time = z_to_t([indstart*0 + halodat[0].z,halodat[istep + 1].z])
     histogramp,haloaccr_time,weight = [reform((gpart[indstart].mass)[0,indstart]),mass[arraybase + istep + 1]],binsize = 0.2,xrange = [0,14],xtitle = 'Time [Gyr]',ytitle = 'Mass (Re)Accreted/Expelled to Halo',cum = debug_cum,title = 'Mass (Re)Accreted/Expelled to Halo'
     stop
  ENDIF

;------------------ Re-outflow of gas from halo
  ind_outflow = where(gloc_diff EQ 3) ;Select for instances particles were first in the halo and then out of it
  istep = fix(ind_outflow/npart)
  ipart = ind_outflow - istep*npart
  gloc_temp = [0]
  gloc_diff = [0]
  reoutflow_iord = gpart[ipart].iord
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,gpart[ipart].iord,'grp' + finalid + '.' + coolext + 'reoutflow_iord.fits',/create ;iord of the reoutflowed particle
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,halodat[istep + 1].z   ,'grp' + finalid + '.' + coolext + 'reoutflow_z.fits',/create ;Redshift at which particle is first out of the halo
  mass = gpart[ipart].mass
  arraybase = lindgen(n_elements(ipart))*n_elements(gpart[0].mass) 
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,mass[arraybase + istep + 1],'grp' + finalid + '.' + coolext + 'mass_at_reoutflow.fits',/create ;Mass of particle when first out of the halo
  IF keyword_set(debug) THEN BEGIN
     outflow_time = z_to_t(halodat[istep + 1].z)
     histogramp,outflow_time, weight = mass[arraybase + istep + 1],binsize = 0.2,color = 100, /overplot,/debug,cum = debug_cum
   ENDIF

;-------------------- Gas (re)accreted to disk, need not have been heated by supernova
  gloc_temp = gloc
  gloc_temp[where(gloc_temp EQ  1 OR gloc_temp EQ 2)] = 0 ;Set all unbound halo gas to zero
  gloc_temp[where(gloc_temp EQ 2)] = 0 ;Set all bound halo gas to zero
  gloc_temp[where(gloc_temp EQ -1)] = 3 ;Set gas that formed stars to disk
  gloc_diff = gloc_temp[*,0:nsteps - 2] - gloc_temp[*,1:nsteps - 1]
  ind_accr = where(gloc_diff EQ -3) ;Select for instances particles were first out of the disk and then in it
  istep_later = [fix(ind_accr/npart)]
  IF (where(gloc_temp[*,0] eq 3))[0] NE -1 THEN BEGIN ;Includes particles that start in the disk
     ipart = [where(gloc_temp[*,0] EQ 3),ind_accr - istep_later*npart]
     istep = [lonarr(n_elements(where(gloc_temp[*,0] eq 3))) - 1,fix(ind_accr/npart)] 
  ENDIF ELSE BEGIN
     ipart = [ind_accr - istep_later*npart]    
     istep = fix(ind_accr/npart) ;If no particles start in the disk
  ENDELSE
  time = z_to_t(halodat[istep + 1].z)
  reaccrdiskall_iord = gpart[ipart].iord
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,gpart[ipart].iord,'grp' + finalid + '.' + coolext + 'reaccrdiskall_iord.fits',/create ;iord of the reaccreted particle
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,halodat[istep + 1].z   ,'grp' + finalid + '.' + coolext + 'reaccrdiskall_z.fits',/create ;Redshift at which particle is first in the disk
  mass = gpart[ipart].mass
  arraybase = lindgen(n_elements(ipart))*n_elements(gpart[0].mass)
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,mass[arraybase + istep + 1],'grp' + finalid + '.' + coolext + 'mass_at_reaccrdiskall.fits',/create ;Mass of particle when first in the disk
  IF keyword_set(debug) THEN BEGIN
     accrdisk_time = z_to_t(halodat[istep + 1].z)
     histogramp,accrdisk_time, weight = mass[arraybase + istep + 1],binsize = 0.2, /overplot,/debug,cum = debug_cum,linestyle = 2
  ENDIF

;-------------------- Re-lost (gas that is heated sufficiently to
;                     longer be part of the disk
;                      -- need not have been heated by supernova)
  gloc_temp = gloc
  gloc_temp[where(gloc_temp EQ  1)] = 0 ;Set disk-unbound halo gas to zero
  gloc_temp[where(gloc_temp EQ  2)] = 0 ;Set disk-bound halo gas to zero
  gloc_temp[where(gloc_temp EQ -1)] = 3 ;Set gas that formed stars to disk
  gloc_diff = gloc_temp[*,0:nsteps - 2] - gloc_temp[*,1:nsteps - 1]
  ind_heat = where(gloc_diff EQ 3) ;Select for instances particles were first in disk and then not it the disk
  istep_heat = fix(ind_heat/npart)
  ipart_heat = ind_heat - istep_heat*npart
  heat_time = z_to_t(halodat[istep_heat + 1].z)
  relost_iord = gpart[ipart_heat].iord
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,gpart[ipart_heat].iord,'grp' + finalid + '.relost_iord.fits',/create ;iord of the reheated particle
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,halodat[istep_heat + 1].z   ,'grp' + finalid + '.relost_z.fits',/create ;Redshift at which particle is first out of the disk
  mass = gpart[ipart_heat].mass
  arraybase = lindgen(n_elements(ipart_heat))*n_elements(gpart[0].mass)
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,mass[arraybase + istep_heat + 1],'grp' + finalid + '.mass_at_relost.fits',/create ;Mass of particle when first out of the disk
  IF keyword_set(debug) THEN BEGIN
     outflow_time = z_to_t(halodat[istep_heat + 1].z)
     histogramp,outflow_time, weight = mass[arraybase + istep_heat + 1],binsize = 0.2, /overplot,/debug,cum = debug_cum,linestyle = 2,color = 100
     legend,['Re-accreted Halo','Re-expelled from Halo'],color = [0,100],linestyle = [0,0]
      stop
  ENDIF

;---------------------- First Accretion (gas that either starts as
;                       disk or becomes disk
  gloc_temp = gloc
  gloc_temp[where(gloc_temp EQ  1 OR gloc_temp EQ  2)] = 0 ;Set all unbound and bound halo gas to zero
  gloc_temp[where(gloc_temp EQ -1)] = 3 ;Set gas that formed stars to disk
  FOR j = 0, nsteps - 2 DO gloc_temp[where((gloc_temp[*,j] EQ 3)),j + 1] = 3
  gloc_diff = gloc_temp[*,0:nsteps - 2] - gloc_temp[*,1:nsteps - 1]
  ind_accrfirst = where(gloc_diff EQ -3) ;Select for instances particles were first out of the disk and then in it
  istepfirst_later = [fix(ind_accrfirst/npart)]
  ipartfirst_later = [ind_accrfirst - istepfirst_later*npart]
;  uniqi_first2 = [where(gloc_temp[*,0] EQ 3),(istepfirst_later + 1)*npart + ipartfirst_later]
  IF (where(gloc_temp[*,0] EQ 3)) NE -1 THEN BEGIN ;Includes particles that start in the disk
     istepfirst = [lonarr(n_elements(where(gloc_temp[*,0] EQ 3))) - 1,istepfirst_later]
     ipartfirst = [where(gloc_temp[*,0] EQ 3),ipartfirst_later]
  ENDIF ELSE BEGIN
     istepfirst = istepfirst_later
     ipartfirst = ipartfirst_later
  ENDELSE
  uniqi_first = [(istepfirst + 1)*npart + ipartfirst]
  gloc_temp = [0]
  gloc_diff = [0]

;-------------------- Re-heated (gas that is heated sufficiently to no
;                     longer be part of the disk
;                      -- must be heated by supernova)
  gloc_temp = gloc;[indcoolon,*]
  IF indcooloff[0] NE -1 THEN gloc_temp[indcooloff,*] = -10
  gloc_temp[where(gloc_temp EQ  1 OR gloc_temp EQ  2)] = 0 ;Set disk-unbound and disk-bound halo gas to zero
  gloc_temp[where(gloc_temp EQ -1)] = 3 ;Set gas that formed stars to disk
  gloc_diff = gloc_temp[*,0:nsteps - 2] - gloc_temp[*,1:nsteps - 1]
  ind_heat = where(gloc_diff EQ 3) ;Select for instances particles were first in disk and then not it the disk
  istep_heat = fix(ind_heat/npart)
  ipart_heat = ind_heat - istep_heat*npart
  heat_time = z_to_t(halodat[istep_heat + 1].z)
;  contour,hist_2d(z_to_t(halodat[istep_heat].z),coolong[indcoolon[ipart_heat]],bin1 = 0.5,bin2 = 0.5),nlevels = 254,/fill,xtitle = 'Time last in disk',ytitle = 'Time Cooling is off'
;  oplot,[0,50],[0,50]
  true = where(z_to_t(halodat[istep_heat].z) LT coolong[ipart_heat])
  IF true[0] NE -1 THEN BEGIN
     gloc_temp = [0]
     gloc_diff = [0]
     reheat_iord = gpart[ipart_heat[true]].iord
     IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
        mwrfits,gpart[ipart_heat[true]].iord,'grp' + finalid + '.' + coolext + 'reheat_iord.fits',/create ;iord of the reheated particle
     IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
        mwrfits,halodat[istep_heat[true] + 1].z   ,'grp' + finalid + '.' + coolext + 'reheat_z.fits',/create ;Redshift at which particle is first out of the disk
     mass = gpart[ipart_heat[true]].mass
     arraybase = lindgen(n_elements(true))*n_elements(gpart[0].mass) 
     IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
        mwrfits,mass[arraybase + istep_heat[true] + 1],'grp' + finalid + '.' + coolext + 'mass_at_reheat.fits',/create ;Mass of particle when first out of the disk
     IF keyword_set(debug) THEN BEGIN
        window,1
        histogramp,heat_time[true], weight = mass[arraybase + istep_heat[true] + 1],binsize = 0.2, xrange = [0,14],cum = debug_cum,xtitle = 'Time [Gyr]',ytitle = 'Mass Re-Accreted/Ejected/Expelled',title = 'Mass Re-Accreted/Heated/Ejected/Expelled',/nodata,yrange = [0,total(mass[arraybase + istep_heat[true] + 1])*2]
        histogramp,heat_time[true], weight = mass[arraybase + istep_heat[true] + 1],binsize = 0.2,color = 100, xrange = [0,14],cum = debug_cum,/overplot
 ;       mass_nottrue = gpart[ipart_heat].mass
 ;       arraybase_nottrue = indgen(n_elements(ipart_heat))*n_elements(gpart[0].mass)
 ;       histogramp,heat_time, weight = mass[arraybase_nottrue + istep_heat + 1],binsize = 0.2,color = 100, xrange = [0,14],/overplot,cum = debug_cum,linestyle = 2
        legend,['Accreted (heated)','Accreted (ejected)','Heated','Ejected','Expelled'],linestyle = [0,0,0,0,0],color = [254,210,100,130,170]
        stop
     ENDIF
  ENDIF
  ipart_heat_true = ipart_heat[true]
  istep_heat_true = istep_heat[true]
  uniq_ipart_heat_true = ipart_heat_true[uniq(ipart_heat_true,sort(ipart_heat_true))]
  uniq_istep_heat_true = lonarr(n_elements(uniq_ipart_heat_true)) - 1
  FOR jstep = 0, nsteps -1 DO BEGIN
     IF (where(istep_heat_true EQ jstep))[0] NE -1 THEN BEGIN
        match,uniq_ipart_heat_true,ipart_heat_true[where(istep_heat_true EQ jstep)],ind1,ind2
        IF ind1[0] NE -1 THEN uniq_istep_heat_true[ind1] = jstep
     ENDIF
  ENDFOR

;--------------------------- Accretion and reaccretion of material heated by supernova
  gloc_temp = gloc
;  gloc_temp[indcooloff,*] = -10
  gloc_temp[where(gloc_temp EQ  1 OR gloc_temp EQ 2)] = 0 ;Set all unbound halo gas to zero
  gloc_temp[where(gloc_temp EQ -1)] = 3 ;Set gas that formed stars to disk
  gloc_diff = gloc_temp[*,0:nsteps - 2] - gloc_temp[*,1:nsteps - 1]
  ind_accr = where(gloc_diff EQ -3) ;Select for instances particles were first out of the disk and then in it
  istep = [fix(ind_accr/npart)]
  ipart = [ind_accr - istep*npart]
  uniqi = [(istep + 1)*npart + ipart]
  match2,uniqi,uniqi_first,ind,indfirst
  uniqi_reaccr = uniqi[where(ind EQ -1)] ;instances when particle is being _re_ accreted
  istep_reaccr = [fix(uniqi_reaccr/npart) - 1]
  ipart_reaccr = [uniqi_reaccr - (istep_reaccr + 1)*npart]
  match2,ipart_reaccr,uniq_ipart_heat_true,ind1,ind2
  ipart_reaccr = ipart_reaccr[where(ind1 NE -1)] ;Gets rid of particles that don't show up as being heated by supernova because cooling turned on again too early
  istep_reaccr = istep_reaccr[where(ind1 NE -1)]
  match2,ipart_reaccr,uniq_ipart_heat_true,ind1,ind2
  istep_reaccr_after = istep_reaccr[where(istep_reaccr GT uniq_istep_heat_true[ind1])] ;reaccreted after last ejection
  istep_reaccr_before = istep_reaccr[where(istep_reaccr LE uniq_istep_heat_true[ind1])]  ;reaccreted before last ejection
  ipart_reaccr_after = ipart_reaccr[where(istep_reaccr GT uniq_istep_heat_true[ind1])]
  ipart_reaccr_before = ipart_reaccr[where(istep_reaccr LE uniq_istep_heat_true[ind1])]
  uniq_ipart_reaccr_after = ipart_reaccr_after[uniq(ipart_reaccr_after,sort(ipart_reaccr_after))]
  uniq_istep_reaccr_after = lonarr(n_elements(uniq_ipart_reaccr_after)) - 1
  FOR jstep = nsteps -1, 0 , -1 DO BEGIN
     IF (where(istep_reaccr_after EQ jstep))[0] NE -1 THEN BEGIN
        match,uniq_ipart_reaccr_after,ipart_reaccr_after[where(istep_reaccr_after EQ jstep)],ind1,ind2
        IF ind1[0] NE -1 THEN uniq_istep_reaccr_after[ind1] = jstep
     ENDIF
  ENDFOR

  gloc_temp = [0]
  gloc_diff = [0]
  ipart_reaccrheat = [ipartfirst,ipart_reaccr_before,uniq_ipart_reaccr_after]
  istep_reaccrheat = [istepfirst,istep_reaccr_before,uniq_istep_reaccr_after]
  reaccrdiskheat_iord = gpart[ipart_reaccrheat].iord
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,gpart[ipart_reaccrheat].iord,'grp' + finalid + '.' + coolext + 'reaccrdiskheat_iord.fits',/create ;iord of the reaccreted particle
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,halodat[istep_reaccrheat + 1].z   ,'grp' + finalid + '.' + coolext + 'reaccrdiskheat_z.fits',/create ;Redshift at which particle is first in the dis
  mass = gpart[ipart_reaccrheat].mass
  arraybase = lindgen(n_elements(ipart_reaccrheat))*n_elements(gpart[0].mass)
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,mass[arraybase + istep_reaccrheat + 1],'grp' + finalid + '.' + coolext + 'mass_at_reaccrdiskheat.fits',/create ;Mass of particle when first in the disk
  IF keyword_set(debug) THEN BEGIN
     accrdisk_t = z_to_t(halodat[istep_reaccrheat + 1].z)
;     window,1
;     histogramp,accrdisk_t,weight = mass[arraybase + istepall + 1],binsize = 0.2,cum = debug_cum,xrange = [0,14],xtitle = 'Time [Gyr]',ytitle = 'Mass Re-Accreted/Ejected/Expelled',title = 'Mass Re-Accreted/Heated/Ejected/Expelled',/nodata
     histogramp,accrdisk_t,weight = mass[arraybase + istep_reaccrheat + 1],binsize = 0.2,cum = debug_cum,/overplot,color = 254
     stop
  ENDIF


;-------------------- Re-ejection (gas that moves from being in disk
;                     to being unbound from it -- NO LONGER must be heated by supernova)
  gloc_temp = gloc;[indcoolon,*]
  IF indcooloff[0] NE -1 THEN gloc_temp[indcooloff,*] = -10
  FOR j = 0,nsteps - 2 DO gloc_temp[where((gloc_temp[*,nsteps - 1 - j] EQ 1 OR gloc_temp[*,nsteps - 1 - j] EQ 0) AND gloc_temp[*,nsteps - 2 - j] EQ 2),nsteps - 2 - j] = 1
  gloc_temp[where(gloc_temp EQ  1)] = 0 ;Set disk-unbound halo gas to zero
;  gloc_temp[where(gloc_temp EQ  2)] = 0 ;Set disk-bound halo gas to zero
  gloc_temp[where(gloc_temp EQ -1)] = 3 ;Set gas that formed stars to disk
  gloc_diff = gloc_temp[*,0:nsteps - 2] - gloc_temp[*,1:nsteps - 1]
  ind_eject = where(gloc_diff EQ 3) ;Select for instances particles were first in disk and then unbound from the disk (marks the last step in the disk)
  istep_eject = fix(ind_eject/npart)
  ipart_eject = ind_eject - istep_eject*npart
  eject_time = z_to_t(halodat[istep_eject + 1].z)
;  contour,hist_2d(z_to_t(halodat[istep_eject].z),coolong[indcoolon[ipart_eject]],bin1 = 0.5,bin2 = 0.5),nlevels = 254,/fill,xtitle = 'Time last in disk',ytitle = 'Time Cooling is off'
;  oplot,[0,50],[0,50]
  true = where(z_to_t(halodat[istep_eject].z) LT coolong[ipart_eject])
  iord_part_true = gpart[ipart_eject[true]].iord
  IF true[0] NE -1 THEN BEGIN
     gloc_temp = [0]
     gloc_diff = [0]
     reeject_iord = gpart[ipart_eject[true]].iord
     IF NOT keyword_set(nowrite) THEN mwrfits,gpart[ipart_eject[true]].iord,'grp' + finalid + '.' + coolext + ejectext + '_iord.fits',/create ;iord of the reejected particle
     IF NOT keyword_set(nowrite) THEN mwrfits,halodat[istep_eject[true] + 1].z   ,'grp' + finalid + '.' + coolext + ejectext + '_z.fits',/create ;Redshift at which particle is first out of the disk
     mass = gpart[ipart_eject[true]].mass
     arraybase = lindgen(n_elements(true))*n_elements(gpart[0].mass) 
     IF NOT keyword_set(nowrite) THEN mwrfits,mass[arraybase + istep_eject[true] + 1],'grp' + finalid + '.' + coolext + 'mass_at_' + ejectext + '.fits',/create ;Mass of particle when first out of the disk
     IF keyword_set(debug) THEN BEGIN
;        mass_nottrue = gpart[indcoolon[ipart_eject]].mass
;        arraybase_nottrue = indgen(n_elements(ipart_eject))*n_elements(gpart[0].mass)
;        histogramp,eject_time, weight = mass[arraybase_nottrue + istep_eject + 1],binsize = 0.2,color = 130, xrange = [0,14],/overplot,cum = debug_cum,linestyle = 2
        histogramp,eject_time[true], weight = mass[arraybase + istep_eject[true] + 1],binsize = 0.2,color = 130, xrange = [0,14],/overplot,cum = debug_cum
        stop
     ENDIF
  ipart_eject_true = ipart_eject[true]
  istep_eject_true = istep_eject[true]
  uniq_ipart_eject_true = ipart_eject_true[uniq(ipart_eject_true,sort(ipart_eject_true))]
  uniq_istep_eject_true = lonarr(n_elements(uniq_ipart_eject_true)) - 1
  FOR jstep = 0, nsteps -1 DO BEGIN
     IF (where(istep_eject_true EQ jstep))[0] NE -1 THEN BEGIN
        match,uniq_ipart_eject_true,ipart_eject_true[where(istep_eject_true EQ jstep)],ind1,ind2
        IF ind1[0] NE -1 THEN uniq_istep_eject_true[ind1] = jstep
     ENDIF
  ENDFOR

;--------------------------- Accretion and reaccretion of ejected material
  gloc_temp = gloc
;  gloc_temp[indcooloff,*] = -10
  gloc_temp[where(gloc_temp EQ  1)] = 0 ;Set all unbound halo gas to zero
  gloc_temp[where(gloc_temp[*,0] EQ 2),0] = 0
  FOR j = 1, nsteps - 1 DO gloc_temp[where((gloc_temp[*,j] EQ 2) AND gloc_temp[*,j - 1] EQ 0),j] = 0
  gloc_temp[where(gloc_temp EQ -1)] = 3 ;Set gas that formed stars to disk
  gloc_diff = gloc_temp[*,0:nsteps - 2] - gloc_temp[*,1:nsteps - 1]
  ind_accr = where(gloc_diff EQ -3) ;Select for instances particles were first out of the disk and then in it
  istep = [fix(ind_accr/npart)]
  ipart = [ind_accr - istep*npart]
  uniqi = [(istep + 1)*npart + ipart]
  match2,uniqi,uniqi_first,ind,indfirst
  uniqi_reaccr = uniqi[where(ind EQ -1)] ;instances when particle is being _re_ accreted
  istep_reaccr = [fix(uniqi_reaccr/npart) - 1]
  ipart_reaccr = [uniqi_reaccr - (istep_reaccr + 1)*npart]
  match2,ipart_reaccr,uniq_ipart_eject_true,ind1,ind2
  ipart_reaccr = ipart_reaccr[where(ind1 NE -1)] ;Gets rid of particles that don't show up as being heated by supernova because cooling turned on again too early
  istep_reaccr = istep_reaccr[where(ind1 NE -1)]
  match2,ipart_reaccr,uniq_ipart_eject_true,ind1,ind2
  istep_reaccr_after = istep_reaccr[where(istep_reaccr GT uniq_istep_eject_true[ind1])] ;reaccreted after last ejection
  istep_reaccr_before = istep_reaccr[where(istep_reaccr LE uniq_istep_eject_true[ind1])]  ;reaccreted before last ejection
  ipart_reaccr_after = ipart_reaccr[where(istep_reaccr GT uniq_istep_eject_true[ind1])]
  ipart_reaccr_before = ipart_reaccr[where(istep_reaccr LE uniq_istep_eject_true[ind1])]
  uniq_ipart_reaccr_after = ipart_reaccr_after[uniq(ipart_reaccr_after,sort(ipart_reaccr_after))]
  uniq_istep_reaccr_after = lonarr(n_elements(uniq_ipart_reaccr_after)) - 1
  FOR jstep = nsteps -1, 0 , -1 DO BEGIN
     IF (where(istep_reaccr_after EQ jstep))[0] NE -1 THEN BEGIN
        match,uniq_ipart_reaccr_after,ipart_reaccr_after[where(istep_reaccr_after EQ jstep)],ind1,ind2
        IF ind1[0] NE -1 THEN uniq_istep_reaccr_after[ind1] = jstep
     ENDIF
  ENDFOR
  gloc_diff = [0]
  ipartall = [ipartfirst,ipart_reaccr_before,uniq_ipart_reaccr_after]
  istepall = [istepfirst,istep_reaccr_before,uniq_istep_reaccr_after]
  
  reaccr_eject_iord = gpart[ipartall].iord
  IF NOT keyword_set(nowrite) THEN mwrfits,gpart[ipartall].iord,'grp' + finalid + '.' + coolext + accrext + '_iord.fits',/create ;iord of the reaccreted particle
  IF NOT keyword_set(nowrite) THEN mwrfits,halodat[istepall + 1].z   ,'grp' + finalid + '.' + coolext + accrext + '_z.fits',/create  ;Redshift at which particle is first in the disk
  mass = gpart[ipartall].mass
  arraybase = lindgen(n_elements(ipartall))*n_elements(gpart[0].mass)
  IF NOT keyword_set(nowrite) THEN mwrfits,mass[arraybase + istepall + 1],'grp' + finalid + '.' + coolext + 'mass_at_' + accrext + '.fits',/create ;Mass of particle when first in the disk
  IF keyword_set(debug) THEN BEGIN
;     debug_cum = 1
     accrdisk_t = z_to_t(halodat[istepall + 1].z)
;     window,1
     histogramp,accrdisk_t,weight = mass[arraybase + istepall + 1],binsize = 0.2,cum = debug_cum,/overplot,color = 210;xrange = [0,14],xtitle = 'Time [Gyr]',ytitle = 'Mass Re-Accreted/Ejected/Expelled',title = 'Mass Re-Accreted/Ejected/Expelled'
     stop
  ENDIF

;------------------- Re-expullsion (gas that moves from being in the
;                    disk to being out of the halo -- must be heated by supernova)  
     gloc_temp = gloc
     IF indcooloff[0] NE -1 THEN gloc_temp[indcooloff,*] = -10
     gloc_temp[ipart_eject_true,istep_eject_true] = 4 ;Mark ejected gas by setting to 4
     gloc_temp[where(gloc_temp EQ -1)] = 3  ;Set gas that formed stars to disk
     FOR j = 1, nsteps - 1 DO gloc_temp[where((gloc_temp[*,j] EQ 1 OR gloc_temp[*,j] EQ 2) AND gloc_temp[*,j - 1] EQ 4),j] = 4
     gloc_diff = gloc_temp[*,0:nsteps - 2] - gloc_temp[*,1:nsteps - 1]
     ind_expell = where(gloc_diff EQ 4) ;Select for instances particles were first in disk and then out of the halo
     istep = fix(ind_expell/npart)
     ipart = ind_expell - istep*npart
     match2,gpart[ipart].iord,iord_part_true,true_expell,indtemp
     ipart = ipart[where(true_expell NE -1)]
     istep = istep[where(true_expell NE -1)]
     gloc_temp = [0]
     gloc_diff = [0]
     reexpell_iord = gpart[ipart].iord
     IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
        mwrfits,gpart[ipart].iord,'grp' + finalid + '.' + coolext + 'reexpell_iord.fits',/create ;iord of the reexpelled particle
     IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
        mwrfits,halodat[istep + 1].z   ,'grp' + finalid + '.' + coolext + 'reexpell_z.fits',/create ;Redshift at which particle is first out of the disk
     mass = gpart[ipart].mass
     arraybase = lindgen(n_elements(ipart))*n_elements(gpart[0].mass)
     IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
        mwrfits,mass[arraybase + istep + 1],'grp' + finalid + '.' + coolext + 'mass_at_reexpell.fits',/create ;Mass of particle when first out of the disk
     IF keyword_set(debug) THEN BEGIN
        expell_time = z_to_t(halodat[istep + 1].z)
;        mass_nottrue = gpart[indcoolon[ipart]].mass
;        arraybase_nottrue = indgen(n_elements(istep))*n_elements(gpart[0].mass)
;        histogramp,expell_time, weight = mass_nottrue[arraybase_nottrue + istep + 1],binsize = 0.2,color = 170, /overplot,cum = debug_cum,linestyle = 2 ;, xrange = [0,10],xtitle = 'Redshift',ytitle = 'Mass (Re)Expelled From Disk'
        histogramp,expell_time, weight = mass[arraybase + istep + 1],binsize = 0.2,color = 170, xrange = [0,14],/overplot,cum = debug_cum
;        legend,['Accreted (heated)','Accreted (ejected)','Heated','Ejected','Expelled'],linestyle = [0,0,0,0,0],color = [254,210,100,130,170]
        stop
     ENDIF
;  stop
ENDIF

;=====================================

;-------------------- Find earlierst step when the particle was in the disk
  everindisk  = where(indiskarr NE 0)
  accrdisk    = intarr(npart) - 1 ;Gas that is ever in the disk
  FOR i = 0L, n_elements(everindisk) - 1 DO BEGIN ;iterate through gas particles were ever in the disk
     indisk = where(gloc[everindisk[i],*] EQ 3)   
     accrdisk[everindisk[i]] = indisk[0]
  ENDFOR
  neverindisk = where(accrdisk EQ -1, complement = everindisk)
  temp = fltarr(n_elements(accrdisk))
  temp[ everindisk] = halodat[accrdisk[everindisk]].z
  temp[neverindisk] = 99
  accrdisk_z = temp
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,temp,'grp' + finalid + '.' + coolext + 'accrdisk_z.fits',/create ;Time when accreted onto the disk

  mass = fltarr(n_elements(accrdisk))
  FOR i = 0L, n_elements(everindisk) - 1 DO mass[everindisk[i]] = gpart[everindisk[i]].mass[accrdisk[everindisk[i]]] ;uneject particles = 0
    accrdisk_mass = mass
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,mass,'grp' + finalid + '.' + coolext + 'mass_at_accrdisk.fits',/create ;Mass when accreted onto the disk
  IF keyword_set(debug) THEN BEGIN
     diskaccr_time = z_to_t(temp[everindisk])
;  histogramp,13.7 - wmap3_lookback(temp[everindisk])/1e9,weight = mass[everindisk], binsize = 0.2,xrange = [1e-5,14];,/cum
;     histogramp,temp[everindisk],weight = mass[everindisk], binsize = 0.2,xrange = [0,10],min = 0, max = 10,xtitle = 'Redshift',ytitle = 'Disk Mass Accreted' ;,/cum
     histogramp,diskaccr_time,weight = mass[everindisk], binsize = 0.2,xrange = [0,14],min = 0, max = 14,xtitle = 'Time [Gyr]',ytitle = 'Disk Mass Accreted',cum = debug_cum,title = 'Disk Mass Accreted'
     stop
  ENDIF

;-------------------- Find earlierst step when the particle was in the halo
  accr = intarr(npart) - 1 ;Gas that is ever in the disk
  FOR i = 0L, n_elements(accr) - 1 DO BEGIN ;iterate through gas particles were ever in the disk
     accreted = where(gloc[i,*] EQ 2 OR gloc[i,*] EQ 1 OR gloc[i,*] EQ 3)   
     accr[i] = accreted[0]
  ENDFOR
  neverinhalo = where(accr EQ -1, complement = everinhalo)
  temp = fltarr(n_elements(accr))
  temp[ everinhalo] = halodat[accr[everinhalo]].z
  IF neverinhalo[0] NE -1 THEN temp[neverinhalo] = 99
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,temp,'grp' + finalid + '.' + coolext + 'accr_rvir_z.fits',/create ;Time when accreted onto the disk

  mass = fltarr(n_elements(accr))
  FOR i = 0L, n_elements(everinhalo) - 1 DO mass[everinhalo[i]] = gpart[everinhalo[i]].mass[accr[everinhalo[i]]] ;uneject particles = 0
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,mass,'grp' + finalid + '.' + coolext + 'mass_at_accr_rvir.fits',/create ;Mass when accreted onto the disk
  IF keyword_set(debug) THEN BEGIN
     haloaccr_time = z_to_t(temp)
;  histogramp,13.7 - wmap3_lookback(temp)/1e9,weight = mass, binsize = 0.2,xrange = [1e-5,14];,/cum
;     histogramp,temp,weight = mass, binsize = 0.2,xrange = [0,10],min = 0, max = 10,xtitle = 'Redshift',ytitle = 'Mass Accreted To Halo' ;,/cum ;,/cum
     histogramp,haloaccr_time,weight = mass, binsize = 0.2,xrange = [0,14],min = 0, max = 14,xtitle = 'Time [Gyr]',ytitle = 'Mass Accreted To Halo',cum = debug_cum,title = 'Mass Accreted To Halo'
     stop
  ENDIF

;-------------------- Gas Loss 
  lost        = intarr(npart) - 1 ;Gas that is in the halo and then out of the halo
  FOR i = 0L, n_elements(lost) - 1 DO BEGIN ;iterate through gas particles that were ever in the disk
     outhalo = where(gloc[i,*] EQ 0)   
     indhalo = where(gloc[i,*] GE 1)
     temp = min(indhalo,firsthalo) ;First timestep that it is in the halo
     firsthalo = indhalo[firsthalo]
;    outhalo = where(gloc[indcoolon[i],*] EQ 0)
     IF (where(outhalo GT firsthalo))[0] NE -1 THEN lost[i] = outhalo[(where(outhalo GT firsthalo))[0]]
  ENDFOR
  neverlost = where(lost EQ -1, complement = everlost)
  temp = fltarr(n_elements(lost))
  temp[ everlost] = halodat[lost[everlost]].z
  temp[neverlost] = 99
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,temp,'grp' + finalid + '.' + coolext + 'outflow_z.fits',/create ;Time when accreted onto the disk

  mass = fltarr(n_elements(lost))
  FOR i = 0L, n_elements(everlost) - 1 DO mass[everlost[i]] = gpart[everlost[i]].mass[lost[everlost[i]]] ;uneject particles = 0
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,mass,'grp' + finalid + '.' + coolext + 'mass_at_outflow.fits',/create ;Mass when accreted onto the disk
  IF keyword_set(debug) THEN BEGIN
     outflow_time = z_to_t(temp[everindisk])
;  histogramp,13.7 - wmap3_lookback(temp[everindisk])/1e9,weight = mass[everindisk], binsize = 0.2,xrange = [1e-5,14];,/cum
;     histogramp,temp[everindisk],weight = mass[everindisk], binsize = 0.2,xrange = [0,10],xtitle = 'Redshift',ytitle = 'Mass Outflowing From Disk' ;,/cum
     histogramp,outflow_time,weight = mass[everindisk], binsize = 0.2,xrange = [0,14],xtitle = 'Redshift',ytitle = 'Mass Outflowing From Disk',cum = debug_cum,title = 'Mass Outflowing From Disk' 
     stop
  ENDIF

;------------------------------- Ejection ------------------------------
  coolong = fltarr(npart)
  match,gpart.iord,iord,indg,inda ;Match the tracked particles for halo 1 to the last gas array
  coolong[indg] = coolon[inda]*units.timeunit/1e9 
;  indcoolon = where(coolong NE 0)
  it = nsteps - 1 ;Index of final step
;  dense = where(gpart.rho[it] GT 100) 
;  print,dense[1] 
;  ip = 0
;  inhalo = where(gloc[*,it] EQ 1)
;  disk = where(gloc[*,it] EQ 2)
;  allejectfinal = where(coolong NE 0 AND gloc[*,it] EQ 0) ;Gas that is heated by SN + out of the halo at the last step
;About 95% of the particles that are ejected from the halo at the final timestep and have been heated by SN have been in the disk

  outflow     = intarr(npart) - 1 ;Gas that is heated by SN and out of the halo
;  reaccret    = intarr(npart) - 1 ;Gas that is heated by SN + in the disk then out of the disk and then in the disk
  eject       = intarr(npart) - 1 ;Gas that is heated by SN + in the disk then out of the disk                     
  expell      = intarr(npart) - 1 ;Gas that is heated by SN + in the disk then out of the halo
  expellfinal = intarr(npart) - 1 ;Gas that is heated by SN + in the disk then out of the halo at the final step 
  indcoolon   = where(coolong NE 0, complement = indcooloff) 

  FOR i = 0L, n_elements(indcoolon) - 1 DO BEGIN ;iterate through gas particles that are heated by SN 
     indisk  = where(gloc[indcoolon[i],*] EQ 3)
     outhalo = where(gloc[indcoolon[i],*] EQ 0)
     IF indisk[0] NE -1 THEN BEGIN
        temp = min(indisk,firstdisk) ;First timestep it is in the disk
        temp = max(indisk,lastdisk)  ;Last  timestep it is in the disk
        firstdisk = indisk[firstdisk] ;Index of step where particle is first in disk
        lastdisk  = indisk[lastdisk] ;Index of step where partcile is last in disk
        outdisk = where(gloc[indcoolon[i],*] EQ 0 OR gloc[indcoolon[i],*] EQ 1); OR gloc[indcoolon[i],*] EQ 2) ;Indices of steps during which particle is out of disk
;        IF (where(outdisk GT firstdisk AND outdisk LT lastdisk))[0] NE -1 THEN reaccret[indcoolon[i]] = (where(outdisk GT firstdisk AND outdisk LT lastdisk))[0]
        IF (where(outdisk GT firstdisk))[0] NE -1 THEN       eject[indcoolon[i]] = outdisk[(where(outdisk GT firstdisk))[0]] ;set the eject array to the timestep when the particle is first out of the disk following being in the disk
        IF (where(outhalo GT firstdisk))[0] NE -1 THEN      expell[indcoolon[i]] = outhalo[(where(outhalo GT firstdisk))[0]] ;set the expell array to the timestep when the particle is first out of the halo following being in the disk
        IF (gloc[indcoolon[i],it] EQ 0)      THEN expellfinal[indcoolon[i]] = outhalo[(where(outhalo GT firstdisk))[0]] ;set the expellfinal array to the timestep when the particle is first out of the halo following being in the disk 
    ENDIF
    indhalo = where(gloc[indcoolon[i],*] GE 1) ;Step indices during which particle is in the halo
    temp = min(indhalo,firsthalo) ;First timestep that it is in the halo
    firsthalo = indhalo[firsthalo]
;    outhalo = where(gloc[indcoolon[i],*] EQ 0)
    IF (where(outhalo GT firsthalo))[0]     NE -1 THEN     outflow[indcoolon[i]] = outhalo[(where(outhalo GT firsthalo))[0]]
  ENDFOR

;  nooutflow = where(outflow EQ -1, complement = yesoutflow)
;  temp = fltarr(n_elements(outflow))
;  temp[yesoutflow] = halodat[outflow[yesoutflow]].z
;  temp[nooutflow]  = 99  
;  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.outflow2_z.fits',/create
;  temp = lonarr(n_elements(outflow))
;  temp[yesoutflow] = gpart[yesoutflow].iord
;  temp[nooutflow]  = 0
;  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.outflow2_iord.fits',/create 

;Output the iord, z, and mass of particles being ejected, when ejected
  temp = gpart.iord
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,temp,'grp' + finalid + '.' + coolext + 'eject_iord.fits',/create ;(number in entropy - early;)

  noeject = where(eject EQ -1, complement = yeseject)

  temp = fltarr(n_elements(eject))
  temp[yeseject] = halodat[eject[yeseject]].z
  temp[noeject]  = 99  
  IF keyword_set(debug) THEN BEGIN
     eject_time = z_to_t(temp[yeseject])
;  histogramp,13.7 - wmap3_lookback(temp[yeseject])/1e9, binsize = 0.2,xrange = [1e-5,14];,/cum
;     histogramp,temp[yeseject], binsize = 0.2,xrange = [0,10],xtitle = 'Redshift',ytitle = 'Mass Ejected From Disk' ;,/cum
     histogramp,eject_time[where(temp[yeseject] ne 99)], binsize = 0.2,xrange = [0,14],xtitle = 'Redshift',ytitle = 'Mass Ejected/Expelled From Disk',cum = debug_cum,title = 'Mass Ejected/Expelled From Disk'
  ENDIF
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,temp,'grp' + finalid + '.' + coolext + 'eject_z.fits',/create         

  temp = fltarr(n_elements(eject))
  FOR i = 0L, n_elements(yeseject) - 1 DO temp[yeseject[i]] = gpart[yeseject[i]].mass[eject[yeseject[i]]] ;uneject particles = 0
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,temp,'grp' + finalid + '.' + coolext + 'mass_at_eject.fits',/create   

  temp = fltarr(n_elements(eject))
  temp[yeseject] = halodat[eject[yeseject] - 1].z
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,temp,'grp' + finalid + '.' + coolext + 'eject_timeindisk.fits',/create ;in redshift (0 = never expelled OR never in disk)

;Output the iord, z, and mass of particles being expelled, when expelled
  noexpell = where(expell EQ -1, complement = yesexpell)

  temp = fltarr(n_elements(expell))
  temp[yesexpell] = halodat[expell[yesexpell]].z
  temp[noexpell]  = 99  
  If Keyword_set(debug) THEN BEGIN
     expell_time = z_to_t(temp[yesexpell])
;  histogramp,13.7 - wmap3_lookback(temp[yesexpell])/1e9, binsize = 0.2,xrange = [1e-5,14],/overplot,linestyle = 2;,/cum
;     histogramp,temp[yesexpell], binsize = 0.2,xrange = [0,10],/overplot,linestyle = 2 ;,/cum
     histogramp,expell_time[where(temp[yesexpell] NE 99)], binsize = 0.2,xrange = [0,14],/overplot,linestyle = 2,cum = debug_cum
     legend,['Ejected','Expelled'],linestyle = [0,2],/left
     stop
  ENDIF

  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,temp,'grp' + finalid + '.' + coolext + 'expell_z.fits',/create         

  temp = fltarr(n_elements(expell))
  FOR i = 0L, n_elements(yesexpell) - 1 DO temp[yesexpell[i]] = gpart[yesexpell[i]].mass[expell[yesexpell[i]]] ;uneject particles = 0
  IF NOT keyword_set(nowrite) AND NOT keyword_set(limitedwrite) THEN $
     mwrfits,temp,'grp' + finalid + '.' + coolext + 'mass_at_expell.fits',/create   

  IF keyword_set(plots) THEN BEGIN
     loadct,39
;----------------------------- Plots -----------------------------
;  ip = (where(reaccret    NE -1))[0]
;  ip = (where(eject       NE -1))[0]
;  ip = (where(expell      NE -1))[0]
  ip = (where(expellfinal NE -1))[0]

  outhalo =  where(gloc[*,it] EQ 0)
  inhalo_unbound = where(gloc[*,it] EQ 1)
  inhalo = where(gloc[*,it] EQ 2)
  disk = where(gloc[*,it] EQ 3)

  disk_early = where(gloc[*,it-2] EQ 3)

  ipdisk = where(gloc[ip,*] EQ 3)
  iphalo_bound = where(gloc [ip,*] EQ 2)
  iphalo_unbound = where(gloc[ip,*] EQ 1)        
  sntime = max(halodat[where(halodat.time LE coolong[ip])].time)
  ipsn = where(halodat.time EQ sntime)
  IF NOT keyword_set(outplot) THEN window,2,xsize = 800, ysize = 500 ELSE  device,filename = 'example.fb.phase.eps',bits_per_pixel= 8,/times,ysize=7,xsize=7,/inch,/color,xoffset = 3, yoffset = 5;,ytitle = 'Temperature [K]',xtitle = 'Density [amu/cc]'
  plot,gpart.rho[it],gpart.temp[it],psym = 3,/xlog,/ylog,xrange = [1e-7,1e4],yrange = [100,1e6],ytitle = 'Temperature [K]',xtitle = 'Density [amu/cc]',xstyle = 1
  oplot,gpart[inhalo].rho[it],gpart[inhalo].temp[it],psym = 3, color = 100
  oplot,gpart[disk].rho[it],gpart[disk].temp[it],psym = 3, color = 50
  oplot,gpart[ip].rho,gpart[ip].temp,psym = -4,color = 190
  IF n_elements(iphalo_unbound) EQ 1 THEN $
     oplot,[gpart[ip].rho[iphalo_unbound],gpart[ip].rho[iphalo_unbound]],[gpart[ip].temp[iphalo_unbound],gpart[ip].temp[iphalo_unbound]],psym = 4,color = 200 ELSE $
     oplot,gpart[ip].rho[iphalo_unbound],gpart[ip].temp[iphalo_unbound],psym = 4,color = 200
  IF n_elements(iphalo_bound) EQ 1 THEN $
     oplot,[gpart[ip].rho[iphalo_bound],gpart[ip].rho[iphalo_bound]],[gpart[ip].temp[iphalo_bound],gpart[ip].temp[iphalo_bound]],psym = 4,color = 215 ELSE $
     oplot,gpart[ip].rho[iphalo_bound],gpart[ip].temp[iphalo_bound],psym = 4,color = 215
  IF n_elements(ipdisk) EQ 1 THEN $
     oplot,[gpart[ip].rho[ipdisk],gpart[ip].rho[ipdisk]],[gpart[ip].temp[ipdisk],gpart[ip].temp[ipdisk]],psym = 4,color = 254 ELSE $
     oplot,gpart[ip].rho[ipdisk],gpart[ip].temp[ipdisk],psym = 4,color = 254
  oplot,[gpart[ip].rho[ipsn],gpart[ip].rho[ipsn]],[gpart[ip].temp[ipsn],gpart[ip].temp[ipsn]],psym = symcat(14),color = 254,symsize = 3
  legend,['Not in halo',' ','In halo',' ','Bound In halo',' ','In Disk',' ','History'],color = [0,190,100,200,100,220,50,254,190],linestyle = [1,0,1,0,1,0,1,0,0],psym = [3,4,3,4,3,4,3,4,3],/bottom,/left

  IF NOT keyword_set(outplot) THEN BEGIN
     stop
     window,6,xsize = 500, ysize = 500
  ENDIF ELSE BEGIN
     device,/close
     device,filename = 'example.fb.fo.eps',bits_per_pixel= 8,/times,ysize=7,xsize=7,/inch,/color,xoffset = 3, yoffset = 5 ;,ytitle = 'Y [kpc]',xtitle = 'X [kpc]'
  ENDELSE
;  verydense = where(gpart.rho[it] GE 100)
  range = 40
  plot,gpart.x[it],gpart.y[it],psym = 3,ytitle = 'Y [kpc]',xtitle = 'X [kpc]',xrange = [mean(gpart[disk].x[it])-range,mean(gpart[disk].x[it])+range], yrange = [mean(gpart[disk].y[it])-range,mean(gpart[disk].y[it])+range],/nodata
;  oplot,gpart[inhalo].x[it],gpart[inhalo].y[it],psym = 3, color = 100
  oplot,gpart[disk].x[it],gpart[disk].y[it],psym = 3, color = 50
  oplot,gpart[disk_early].x[it],gpart[disk_early].y[it],psym = 3, color = 250
;  oplot,gpart[verydense].x[it],gpart[verydense].y[it],psym = 3, color = 20
  oplot,gpart[ip].x,gpart[ip].y,psym = -4,color = 190
  IF n_elements(iphalo_unbound) EQ 1 THEN $
     oplot,[gpart[ip].x[iphalo_unbound],gpart[ip].x[iphalo_unbound]],[gpart[ip].y[iphalo_unbound],gpart[ip].y[iphalo_unbound]],psym = 4,color = 200 ELSE $
     oplot,gpart[ip].x[iphalo_unbound],gpart[ip].y[iphalo_unbound],psym = 4,color = 200
  IF n_elements(iphalo_bound) EQ 1 THEN $
     oplot,[gpart[ip].x[iphalo_bound],gpart[ip].x[iphalo_bound]],[gpart[ip].y[iphalo_bound],gpart[ip].y[iphalo_bound]],psym = 4,color = 220 ELSE $
     oplot,gpart[ip].x[iphalo_bound],gpart[ip].y[iphalo_bound],psym = 4,color = 220
  IF n_elements(ipdisk) EQ 1 THEN $
     oplot,[gpart[ip].x[ipdisk],gpart[ip].x[ipdisk]],[gpart[ip].y[ipdisk],gpart[ip].y[ipdisk]],psym = 4,color = 254 ELSE $
     oplot,gpart[ip].x[ipdisk],gpart[ip].y[ipdisk],psym = 4,color = 254
  oplot,[gpart[ip].x[ipsn],gpart[ip].x[ipsn]],[gpart[ip].y[ipsn],gpart[ip].y[ipsn]],psym = symcat(14),color = 254,symsize = 3

  IF NOT keyword_set(outplot) THEN stop ELSE BEGIN
     device,/close
     device,filename = 'example.fb.eo.eps',bits_per_pixel= 8,/times,ysize=7,xsize=7,/inch,/color,xoffset = 3, yoffset = 5
  ENDELSE
  plot,gpart.x[it],gpart.z[it],psym = 3,xrange = [-1*range,range], yrange = [-1*range,range],ytitle = 'Z [kpc]',xtitle = 'X [kpc]'
  oplot,gpart[inhalo].x[it],gpart[inhalo].z[it],psym = 3, color = 100
  oplot,gpart[disk].x[it],gpart[disk].z[it],psym = 3, color = 50
  oplot,gpart[ip].x,gpart[ip].z,psym = -4,color = 190
  IF n_elements(iphalo_unbound) EQ 1 THEN $
     oplot,[gpart[ip].x[iphalo_unbound],gpart[ip].x[iphalo_unbound]],[gpart[ip].z[iphalo_unbound],gpart[ip].z[iphalo_unbound]],psym = 4,color = 200 ELSE $
     oplot,gpart[ip].x[iphalo_unbound],gpart[ip].z[iphalo_unbound],psym = 4,color = 200
  IF n_elements(iphalo_bound) EQ 1 THEN $
     oplot,[gpart[ip].x[iphalo_bound],gpart[ip].x[iphalo_bound]],[gpart[ip].z[iphalo_bound],gpart[ip].z[iphalo_bound]],psym = 4,color = 220 ELSE $
     oplot,gpart[ip].x[iphalo_bound],gpart[ip].z[iphalo_bound],psym = 4,color = 220
  IF n_elements(ipdisk) EQ 1 THEN $
     oplot,[gpart[ip].x[ipdisk],gpart[ip].x[ipdisk]],[gpart[ip].z[ipdisk],gpart[ip].z[ipdisk]],psym = 4,color = 254 ELSE $
     oplot,gpart[ip].x[ipdisk],gpart[ip].z[ipdisk],psym = 4,color = 254
  oplot,[gpart[ip].x[ipsn],gpart[ip].x[ipsn]],[gpart[ip].z[ipsn],gpart[ip].z[ipsn]],psym = symcat(14),color = 254,symsize = 3
  IF keyword_set(outplot) THEN device,/close
  ENDIF
END


PRO master_gas_gal_disk_history

  dir = '/nobackupp2/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
  finalid = '2'
  gas_gal_disk_history,dir = dir,finalid = finalid
  finalid = '1'
  gas_gal_disk_history,dir = dir,finalid = finalid  

  dir = '/nobackupp2/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
  finalid = '2'
  gas_gal_disk_history,dir = dir,finalid = finalid
  finalid = '1'
  gas_gal_disk_history,dir = dir,finalid = finalid

  dir = '/nobackupp2/crchrist/MolecH/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
  finalid = '3'
  gas_gal_disk_history,dir = dir,finalid = finalid
  finalid = '2'
  gas_gal_disk_history,dir = dir,finalid = finalid
  finalid = '1'
  gas_gal_disk_history,dir = dir,finalid = finalid

  dir = '/nobackupp2/crchrist/MolecH/h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
  finalid = '2'
  gas_gal_disk_history,dir = dir,finalid = finalid
  finalid = '1'
  gas_gal_disk_history,dir = dir,finalid = finalid
END
