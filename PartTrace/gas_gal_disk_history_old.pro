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

PRO gas_gal_disk_history,dir = dir,finalid = finalid,laststep = laststep,centralhalo = centralhalo,plots = plots,nowrite = nowrite,outplot = outplot,debug = debug
  formatplot,outplot = outplot

;-------------------- Read in simulation data for final step
  IF NOT keyword_set(finalid) THEN finalid = '1' 
  IF NOT keyword_set(laststep) THEN laststep = '512'
  IF NOT keyword_set(dir) THEN dir = '.'
  cd,dir
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.coolontime',file_coolon
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file
  rtipsy,file[0],h,g,d,s,/justhead
  readarr,file_coolon[0],h,coolon,part = 'gas',/ascii
  readarr,file_iord[0],h,iord,part = 'gas',/ascii,type = 'long'  
  spawn,'ls ' + dir + '/h*.param',pfile
  units = tipsyunits(pfile[0])

;-------------------- Read in the information about the main halo at each timestep
;  readcol,'alignment.txt',files,haloid,time,z,xc,yc,zc,xa,ya,za,format = '(A,I,F,F,F,F,F,F,F,F)'
;  center = [[time],[z],[xc],[yc],[zc]]
;  center[*,2] = center[*,2]*units.lengthunit ;"center" is in tipsy units.  I need to change it to comoving units
;  center[*,3] = center[*,3]*units.lengthunit
;  center[*,4] = center[*,4]*units.lengthunit
;  center[*,2] = center[*,2]*1000.0 - units.lengthunit/2.0 ;"center" is in Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
;  center[*,3] = center[*,3]*1000.0 - units.lengthunit/2.0
;  center[*,4] = center[*,4]*1000.0 - units.lengthunit/2.0
;  a = [[xa],[ya],[za]]

;-------------------- Read in the information on the halos
  halodat = mrdfits('grp' + finalid + '.alignment.fits',1)
  halodat.xc = halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.yc = halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.zc = halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  center = [[halodat.time],[halodat.z],[halodat.xc],[halodat.yc],[halodat.zc]]
  center[*,2] = center[*,2]*1000.0 - units.lengthunit/2.0 ;"center" is in Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
  center[*,3] = center[*,3]*1000.0 - units.lengthunit/2.0
  center[*,4] = center[*,4]*1000.0 - units.lengthunit/2.0
  a = [[halodat.xa],[halodat.ya],[halodat.za]]

;-------------------- Read in the information about the gas particles at each timestep
  gpart = mrdfits('grp' + finalid + '.allgas.entropy.fits',1)

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
  IF NOT n_elements(halodat.file) EQ n_elements(gpart[0].mass) THEN print,'Files are for different time steps'

  indiskarr = intarr(n_elements(gpart)) ; List of gas particles that are ever in the disk
  gloc = {loc: intarr(n_elements(halodat.file))}
  gloc = replicate(gloc, n_elements(gpart)) ;fill this matrix with numbers representing whether the gas is gone (-1), not in the halo (0), in the halo but not the disk (1), or in the disk (2)
  
;-------------------- Step through outputs to find the location of the particles
  FOR i = 0, n_elements(halodat.file) - 1 DO BEGIN
     stat = read_stat_struc_amiga(halodat[i].file + '.amiga.stat')
     main = where(halodat[i].haloid EQ stat.group)
     satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir AND stat.ngas gt 0)
     satellites2 = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
     IF keyword_set(debug) THEN print,halodat[i].file,satellites,satellites2
     z = center[i,1] 
     scale = 1.0/(1.0 + z)

;-------------------- Convert the position information gpart to the frame of the halo
     az = reform(a[i,*])
     az0 = az[0]
     az1 = az[1]
     az2 = az[2]
     ax = [az2/sqrt(az0*az0 + az2*az2),0,-1.0*az0/sqrt(az0*az0 + az2*az2)]
     ay = crossp(az,ax)         ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
     basis = [[ax],[ay],[az]]
     gpos = [[gpart.x[i] - center[i,2]*scale],[gpart.y[i] - center[i,3]*scale],[gpart.z[i] - center[i,4]*scale]];*scale
     gpos = transpose(transpose(basis)#transpose(gpos))
     gpart.x[i] = gpos[*,0]
     gpart.y[i] = gpos[*,1]
     gpart.z[i] = gpos[*,2]

;-------------------- Gas particles that have been deleted are set to -1
     gone = where(gpart.grp[i] EQ -1)
     IF gone[0] NE -1 THEN gloc[gone].loc[i] = -1

;-------------------- Gas particles that are in the halo/satellites of the halo are set to 0 
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
        gloc[inhalo].loc[i] = 1
     ENDIF ELSE BEGIN
        print,'PROBLEM'
        stop
     ENDELSE

;-------------------- Gas particles that are part of the disk (in the galaxy, with density GE 0.1, temperature LE 2e4, and abs(z) LE 10kpc) are set to 2
     disk = where(gloc.loc[i] EQ 1 AND gpart.rho[i] GE 0.1 AND gpart.temp[i] LE 2e4 AND abs(gpart.z[i]) LE 10.0)
;     disk = where(gpart.grp[i] EQ halodat[i].haloid AND gpart.rho[i] GE 0.1 AND gpart.temp[i] LE 2e4 AND abs(gpart.z[i]) LE 10.0)
     gloc[disk].loc[i] = 2
     indiskarr[disk] = 1

     IF keyword_set(debug) THEN BEGIN
        print,halodat[i].file,minmax(gpart[inhalo].grp[i])
        window,2,xsize = 400,ysize = 400
        plot,gpart[inhalo].x[i],gpart[inhalo].z[i],psym = 3,title = halodat[i].file,xtitle = 'X',ytitle = 'Y',yrange = [-1*stat[main].rvir,stat[main].rvir],xrange = [-1*stat[main].rvir,stat[main].rvir]
        oplot,gpart[disk].x[i],gpart[disk].z[i],psym = 3,color = 245
        window,3
        plot,gpart[inhalo].rho[i],gpart[inhalo].temp[i],psym = 3,/xlog,/ylog,xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]',xrange = [1e-7,1e4],yrange = [10,1e8],title = halodat[i].file,xstyle = 1
        oplot,gpart[disk].rho[i],gpart[disk].temp[i],psym = 3,color = 245
;        stop
     ENDIF
  ENDFOR
  
;-------------------- Find earlierst step when the particle was in the disk
  everindisk  = where(indiskarr NE 0)
  accrdisk    = intarr(n_elements(gloc)) - 1 ;Gas that is ever in the disk
  FOR i = 0L, n_elements(everindisk) - 1 DO BEGIN ;iterate through gas particles were ever in the disk
     indisk = where(gloc[everindisk[i]].loc EQ 2)   
     accrdisk[everindisk[i]] = indisk[0]
  ENDFOR
  neverindisk = where(accrdisk EQ -1, complement = everindisk)
  temp = fltarr(n_elements(accrdisk))
  temp[ everindisk] = halodat[accrdisk[everindisk]].z
  temp[neverindisk] = 99
  accrdisk_z = temp
  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.accrdisk_z.fits',/create ;Time when accreted onto the disk

  mass = fltarr(n_elements(accrdisk))
  FOR i = 0L, n_elements(everindisk) - 1 DO mass[everindisk[i]] = gpart[everindisk[i]].mass[accrdisk[everindisk[i]]] ;uneject particles = 0
    accrdisk_mass = mass
  IF NOT keyword_set(nowrite) THEN mwrfits,mass,'grp' + finalid + '.mass_at_accrdisk.fits',/create  ;Mass when accreted onto the disk
  IF keyword_set(debug) THEN BEGIN
;  histogramp,13.7 - wmap3_lookback(temp[everindisk])/1e9,weight = mass[everindisk], binsize = 0.2,xrange = [1e-5,14];,/cum
     histogramp,temp[everindisk],weight = mass[everindisk], binsize = 0.2,xrange = [0,10],min = 0, max = 10,xtitle = 'Redshift',ytitle = 'Disk Mass Accreted' ;,/cum
     stop
  ENDIF

;-------------------- Find earlierst step when the particle was in the halo
  accr = intarr(n_elements(gloc)) - 1 ;Gas that is ever in the disk
  FOR i = 0L, n_elements(accr) - 1 DO BEGIN ;iterate through gas particles were ever in the disk
     accreted = where(gloc[i].loc EQ 2 OR gloc[i].loc EQ 1)   
     accr[i] = accreted[0]
  ENDFOR
  neverinhalo = where(accr EQ -1, complement = everinhalo)
  temp = fltarr(n_elements(accr))
  temp[ everinhalo] = halodat[accr[everinhalo]].z
  IF neverinhalo[0] NE -1 THEN temp[neverinhalo] = 99
  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.accr_rvir_z.fits',/create ;Time when accreted onto the disk

  mass = fltarr(n_elements(accr))
  FOR i = 0L, n_elements(everinhalo) - 1 DO mass[everinhalo[i]] = gpart[everinhalo[i]].mass[accr[everinhalo[i]]] ;uneject particles = 0
  IF NOT keyword_set(nowrite) THEN mwrfits,mass,'grp' + finalid + '.mass_at_accr_rvir.fits',/create  ;Mass when accreted onto the disk
  IF keyword_set(debug) THEN BEGIN
;  histogramp,13.7 - wmap3_lookback(temp)/1e9,weight = mass, binsize = 0.2,xrange = [1e-5,14];,/cum
     histogramp,temp,weight = mass, binsize = 0.2,xrange = [0,10],min = 0, max = 10,xtitle = 'Redshift',ytitle = 'Mass Accreted To Halo' ;,/cum ;,/cum
     stop
  ENDIF

;-------------------- Gas Loss 
  lost        = intarr(n_elements(gloc)) - 1 ;Gas that is in the halo and then out of the halo
  FOR i = 0L, n_elements(lost) - 1 DO BEGIN ;iterate through gas particles that were ever in the disk
     outhalo = where(gloc[i].loc EQ 0)   
     indhalo = where(gloc[i].loc GE 1)
     temp = min(indhalo,firsthalo) ;First timestep that it is in the halo
     firsthalo = indhalo[firsthalo]
;    outhalo = where(gloc[indcoolon[i]].loc EQ 0)
     IF (where(outhalo GT firsthalo))[0] NE -1 THEN lost[i] = outhalo[(where(outhalo GT firsthalo))[0]]
  ENDFOR
  neverlost = where(lost EQ -1, complement = everlost)
  temp = fltarr(n_elements(lost))
  temp[ everlost] = halodat[lost[everlost]].z
  temp[neverlost] = 99
  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.outflow_z.fits',/create ;Time when accreted onto the disk

  mass = fltarr(n_elements(lost))
  FOR i = 0L, n_elements(everlost) - 1 DO mass[everlost[i]] = gpart[everlost[i]].mass[lost[everlost[i]]] ;uneject particles = 0
  IF NOT keyword_set(nowrite) THEN mwrfits,mass,'grp' + finalid + '.mass_at_outflow.fits',/create  ;Mass when accreted onto the disk
  IF keyword_set(debug) THEN BEGIN
;  histogramp,13.7 - wmap3_lookback(temp[everindisk])/1e9,weight = mass[everindisk], binsize = 0.2,xrange = [1e-5,14];,/cum
     histogramp,temp[everindisk],weight = mass[everindisk], binsize = 0.2,xrange = [0,10],xtitle = 'Redshift',ytitle = 'Mass Outflowing From Disk' ;,/cum
     stop
  ENDIF

  reoutflow_iord = [0]           ;Gas that is heated by SN + in the disk then out of the halo
  reoutflow_mass = [0.0]
  reoutflow_z    = [0.0]
  reaccr_iord    = [0]           ;Gas that is heated by SN + in the disk then out of the halo
  reaccr_mass    = [0.0]
  reaccr_z       = [0.0]

  IF 0 THEN BEGIN ;too costly and I don't use this information much
  FOR i = 0L, n_elements(gloc.loc[0]) - 1 DO BEGIN ;iterate through all gas particles
     inhalo  = where(gloc[i].loc EQ 1 OR gloc[i].loc EQ 2)
     outhalo = where(gloc[i].loc EQ 0)
     IF indisk[0] NE -1 THEN BEGIN
        lastoutflow = -1
 ;       count = 0
        WHILE 1 DO BEGIN
           IF (where(inhalo GT lastoutflow))[0] EQ -1 THEN BREAK
           firsthalo = inhalo[(where(inhalo GT lastoutflow))[0]]
           reaccr_iord = [reaccr_iord,gpart[i].iord]
           reaccr_mass = [reaccr_mass,gpart[i].mass[firsthalo]]
           reaccr_z    = [reaccr_z,halodat[firsthalo].z]

           IF (where(outhalo GT firsthalo))[0] EQ -1 THEN BREAK
           lastoutflow = outhalo[(where(outhalo GT firsthalo))[0]]
           reoutflow_iord = [reoutflow_iord,gpart[i].iord]
           reoutflow_mass = [reoutflow_mass,gpart[i].mass[lastoutflow]]
           reoutflow_z    = [reoutflow_z,halodat[lastoutflow].z]
 ;          count = count + 1
 ;          IF count GT 1 THEN stop
        ENDWHILE
    ENDIF
  ENDFOR
  IF n_elements(reaccr_iord) NE 1 THEN BEGIN
     reaccr_iord = reaccr_iord[1:n_elements(reaccr_iord) - 1]
     reaccr_mass = reaccr_mass[1:n_elements(reaccr_mass) - 1]
     reaccr_z    = reaccr_z   [1:n_elements(reaccr_z   ) - 1]
  ENDIF
  IF n_elements(reoutflow_iord) NE 1 THEN BEGIN
     reoutflow_iord = reoutflow_iord[1:n_elements(reoutflow_iord) - 1]
     reoutflow_mass = reoutflow_mass[1:n_elements(reoutflow_mass) - 1]
     reoutflow_z    = reoutflow_z   [1:n_elements(reoutflow_z   ) - 1]
  ENDIF

  IF NOT keyword_set(nowrite) THEN mwrfits,reaccr_iord,'grp' + finalid + '.reaccr_iord.fits',/create
  IF NOT keyword_set(nowrite) THEN mwrfits,reaccr_z   ,'grp' + finalid + '.reaccr_z.fits',/create         
  IF NOT keyword_set(nowrite) THEN mwrfits,reaccr_mass,'grp' + finalid + '.mass_at_reaccr.fits',/create
  IF NOT keyword_set(nowrite) THEN mwrfits,reoutflow_iord,'grp' + finalid + '.reoutflow_iord.fits',/create
  IF NOT keyword_set(nowrite) THEN mwrfits,reoutflow_z   ,'grp' + finalid + '.reoutflow_z.fits',/create         
  IF NOT keyword_set(nowrite) THEN mwrfits,reoutflow_mass,'grp' + finalid + '.mass_at_reoutflow.fits',/create   
  IF keyword_set(debug) THEN BEGIN
     histogramp,reoutflow_z, weight = reoutflow_mass,binsize = 0.2,xrange = [0,10],xtitle = 'Redshift',ytitle = 'Mass (Re)Outflow From Disk' ;,/cum
;  histogramp,13.7 - wmap3_lookback(temp[yesoutflow])/1e9, binsize = 0.2,xrange = [1e-5,14];,/cum
     stop
  ENDIF
  ENDIF

;------------------------------- Ejection ------------------------------
  coolong = fltarr(n_elements(gpart.iord))
  match,gpart.iord,iord,indg,inda ;Match the tracked particles for halo 1 to the last gas array
  coolong[indg] = coolon[inda]*units.timeunit/1e9 
;  indcoolon = where(coolong NE 0)
  it = n_elements(gpart[0].rho) - 1 ;Index of final step
;  dense = where(gpart.rho[it] GT 100) 
;  print,dense[1] 
;  ip = 0
;  inhalo = where(gloc.loc[it] EQ 1)
;  disk = where(gloc.loc[it] EQ 2)
;  allejectfinal = where(coolong NE 0 AND gloc.loc[it] EQ 0) ;Gas that is heated by SN + out of the halo at the last step
;About 95% of the particles that are ejected from the halo at the final timestep and have been heated by SN have been in the disk

  outflow     = intarr(n_elements(gloc)) - 1 ;Gas that is heated by SN and out of the halo
  reaccret    = intarr(n_elements(gloc)) - 1 ;Gas that is heated by SN + in the disk then out of the disk and then in the disk
  eject       = intarr(n_elements(gloc)) - 1 ;Gas that is heated by SN + in the disk then out of the disk                     
  expell      = intarr(n_elements(gloc)) - 1 ;Gas that is heated by SN + in the disk then out of the halo
  expellfinal = intarr(n_elements(gloc)) - 1 ;Gas that is heated by SN + in the disk then out of the halo at the final step 
  indcoolon   = where(coolong NE 0) 

  FOR i = 0L, n_elements(indcoolon) - 1 DO BEGIN ;iterate through gas particles that are heated by SN 
     indisk = where(gloc[indcoolon[i]].loc EQ 2)
     outhalo = where(gloc[indcoolon[i]].loc EQ 0)
     IF indisk[0] NE -1 THEN BEGIN
        temp = min(indisk,firstdisk) ;First timestep it is in the disk
        temp = max(indisk,lastdisk)  ;Last  timestep it is in the disk
        firstdisk = indisk[firstdisk]
        lastdisk  = indisk[lastdisk]
        outdisk = where(gloc[indcoolon[i]].loc EQ 0 OR gloc[indcoolon[i]].loc EQ 1)
        IF (where(outdisk GT firstdisk AND outdisk LT lastdisk))[0] NE -1 THEN reaccret[indcoolon[i]] = (where(outdisk GT firstdisk AND outdisk LT lastdisk))[0]
        IF (where(outdisk GT firstdisk))[0] NE -1 THEN       eject[indcoolon[i]] = outdisk[(where(outdisk GT firstdisk))[0]] ;set the eject array to the timestep when the particle is first out of the disk following being in the disk
        IF (where(outhalo GT firstdisk))[0] NE -1 THEN      expell[indcoolon[i]] = outhalo[(where(outhalo GT firstdisk))[0]] ;set the expell array to the timestep when the particle is first out of the halo following being in the disk
        IF (gloc[indcoolon[i]].loc[it] EQ 0)      THEN expellfinal[indcoolon[i]] = outhalo[(where(outhalo GT firstdisk))[0]] ;set the expellfinal array to the timestep when the particle is first out of the halo following being in the disk 
    ENDIF
    indhalo = where(gloc[indcoolon[i]].loc GE 1)
    temp = min(indhalo,firsthalo) ;First timestep that it is in the halo
    firsthalo = indhalo[firsthalo]
;    outhalo = where(gloc[indcoolon[i]].loc EQ 0)
    IF (where(outhalo GT firsthalo))[0]     NE -1 THEN     outflow[indcoolon[i]] = outhalo[(where(outhalo GT firsthalo))[0]]
  ENDFOR

  nooutflow = where(outflow EQ -1, complement = yesoutflow)
  temp = fltarr(n_elements(outflow))
  temp[yesoutflow] = halodat[outflow[yesoutflow]].z
  temp[nooutflow]  = 99  
;  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.outflow2_z.fits',/create
  temp = lonarr(n_elements(outflow))
  temp[yesoutflow] = gpart[yesoutflow].iord
  temp[nooutflow]  = 0
;  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.outflow2_iord.fits',/create 

;Output the iord, z, and mass of particles being ejected, when ejected
  temp = gpart.iord
  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.eject_iord.fits',/create ;(number in entropy - early;)

  noeject = where(eject EQ -1, complement = yeseject)

  temp = fltarr(n_elements(eject))
  temp[yeseject] = halodat[eject[yeseject]].z
  temp[noeject]  = 99  
  IF keyword_set(debug) THEN BEGIN
;  histogramp,13.7 - wmap3_lookback(temp[yeseject])/1e9, binsize = 0.2,xrange = [1e-5,14];,/cum
     histogramp,temp[yeseject], binsize = 0.2,xrange = [0,10],xtitle = 'Redshift',ytitle = 'Mass Ejected From Disk' ;,/cum
  ENDIF
  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.eject_z.fits',/create         

  temp = fltarr(n_elements(eject))
  FOR i = 0L, n_elements(yeseject) - 1 DO temp[yeseject[i]] = gpart[yeseject[i]].mass[eject[yeseject[i]]] ;uneject particles = 0
  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.mass_at_eject.fits',/create   

  temp = fltarr(n_elements(eject))
  temp[yeseject] = halodat[eject[yeseject] - 1].z
  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.eject_timeindisk.fits',/create ;in redshift (0 = never expelled OR never in disk)

;Output the iord, z, and mass of particles being expelled, when expelled
  noexpell = where(expell EQ -1, complement = yesexpell)

  temp = fltarr(n_elements(expell))
  temp[yesexpell] = halodat[expell[yesexpell]].z
  temp[noexpell]  = 99  
  IF keyword_set(debug) THEN BEGIN
;  histogramp,13.7 - wmap3_lookback(temp[yesexpell])/1e9, binsize = 0.2,xrange = [1e-5,14],/overplot,linestyle = 2;,/cum
     histogramp,temp[yesexpell], binsize = 0.2,xrange = [0,10],/overplot,linestyle = 2 ;,/cum
     legend,['Ejected','Expelled'],linestyle = [0,2]
     stop
  ENDIF

  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.expell_z.fits',/create         

  temp = fltarr(n_elements(expell))
  FOR i = 0L, n_elements(yesexpell) - 1 DO temp[yesexpell[i]] = gpart[yesexpell[i]].mass[expell[yesexpell[i]]] ;uneject particles = 0
  IF NOT keyword_set(nowrite) THEN mwrfits,temp,'grp' + finalid + '.mass_at_expell.fits',/create   

;---------------------- Re-expullsion, and Re-ejection
  coolong = fltarr(n_elements(gpart.iord))
  match,gpart.iord,iord,indg,inda ;Match the tracked particles to the last gas array
  coolong[indg] = coolon[inda]*units.timeunit/1e9 
  it = n_elements(gpart[0].rho) - 1 ;Index of final step

  reeject_iord = [0] ;Gas that is heated by SN + in the disk then out of the disk
  reeject_mass = [0.0]
  reeject_z    = [0.0]
  reexpell_iord = [0] ;Gas that is heated by SN + in the disk then out of the halo
  reexpell_mass = [0.0]
  reexpell_z    = [0.0]
  indcoolon   = where(coolong NE 0,complement = indcooloff) 
  reaccrdisk_iord = [gpart[indcooloff].iord] ;Gas particles that will never be ejected
  reaccrdisk_mass = [accrdisk_mass[indcooloff]]
  reaccrdisk_z    = [accrdisk_z[indcooloff]] 
  indisk = where(reaccrdisk_z NE 99) ;Those of the non-ejected gas particles that are ever in the disk
  reaccrdisk_iord = reaccrdisk_iord[indisk]
  reaccrdisk_mass = reaccrdisk_mass[indisk]
  reaccrdisk_z    = reaccrdisk_z[indisk]


  FOR i = 0L, n_elements(indcoolon) - 1 DO BEGIN ;iterate through gas particles that are heated by SN 
     indisk = where(gloc[indcoolon[i]].loc EQ 2)
 ;    reaccrdisk_iord = [reaccrdisk_iord,gpart[indcoolon[i]].iord]
 ;    reaccrdisk_mass = [reaccrdisk_mass,gpart[indcoolon[i]].mass[indisk[0]]]
 ;    reaccrdisk_z    = [reaccrdisk_z,halodat[indisk[0]].z]
     outdisk = where(gloc[indcoolon[i]].loc EQ 1 OR gloc[indcoolon[i]].loc EQ 0)
     IF indisk[0] NE -1 THEN BEGIN
        lasteject = -1
;        count = 0
        WHILE 1 DO BEGIN
           IF (where(indisk GT lasteject))[0] EQ -1 THEN BREAK ;If the particle is never again in the disk after being ejected
           firstdisk = indisk[(where(indisk GT lasteject))[0]]
           reaccrdisk_iord = [reaccrdisk_iord,gpart[indcoolon[i]].iord]
           reaccrdisk_mass = [reaccrdisk_mass,gpart[indcoolon[i]].mass[firstdisk]]
           reaccrdisk_z    = [reaccrdisk_z,halodat[firstdisk].z]
           IF (where(outdisk GT firstdisk))[0] EQ -1 THEN BREAK ;If the particle is never again ejected
           lasteject = outdisk[(where(outdisk GT firstdisk))[0]]
           reeject_iord = [reeject_iord,gpart[indcoolon[i]].iord]
           reeject_mass = [reeject_mass,gpart[indcoolon[i]].mass[lasteject]]
           reeject_z    = [reeject_z,halodat[lasteject].z]
;           count = count + 1
;           IF count GT 1 THEN stop
        ENDWHILE
    ENDIF
  ENDFOR
  IF n_elements(reeject_iord) NE 1 THEN BEGIN
     reeject_iord = reeject_iord[1:n_elements(reeject_iord) - 1]
     reeject_mass = reeject_mass[1:n_elements(reeject_mass) - 1]
     reeject_z    = reeject_z   [1:n_elements(reeject_z   ) - 1]
  ENDIF

  IF NOT keyword_set(nowrite) THEN mwrfits,reeject_iord,'grp' + finalid + '.reeject_iord.fits',/create
  IF NOT keyword_set(nowrite) THEN mwrfits,reeject_z   ,'grp' + finalid + '.reeject_z.fits',/create         
  IF NOT keyword_set(nowrite) THEN mwrfits,reeject_mass,'grp' + finalid + '.mass_at_reeject.fits',/create   
  IF NOT keyword_set(nowrite) THEN mwrfits,reaccrdisk_iord,'grp' + finalid + '.reaccrdisk_iord.fits',/create
  IF NOT keyword_set(nowrite) THEN mwrfits,reaccrdisk_z   ,'grp' + finalid + '.reaccrdisk_z.fits',/create         
  IF NOT keyword_set(nowrite) THEN mwrfits,reaccrdisk_mass,'grp' + finalid + '.mass_at_reaccrdisk.fits',/create 
  IF keyword_set(debug) THEN BEGIN
;  histogramp,13.7 - wmap3_lookback(temp[yeseject])/1e9, binsize = 0.2,xrange = [1e-5,14];,/cum
     histogramp,reeject_z, weight = reeject_mass,binsize = 0.2,xrange = [0,10],xtitle = 'Redshift',ytitle = 'Mass (Re)Ejected From Disk' ;,/cum
     histogramp,reaccrdisk_z,weight = reaccrdisk_mass,binsize = 0.2,/overplot,color = 100
     histogramp,accrdisk_z,weight = accrdisk_mass,binsize = 0.2,/overplot,color = 50 
     stop
  ENDIF

  FOR i = 0L, n_elements(indcoolon) - 1 DO BEGIN ;iterate through gas particles that are heated by SN 
     indisk = where(gloc[indcoolon[i]].loc EQ 2)
     outdisk = where(gloc[indcoolon[i]].loc EQ 1)
     outhalo = where(gloc[indcoolon[i]].loc EQ 0)
     IF indisk[0] NE -1 THEN BEGIN
        lastexpell = -1
        WHILE 1 DO BEGIN
           IF (where(indisk GT lastexpell))[0] EQ -1 THEN BREAK
           firstdisk = indisk[(where(indisk GT lastexpell))[0]]
           IF (where(outhalo GT firstdisk))[0] EQ -1 THEN BREAK
           lastexpell = outhalo[(where(outhalo GT firstdisk))[0]]
           reexpell_iord = [reexpell_iord,gpart[indcoolon[i]].iord]
           reexpell_mass = [reexpell_mass,gpart[indcoolon[i]].mass[lastexpell]]
           reexpell_z    = [reexpell_z,halodat[lastexpell].z]
;           stop
        ENDWHILE
    ENDIF
  ENDFOR
  IF n_elements(reexpell_iord) NE 1 THEN BEGIN
     reexpell_iord = reexpell_iord[1:n_elements(reexpell_iord) - 1]
     reexpell_mass = reexpell_mass[1:n_elements(reexpell_mass) - 1]
     reexpell_z    = reexpell_z   [1:n_elements(reexpell_z   ) - 1]
  ENDIF

  IF NOT keyword_set(nowrite) THEN mwrfits,reexpell_iord,'grp' + finalid + '.reexpell_iord.fits',/create
  IF NOT keyword_set(nowrite) THEN mwrfits,reexpell_z   ,'grp' + finalid + '.reexpell_z.fits',/create         
  IF NOT keyword_set(nowrite) THEN mwrfits,reexpell_mass,'grp' + finalid + '.mass_at_reexpell.fits',/create   
  IF keyword_set(debug) THEN BEGIN
;  histogramp,13.7 - wmap3_lookback(temp[yeseject])/1e9, binsize = 0.2,xrange = [1e-5,14];,/cum
     histogramp,reexpell_z, weight = reexpell_mass,binsize = 0.2,xrange = [0,10],xtitle = 'Redshift',ytitle = 'Mass (Re)Expelled From Disk' ;,/cum
     stop
  ENDIF

  IF keyword_set(plots) THEN BEGIN
     loadct,39
;----------------------------- Plots -----------------------------
;  ip = (where(reaccret    NE -1))[0]
;  ip = (where(eject       NE -1))[0]
;  ip = (where(expell      NE -1))[0]
  ip = (where(expellfinal NE -1))[0]

  ipdisk = where(gloc[ip].loc EQ 2)
  iphalo = where(gloc[ip].loc EQ 1 OR gloc[ip].loc EQ 2)
  sntime = max(halodat[where(halodat.time LE coolong[ip])].time)
  ipsn = where(halodat.time EQ sntime)
  IF NOT keyword_set(outplot) THEN window,2,xsize = 800, ysize = 500 ELSE  device,filename = 'example.fb.phase.eps',bits_per_pixel= 8,/times,ysize=7,xsize=7,/inch,/color,xoffset = 3, yoffset = 5;,ytitle = 'Temperature [K]',xtitle = 'Density [amu/cc]'
  plot,gpart.rho[it],gpart.temp[it],psym = 3,/xlog,/ylog,xrange = [1e-7,1e4],yrange = [100,1e6],ytitle = 'Temperature [K]',xtitle = 'Density [amu/cc]',xstyle = 1
  oplot,gpart[inhalo].rho[it],gpart[inhalo].temp[it],psym = 3, color = 100
  oplot,gpart[disk].rho[it],gpart[disk].temp[it],psym = 3, color = 50
  oplot,gpart[ip].rho,gpart[ip].temp,psym = -4,color = 190
  IF n_elements(iphalo) EQ 1 THEN $
     oplot,[gpart[ip].rho[iphalo],gpart[ip].rho[iphalo]],[gpart[ip].temp[iphalo],gpart[ip].temp[iphalo]],psym = 4,color = 210 ELSE $
     oplot,gpart[ip].rho[iphalo],gpart[ip].temp[iphalo],psym = 4,color = 210
  IF n_elements(ipdisk) EQ 1 THEN $
     oplot,[gpart[ip].rho[ipdisk],gpart[ip].rho[ipdisk]],[gpart[ip].temp[ipdisk],gpart[ip].temp[ipdisk]],psym = 4,color = 254 ELSE $
     oplot,gpart[ip].rho[ipdisk],gpart[ip].temp[ipdisk],psym = 4,color = 254
  oplot,[gpart[ip].rho[ipsn],gpart[ip].rho[ipsn]],[gpart[ip].temp[ipsn],gpart[ip].temp[ipsn]],psym = symcat(14),color = 254,symsize = 3
  legend,['Not in halo',' ','In halo',' ','In Disk',' ','History'],color = [0,190,100,210,50,254,190],linestyle = [1,0,1,0,1,0,0],psym = [-3,4,-3,4,-3,4,-3],/bottom,/left

  IF NOT keyword_set(outplot) THEN BEGIN
     stop
     window,6,xsize = 500, ysize = 500
  ENDIF ELSE BEGIN
     device,/close
     device,filename = 'example.fb.fo.eps',bits_per_pixel= 8,/times,ysize=7,xsize=7,/inch,/color,xoffset = 3, yoffset = 5 ;,ytitle = 'Y [kpc]',xtitle = 'X [kpc]'
  ENDELSE
;  verydense = where(gpart.rho[it] GE 100)
  range = 40
  plot,gpart.x[it],gpart.y[it],psym = 3,xrange = [-1*range,range], yrange = [-1*range,range],ytitle = 'Y [kpc]',xtitle = 'X [kpc]'
  oplot,gpart[inhalo].x[it],gpart[inhalo].y[it],psym = 3, color = 100
  oplot,gpart[disk].x[it],gpart[disk].y[it],psym = 3, color = 50
;  oplot,gpart[verydense].x[it],gpart[verydense].y[it],psym = 3, color = 20
  oplot,gpart[ip].x,gpart[ip].y,psym = -4,color = 190
  IF n_elements(iphalo) EQ 1 THEN $
     oplot,[gpart[ip].x[iphalo],gpart[ip].x[iphalo]],[gpart[ip].y[iphalo],gpart[ip].y[iphalo]],psym = 4,color = 210 ELSE $
     oplot,gpart[ip].x[iphalo],gpart[ip].y[iphalo],psym = 4,color = 210
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
  IF n_elements(iphalo) EQ 1 THEN $
     oplot,[gpart[ip].x[iphalo],gpart[ip].x[iphalo]],[gpart[ip].z[iphalo],gpart[ip].z[iphalo]],psym = 4,color = 210 ELSE $
     oplot,gpart[ip].x[iphalo],gpart[ip].z[iphalo],psym = 4,color = 210
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
