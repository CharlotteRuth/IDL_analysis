
PRO ejectz_v_mass_data,dir,file,halo,z_bins,zmax,mlbin
  zmax = 3.130
  print,file + '.' + halo
  nz = n_elements(z_bins)
  relostmass = fltarr(nz) ;Mass of gas lost from disk over entire history
  relostmassmet = fltarr(nz) ;Mass of metals lost from disk over entire history
  relostmassr = fltarr(nz) ;Mass of gas lost from disk over entire history
  relostmassmetr = fltarr(nz) ;Mass of metals lost from disk over entire history
  reejectmass  = fltarr(nz) ;Mass of gas (re)ejected over entire history
  reejectmassmet = fltarr(nz) ;Mass of metals (re)ejected over entire history 
  reejectmassr  = fltarr(nz) ;Mass of gas (re)ejected in redshift bins
  reejectmassmetr = fltarr(nz) ;Mass of metals (re)ejected in redshift bins
  reexpellmass = fltarr(nz) ;Mass of gas (re)expelled over entire history
  reexpellmassmet = fltarr(nz) ;Mass of metals (re)expelled over entire history
  reexpellmassr = fltarr(nz) ;Mass of gas (re)expelled in redshift bins
  reexpellmassmetr = fltarr(nz) ;Mass of metals (re)expelled in redshift bins
  diskgmass = fltarr(nz) ;Mass reaccreted onto disk after ejection
  diskgmassmet = fltarr(nz) ;Mass of metals reaccreted onto disk after ejection
  diskgmassr = fltarr(nz) ;Mass reaccreted onto disk after ejection
  diskgmassmetr = fltarr(nz) ;Mass of metals reaccreted onto disk after ejection
  diskgmass_lost = fltarr(nz) ;Mass of accretions onto the disk
  diskgmass_lostmet = fltarr(nz) ;Metal mass of accretions onto the disk
  diskgmass_lostr = fltarr(nz) ;Mass of accretions onto the disk
  diskgmass_lostmetr = fltarr(nz) ;Metal mass of accretions onto the disk
  diskgmass_expell = fltarr(nz) ;Mass reaccreted onto disk after expulsion
  diskgmass_expellmet = fltarr(nz);Mass of metals reaccreted onto disk after expulsion
  diskgmass_expellr = fltarr(nz) ;Mass reaccreted onto disk after expulsion
  diskgmass_expellmetr = fltarr(nz);Mass of metals reaccreted onto disk after expulsion
  halogmass = fltarr(nz) ;Mass accreted onto the halo
  halogmassmet = fltarr(nz) ;Metal mass accreted onto the halo
  sfmassr = fltarr(nz)

  ageUniverse = wmap3_lookback(1000)
  t_bins = (ageUniverse - wmap3_lookback(z_bins))/1e9
  tmin = (z_to_t(zmax))[0]

  units = tipsyunits(dir + file + '.param')
  sldata = mrdfits(dir + 'starlog.' + halo + '.fits',1,/silent)
  sldata.timeform = sldata.timeform*units.timeunit/1e9
  sfmass_total = total(sldata[where(sldata.timeform GT tmin)].massform)*units.massunit
;       match,iordstar,sldata.iorderstar,ind1,ind2
;       sfmass_z0 = total(s[ind1].mass)*units.massunit
;       smass_accr[i] = smass[i] - sfmass_z0
  readcol,dir + '/grp' + halo + '.sfmetals.txt',zstars,oxstars,festars
  sfmassmet_total = total(zstars)

  accrm = mrdfits(dir + '/grp' + halo + '.mass_at_accr_rvir.fits',0,/silent) 
  allgmass_total = total(accrm) ;Total mass accreted onto the virial radius

;  diskearlym = mrdfits(dir + '/grp' + halo + '.earlydisk_mass.fits',0,/silent)
;  diskearlyi = mrdfits(dir + '/grp' + halo + '.earlydisk_iord.fits',0,/silent)
;;  diskearlymet = mrdfits(dir + '/grp' + halo + '.earlydisk_zmetals.fits',0,/silent) ;All these values seem to be zero?
;  diskearlymet = diskearlym * 0
;  diskearlyz = fltarr(n_elements(diskearlyi)) + 99

;  haloearlym = mrdfits(dir + '/grp' + halo + '.earlyhalo_mass.fits',0,/silent)
;  haloearlyi = mrdfits(dir + '/grp' + halo + '.earlyhalo_iord.fits',0,/silent)
;  haloearlymet = mrdfits(dir + '/grp' + halo + '.earlyhalo_zmetals.fits',0,/silent)
;  haloearlyz = fltarr(n_elements(haloearlyi)) + 99
  
; Gas that is accreted for the first time or after being "lost"
  diskaccrm_lost = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdiskall.fits',0,/silent)
  diskaccrm_lost[where(diskaccrm_lost LE 0)] = units.istarmass
  diskaccrz_lost = mrdfits(dir + '/grp' + halo + '.reaccrdiskall_z.fits',0,/silent)
  diskaccri_lost = mrdfits(dir + '/grp' + halo + '.reaccrdiskall_iord.fits',0,/silent)
  diskaccrhist_lost = mrdfits(dir + '/grp' + halo + '.reaccrdiskall_history.fits',1,/silent)
  diskaccrt_lost = z_to_t(diskaccrz_lost)
;  diskregmassi_lost = [diskearlyi,diskaccri_lost];[where(diskaccrz_lost LT zmax)]]
;  diskregmassm_lost = [diskearlym,diskaccrm_lost];[where(diskaccrz_lost LT zmax)]]
;  diskregmassz_lost = [diskearlyz,diskaccrz_lost];[where(diskaccrz_lost LT zmax)]]
;  diskregmassmet_lost = [diskearlymet,diskaccrhist_lost.metallicity];[where(diskaccrz_lost LT zmax)].metallicity]
  diskregmassi_lost = [diskaccri_lost[where(diskaccrz_lost LT zmax)]]
  diskregmassm_lost = [diskaccrm_lost[where(diskaccrz_lost LT zmax)]]
  diskregmassz_lost = [diskaccrz_lost[where(diskaccrz_lost LT zmax)]]
  diskregmassmet_lost = [diskaccrhist_lost[where(diskaccrz_lost LT zmax)].metallicity]
  diskregmasst_lost = z_to_t(diskregmassz_lost)
  ind_diskaccri_lost_uniq = rem_dup(diskregmassi_lost);
  diskaccri_lost_uniq = diskregmassi_lost[ind_diskaccri_lost_uniq]
  diskaccrz_lost_uniq = diskregmassz_lost[ind_diskaccri_lost_uniq]
;  diskaccrm_lost_uniq = diskregmassm_lost[ind_diskaccri_lost_uniq]

  inflowm = mrdfits(dir + '/grp' + halo + '.mass_at_reaccr.fits',/silent) ;'mass_at_accr_rvir';'.mass_at_accr.fits'
  inflowi = mrdfits(dir + '/grp' + halo + '.reaccr_iord.fits',/silent)
  inflow_data = mrdfits(dir + '/grp' + halo + '.reaccr_history.fits',1,/silent)
  ind_inflow_uniq = rem_dup(inflowi) ;uniq(inflowi,sort(inflowi)) 
  inflowm_uniq = inflowm[ind_inflow_uniq]
  inflowi_uniq = inflowi[ind_inflow_uniq]
  inflow_data_uniq = inflow_data[ind_inflow_uniq]

  match2,inflowi_uniq,diskaccri_lost_uniq,ind_halo,ind_disk  
  IF n_elements(diskaccri_lost_uniq) NE n_elements(where(ind_halo NE -1)) THEN stop ;Sanity check that all particle accreted to disk are also accreted to the halo
  diskaccrm_lost_uniq = inflowm_uniq[ind_halo[where(ind_halo NE -1)]]
  diskaccrmet_lost_uniq = inflow_data_uniq[ind_halo[where(ind_halo NE -1)]].metallicity
  diskregmassm_lost[ind_diskaccri_lost_uniq] = diskregmassm_lost[ind_diskaccri_lost_uniq] - diskaccrm_lost_uniq
  diskregmassmet_lost[ind_diskaccri_lost_uniq] = diskregmassmet_lost[ind_diskaccri_lost_uniq] - diskaccrmet_lost_uniq 

;  diskregmassm_lost[ind_diskaccri_lost_uniq] = 0 ;Set mass and metallicity of first accretion to zero. That way, diskregmassm_lost is only for stuff that is reaccreted after being lost
;  diskregmassmet_lost[ind_diskaccri_lost_uniq] = 0

;--------------- Note that these include the initial accretion as well
;                as the accretion after ejection ------------
;  diskaccrm = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdisk.fits',0,/silent)
;  diskaccrm[where(diskaccrm LE 0)] = units.istarmass
;  diskaccrz = mrdfits(dir + '/grp' + halo + '.reaccrdisk_z.fits',0,/silent)
;  diskaccri = mrdfits(dir + '/grp' + halo + '.reaccrdisk_iord.fits',0,/silent)
;---------------------------------

;Gas that is reaccreted after ejection
  diskaccrhist = mrdfits(dir + '/grp' + halo + '.reaccrdisk_history.fits',1,/silent) ;Only reaccretion after ejection
  diskaccri = [diskaccrhist.iord]
  diskaccrm = [diskaccrhist.mass]
  diskaccrz = [diskaccrhist.red]
  diskaccrmet = [diskaccrhist.metallicity]
  diskaccrt = z_to_t(diskaccrz)
 ; diskregmassi = [diskearlyi,diskaccri];[where(diskaccrz LT zmax)]] 
 ; diskregmassm = [diskearlym,diskaccrm];[where(diskaccrz LT zmax)]] 
 ; diskregmassmet = [diskearlymet,diskaccrmet];[where(diskaccrz LT zmax)]]
 ; diskregmassz = [diskearlyz,diskaccrz];[where(diskaccrz LT zmax)]]
  diskregmassi = [diskaccri[where(diskaccrz LT zmax)]] 
  diskregmassm = [diskaccrm[where(diskaccrz LT zmax)]] 
  diskregmassmet = [diskaccrmet[where(diskaccrz LT zmax)]]
  diskregmassz = [diskaccrz[where(diskaccrz LT zmax)]]
  diskregmasst = [diskaccrt[where(diskaccrz LT zmax)]]

;Gas that is (re)accreted onto the halo  
  haloaccrm = mrdfits(dir + '/grp' + halo + '.mass_at_reaccr.fits',0,/silent)
  haloaccrz = mrdfits(dir + '/grp' + halo + '.reaccr_z.fits',0,/silent)
  haloaccri = mrdfits(dir + '/grp' + halo + '.reaccr_iord.fits',0,/silent)
  haloaccrhist = mrdfits(dir + '/grp' + halo + '.reaccr_history.fits',1,/silent)
  haloaccrt = z_to_t(diskaccrz)  

;Unique accretions onto the halo
  ind_haloaccri_uniq = rem_dup(haloaccri) ;
  haloaccri_uniq = haloaccri[ind_haloaccri_uniq]
  haloaccrz_uniq = haloaccrz[ind_haloaccri_uniq]
  haloaccrm_uniq = haloaccrm[ind_haloaccri_uniq]
  haloaccrmet_uniq = haloaccrhist[ind_haloaccri_uniq].metallicity

;  haloaccrm_clip = haloaccrm[where(haloaccrz LT zmax)]
;  haloaccri_clip = haloaccri[where(haloaccrz LT zmax)]
;  haloaccri = [haloearlyi,haloaccri];[where(haloaccrz LT zmax)]]
;  haloaccrm = [haloearlym,haloaccrm];[where(haloaccrz LT zmax)]]
;  haloaccrmet = [haloearlymet,haloaccrhist.metallicity];[where(haloaccrz LT zmax)].metallicity]
;  haloaccrz = [fltarr(n_elements(haloearlymet))+max(haloaccrz),haloaccrz];[where(haloaccrz LT zmax)]
  haloaccri = [haloaccri[where(haloaccrz LT zmax)]]
  haloaccrm = [haloaccrm[where(haloaccrz LT zmax)]]
  haloaccrmet = [haloaccrhist[where(haloaccrz LT zmax)].metallicity]
  haloaccrz = [haloaccrz[where(haloaccrz LT zmax)]]
 
;Unique accretion onto disk
;Select for those particles the first time they are accreted onto the disk
;  ind_diskaccri_lost_uniq = rem_dup(diskregmassi_lost);
;  diskaccri_lost_uniq = diskregmassi_lost[ind_diskaccri_lost_uniq]
;  diskaccrz_lost_uniq = diskregmassz_lost[ind_diskaccri_lost_uniq]
;  diskaccrm_lost_uniq = diskregmassm_lost[ind_diskaccri_lost_uniq]
; Metallicity at accretion onto the disk
;  diskaccrmet_lost_uniq = diskregmassmet_lost[ind_diskaccri_lost_uniq]
; Metallicity when accreted to the halo
  match2,haloaccri_uniq,diskaccri_lost_uniq,ind_halo,ind_disk  
  IF n_elements(diskaccri_lost_uniq) NE n_elements(where(ind_halo NE -1)) THEN stop ;Sanity check that all particle accreted to disk are also accreted to the halo
  diskaccrmet_lost_uniq = haloaccrmet_uniq[ind_halo[where(ind_halo NE -1)]] 
  
;Only select unique accretions since z = zmax
  diskgaccri_lost_uniq = diskaccri_lost_uniq[where(diskaccrz_lost_uniq LT zmax)]
  diskgaccrm_lost_uniq = diskaccrm_lost_uniq[where(diskaccrz_lost_uniq LT zmax)]
  diskgaccrmet_lost_uniq = diskaccrmet_lost_uniq[where(diskaccrz_lost_uniq LT zmax)]
  diskgaccrz_lost_uniq = diskaccrz_lost_uniq[where(diskaccrz_lost_uniq LT zmax)]
  diskgmass_total = total(diskaccrm_lost_uniq) ;Total unique mass accreted
  diskgmet_total = total(diskaccrmet_lost_uniq*diskaccrm_lost_uniq)

  haloaccri_uniq = haloaccri_uniq[where(haloaccrz_uniq LT zmax)]
  haloaccrm_uniq = haloaccrm_uniq[where(haloaccrz_uniq LT zmax)]
  haloaccrmet_uniq = haloaccrmet_uniq[where(haloaccrz_uniq LT zmax)]
  haloaccrz_uniq = haloaccrz_uniq[where(haloaccrz_uniq LT zmax)]
  halogmass_total = total(haloaccrm);Total unique accretions onto the halo

;Gas that leaves the disk
  relostm = mrdfits(dir + '/grp' + halo + '.mass_at_relost.fits',0,/silent) 
  relosti = mrdfits(dir + '/grp' + halo + '.relost_iord.fits',0,/silent)
  relostz = mrdfits(dir + '/grp' + halo + '.relost_z.fits',0,/silent)
  relosthist = mrdfits(dir + '/grp' + halo + '.relost_history.fits',1,/silent)
  relostmet = relosthist.metallicity
  relosti = [relosti[where(relostz LT zmax)]]
  relostm = [relostm[where(relostz LT zmax)]]
  relostmet = [relostmet[where(relostz LT zmax)]]
  relostz = [relostz[where(relostz LT zmax)]]
  relostt = z_to_t(relostz)
;  relost_total = total(relostm[uniq(relosti,sort(relosti))])
;  relostmet_total=total(relosthist[uniq(relosti,sort(relosti))].metallicity*relostm[uniq(relosti,sort(relosti))])
  relost_total = total(relostm[rem_dup(relosti)])
  relostmet_total=total(relostmet[rem_dup(relosti)]*relostm[rem_dup(relosti)])

;Gas that is (re)ejected from the disk
  reejecti = mrdfits(dir + '/grp' + halo + '.reeject_iord.fits',0,/silent)
  reejectm = mrdfits(dir + '/grp' + halo + '.mass_at_reeject.fits',0,/silent)
  reejectz = mrdfits(dir + '/grp' + halo + '.reeject_z.fits',0,/silent)
  reejecthist = mrdfits(dir + '/grp' + halo + '.reeject_halo.fits',1,/silent)
  reejectmet = reejecthist.metallicity
;  reeject_total = total(reejectm[uniq(reejecti,sort(reejecti))])
;  reejectmet_total = total(reejecthist[uniq(reejecti,sort(reejecti))].metallicity*reejectm[uniq(reejecti,sort(reejecti))])
  reejecti = [reejecti[where(reejectz LT zmax)]]
  reejectm = [reejectm[where(reejectz LT zmax)]]
  reejectmet = [reejectmet[where(reejectz LT zmax)]]
  reejectz = [reejectz[where(reejectz LT zmax)]]
  reeject_total = total(reejectm[rem_dup(reejecti)])
  reejectmet_total = total(reejectmet[rem_dup(reejecti)]*reejectm[rem_dup(reejecti)])
  reejectt = z_to_t(reejectz)

;Gas that is (re)expelled from the disk
  reexpelli = mrdfits(dir + '/grp' + halo + '.reexpell_iord.fits',0,/silent)
  reexpellm = mrdfits(dir + '/grp' + halo + '.mass_at_reexpell.fits',0,/silent) 
  reexpellz = mrdfits(dir + '/grp' + halo + '.reexpell_z.fits',0,/silent)
  reexpellt = z_to_t(reexpellz)
  eject_expell_iord = mrdfits(dir+'grp'+halo+'.eject_expell_iord.fits')
  reexpellhist = reejecthist[eject_expell_iord]
  reexpellmet = reexpellhist.metallicity;[where(reexpellz LT zmax)].metallicity
  reexpelli = [reexpelli[where(reexpellz LT zmax)]]
  reexpellm = [reexpellm[where(reexpellz LT zmax)]]
  reexpellz = [reexpellz[where(reexpellz LT zmax)]]
  reexpell_total = total(reexpellm[rem_dup(reexpelli)]) ;total(reexpellm[uniq(reexpelli,sort(reexpelli))])  
;  reexpellmet_total = total(reexpellhist[uniq(reexpelli,sort(reexpelli))].metallicity*reexpellm[uniq(reexpelli,sort(reexpelli))])
  reexpellmet_total = total(reexpellmet[rem_dup(reexpelli)]*reexpellm[rem_dup(reexpelli)])

;Gas that is reaccreted after expullsion
  match,diskaccrhist.iord,reexpelli,ind3,ind4 ;Find iords of all gas particles that were expelled and accreted a second time (i.e., were potentiall reaccreted to the halo after being expelled
  ind_reinflow_expell=[0]
  FOR j = 0,n_elements(reexpelli[ind4]) -1 DO BEGIN
     IF (where(diskaccrhist.red LT reexpellz[ind4[j]] AND diskaccrhist.iord EQ reexpelli[ind4[j]]))[0] NE -1 THEN $ ;If there is an instance of reaccretion after expulsion
        ind_reinflow_expell = [ind_reinflow_expell,(where(diskaccrhist.red LT reexpellz[ind4[j]] AND diskaccrhist.iord EQ reexpelli[ind4[j]]))[0]] ;Save the index to the first one
  ENDFOR
  IF n_elements(ind_reinflow_expell) NE 1 THEN BEGIN
     ind_reinflow_expell = ind_reinflow_expell[1:n_elements(ind_reinflow_expell)-1]
;     diskaccrm_expell = [diskearlym,diskaccrhist[ind_reinflow_expell].mass]
;     diskaccrmet_expell = [diskearlymet,diskaccrhist[ind_reinflow_expell].metallicity]
 ;    diskaccrz_expell = [diskearlyz,diskaccrhist[ind_reinflow_expell].red]
     diskaccrm_expell = [diskaccrhist[ind_reinflow_expell].mass]
     diskaccrmet_expell = [diskaccrhist[ind_reinflow_expell].metallicity]
     diskaccrz_expell = [diskaccrhist[ind_reinflow_expell].red]
     diskaccrt_expell = z_to_t(diskaccrz_expell)
  ENDIF ELSE BEGIN
     diskaccrm_expell = [0]
     diskaccrmet_expell = [0]
     diskaccrz_expell = [0]
     diskaccrt_expell = [0]
  ENDELSE

;Accreted stars
  smass_accrm = mrdfits(dir + '/grp' + halo + '.accrstars_m.fits',/silent)
  smass_accrz = mrdfits(dir + '/grp' + halo + '.accrstars_z.fits',/silent)
  smass_accrm = smass_accrm[where(smass_accrz NE 99)]
  smass_accrz = smass_accrz[where(smass_accrz NE 99)]
;  smass_accrt = fltarr(nz) ;z_to_t(smass_accrz)
  smass_accr_total = total(smass_accrm)

;********* Totals of unique accretions since the start ***********
;Total mass accreted onto the virial radius (not necessarily unique)
;Total pristine accretion onto the halo
;Total pristine accretion onto the disk since z = zmax
;Total unique mass loss from the disk since z = zmax
;Total unique ejection since z = zmax
;Total mass of stars formed in the main halo
;Total metals formed by stars in main halo
  openw,lun,dir + '/grp' + halo + '.ejectz_quant.txt',/get_lun
  print,     allgmass_total,halogmass_total,float(diskgmass_total),float(diskgmet_total),float(relost_total),float(relostmet_total),reeject_total,float(reejectmet_total),reexpell_total,float(reexpellmet_total),float(sfmass_total),sfmassmet_total,smass_accr_total,format='(13E14)'
  printf,lun,allgmass_total,halogmass_total,float(diskgmass_total),float(diskgmet_total),float(relost_total),float(relostmet_total),reeject_total,float(reejectmet_total),reexpell_total,float(reexpellmet_total),float(sfmass_total),sfmassmet_total,smass_accr_total,format='(13E14)'  
  FOR iz = 0, n_elements(z_bins) - 1 DO BEGIN
     IF z_bins[iz] LT 0.1 THEN mlbin2 = mlbin ELSE mlbin2 = mlbin/2
;     diskaccr_prist_mass[iz]=total([    where(diskregmassz_lost     GE z_bins[iz] AND diskregmassz_lost LT zmax)])
     relostmass[iz]       = total(relostm[      where(relostz      GE z_bins[iz] AND relostz   LT zmax)])
     relostmassmet[iz]    = total(relostm[      where(relostz      GE z_bins[iz] AND relostz   LT zmax)]*relostmet[      where(relostz      GE z_bins[iz] AND relostz   LT zmax)])
     relostmassr[iz]       = total(relostm[     where(relostt      GE t_bins[iz] - mlbin2 AND relostt      LE t_bins[iz] + mlbin2 AND relostz LT zmax)])/mlbin
     relostmassmetr[iz]    = total(relostm[     where(relostt      GE t_bins[iz] - mlbin2 AND relostt      LE t_bins[iz] + mlbin2 AND relostz LT zmax)]*relostmet[     where(relostt     GE t_bins[iz] - mlbin2 AND relostt      LE t_bins[iz] + mlbin2 AND relostz LT zmax)])/mlbin    
     reejectmass[iz]      = total(reejectm[     where(reejectz     GE z_bins[iz] AND reejectz  LT zmax)])
     reejectmassmet[iz]   = total(reejectm[     where(reejectz     GE z_bins[iz] AND reejectz  LT zmax)]*reejectmet[     where(reejectz     GE z_bins[iz] AND reejectz  LT zmax)])
     reejectmassr[iz]     = total(reejectm[     where(reejectt     GE t_bins[iz] - mlbin2 AND reejectt      LE t_bins[iz] + mlbin2 AND reejectz LT zmax)])/mlbin
     reejectmassmetr[iz]  = total(reejectm[     where(reejectt     GE t_bins[iz] - mlbin2 AND reejectt      LE t_bins[iz] + mlbin2 AND reejectz LT zmax)]*reejectmet[     where(reejectt     GE t_bins[iz] - mlbin2 AND reejectt      LE t_bins[iz] + mlbin2 AND reejectz LT zmax)])/mlbin
     reexpellmass[iz]     = total(reexpellm[    where(reexpellz    GE z_bins[iz] AND reexpellz LT zmax)])
     reexpellmassmet[iz]  = total(reexpellm[    where(reexpellz    GE z_bins[iz] AND reexpellz LT zmax)]*reexpellmet[    where(reexpellz    GE z_bins[iz] AND reexpellz LT zmax)])
     reexpellmassr[iz]    = total(reexpellm[    where(reexpellt    GE t_bins[iz] - mlbin2 AND reexpellt     LE t_bins[iz] + mlbin2 AND reexpellz LT zmax)])/mlbin
     reexpellmassmetr[iz] = total(reexpellm[    where(reexpellt    GE t_bins[iz] - mlbin2 AND reexpellt     LE t_bins[iz] + mlbin2 AND reexpellz LT zmax)]*reexpellmet[    where(reexpellt    GE t_bins[iz] - mlbin2 AND reexpellt     LE t_bins[iz] + mlbin2 AND reexpellz LT zmax)])/mlbin
;Gas reaccreted after being lost
     diskgmass_lost[iz]   = total(diskregmassm_lost[    where(diskregmassz_lost     GE z_bins[iz] AND diskregmassz_lost LT zmax)])
     diskgmass_lostmet[iz]= total(diskregmassm_lost[    where(diskregmassz_lost     GE z_bins[iz] AND diskregmassz_lost LT zmax)]*diskregmassmet_lost[    where(diskregmassz_lost     GE z_bins[iz] AND diskregmassz_lost LT zmax)])
     diskgmass_lostr[iz]   = total(diskregmassm_lost[    where(diskregmasst_lost    GE t_bins[iz] - mlbin2 AND diskregmasst_lost     LE t_bins[iz] + mlbin2 AND diskregmassz_lost LT zmax)])/mlbin
     diskgmass_lostmetr[iz]= total(diskregmassm_lost[    where(diskregmasst_lost    GE t_bins[iz] - mlbin2 AND diskregmasst_lost     LE t_bins[iz] + mlbin2 AND diskregmassz_lost LT zmax)]*diskregmassmet_lost[    where(diskregmasst_lost    GE t_bins[iz] - mlbin2 AND diskregmasst_lost     LE t_bins[iz] + mlbin2 AND diskregmassz_lost LT zmax)])/mlbin
;Gas either accreted for the first time or after being ejected
     diskgmass[iz]        = total(diskregmassm[    where(diskregmassz     GE z_bins[iz] AND diskregmassz LT zmax)])
     diskgmassmet[iz]     = total(diskregmassm[    where(diskregmassz     GE z_bins[iz] AND diskregmassz LT zmax)]*diskregmassmet[    where(diskregmassz     GE z_bins[iz] AND diskregmassz LT zmax)])
     diskgmassr[iz]        = total(diskregmassm[    where(diskregmasst    GE t_bins[iz] - mlbin2 AND diskregmasst     LE t_bins[iz] + mlbin2 AND diskregmassz LT zmax)])/mlbin
     diskgmassmetr[iz]     = total(diskregmassm[    where(diskregmasst    GE t_bins[iz] - mlbin2 AND diskregmasst     LE t_bins[iz] + mlbin2 AND diskregmassz LT zmax)]*reexpellmet[    where(diskregmasst    GE t_bins[iz] - mlbin2 AND diskregmasst     LE t_bins[iz] + mlbin2 AND diskregmassz LT zmax)])/mlbin
;Gas accreted to the DISK for the first time or after being expelled
     diskgmass_expell[iz] = total(diskaccrm_expell[    where(diskaccrz_expell     GE z_bins[iz] AND diskaccrz_expell LT zmax)])
     diskgmass_expellmet[iz] = total(diskaccrm_expell[    where(diskaccrz_expell     GE z_bins[iz] AND diskaccrz_expell LT zmax)]*diskaccrmet_expell[    where(diskaccrz_expell     GE z_bins[iz] AND diskaccrz_expell LT zmax)])
     diskgmass_expellr[iz] = total(diskaccrm_expell[    where(diskaccrt_expell    GE t_bins[iz] - mlbin2 AND diskaccrt_expell     LE t_bins[iz] + mlbin2 AND diskaccrz_expell LT zmax)])/mlbin
     diskgmass_expellmetr[iz] = total(diskaccrm_expell[    where(diskaccrt_expell    GE t_bins[iz] - mlbin2 AND diskaccrt_expell     LE t_bins[iz] + mlbin2 AND diskaccrz_expell LT zmax)]*diskaccrmet_expell[    where(diskaccrt_expell    GE t_bins[iz] - mlbin2 AND diskaccrt_expell     LE t_bins[iz] + mlbin2 AND diskaccrz_expell LT zmax)])/mlbin

     halogmass[iz]        = total(haloaccrm[    where(haloaccrz    GE z_bins[iz] AND haloaccrz LT zmax)])
     halogmassmet[iz]     = total(haloaccrm[    where(haloaccrz    GE z_bins[iz] AND haloaccrz LT zmax)]*haloaccrmet[    where(haloaccrz    GE z_bins[iz] AND haloaccrz LT zmax)])
     sfmassr[iz] = total(sldata[where(sldata.timeform GE t_bins[iz] - mlbin2 AND sldata.timeform LE t_bins[iz] + mlbin2 AND sldata.timeform GT tmin)].massform)*units.massunit/mlbin
     printf,lun,z_bins[iz],relostmass[iz],relostmassmet[iz],relostmassr[iz],relostmassmetr[iz],reejectmass[iz],reejectmassmet[iz],reejectmassr[iz],reejectmassmetr[iz],reexpellmass[iz],reexpellmassmet[iz],reexpellmassr[iz],reexpellmassmetr[iz],float(diskgmass_lost[iz]),diskgmass_lostmet[iz],float(diskgmass_lostr[iz]),diskgmass_lostmetr[iz],diskgmass[iz],diskgmassmet[iz],diskgmassr[iz],diskgmassmetr[iz],diskgmass_expell[iz],diskgmass_expellmet[iz],diskgmass_expellr[iz],diskgmass_expellmetr[iz],halogmass[iz],halogmassmet[iz],sfmassr[iz],format='(F,28E18)'
     print,     z_bins[iz],relostmass[iz],relostmassmet[iz],relostmassr[iz],relostmassmetr[iz],reejectmass[iz],reejectmassmet[iz],reejectmassr[iz],reejectmassmetr[iz],reexpellmass[iz],reexpellmassmet[iz],reexpellmassr[iz],reexpellmassmetr[iz],float(diskgmass_lost[iz]),diskgmass_lostmet[iz],float(diskgmass_lostr[iz]),diskgmass_lostmetr[iz],diskgmass[iz],diskgmassmet[iz],diskgmassr[iz],diskgmassmetr[iz],diskgmass_expell[iz],diskgmass_expellmet[iz],diskgmass_expellr[iz],diskgmass_expellmetr[iz],halogmass[iz],halogmassmet[iz],sfmassr[iz],format='(F,28E18)'
;stop
;allgmass_total_temp,halogmass_total_temp,diskgmass_total_temp,relost_total_temp,reeject_total,reexpell_total_temp,sfmass_total_temp,smass_accr_total_temp
;z_bins_temp,relostmass_temp,reejectmass_temp,reejectmassr_temp,reexpellmass_temp,reexpellmassr_temp,diskgmass_lost_temp,diskgmass_temp,halogmass_temp,sfmassr_temp
;     stop
  ENDFOR
  close,lun
  free_lun,lun
END
