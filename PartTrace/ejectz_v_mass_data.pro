
PRO ejectz_v_mass_data,dir,file,halo,z_bins,zmax,mlbin
  print,file + '.' + halo
  nz = n_elements(z_bins)
  relostmass = fltarr(nz)
  reejectmass  = fltarr(nz)
  reejectmassr  = fltarr(nz)
  reexpellmass = fltarr(nz)
  reexpellmassr = fltarr(nz)
  diskgmass = fltarr(nz)
  diskgmass_h = fltarr(nz)
  diskgmass_lost = fltarr(nz)
  halogmass = fltarr(nz)
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
  
  accrm = mrdfits(dir + '/grp' + halo + '.mass_at_accr_rvir.fits',0,/silent) 
  allgmass_total = total(accrm)

  diskearlym = mrdfits(dir + '/grp' + halo + '.earlydisk_mass.fits',0,/silent)
  diskearlyi = mrdfits(dir + '/grp' + halo + '.earlydisk_iord.fits',0,/silent)
;  diskearlymet = mrdfits(dir + '/grp' + halo + '.earlydisk_zmetals.fits',0,/silent)
  diskearlymet = diskearlym * 0
  diskearlyz = fltarr(n_elements(diskearlyi)) + 99

  haloearlym = mrdfits(dir + '/grp' + halo + '.earlyhalo_mass.fits',0,/silent)
  haloearlyi = mrdfits(dir + '/grp' + halo + '.earlyhalo_iord.fits',0,/silent)
  haloearlymet = mrdfits(dir + '/grp' + halo + '.earlyhalo_zmetals.fits',0,/silent)
  haloearlyz = fltarr(n_elements(haloearlyi)) + 99
  
  diskaccrm_lost = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdiskall.fits',0,/silent)
  diskaccrm_lost[where(diskaccrm_lost LE 0)] = units.istarmass
  diskaccrz_lost = mrdfits(dir + '/grp' + halo + '.reaccrdiskall_z.fits',0,/silent)
  diskaccri_lost = mrdfits(dir + '/grp' + halo + '.reaccrdiskall_iord.fits',0,/silent)
  diskaccrhist_lost = mrdfits(dir + '/grp' + halo + '.reaccrdiskall_history.fits',1,/silent)
  diskaccrt_lost = z_to_t(diskaccrz_lost)
  diskregmassi_lost = [diskearlyi,diskaccri_lost[where(diskaccrz_lost LT zmax)]]
  diskregmassm_lost = [diskearlym,diskaccrm_lost[where(diskaccrz_lost LT zmax)]]
  diskregmassz_lost = [diskearlyz,diskaccrz_lost[where(diskaccrz_lost LT zmax)]]
  diskregmassmet_lost = [diskearlymet,diskaccrhist_lost[where(diskaccrz_lost LT zmax)].metallicity]

  ind_diskaccri_lost_uniq = uniq(diskaccri_lost,sort(diskaccri_lost))
;Unique accretion onto disk
;Select for those particles the first time they are accreted onto the disk
  diskaccri_lost_uniq = diskaccri_lost[ind_diskaccri_lost_uniq]
  diskaccrz_lost_uniq = fltarr(n_elements(diskaccri_lost_uniq))
  diskaccrm_lost_uniq = fltarr(n_elements(diskaccri_lost_uniq)) 
  diskaccrmet_lost_uniq = fltarr(n_elements(diskaccri_lost_uniq))
  FOR ind = 0, n_elements(diskaccri_lost_uniq) - 1 DO BEGIN
 ;Iterate through, saving the values that correspond to the highest redshift
     indiord = where(diskaccri_lost EQ diskaccri_lost_uniq[ind])
     temp = max(diskaccrz_lost[indiord],maxz_ind)
     diskaccrz_lost_uniq[ind] = diskaccrz_lost[indiord[maxz_ind]]
     diskaccrm_lost_uniq[ind] = diskaccrm_lost[indiord[maxz_ind]]
     diskaccrmet_lost_uniq[ind] = diskaccrhist_lost[indiord[maxz_ind]].metallicity
  ENDFOR
  diskgmass_total = total(diskaccrm_lost_uniq)

;--------------- Note that these include the initial accretion as well
;                as the accretion after ejection ------------
;  diskaccrm = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdisk.fits',0,/silent)
;  diskaccrm[where(diskaccrm LE 0)] = units.istarmass
;  diskaccrz = mrdfits(dir + '/grp' + halo + '.reaccrdisk_z.fits',0,/silent)
;  diskaccri = mrdfits(dir + '/grp' + halo + '.reaccrdisk_iord.fits',0,/silent)
;---------------------------------
  diskaccrhist = mrdfits(dir + '/grp' + halo + '.reaccrdisk_history.fits',1,/silent) ;Only accretion after ejection
  diskaccri = [diskaccri_lost_uniq,diskaccrhist.iord]
  diskaccrm = [diskaccrm_lost_uniq,diskaccrhist.mass]
  diskaccrz = [diskaccrz_lost_uniq,diskaccrhist.red]
  diskaccrmet = [diskaccrmet_lost_uniq,diskaccrhist.metallicity]
  diskaccrt = z_to_t(diskaccrz)
;  diskaccrm_clip = diskaccrm[where(diskaccrz LT zmax)]
;  diskaccri_clip = diskaccri[where(diskaccrz LT zmax)]
  diskregmassi = [diskearlyi,diskaccri[where(diskaccrz LT zmax)]] 
  diskregmassm = [diskearlym,diskaccrm[where(diskaccrz LT zmax)]] 
  diskregmassmet= [diskearlymet,diskaccrmet]
  diskregmassz = [diskearlyz,diskaccrz[where(diskaccrz LT zmax)]]
  
  haloaccrm = mrdfits(dir + '/grp' + halo + '.mass_at_reaccr.fits',0,/silent)
  haloaccrz = mrdfits(dir + '/grp' + halo + '.reaccr_z.fits',0,/silent)
  haloaccri = mrdfits(dir + '/grp' + halo + '.reaccr_iord.fits',0,/silent)
  haloaccrhist = mrdfits(dir + '/grp' + halo + '.reaccr_history.fits',1,/silent)
  halogmass_total = total(haloaccrm[uniq(haloaccri,sort(haloaccri))])
  haloaccrt = z_to_t(diskaccrz)
;  haloaccrm_clip = haloaccrm[where(haloaccrz LT zmax)]
;  haloaccri_clip = haloaccri[where(haloaccrz LT zmax)]
  haloaccri = [haloearlyi,haloaccri[where(haloaccrz LT zmax)]]
  haloaccrm = [haloearlym,haloaccrm[where(haloaccrz LT zmax)]]
  haloaccrmet = [haloearlymet,haloaccrhist[where(haloaccrz LT zmax)].metallicity]
  haloaccrz = haloaccrz[where(haloaccrz LT zmax)]

  relostm = mrdfits(dir + '/grp' + halo + '.mass_at_relost.fits',0,/silent) 
  relosti = mrdfits(dir + '/grp' + halo + '.relost_iord.fits',0,/silent)
  relostz = mrdfits(dir + '/grp' + halo + '.relost_z.fits',0,/silent)
  relosthist = mrdfits(dir + '/grp' + halo + '.relost_history.fits',1,/silent)
  relost_total = total(relostm[uniq(relosti,sort(relosti))])
  relostm = relostm[where(relostz LT zmax)]
  relosti = relosti[where(relostz LT zmax)]
  relostmet = relosthist[where(relostz LT zmax)].metallicity
  relostz = relostz[where(relostz LT zmax)]
  relostt = z_to_t(relostz)

  reejecti = mrdfits(dir + '/grp' + halo + '.reeject_iord.fits',0,/silent)
  reejectm = mrdfits(dir + '/grp' + halo + '.mass_at_reeject.fits',0,/silent)
  reejectz = mrdfits(dir + '/grp' + halo + '.reeject_z.fits',0,/silent)
  reejecthist = mrdfits(dir + '/grp' + halo + '.reeject_halo.fits',1,/silent)
  reeject_total = total(reejectm[uniq(reejecti,sort(reejecti))])
  reejectm = reejectm[where(reejectz LT zmax)]
  reejecti = reejecti[where(reejectz LT zmax)]
  reejectz = reejectz[where(reejectz LT zmax)]
  reejectmet = reejecthist[where(reejectz LT zmax)].metallicity
  reejectt = z_to_t(reejectz)
;    reejecti_uniq = reejecti[uniq(reejecti,sort(reejecti))]
;    reejectmcum  = weighted_histogram(reejectt, weight =  reejectm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)

  reexpelli = mrdfits(dir + '/grp' + halo + '.reexpell_iord.fits',0,/silent)
  reexpellm = mrdfits(dir + '/grp' + halo + '.mass_at_reexpell.fits',0,/silent) 
  reexpellz = mrdfits(dir + '/grp' + halo + '.reexpell_z.fits',0,/silent)
  reexpell_total = total(reexpellm[uniq(reexpelli,sort(reexpelli))])  
  eject_expell_iord = mrdfits(dir+'grp'+halo+'.eject_expell_iord.fits')
  reexpellhist = reejecthist[eject_expell_iord]
  reexpellm = reexpellm[where(reexpellz LT zmax)]
  reexpellz = reexpellz[where(reexpellz LT zmax)]
  reexpellt = z_to_t(reexpellz)
  reexpellmet = reexpellhist[where(reexpellz LT zmax)].metallicity

  smass_accrm = mrdfits(dir + '/grp' + halo + '.accrstars_m.fits',/silent)
  smass_accrz = mrdfits(dir + '/grp' + halo + '.accrstars_z.fits',/silent)
  smass_accrm = smass_accrm[where(smass_accrz NE 99)]
  smass_accrz = smass_accrz[where(smass_accrz NE 99)]
;  smass_accrt = fltarr(nz) ;z_to_t(smass_accrz)
  smass_accr_total = total(smass_accrm)

  openw,lun,dir + '/grp' + halo + '.ejectz_quant.txt',/get_lun
  printf,lun,allgmass_total,halogmass_total,diskgmass_total,relost_total,reeject_total,reexpell_total,float(sfmass_total),smass_accr_total,format='(8E14)'  
  FOR iz = 0, n_elements(z_bins) - 1 DO BEGIN
     IF z_bins[iz] LT 0.1 THEN mlbin2 = mlbin ELSE mlbin2 = mlbin/2
     relostmass[iz]       = total(relostmet[      where(relostz      GE z_bins[iz] AND relostz   LT zmax)]*relostm[      where(relostz      GE z_bins[iz] AND relostz   LT zmax)])
     reejectmass[iz]      = total(reejectm[     where(reejectz     GE z_bins[iz] AND reejectz  LT zmax)])
     reejectmassr[iz]     = total(reejectm[     where(reejectt     GE t_bins[iz] - mlbin2 AND reejectt      LE t_bins[iz] + mlbin2 AND reejectz LT zmax)])/mlbin
     reexpellmass[iz]     = total(reexpellm[    where(reexpellz    GE z_bins[iz] AND reexpellz LT zmax)])
     reexpellmassr[iz]    = total(reexpellm[    where(reexpellt    GE t_bins[iz] - mlbin2 AND reexpellt     LE t_bins[iz] + mlbin2 AND reexpellz LT zmax)])/mlbin
     diskgmass_lost[iz]        = total(diskregmassm_lost[    where(diskregmassz_lost     GE z_bins[iz])])
     diskgmass[iz]        = total(diskregmassm[    where(diskregmassz     GE z_bins[iz])])
     halogmass[iz]        = total(haloaccrm[    where(haloaccrz    GE z_bins[iz])])
     sfmassr[iz] = total(sldata[where(sldata.timeform GE t_bins[iz] - mlbin2 AND sldata.timeform LE t_bins[iz] + mlbin2 AND sldata.timeform GT tmin)].massform)*units.massunit/mlbin
     printf,lun,z_bins[iz],relostmass[iz],reejectmass[iz],reejectmassr[iz],reexpellmass[iz],reexpellmassr[iz],diskgmass_lost[iz],diskgmass[iz],halogmass[iz],sfmassr[iz],format='(F,10E14)'
  ENDFOR
  close,lun
  free_lun,lun
END
