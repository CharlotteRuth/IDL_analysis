
PRO eject_v_mass_data,dir,file,halo,z_bins,zmax,mlbin
  print,file + '.' + halo
  nz = n_elements(z_bins)
  relostmass = fltarr(nz)
  reheatmass = fltarr(nz)
  reheatmassr = fltarr(nz)
  reejectmass_rv2  = fltarr(nz)
  reejectmassr_rv2  = fltarr(nz)
  reejectmass_rv5  = fltarr(nz)
  reejectmassr_rv5  = fltarr(nz)
  reejectmass  = fltarr(nz)
  reejectmassr  = fltarr(nz)
  reejectmass_cool  = fltarr(nz)
  reejectmassr_cool  = fltarr(nz)
  reexpellmass = fltarr(nz)
  reexpellmassr = fltarr(nz)
  reexpellmass_cool = fltarr(nz)
  reexpellmassr_cool = fltarr(nz)
  diskgmass = fltarr(nz)
  diskgmass_h = fltarr(nz)
  diskgmass_rv2 = fltarr(nz)
  diskgmass_rv5 = fltarr(nz)
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
  diskearlyz = fltarr(n_elements(diskearlyi)) + 99

  haloearlym = mrdfits(dir + '/grp' + halo + '.earlyhalo_mass.fits',0,/silent)
  haloearlyi = mrdfits(dir + '/grp' + halo + '.earlyhalo_iord.fits',0,/silent)
  haloearlyz = fltarr(n_elements(haloearlyi)) + 99

  diskaccrm_lost = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdiskall.fits',0,/silent)
  diskaccrm_lost[where(diskaccrm_lost LE 0)] = units.istarmass
  diskaccrz_lost = mrdfits(dir + '/grp' + halo + '.reaccrdiskall_z.fits',0,/silent)
  diskaccri_lost = mrdfits(dir + '/grp' + halo + '.reaccrdiskall_iord.fits',0,/silent)
  diskaccrt_lost = z_to_t(diskaccrz_lost)
  diskregmassi_lost = [diskearlyi,diskaccri_lost[where(diskaccrz_lost LT zmax)]]
  diskregmassm_lost = [diskearlym,diskaccrm_lost[where(diskaccrz_lost LT zmax)]]
  diskregmassz_lost = [diskearlyz,diskaccrz_lost[where(diskaccrz_lost LT zmax)]]

  diskaccrm = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdisk.fits',0,/silent)
  diskaccrm[where(diskaccrm LE 0)] = units.istarmass
  diskaccrz = mrdfits(dir + '/grp' + halo + '.reaccrdisk_z.fits',0,/silent)
  diskaccri = mrdfits(dir + '/grp' + halo + '.reaccrdisk_iord.fits',0,/silent)
  diskgmass_total = total(diskaccrm[uniq(diskaccri,sort(diskaccri))])
  diskaccrt = z_to_t(diskaccrz)
;  diskaccrm_clip = diskaccrm[where(diskaccrz LT zmax)]
;  diskaccri_clip = diskaccri[where(diskaccrz LT zmax)]
  diskregmassi = [diskearlyi,diskaccri[where(diskaccrz LT zmax)]] 
  diskregmassm = [diskearlym,diskaccrm[where(diskaccrz LT zmax)]] 
  diskregmassz = [diskearlyz,diskaccrz[where(diskaccrz LT zmax)]]

  diskaccrm_h = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdiskheat.fits',0,/silent)
  diskaccrz_h = mrdfits(dir + '/grp' + halo + '.reaccrdiskheat_z.fits',0,/silent)
  diskaccri_h = mrdfits(dir + '/grp' + halo + '.reaccrdiskheat_iord.fits',0,/silent)
  diskregmassi_h = [diskearlyi,diskaccri_h[where(diskaccrz_h LT zmax)]]
  diskregmassm_h = [diskearlym,diskaccrm_h[where(diskaccrz_h LT zmax)]]
  diskregmassz_h = [diskearlyz,diskaccrz_h[where(diskaccrz_h LT zmax)]]

  diskaccrm_rv2 = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdisk_rvir0.2.fits',0,/silent)
  diskaccrz_rv2 = mrdfits(dir + '/grp' + halo + '.reaccrdisk_rvir0.2_z.fits',0,/silent)
  diskaccri_rv2 = mrdfits(dir + '/grp' + halo + '.reaccrdisk_rvir0.2_iord.fits',0,/silent)
  diskregmassi_rv2 = [diskearlyi,diskaccri_rv2[where(diskaccrz_rv2 LT zmax)]]
  diskregmassm_rv2 = [diskearlym,diskaccrm_rv2[where(diskaccrz_rv2 LT zmax)]]
  diskregmassz_rv2 = [diskearlyz,diskaccrz_rv2[where(diskaccrz_rv2 LT zmax)]]

  diskaccrm_rv5 = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdisk_rvir0.5.fits',0,/silent)
  diskaccrz_rv5 = mrdfits(dir + '/grp' + halo + '.reaccrdisk_rvir0.5_z.fits',0,/silent)
  diskaccri_rv5 = mrdfits(dir + '/grp' + halo + '.reaccrdisk_rvir0.5_iord.fits',0,/silent)
  diskregmassi_rv5 = [diskearlyi,diskaccri_rv5[where(diskaccrz_rv5 LT zmax)]]
  diskregmassm_rv5 = [diskearlym,diskaccrm_rv5[where(diskaccrz_rv5 LT zmax)]]
  diskregmassz_rv5 = [diskearlyz,diskaccrz_rv5[where(diskaccrz_rv5 LT zmax)]]
  
  haloaccrm = mrdfits(dir + '/grp' + halo + '.mass_at_reaccr.fits',0,/silent)
  haloaccrz = mrdfits(dir + '/grp' + halo + '.reaccr_z.fits',0,/silent)
  haloaccri = mrdfits(dir + '/grp' + halo + '.reaccr_iord.fits',0,/silent)
  halogmass_total = total(haloaccrm[uniq(haloaccri,sort(haloaccri))])
  haloaccrt = z_to_t(diskaccrz)
;  haloaccrm_clip = haloaccrm[where(haloaccrz LT zmax)]
;  haloaccri_clip = haloaccri[where(haloaccrz LT zmax)]
  haloregmass1i = [haloearlyi,haloaccri[where(haloaccrz LT zmax)]]
  haloregmass1m = [haloearlym,haloaccrm[where(haloaccrz LT zmax)]]

  relostm = mrdfits(dir + '/grp' + halo + '.mass_at_relost.fits',0,/silent) 
  relosti = mrdfits(dir + '/grp' + halo + '.relost_iord.fits',0,/silent)
  relostz = mrdfits(dir + '/grp' + halo + '.relost_z.fits',0,/silent)
  relostm = relostm[where(relostz LT zmax)]
  relosti = relosti[where(relostz LT zmax)]
  relostz = relostz[where(relostz LT zmax)]
  relostt = z_to_t(relostz)

  reheatm = mrdfits(dir + '/grp' + halo + '.mass_at_reheat.fits',0,/silent)
  reheatz = mrdfits(dir + '/grp' + halo + '.reheat_z.fits',0,/silent)
  reheati = mrdfits(dir + '/grp' + halo + '.reheat_iord.fits',0,/silent)
  reheat_total = total(reheatm[uniq(reheati,sort(reheati))])
  reheatm = reheatm[where(reheatz LT zmax)]
  reheati = reheati[where(reheatz LT zmax)]
  reheatz = reheatz[where(reheatz LT zmax)]
  reheatt = z_to_t(reheatz)

  reeject_rv2m = mrdfits(dir + '/grp' + halo + '.mass_at_reeject_rvir0.2.fits',0,/silent)
  reeject_rv2z = mrdfits(dir + '/grp' + halo + '.reeject_rvir0.2_z.fits',0,/silent)
  reeject_rv2i = mrdfits(dir + '/grp' + halo + '.reeject_rvir0.2_iord.fits',0,/silent)
  reeject_rv2_total = total(reeject_rv2m[uniq(reeject_rv2i,sort(reeject_rv2i))])
  reeject_rv2m = reeject_rv2m[where(reeject_rv2z LT zmax)]
  reeject_rv2i = reeject_rv2i[where(reeject_rv2z LT zmax)]
  reeject_rv2z = reeject_rv2z[where(reeject_rv2z LT zmax)]
  reeject_rv2t = z_to_t(reeject_rv2z)

  reeject_rv5m = mrdfits(dir + '/grp' + halo + '.mass_at_reeject_rvir0.5.fits',0,/silent)
  reeject_rv5z = mrdfits(dir + '/grp' + halo + '.reeject_rvir0.5_z.fits',0,/silent)
  reeject_rv5i = mrdfits(dir + '/grp' + halo + '.reeject_rvir0.5_iord.fits',0,/silent)
  reeject_rv5_total = total(reeject_rv5m[uniq(reeject_rv5i,sort(reeject_rv5i))])
  reeject_rv5m = reeject_rv5m[where(reeject_rv5z LT zmax)]
  reeject_rv5i = reeject_rv5i[where(reeject_rv5z LT zmax)]
  reeject_rv5z = reeject_rv5z[where(reeject_rv5z LT zmax)]
  reeject_rv5t = z_to_t(reeject_rv5z)

  reejecti = mrdfits(dir + '/grp' + halo + '.reeject_iord.fits',0,/silent)
  reejectm = mrdfits(dir + '/grp' + halo + '.mass_at_reeject.fits',0,/silent)
  reejectz = mrdfits(dir + '/grp' + halo + '.reeject_z.fits',0,/silent)
  reeject_total = total(reejectm[uniq(reejecti,sort(reejecti))])
  reejectm = reejectm[where(reejectz LT zmax)]
  reejecti = reejecti[where(reejectz LT zmax)]
  reejectz = reejectz[where(reejectz LT zmax)]
  reejectt = z_to_t(reejectz)
;    reejecti_uniq = reejecti[uniq(reejecti,sort(reejecti))]
;    reejectmcum  = weighted_histogram(reejectt, weight =  reejectm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)

  reeject_cooli = mrdfits(dir + '/grp' + halo + '.coolon_reeject_iord.fits',0,/silent)
  reeject_coolm = mrdfits(dir + '/grp' + halo + '.coolon_mass_at_reeject.fits',0,/silent)
  reeject_coolz = mrdfits(dir + '/grp' + halo + '.coolon_reeject_z.fits',0,/silent)
  reeject_cooltotal = total(reeject_coolm[uniq(reeject_cooli,sort(reeject_cooli))])
  reeject_coolm = reeject_coolm[where(reeject_coolz LT zmax)]
  reeject_cooli = reeject_cooli[where(reeject_coolz LT zmax)]
  reeject_coolz = reeject_coolz[where(reeject_coolz LT zmax)]
  reeject_coolt = z_to_t(reeject_coolz)

;    reexpelli = mrdfits(dir + '/grp' + halo + '.reexpell_iord.fits',0,/silent)
  reexpellm = mrdfits(dir + '/grp' + halo + '.mass_at_reexpell.fits',0,/silent) 
  reexpellz = mrdfits(dir + '/grp' + halo + '.reexpell_z.fits',0,/silent)
  reexpellt = z_to_t(reexpellz)
  reexpellm = reexpellm[where(reexpellz LT zmax)]
  reexpellz = reexpellz[where(reexpellz LT zmax)]
;    reexpelli_uniq = reexpelli[uniq(reexpelli,sort(reexpelli))]
;    reexpellmcum  = weighted_histogram(reexpellt, weight =  reexpellm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)
;    reexpellm1 = fltarr(n_elements(reexpellm))
;    reexpellz1 = fltarr(n_elements(reexpellz))
;    exi = 0L
;    FOR eji = 0, n_elements(reejecti) -1 DO BEGIN
;       IF exi GE n_elements(reexpelli) THEN BREAK
;       IF reejecti[eji] EQ reexpelli[exi] THEN BEGIN
;          reexpellm1[exi] = reejectm[eji]
;          reexpellz1[exi] = reejectz[eji]
;          exi = exi + 1
;       ENDIF
;    ENDFOR
;    reexpellt1 = z_to_t(uniqz,uniqt,reexpellz1)

  reexpell_coolm = mrdfits(dir + '/grp' + halo + '.coolon_mass_at_reexpell.fits',0,/silent) 
  reexpell_coolz = mrdfits(dir + '/grp' + halo + '.coolon_reexpell_z.fits',0,/silent)
  reexpell_coolt = z_to_t(reexpell_coolz)
  reexpell_coolm = reexpell_coolm[where(reexpell_coolz LT zmax)]
  reexpell_coolz = reexpell_coolz[where(reexpell_coolz LT zmax)]

  smass_accrm = mrdfits(dir + '/grp' + halo + '.accrstars_m.fits',/silent)
  smass_accrz = mrdfits(dir + '/grp' + halo + '.accrstars_z.fits',/silent)
  smass_accrm = smass_accrm[where(smass_accrz NE 99)]
  smass_accrz = smass_accrz[where(smass_accrz NE 99)]
;  smass_accrt = fltarr(nz) ;z_to_t(smass_accrz)
  smass_accr_total = total(smass_accrm)

  openw,lun,dir + '/grp' + halo + '.eject_quant.txt',/get_lun
  printf,lun,allgmass_total,halogmass_total,diskgmass_total,reheat_total,reeject_rv2_total,reeject_rv2_total,reeject_total,sfmass_total,smass_accr_total  
  FOR iz = 0, n_elements(z_bins) - 1 DO BEGIN
     IF z_bins[iz] LT 0.1 THEN mlbin2 = mlbin ELSE mlbin2 = mlbin/2
     relostmass[iz]     = total(relostm[    where(relostz    GE z_bins[iz] AND relostz LT zmax)])
     diskgmass_lost[iz]        = total(diskregmassm_lost[    where(diskregmassz_lost     GE z_bins[iz])])
     reheatmass[iz]       = total(reheatm[      where(reheatz      GE z_bins[iz] AND reheatz LT zmax)])
     reheatmassr[iz]      = total(reheatm[      where(reheatt      GE t_bins[iz] - mlbin2 AND reheatt       LE t_bins[iz] + mlbin2 AND reheatz LT zmax)])/mlbin
     reejectmass_rv2[iz]  = total(reeject_rv2m[ where(reeject_rv2z GE z_bins[iz] AND reeject_rv2z LT zmax)])
     reejectmassr_rv2[iz] = total(reeject_rv2m[ where(reeject_rv2t GE t_bins[iz] - mlbin2 AND reeject_rv2t  LE t_bins[iz] + mlbin2 AND reeject_rv2z LT zmax)])/mlbin
     reejectmass_rv5[iz]  = total(reeject_rv5m[ where(reeject_rv5z GE z_bins[iz] AND reeject_rv5z LT zmax)])
     reejectmassr_rv5[iz] = total(reeject_rv5m[ where(reeject_rv5t GE t_bins[iz] - mlbin2 AND reeject_rv5t  LE t_bins[iz] + mlbin2 AND reeject_rv5z LT zmax)])/mlbin
     reejectmass[iz]      = total(reejectm[     where(reejectz     GE z_bins[iz] AND reejectz LT zmax)])
     reejectmassr[iz]     = total(reejectm[     where(reejectt     GE t_bins[iz] - mlbin2 AND reejectt      LE t_bins[iz] + mlbin2 AND reejectz LT zmax)])/mlbin
     reejectmass_cool[iz] = total(reeject_coolm[where(reeject_coolz GE z_bins[iz] AND reeject_coolz LT zmax)])
     reejectmassr_cool[iz]= total(reeject_coolm[where(reeject_coolt GE t_bins[iz] - mlbin2 AND reeject_coolt LE t_bins[iz] + mlbin2 AND reeject_coolz LT zmax)])/mlbin
     reexpellmass[iz]     = total(reexpellm[    where(reexpellz    GE z_bins[iz] AND reexpellz LT zmax)])
     reexpellmassr[iz]    = total(reexpellm[    where(reexpellt    GE t_bins[iz] - mlbin2 AND reexpellt     LE t_bins[iz] + mlbin2 AND reexpellz LT zmax)])/mlbin
     reexpellmass_cool[iz]= total(reexpell_coolm[where(reexpell_coolz GE z_bins[iz] AND reexpell_coolz LT zmax)])
     reexpellmassr_cool[iz]= total(reexpell_coolm[where(reexpell_coolt GE t_bins[iz] - mlbin2 AND reexpell_coolt     LE t_bins[iz] + mlbin2 AND reexpell_coolz LT zmax)])/mlbin
     diskgmass[iz]        = total(diskregmassm[    where(diskregmassz     GE z_bins[iz])])
     diskgmass_h[iz]      = total(diskregmassm_h[  where(diskregmassz_h   GE z_bins[iz])])
     diskgmass_rv2[iz]    = total(diskregmassm_rv2[where(diskregmassz_rv2 GE z_bins[iz])])
     diskgmass_rv5[iz]    = total(diskregmassm_rv5[where(diskregmassz_rv5 GE z_bins[iz])])
     halogmass[iz]        = total(haloaccrm[    where(haloaccrz    GE z_bins[iz])])
     sfmassr[iz] = total(sldata[where(sldata.timeform GE t_bins[iz] - mlbin2 AND sldata.timeform LE t_bins[iz] + mlbin2 AND sldata.timeform GT tmin)].massform)*units.massunit/mlbin
     printf,lun,z_bins[iz],relostmass[iz],reheatmass[iz],reheatmassr[iz],reejectmass_rv2[iz],reejectmassr_rv2[iz],reejectmass_rv5[iz],reejectmassr_rv5[iz],reejectmass[iz],reejectmassr[iz],reejectmass_cool[iz],reejectmassr_cool[iz],reexpellmass[iz],reexpellmassr[iz],reexpellmass_cool[iz],reexpellmassr_cool[iz],diskgmass_lost[iz],diskgmass[iz],diskgmass_h[iz],diskgmass_rv2[iz],diskgmass_rv5[iz],halogmass[iz],sfmassr[iz],format='(F,22E14)'
  ENDFOR
  close,lun
  free_lun,lun
END
