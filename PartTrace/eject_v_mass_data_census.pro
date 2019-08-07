;Same as eject_v_mass_data but generates only that data used for the baryon census
;eject_v_mass_data_census,'/home/christenc/Data/Sims/cptmarvel.cosmo25cmb/cptmarvel.cosmo25cmb.4096g5HbwK1BH/','cptmarvel.cosmo25cmb.4096g5HbwK1BH','1'

;starlog.*.fits *reaccr_* *reaccrdisk_* *at_reaccr.* *at_reaccrdisk.* grp*.accrstars_?.fits
;
;

PRO eject_v_mass_data_census,dir,file,halo
  print,file + '.' + halo

  ageUniverse = wmap3_lookback(1000)

  units = tipsyunits(dir + file + '.param')
  sldata = mrdfits(dir + 'starlog.' + halo + '.fits',1,/silent)
  sldata.timeform = sldata.timeform*units.timeunit/1e9
  sfmass_total = total(sldata.massform)*units.massunit

  diskaccrm = mrdfits(dir + '/grp' + halo + '.mass_at_reaccrdisk.fits',0,/silent)
  diskaccrm[where(diskaccrm LE 0)] = units.istarmass
  diskaccrz = mrdfits(dir + '/grp' + halo + '.reaccrdisk_z.fits',0,/silent)
  diskaccri = mrdfits(dir + '/grp' + halo + '.reaccrdisk_iord.fits',0,/silent)
  diskgmass_total = total(diskaccrm[uniq(diskaccri,sort(diskaccri))])
  
  haloaccrm = mrdfits(dir + '/grp' + halo + '.mass_at_reaccr.fits',0,/silent)
  haloaccrz = mrdfits(dir + '/grp' + halo + '.reaccr_z.fits',0,/silent)
  haloaccri = mrdfits(dir + '/grp' + halo + '.reaccr_iord.fits',0,/silent)
  halogmass_total = total(haloaccrm[uniq(haloaccri,sort(haloaccri))])

  smass_accrm = mrdfits(dir + '/grp' + halo + '.accrstars_m.fits',/silent)
  smass_accrz = mrdfits(dir + '/grp' + halo + '.accrstars_z.fits',/silent)
  smass_accrm = smass_accrm[where(smass_accrz NE 99)]
  smass_accrz = smass_accrz[where(smass_accrz NE 99)]
;  smass_accrt = fltarr(nz) ;z_to_t(smass_accrz)
  smass_accr_total = total(smass_accrm)
  stop
  
  openw,lun,dir + '/grp' + halo + '.eject_quant_census.txt',/get_lun
  printf,lun,halogmass_total,diskgmass_total,smass_accr_total  
  close,lun
  free_lun,lun
END
