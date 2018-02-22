PRO ejectz_v_mass,dirs,files,halo = halo,finalstep = finalstep,colors = colors,outplot = outplot,stellarmass = stellarmass,z_cut = z_cut,symbols = symbols,formatthick = formatthick,z_colors = z_colors,gmass = gmass,accr_star = accr_star,rewrite = rewrite

formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
    IF keyword_set(stellarmass) THEN outplot = outplot + '_sm' ELSE outplot = outplot + '_vm'
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 600
    ysize = 400
ENDELSE
device,decomposed = 0

n = n_elements(dirs)

zmax = 3.130 ;100;3.130 ;Instituted because I am frequently missing steps below 80.
tmin = (z_to_t(zmax))[0]
mlbin = 1.0 ;in billion
f_bar = 0.16510
IF NOT keyword_set(halo) THEN halo = strarr(n) + '1'
IF NOT keyword_set(finalstep) THEN finalstep = '00512'
IF keyword_set(z_cut) THEN $
 IF n_elements(z_cut) EQ 1 THEN z_bins = [5,4,3,2,1,0.5,-1e-10] ELSE z_bins = z_cut $
 ELSE z_bins = [0]
z_bins_legend = strtrim(fix(z_bins*10.0)/10.0,2)
nz = n_elements(z_bins)
vmass = fltarr(n)
smass = fltarr(n)
vcirc = fltarr(n)
gmass_hot = fltarr(n)
vmasst = fltarr(n,nz)
smasst = fltarr(n,nz)
vcirct = fltarr(n,nz)

diskgmass_uniq_total  = fltarr(n)
diskgmassmet_uniq_total  = fltarr(n)
halogmass_total  = fltarr(n)
allgmass_total   = fltarr(n)
sfmass_total     = fltarr(n)
smass_accr_total = fltarr(n)
uniqaccrgmass_total = fltarr(n)
uniqaccrzmass_total = fltarr(n)  
zproduced_total = fltarr(n)
feproduced_total = fltarr(n)
oxproduced_total = fltarr(n)
zavail_zmax_total = fltarr(n)

relostmass = fltarr(n,nz)
relostmassmet = fltarr(n,nz)
relostmassr = fltarr(n,nz)
relostmassmetr = fltarr(n,nz)
;reheatmass = fltarr(n,nz)
;reheatmassr = fltarr(n,nz)
;reejectmass_rv2  = fltarr(n,nz)
;reejectmassr_rv2  = fltarr(n,nz)
;reejectmass_rv5  = fltarr(n,nz)
;reejectmassr_rv5  = fltarr(n,nz)
reejectmass  = fltarr(n,nz)
reejectmassmet  = fltarr(n,nz)
reejectmassr  = fltarr(n,nz)
reejectmassmetr  = fltarr(n,nz)
;reejectmass_cool  = fltarr(n,nz)
;reejectmass_coolr  = fltarr(n,nz)
reexpellmass = fltarr(n,nz)
reexpellmassmet = fltarr(n,nz)
reexpellmassr = fltarr(n,nz)
reexpellmassmetr = fltarr(n,nz)
;reexpellmass_cool = fltarr(n,nz)
;reexpellmass_coolr = fltarr(n,nz)
diskgmass = fltarr(n,nz)
diskgmassmet = fltarr(n,nz)
diskgmassr = fltarr(n,nz)
diskgmassmetr = fltarr(n,nz)
;diskgmass_h = fltarr(n,nz)
;diskgmass_rv2 = fltarr(n,nz)
;diskgmass_rv5 = fltarr(n,nz)
diskgmass_lost = fltarr(n,nz)
diskgmass_lostmet = fltarr(n,nz)
diskgmass_lostr = fltarr(n,nz)
diskgmass_lostmetr = fltarr(n,nz)
diskgmass_expell = fltarr(n,nz)
diskgmass_expellmet = fltarr(n,nz)
diskgmass_expellr = fltarr(n,nz)
diskgmass_expellmetr = fltarr(n,nz)
halogmass = fltarr(n,nz)
halogmassmet = fltarr(n,nz)
sfmassr = fltarr(n,nz)
reoutflowmass = fltarr(n,nz)
reoutflowmassr = fltarr(n,nz)
reoutflowmass_uniq = fltarr(n,nz)
metallicity_ISMr = fltarr(n,nz)

FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   print,files[i] + '.' + finalstep + '.halo.' + halo[i]
   stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group EQ halo[i]))[0]
   vmass[i] = stat[ind].m_tot
   smass[i] = stat[ind].m_star
   gmass_hot[i] = stat[ind].m_gas
;   vcirc[i] = stat[ind].vc
   vcirc[i] = sqrt(6.67e-11*stat[ind].m_tot*1.9891e30/(3.08567758e19*stat[ind].rvir))/1000.
   print,vcirc[i],stat[ind].vc
   align = mrdfits(dirs[i] + '/grp' + halo[i] + '.alignment.fits',1,/silent)
   readcol,dirs[i] + '/grp' + halo[i] + '.mass_metal_track.dat',halo_t,time_t,z_t,mtot_t,mgas_t,mstar_t,mdark_t,metal_t,ox_t,fe_t,HI_t,H2_t,coldg_t,/silent
   ageUniverse = wmap3_lookback(1000)
   t_bins = (ageUniverse - wmap3_lookback(z_bins))/1e9
   readcol,dirs[i] + '/grp' + halo[i] + '.metals.txt',zmetals,ox,fe,coldg,zmetals_H,ox_H,fe_H,mHI_m,mH2_m,Hgas,/silent
   FOR iz = 0, n_elements(z_bins) - 1 DO BEGIN
      temp = min(abs(time_t - t_bins[iz]),ind_t)
      vmasst[i,iz] = mtot_t[ind_t]
      smasst[i,iz] = mstar_t[ind_t] 
      temp = min(abs(align.z - z_bins[iz]),ind_t)
      stat = read_stat_struc_amiga(dirs[i] + '/' + align[ind_t].file + '.amiga.stat')
      ind = (where(stat.group EQ align[ind_t].haloid))[0]
;      vcirct[i,iz] = stat[ind].vc
      vcirct[i,iz] = sqrt(6.67e-11*stat[ind].m_tot*1.9891e30/(3.08567758e19*stat[ind].rvir))/1000.
       metallicity_ISMr[i,iz] = zmetals[ind_t]/coldg[ind_t]
   ENDFOR
   readcol,dirs[i] + '/grp' + halo[i] + '.sfmetals.txt',zstars,oxstars,festars
   ind_zrange = where(z_t LE zmax) ;Only look at metals produced by stars since the max redshift
   zproduced_total[i] = total(zstars[ind_zrange])
   oxproduced_total[i] = total(oxstars[ind_zrange])
   feproduced_total[i] = total(festars[ind_zrange])
   zavail_zmax_total[i] = zproduced_total[i]+zmetals[ind_zrange[0]]
   print,'Metals since z=zmax, produced: ',zproduced_total[i],', available: ',zproduced_total[i]+zmetals[ind_zrange[0]]
   IF keyword_set(rewrite) OR NOT file_test(dirs[i] + '/grp' + halo[i] + '.ejectz_quant.txt') THEN $
      ejectz_v_mass_data,dirs[i],files[i],halo[i],z_bins,zmax,mlbin
   openr,lun,dirs[i] + '/grp' + halo[i] + '.ejectz_quant.txt',/get_lun
   readf,lun,allgmass_total_temp,halogmass_total_temp,diskgmass_uniq_total_temp,diskgmassmet_uniq_total_temp,relost_total_temp,relostmet_total_temp,reeject_total,reejectmet_total_temp,reexpell_total_temp,reexpellmet_total_temp,sfmass_total_temp,sfmassmet_total_temp,smass_accr_total_temp
   close,lun
   free_lun,lun
   allgmass_total[i] = allgmass_total_temp
   halogmass_total[i] = halogmass_total_temp
   diskgmass_uniq_total[i] = diskgmass_uniq_total_temp
   diskgmassmet_uniq_total[i] = diskgmassmet_uniq_total_temp
   sfmass_total[i] = sfmass_total_temp
   smass_accr_total[i] = smass_accr_total_temp
   readcol,dirs[i] + '/grp' + halo[i] + '.ejectz_quant.txt',z_bins_temp,relostmass_temp,relostmassmet_temp,relostmassr_temp,relostmassmetr_temp,reejectmass_temp,reejectmassmet_temp,reejectmassr_temp,reejectmassmetr_temp,reexpellmass_temp,reexpellmassmet_temp,reexpellmassr_temp,reexpellmassmetr_temp,diskgmass_lost_temp,diskgmass_lostmet_temp,diskgmass_lostr_temp,diskgmass_lostmetr_temp,diskgmass_temp,diskgmassmet_temp,diskgmassr_temp,diskgmassmetr_temp,diskgmass_expell_temp,diskgmass_expellmet_temp,diskgmass_expellr_temp,diskgmass_expellmetr_temp,halogmass_temp,halogmassmet_temp,sfmassr_temp
;   stop

   relostmass[i,*] = relostmass_temp
   relostmassmet[i,*] = relostmassmet_temp
   relostmassr[i,*] = relostmassr_temp
   relostmassmetr[i,*] = relostmassmetr_temp
;   reheatmass[i,*] = reheatmass_temp
;   reheatmassr[i,*] = reheatmassr_temp
;   reejectmass_rv2[i,*] = reejectmass_rv2_temp
;   reejectmassr_rv2[i,*] = reejectmassr_rv2_temp
;   reejectmass_rv5[i,*] = reejectmass_rv5_temp
;   reejectmassr_rv5[i,*] = reejectmassr_rv5_temp
   reejectmass[i,*] = reejectmass_temp
   reejectmassmet[i,*] = reejectmassmet_temp
   reejectmassr[i,*] = reejectmassr_temp
   reejectmassmetr[i,*] = reejectmassmetr_temp
;   reejectmass_cool[i,*] = reejectmass_cool_temp
;   reejectmass_coolr[i,*] = reejectmass_coolr_temp
   reexpellmass[i,*] = reexpellmass_temp
   reexpellmassmet[i,*] = reexpellmassmet_temp
   reexpellmassr[i,*] = reexpellmassr_temp
   reexpellmassmetr[i,*] = reexpellmassmetr_temp
;   reexpellmass_cool[i,*] = reexpellmass_cool_temp
;   reexpellmass_coolr[i,*] = reexpellmass_coolr_temp
   diskgmass_lost[i,*] = diskgmass_lost_temp ;Gas accreted for the 1st time or after being lost
   diskgmass_lostmet[i,*] = diskgmass_lostmet_temp
   diskgmass_lostr[i,*] = diskgmass_lostr_temp ;Gas accreted for the 1st time or after being lost
   diskgmass_lostmetr[i,*] = diskgmass_lostmetr_temp
   diskgmass[i,*] = diskgmass_temp ;Gas accreted after being ejected
   diskgmassmet[i,*] = diskgmassmet_temp
   diskgmassr[i,*] = diskgmassr_temp ;Gas accreted after being ejected
   diskgmassmetr[i,*] = diskgmassmetr_temp
   diskgmass_expell[i,*] = diskgmass_expell_temp ;Gas accreted after being ejected
   diskgmass_expellmet[i,*] = diskgmass_expellmet_temp
   diskgmass_expellr[i,*] = diskgmass_expellr_temp ;Gas accreted after being ejected
   diskgmass_expellmetr[i,*] = diskgmass_expellmetr_temp
;   diskgmass_h[i,*] = diskgmass_h_temp
;   diskgmass_rv2[i,*] = diskgmass_rv2_temp
;   diskgmass_rv5[i,*] = diskgmass_rv5_temp
   halogmass[i,*] = halogmass_temp
   halogmassmet[i,*] = halogmassmet_temp
   sfmassr[i,*] = sfmassr_temp

;For Rachel
   IF 0 THEN BEGIN
      units = tipsyunits(dirs[i] + '/' + files[i] + '.param')
      reoutflowm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reoutflow.fits',0,/silent) 
      reoutflowi = mrdfits(dirs[i] + '/grp' + halo[i] + '.reoutflow_iord.fits',0,/silent)
      reoutflowz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reoutflow_z.fits',0,/silent)
      reoutflowt = z_to_t(reoutflowz)
      FOR iz = 0, n_elements(z_bins) - 1 DO $
         reoutflowmass[i,iz] = total(reoutflowm[where(reoutflowz GE z_bins[iz] AND reoutflowz LT zmax)])
      reoutflowi_uniq_ind = uniq(reoutflowi,sort(reoutflowi))
      FOR iz = 0, n_elements(z_bins) - 1 DO BEGIN
         IF z_bins[iz] LT 0.1 THEN mlbin2 = mlbin ELSE mlbin2 = mlbin/2
         reoutflowmassr[i,iz]= total(reoutflowm[where(reoutflowt GE t_bins[iz] - mlbin2 AND reoutflowt LE t_bins[iz] + mlbin2 AND reoutflowz LT zmax)])/mlbin
      ENDFOR
      reoutflowm_uniq = reoutflowm[ reoutflowi_uniq_ind]
      reoutflowz_uniq = reoutflowz[ reoutflowi_uniq_ind]
      reoutflowt_uniq = reoutflowt[ reoutflowi_uniq_ind]
      FOR iz = 0, n_elements(z_bins) - 1 DO $
         reoutflowmass_uniq[i,iz] = total(reoutflowm_uniq[where(reoutflowz_uniq GE z_bins[iz] AND reoutflowz_uniq LT zmax)])
   ENDIF

;plot,xmass,((halogmass_total + smass_accr_total) - (smass + gmass_hot))/reexpellmass[*,3],psym = 4,/xlog
;plot,xmass,(halogmass_total + smass_accr_total) - (smass + gmass_hot),psym = 4,/xlog,/ylog
;oplot,xmass,reexpellmass[*,3],psym = 5

   IF 0 THEN BEGIN
      units = tipsyunits(dirs[i] + '/' + files[i] + '.param')
      relostm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_relost.fits',0,/silent) 
      relosti = mrdfits(dirs[i] + '/grp' + halo[i] + '.relost_iord.fits',0,/silent)
      relostz = mrdfits(dirs[i] + '/grp' + halo[i] + '.relost_z.fits',0,/silent)
      relostt = z_to_t(relostz)
      
      diskaccrm_lost = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdiskall.fits',0,/silent)
      diskaccrm_lost[where(diskaccrm_lost LE 0)] = units.istarmass
      diskaccrz_lost = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdiskall_z.fits',0,/silent)
      diskaccri_lost = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdiskall_iord.fits',0,/silent)
      diskaccrt_lost = z_to_t(diskaccrz_lost)
      rtipsy,dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' + finalstep + '.halo.' + halo[i] + '.std',h,g,d,s
         
      IF keyword_set(Hcut) THEN BEGIN 
         read_tipsy_arr,dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' + finalstep + '.halo.' + halo[i] + '.HI',h,HI,part = 'gas'
         read_tipsy_arr,dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' + finalstep + '.halo.' + halo[i] + '.H2',h,H2,part = 'gas'
         Hall = (HI + 2*H2)/max(HI + 2*H2)
         diskind = where(Hall GE 0.5 AND abs(g.z*units.lengthunit) LE 3)
      ENDIF ELSE BEGIN
         diskind = where(g.dens*units.rhounit GE 0.1 AND g.tempg LE 1.2e4 AND abs(g.z*units.lengthunit) LE 3)
      ENDELSE
;      gmass[i] = total(g[diskind].mass*units.massunit)
;   print,'Cold Gas Mass',coldg_t[n_elements(coldg_t) - 1],gmass[i],total(g[diskind].mass*units.massunit)
      
      read_tipsy_arr,dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' + finalstep + '.halo.' + halo[i] + '.iord',h,iord,part = 'gas',type = 'long'
      read_tipsy_arr,dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' + finalstep + '.halo.' + halo[i] + '.iord',h,iordstar,part = 'star',type = 'long'
      spawn,'ls ' + dirs[i] + '/' + files[i] + '.starlog.fits',files_sl
      IF file_test(files_sl[0]) THEN BEGIN
         sl = mrdfits(files_sl[0],1,/silent)
         print,files_sl[0]
      ENDIF ELSE BEGIN
         spawn,'ls ' + dirs[i] + '/*_2merge.starlog.fits',files_sl
         IF file_test(files_sl[0]) THEN BEGIN
            sl = mrdfits(files_sl[0],1,/silent)
            print,files_sl[0]
         ENDIF ELSE BEGIN
            spawn,'ls ' + dirs[i] + '/*.starlog',files_sl
            sl = rstarlog(files_sl[0],/molecularH)
            print,files_sl[0]
         ENDELSE
      ENDELSE
      match,sl.iorderstar,iordstar,ind4,ind5 ;stars in the halo
      sl = sl[ind4]
      
      relosti_uniq = relosti[uniq(relosti,sort(relosti))]
      relostt_max_uniq = relostt[uniq(relosti,sort(relosti))]
      relostm_max_uniq = relostm[uniq(relosti,sort(relosti))]
      
      diskaccri_uniq_lost = diskaccri_lost[uniq(diskaccri_lost,sort(diskaccri_lost))]
      diskaccrt_max_uniq_lost = diskaccrt_lost[uniq(diskaccri_lost,sort(diskaccri_lost))]
      diskaccrm_max_uniq_lost = diskaccrm_lost[uniq(diskaccri_lost,sort(diskaccri_lost))]
      
      match2,relosti_uniq, diskaccri_uniq_lost,ind0,ind1
      neverlefti = diskaccri_uniq_lost[where(ind1 EQ -1)]
      neverleftm = diskaccrm_max_uniq_lost[where(ind1 EQ -1)]
      match,relosti_uniq, diskaccri_uniq_lost,ind0,ind1
      relosti_uniq2 = relosti_uniq[ind0]
      diskaccri_uniq_lost2 = diskaccri_uniq_lost[ind1]
      relostt_max_uniq2 = relostt_max_uniq[ind0]
      diskaccrt_max_uniq_lost2 = diskaccrt_max_uniq_lost[ind1]
      relostm_max_uniq2 = relostm_max_uniq[ind0]
      diskaccrm_max_uniq_lost2 = diskaccrm_max_uniq_lost[ind1]
      stilltherei = diskaccri_uniq_lost2[where(diskaccrt_max_uniq_lost2 GT relostt_max_uniq2)] ;Where they were accreted more recently than lost
      stilltherei = [neverlefti,stilltherei]
      stilltherem = diskaccrm_max_uniq_lost2[where(diskaccrt_max_uniq_lost2 GT relostt_max_uniq2)] ;Where they were accreted more recently than lost
      stilltherem = [neverleftm,stilltherem]
      
      smass_accri = mrdfits(dirs[i] + '/grp' + halo[i] + '.accrstars_i.fits',/silent)
      smass_accrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.accrstars_z.fits',/silent)
      smass_accri = smass_accri[where(smass_accrz NE 99)]
      match2,smass_accri,sl.iorderstar,ind0,ind1
      sl_insitu = sl[where(ind1 EQ -1)] ;stars that weren't accreted
      
      match2,sl_insitu.iordergas,[diskaccri_lost],ind0,ind1 
      help,where(ind0 EQ -1) ;Number of stars that are in the final halo but are neither accreted nor form out of accreted gas
      
      match2,iord[diskind],diskaccri_lost,ind0,ind1
      help,where(ind0 EQ -1) ;the number of gparticles identified in final gal as being part of the disk that were never identified as being accreted.  Ideally, should be zero but it is probably fine if it is small
      
      match2,[sl_insitu.iordergas,iord[diskind]],diskaccri_lost,ind0,ind1
      
      match2,iord[diskind],stilltherei,ind2,ind3
      IF (where(ind2 EQ -1))[0] NE -1 THEN indmissing = diskind[where(ind2 eq -1)]
   ENDIF 

   IF 0 THEN BEGIN
; Unique accretion onto halo here
      inflowz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccr_z.fits') 
      inflowm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccr.fits') 
      inflowi = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccr_iord.fits')
      inflow_data = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccr_history.fits',1)
      inflowm[where(inflowm LE 0)] = 0 ;Set mass to zero when it comes in as a star
      inflowt = z_to_t(inflowz)
      inflowmcum  = weighted_histogram(inflowt, weight =  inflowm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
      inflowzcum  = weighted_histogram(inflowt, weight =  inflowm*inflow_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
      nbins = n_elements(inflowmcum)
      
      ind_inflow_uniq = uniq(inflowi,sort(inflowi))
      inflowz_uniq = inflowz[ind_inflow_uniq]
      inflowm_uniq = inflowm[ind_inflow_uniq]
      inflowi_uniq = inflowi[ind_inflow_uniq]
      inflow_data_uniq = inflow_data[ind_inflow_uniq]
      inflowt_uniq = z_to_t(inflowz_uniq)    
      inflowmcum_uniq = weighted_histogram(inflowt_uniq, weight =  inflowm_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
      inflowzcum_uniq = weighted_histogram(inflowt_uniq, weight =  inflowm_uniq*inflow_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
      uniqaccrgmass_total[i] = inflowmcum_uniq[n_elements(inflowzcum_uniq) - 1]
      uniqaccrzmass_total[i] = inflowzcum_uniq[n_elements(inflowzcum_uniq) - 1]
   ENDIF
ENDFOR


;------------ Expelled and Ejected Metal Frac as a Function of Halo Mass
readcol,'~/Datafiles/HIcubes/moster.stars.z0',mhalo,mstar_mhalo,mstar,logmstar
readcol,'~/Zproduced_z0.txt',files_readz,halo_readz,currenthalo_stars_zmass,primehalo_stars_zmass,accrh_stars_zmass,s_zmass,ism_zmass,accr_stars_halo_zmass,format = '(A,I,F,F,F,F,F,F)' ;Produced by metals_v_mass.pro
IF n_elements(files_readz) NE n_elements(files) THEN BEGIN
   print,'The file "Zproduced_z0.txt" is not compatible with this program'
   stop
ENDIF
;Metals produced since the start of the simulation
match2,files+halo,files_readz+strtrim(string(halo_readz),2),indread,indorig
files_readz = files_readz[indread]
halo_readz = halo_readz[indread]
currenthalo_stars_zmass = currenthalo_stars_zmass[indread]
primehalo_stars_zmass = primehalo_stars_zmass[indread];compare with zproduced_total
;accrh_stars_zmass = accrh_stars_zmass[indread]
;s_zmass = s_zmass[indread]
;ism_zmass = ism_zmass[indread]
;accr_stars_halo_zmass = accr_stars_halo_zmass[indread]

cheat = 7 ;150 ;yellow
cexpell = 2 ;60 ;blue
ceject = 10;254 ;red
cpristine = 11 ;20 ;purple
cstar = 5;210 ;green
IF keyword_set(colors) THEN BEGIN
    loadct,39
    distinct_colors,n_colors = 12
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*254/n else colors = colors
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(symbols) THEN symbols = [14,15,17,18,16,34];[4,5]
    colore = [cexpell,ceject,180,200,cheat,30]
    ;expelled, ejected, rvir/2,rvir/5,heated,lost
    IF keyword_set(z_cut) THEN BEGIN
        IF NOT keyword_set(z_colors) THEN z_colors = [7,2,11,0] ;(findgen(nz) + 1)/n_elements(z_bins)*254
        IF NOT keyword_set(z_psym)   THEN z_psym   =  fltarr( nz) + 16;4
    ENDIF
    color_s = cstar;fgcolor
    symbol_s = 46
    color_a = fgcolor
    symbol_a = 16
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(symbols) THEN symbols = [14,15,17,18,16,34];[4,5]
    colore =  [50,254,210,30,110,30]
    IF keyword_set(z_cut) THEN BEGIN
        IF NOT keyword_set(z_colors) THEN z_colors = fltarr( nz) + fgcolor
        IF NOT keyword_set(z_psym)   THEN z_psym   = findgen(nz) + 16;1
    ENDIF
    color_s = 120;fgcolor
    symbol_s = 46
    color_a = fgcolor
    symbol_a = 16
ENDELSE
IF NOT keyword_set(symsize) THEN $
   IF keyword_set(formatthick) THEN symsize = 1.5 ELSE symsize = 2
;IF NOT keyword_set(z_psym) THEN BEGIN
   ej_psym = fltarr(nz) + symbols[1]
   ex_psym = fltarr(nz) + symbols[0]
;ENDIF ELSE BEGIN
;   ej_psym = z_psym
;   ex_psym = z_psym
;ENDELSE

IF keyword_set(stellarmass) THEN BEGIN
   xtitle = 'Stellar Mass [M' + sunsymbol() + ']'
   xrange = [1e6,1e11]
ENDIF ELSE BEGIN
   xtitle = 'Virial Mass [M' + sunsymbol() + ']'
   xrange = [3e9,1e12]
ENDELSE
IF keyword_set(stellarmass) THEN BEGIN
   xmass = smass 
   xmasst = smasst
ENDIF ELSE BEGIN
   xmass = vmass
   xmasst = vmasst
ENDELSE

;------------Net Expelled and Ejected Metal Mass as a Function of Halo Mass
IF keyword_set(outplot) THEN  device,filename = outplot + '_neteject_expell_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,primehalo_stars_zmass,psym = symcat(symbols[0]),/ylog,/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Metal Mass [M' + sunsymbol() + ']',symsize = symsize,/nodata,yrange = [1e4,8e8];Metals produced by stars since start of time
oplot,xmass,(relostmassmet[*,nz - 1] - diskgmass_lostmet[*,nz - 1]),psym = symcat(symbols[4]),color = colore[4],symsize = symsize;Net Metals lost
oplot,xmass,(relostmassmet[*,nz - 1] - diskgmass_lostmet[*,nz - 1]),psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),color = fgcolor,symsize = symsize
;oplot,xmass,(reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1]),psym = symcat(symbols[1]),color = colore[1],symsize = symsize;Net Metals ejected
;oplot,xmass,(reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1]),psym = symcat(sym_outline(symbols[1])),color = fgcolor,symsize = symsize
oplot,xmass,reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1],psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),color = fgcolor,symsize = symsize
;oplot,xmass,(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(17),color = fgcolor,symsize = symsize
;oplot,xmass,(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[4])),color = fgcolor,symsize = symsize
oplot,xmass, primehalo_stars_zmass + diskgmassmet_uniq_total, psym = symcat(symbol_s),color = color_s,symsize = symsize ;Metals produced by the main progenitor (consider using currenthalo_stars_zmass if I want metals produced by stars in main halo at z = 0)
oplot,xmass, primehalo_stars_zmass + diskgmassmet_uniq_total, psym = symcat(sym_outline(symbol_s),thick = thicks[0]),color = fgcolor,symsize = symsize
;legend,['Net metals expelled','Net metals ejected','Net metals lost','Metals produced','Metals produced + accreted'],psym = [symbols[0:1],symbols[4],symbol_s,17],color = [colore[0:1],colore[4],color_s,fgcolor],/bottom,/right,box = 0
;legend,['Net metals expelled','Net metals lost','Metals produced','Metals produced + accreted'],psym = [symbols[0],symbols[4],symbol_s,17],color = [colore[0],colore[4],color_s,fgcolor],/bottom,/right,box = 0
legend,['Net metals expelled','Net metals removed','Metals produced + accreted'],psym = [symbols[0],symbols[4],symbol_s],color = [colore[0],colore[4],color_s],/bottom,/right,box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;Total Metals
distinct_colors,n_colors = 12
IF keyword_set(outplot) THEN  device,filename = outplot + '_eject_expell_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,reexpellmassmet[*,nz - 1],psym = symcat(symbols[0]),/ylog,/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Metal Mass [M' + sunsymbol() + ']',symsize = symsize,/nodata,yrange = [1e4,2e9]
oplot,xmass,reexpellmassmet[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),color = fgcolor,symsize = symsize
oplot,xmass,reejectmassmet[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize;Metals ejected
oplot,xmass,reejectmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),color = fgcolor,symsize = symsize
oplot,xmass,relostmassmet[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize;Metals lost
oplot,xmass,relostmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),color = fgcolor,symsize = symsize
;oplot,xmass,(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(17),color = fgcolor,symsize = symsize
;oplot,xmass,(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[4])),color = fgcolor,symsize = symsize
oplot,xmass, primehalo_stars_zmass, psym = symcat(symbol_s),color = color_s,symsize = symsize ;Metals produced by the main progenitor (consider using currenthalo_stars_zmass if I want metals produced by stars in main halo at z = 0)
oplot,xmass, primehalo_stars_zmass, psym = symcat(sym_outline(symbol_s),thick = thicks[0]),color = fgcolor,symsize = symsize
;legend,['Metals expelled','Metals ejected','Metals lost','Metals produced','Metals produced + accreted'],psym = [symbols[0:1],symbols[4],symbol_s,17],color = [colore[0:1],colore[4],color_s,fgcolor],/bottom,/right,box = 0
legend,['Metals expelled','Metals ejected','Metals removed','Metals produced'],psym = [symbols[0:1],symbols[4],symbol_s],color = [colore[0:1],colore[4],color_s],/bottom,/right,box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;What fraction of metals were produced by the stars in the main halo?
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_produced_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
;plot,xmass,primehalo_stars_zmass/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(4,thick = 8),/xlog,xrange = xrange,xtitle = xtitle,ytitle = textoidl('Z_{produced}/(Z_{produced} + Z_{externally accreted})'),symsize = symsize,yrange = [0.85,1.01]
;oplot,[1e6,1e12],[1,1]
binsize = 0.1
plot,alog10([vmass[0],vmass[0]]),[diskgmassmet_uniq_total[0],diskgmassmet_uniq_total[0]]/primehalo_stars_zmass[0],xrange = [min(alog10(xrange)) - binsize/2,max(alog10(xrange)) + binsize/2],xtitle = textoidl('Log(M_{vir}/M')+sunsymbol()+')',ytitle = textoidl('Z_{accreted as gas}/Z_{produced}'),thick = 14,/nodata,yrange = [0,0.15]
distinct_colors,n_colors = 12
FOR i = 0, n_elements(dirs) - 1 DO polyfill,[alog10(xmass[i]) - binsize/2,alog10(xmass[i]) - binsize/2,alog10(xmass[i]) + binsize/2,alog10(xmass[i]) + binsize/2],[0,diskgmassmet_uniq_total[i],diskgmassmet_uniq_total[i],0]/primehalo_stars_zmass[i],/line_fill,orientation = 45,color = cpristine
FOR i = 0, n_elements(dirs) - 1 DO oplot,[alog10(xmass[i]) - binsize/2,alog10(xmass[i]) - binsize/2,alog10(xmass[i]) + binsize/2,alog10(xmass[i]) + binsize/2],[0,diskgmassmet_uniq_total[i],diskgmassmet_uniq_total[i],0]/primehalo_stars_zmass[i]
;,thick = 14,linestyle = 0;,color = 
;Metals accreted as gas/metals produced by stars in main progenitor
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Expelled Metal Mass as a Function of Halo Mass and Time
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_zexpell_mass_time.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xmass,reexpellmassmet[*,0],/ylog,/xlog,xrange = xrange,yrange = [1e3,2e8],xtitle = xtitle,ytitle = 'Mass of Metals Expelled [M' + sunsymbol() + ']',/nodata
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xmass,reexpellmassmet[*,iz],psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize  
       oplot,xmass,reexpellmassmet[*,iz],psym = symcat(sym_outline(ex_psym[iz]),thick = thicks[0]),symsize = symsize
    ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = ex_psym,color  = z_colors,/right,/bottom,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
ENDIF

;------------ Ejected Metal Mass as a Function of Halo Mass and Time
distinct_colors,n_colors = 12
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_zeject_mass_time.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xmass,reejectmassmet[*,0],/ylog,/xlog,xrange = xrange,yrange = [1e3,2e8],xtitle = xtitle,ytitle = 'Mass of Metals Ejected [M' + sunsymbol() + ']',/nodata
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xmass,reejectmassmet[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz] ,symsize = symsize 
       oplot,xmass,reejectmassmet[*,iz],psym = symcat(sym_outline(ej_psym[iz]),thick = thicks[0]),symsize = symsize
    ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = ej_psym,color  = z_colors,/right,/bottom,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
ENDIF

;------------ Fraction of Metal Mass Lost that is Expelled & Ejected
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_lostgas_eject_expell_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,reexpellmassmet[*,nz - 1]/relostmassmet[*,nz - 1],/xlog,xrange = xrange,yrange = [0,0.6],xtitle = xtitle,ytitle = 'Metal Outflow Mass/Metal Mass Removed',/nodata
distinct_colors,n_colors = 12
oplot,xmass,reexpellmassmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmassmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xmass, reejectmassmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, reejectmassmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xmass,(relostmassmet[*,nz - 1] - diskgmass_lostmet[*,nz - 1])/relostmassmet[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize;Net Metals lost
oplot,xmass,(relostmassmet[*,nz - 1] - diskgmass_lostmet[*,nz - 1])/relostmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),color = fgcolor,symsize = symsize
legend,['Metal Mass Expelled','Metal Mass Ejected','Net Metal Mass Removed'],psym = [symbols[0:1],symbols[4]],color = [colore[0:1],colore[4]],box = 0,/top,/right;position = [2.2e11,0.4];/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop


IF 1 THEN BEGIN
;------------ Fraction of Net Metals in produced+accreted that are Expelled & Ejected
;XXX Change this to net metals ejected/expelled
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_netgas_eject_expell_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
distinct_colors,n_colors = 12
plot,[1,1e12],[1,1],/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Fraction of Net Metals Removed/Ejected/Expelled',yrange = [0,1.2]
;zavail_zmax_total + diskgmassmet_uniq_total
oplot,xmass,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xmass, (reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, (reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xmass, (relostmassmet[*,nz - 1] - diskgmass_lostmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass, (relostmassmet[*,nz - 1]- diskgmass_lostmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
legend,['Net Metal Mass Expelled','Net Metal Mass Ejected','Net Metal Mass Removed'],psym = [symbols[0:1],symbols[4]],color = [colore[0:1],colore[4]],box = 0,/top,/right;position = [2.2e11,0.4];/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

IF 1 THEN BEGIN
;------------ Fraction of Metals in produced+accreted that are Expelled & Ejected
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_eject_expell_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
distinct_colors,n_colors = 12
plot,[1,1e12],[1,1],/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Fraction of Metals Removed/Ejected/Expelled',yrange = [0,3.5]
;zavail_zmax_total + diskgmassmet_uniq_total
oplot,xmass,(reexpellmassmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,(reexpellmassmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xmass, (reejectmassmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, (reejectmassmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xmass, (relostmassmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass, (relostmassmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
;oplot,xmass,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[0])),symsize = symsize,color = colore[0]
;oplot,xmass, (reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[1])),color = colore[1],symsize = symsize
;oplot,xmass, (relostmassmet[*,nz - 1]- diskgmass_lostmet[*,nz - 1])/(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[4])),symsize = symsize,color = colore[4]
legend,['Metal Mass Expelled','Metal Mass Ejected','Metal Mass Removed'],psym = [symbols[0:1],symbols[4]],color = [colore[0:1],colore[4]],box = 0,/top,/left;position = [2.2e11,0.4];/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

;Fraction of metals reaccreted
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_reaccreted_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
distinct_colors,n_colors = 12
plot,xmass,diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1],/nodata,/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Fraction of Metals and Mass Reaccreted',yrange = [0,1.1],title = 'z = 0'
oplot,[1e9,1e12],[1,1],linestyle = 1
oplot,xmass,diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass,diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
;oplot,xmass,diskgmass_expellmet[*,nz - 1]/reexpellmassmet[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
;oplot,xmass,diskgmass_expellmet[*,nz - 1]/reexpellmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[0])),symsize = symsize
;oplot,xmass,diskgmassmet[*,nz - 1]/reejectmassmet[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
;oplot,xmass,diskgmassmet[*,nz - 1]/reejectmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[1])),symsize = symsize

;oplot,xmass,diskgmass_lost[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass,diskgmass_lost[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[4]),thick = 8),symsize = symsize,color = colore[4]
;oplot,xmass,diskgmass_expell[*,nz - 1]/reexpellmass[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,diskgmass_expell[*,nz - 1]/reexpellmass[*,nz - 1],psym = symcat(sym_outline(symbols[0]),thick = 8),symsize = symsize,color = colore[0]
;oplot,xmass,diskgmass[*,nz - 1]/reejectmass[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass,diskgmass[*,nz - 1]/reejectmass[*,nz - 1],psym = symcat(sym_outline(symbols[1]),thick = 8),symsize = symsize,color = colore[1]
legend,['Reaccretion after removal from disk','Reaccretion after ejection',textoidl('Reaccretion after expullsion beyond R_{vir}')],psym = [sym_outline(symbols[4]),sym_outline(symbols[1]),sym_outline(symbols[0])],color = [colore[4],colore[1],colore[0]],box = 0,/top,/left

IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of metals reaccreted
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_reaccreted_zmass_z0_z2.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*1.6,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize*1.6
distinct_colors,n_colors = 12
multiplot,[1,2],mytitle = 'Fraction of Metals and Mass Reaccreted',mxtitle = xtitle,mxTitSize = !x.charsize,myTitSize = !y.charsize
plot,xmass,diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1],/nodata,/xlog,xrange = xrange,yrange = [0,1]
oplot,[1e9,1e12],[1,1],linestyle = 1
oplot,xmass,diskgmass_lost[*,0]/relostmass[*,0],psym = symcat(sym_outline(symbols[4]),thick = 8),symsize = symsize,color = colore[4]
oplot,xmass,diskgmass_expell[*,0]/reexpellmass[*,0],psym = symcat(sym_outline(symbols[0]),thick = 8),symsize = symsize,color = colore[0]
oplot,xmass,diskgmass[*,0]/reejectmass[*,0],psym = symcat(sym_outline(symbols[1]),thick = 8),symsize = symsize,color = colore[1]
oplot,xmass,diskgmass_lostmet[*,0]/relostmassmet[*,0],psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass,diskgmass_lostmet[*,0]/relostmassmet[*,0],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
legend,['Reaccretion after removal from disk','Reaccretion after ejection',textoidl('Reaccretion after expullsion beyond R_{vir}')],psym = [sym_outline(symbols[4]),sym_outline(symbols[1]),sym_outline(symbols[0])],color = [colore[4],colore[1],colore[0]],box = 0,/top,/right
legend,['z = 2'],box = 0,/top,/left
multiplot

plot,xmass,diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1],/nodata,/xlog,xrange = xrange,yrange = [0,1]
distinct_colors,n_colors = 12
oplot,[1e9,1e12],[1,1],linestyle = 1
oplot,xmass,diskgmass_lost[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[4]),thick = 8),symsize = symsize,color = colore[4]
oplot,xmass,diskgmass_expell[*,nz - 1]/reexpellmass[*,nz - 1],psym = symcat(sym_outline(symbols[0]),thick = 8),symsize = symsize,color = colore[0]
oplot,xmass,diskgmass[*,nz - 1]/reejectmass[*,nz - 1],psym = symcat(sym_outline(symbols[1]),thick = 8),symsize = symsize,color = colore[1]
oplot,xmass,diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass,diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
legend,['z = 0'],box = 0,/top,/left
multiplot,/reset
IF keyword_set(outplot) THEN device, /close ELSE stop


;Fraction of mass reaccreted
;IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_reaccreted_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
;plot,xmass,diskgmass_lost[*,nz - 1]/relostmass[*,nz - 1],/nodata,/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Fraction of Mass Reaccreted',yrange = [0,1]
;oplot,xmass,diskgmass_lost[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize
;oplot,xmass,diskgmass_lost[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[4])),symsize = symsize
;oplot,xmass,diskgmass_expell[*,nz - 1]/reexpellmass[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
;oplot,xmass,diskgmass_expell[*,nz - 1]/reexpellmass[*,nz - 1],psym = symcat(sym_outline(symbols[0])),symsize = symsize
;oplot,xmass,diskgmass[*,nz - 1]/reejectmass[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
;oplot,xmass,diskgmass[*,nz - 1]/reejectmass[*,nz - 1],psym = symcat(sym_outline(symbols[1])),symsize = symsize
;IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_reaccreted_metallicity.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
distinct_colors,n_colors = 12
plot,xmass,diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1],/nodata,/xlog,xrange = xrange,xtitle = xtitle,ytitle = 'Fraction of Metals Reaccreted/Fraction Mass Reaccreted',yrange = [0,1]
;oplot,[1e9,1e12],[1,1],linestyle = 1
oplot,xmass,(diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1])/(diskgmass_lost[*,nz - 1]/relostmass[*,nz - 1]),psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass,(diskgmass_lostmet[*,nz - 1]/relostmassmet[*,nz - 1])/(diskgmass_lost[*,nz - 1]/relostmass[*,nz - 1]),psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
oplot,xmass,(diskgmass_expellmet[*,nz - 1]/reexpellmassmet[*,nz - 1])/(diskgmass_expell[*,nz - 1]/reexpellmass[*,nz - 1]),psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,(diskgmass_expellmet[*,nz - 1]/reexpellmassmet[*,nz - 1])/(diskgmass_expell[*,nz - 1]/reexpellmass[*,nz - 1]),psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xmass,(diskgmassmet[*,nz - 1]/reejectmassmet[*,nz - 1])/(diskgmass[*,nz - 1]/reejectmass[*,nz - 1]),psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass,(diskgmassmet[*,nz - 1]/reejectmassmet[*,nz - 1])/(diskgmass[*,nz - 1]/reejectmass[*,nz - 1]),psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
legend,['Reaccretion after removal from disk','Reaccretion after ejection',textoidl('Reaccretion after expullsion beyond R_{vir}')],psym = [symbols[4],symbols[1],symbols[0]],color = [colore[4],colore[1],colore[0]],box = 0,/bottom,/left

IF keyword_set(outplot) THEN device, /close ELSE stop

IF 1 THEN BEGIN
;------------ Fraction of Metals Expelled
;XXX Need total metals in the disk
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_expelled_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
distinct_colors,n_colors = 12
plot,xmass,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),/xlog,xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = 'Net Fraction of Metals Expelled',/nodata;yrange = [0,1]
oplot,xmass,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

IF 1 THEN BEGIN
;------------ Fraction of Gas Mass Expelled Over Time
;Fraction of metals expelled over time
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_expelled_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    distinct_colors,n_colors = 12
    plot,xmass,(reexpellmassmet - diskgmass_expellmet)/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = 'Net Fraction of Metals Expelled',/nodata,symsize = symsize; yrange = [0,1]
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xmass,(reexpellmassmet[*,iz] - diskgmass_expellmet[*,iz])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize
       oplot,xmass,(reexpellmassmet[*,iz] - diskgmass_expellmet[*,iz])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(sym_outline(ex_psym[iz]),thick = thicks[0]),symsize = symsize
    ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = ex_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF
ENDIF

IF 1 THEN BEGIN
;------------ Fraction of Metals produced+accreted that is Ejected
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_ejected_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
distinct_colors,n_colors = 12
plot,xmass, (reejectmassmet[*,nz - 1]- diskgmassmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),/xlog,xrange = xrange,yrange = [0,1.0],xtitle = xtitle,ytitle = 'Net Fraction of Metals Ejected',/nodata;yrange = [0,1]
oplot,xmass, (reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, (reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

IF 1 THEN BEGIN
;------------ Fraction of Metal Mass Ejected Over Time
;XXX Need total metals in the disk
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_ejected_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    distinct_colors,n_colors = 12
    plot,xmass,(reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange = [0,1.3],xtitle = xtitle,symsize = symsize,ytitle = 'Net Fraction of Metals Ejected',/nodata; yrange = [0,1]
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xmass,(reejectmassmet[*,iz] - diskgmassmet[*,iz])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(ej_psym[iz]),color  = z_colors[iz],symsize = symsize
       oplot,xmass,(reejectmassmet[*,iz] - diskgmassmet[*,iz])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(sym_outline(ej_psym[iz]),thick = thicks[0]),symsize = symsize
    ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF
ENDIF 

;------------Fraction of Metal lost?

;------------ Fraction of Lost Metals Ejected or Expelled
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_expelled_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
distinct_colors,n_colors = 12
plot,xmass,reejectmassmet[*,nz - 1]/relostmassmet[*,nz - 1],ytitle = textoidl('M_Z/M_{Z,lost}'),xtitle = xtitle,yrange = [0,0.6],xrange = xrange,psym = symcat(18),/nodata,/xlog
;loadct,0
IF keyword_set(outplot) THEN bar_thick = 14 ELSE bar_thick = 4
FOR i = 0,n_elements(xmass) - 1 DO oplot,[xmass[i],xmass[i]],[0,1.0],color = 100,thick = bar_thick
IF keyword_set(colors) THEN distinct_colors,n_colors = 12 ;loadct,39
FOR i = 0,n_elements(xmass) - 1 DO oplot,[xmass[i],xmass[i]],[0,reejectmassmet[i,nz - 1]/relostmassmet[i,nz - 1]],color = colore[1],thick = bar_thick
FOR i = 0,n_elements(xmass) - 1 DO oplot,[xmass[i],xmass[i]],[0,reexpellmassmet[i,nz - 1]/relostmassmet[i,nz - 1]],color = colore[0],thick = bar_thick
oplot,xmass,reejectmassmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
oplot,xmass,reejectmassmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xmass,reexpellmassmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
oplot,xmass,reexpellmassmet[*,nz - 1]/relostmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
legend,['Metals Expelled','Metals Ejected'],psym = symbols[0:1],color = colore[0:1],/top,/right,box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Fraction of Ejected Metals Expelled Over Time
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_expelled_zmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    distinct_colors,n_colors = 12
    plot,xmass,reexpellmassmet[*,nz - 1]/reejectmassmet[*,nz - 1],/xlog,xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = 'Fraction of Metal Mass Ejected that is Expelled',psym = 4,symsize = symsize
;    FOR iz = 0, nz - 1 DO BEGIN
;       oplot,xmass,reexpellmassmet[*,iz]/reejectmassmet[*,iz],psym = symcat(z_psym[iz]),color  = z_colors[iz],symsize = symsize
;       oplot,xmass,reexpellmassmet[*,iz]/reejectmassmet[*,iz],psym = symcat(sym_outline(z_psym[iz]),thick = thicks[0]),symsize = symsize
;    ENDFOR
;    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = z_psym,color  = z_colors,/right,/top,box=0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
ENDIF

mlrange = [1e9,1e12]
mlrange = [10,1000]
;------------ Mass Loading
xarr = vcirc ;xmass
xarrt = vcirct
yield = 0.01 ;0.015
yield_str = '0.01' ;'0.015'
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,reexpellmassmet[*,nz - 1]/sfmass_total/yield,yrange = [0.001/yield,0.05/yield],ytitle = textoidl('Effective metal mass loading factor/' + yield_str),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [16,140],/xlog,xstyle = 1
distinct_colors,n_colors = 12
oplot,xarr,reexpellmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xarr,reexpellmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xarr,reejectmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xarr,reejectmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xarr,relostmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xarr,relostmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
fits_expell = robust_linefit( alog10(xarr), alog10(reexpellmass[*,nz - 1]/sfmass_total), reexpellmass_fit, sigma_expell )
fits_eject = robust_linefit( alog10(xarr), alog10(reejectmass[*,nz - 1]/sfmass_total/yield), reejectmass_fit, sigma_eject )
print,'Expell (total): ',fits_expell
print,'Eject (total): ',fits_eject
;oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),color = colore[0],linestyle = 0,psym = -3,thick = thicks[0]
;oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),color = colore[1],linestyle = 0,psym = -3,thick = thicks[0]
;legend,['Gas Mass Expelled/Stellar Mass Formed','Gas Mass Ejected/Stellar Mass Formed'],psym = symbols,color = colore,/bottom,/left,box=0
;legend,['Metals Expelled','Metals Ejected'],psym = symbols[0:1],color = colore[0:1],/bottom,/left,box=0
legend,['Metals Expelled','Metals Ejected','Metals Removed from Disk'],psym = [symbols[0],symbols[1],symcat(symbols[4])],color = [colore[0],colore[1],colore[4]],/bottom,/left,box=0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Mass Loading using Net metal loss rather than total
;             metal loss
xarr = vcirc ;xmass
xarrt = vcirct
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_netmass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/sfmass_total/yield,yrange = [0.001/yield,0.05/yield],ytitle = textoidl('Net metal loss/stellar mass formed/' + yield_str),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [16,140],/xlog,xstyle = 1
distinct_colors,n_colors = 12
oplot,xarr,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/sfmass_total/yield,psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xarr,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/sfmass_total/yield,psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xarr,(reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1])/sfmass_total/yield,psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xarr,(reejectmassmet[*,nz - 1] - diskgmassmet[*,nz - 1])/sfmass_total/yield,psym = symcat(sym_outline(symbols[1])),thick = thicks[0],symsize = symsize
oplot,xarr,(relostmassmet[*,nz - 1] - diskgmass_lostmet[*,nz - 1])/sfmass_total/yield,psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xarr,(relostmassmet[*,nz - 1] - diskgmass_lostmet[*,nz - 1])/sfmass_total/yield,psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
fits_expell = robust_linefit( alog10(xarr), alog10((reexpellmassmet[*,nz - 1]-diskgmass_expellmet[*,nz - 1])/sfmass_total/yield), reexpellmass_fit, sigma_expell )
fits_eject = robust_linefit( alog10(xarr), alog10((reejectmassmet[*,nz - 1]-diskgmassmet[*,nz - 1])/sfmass_total/yield), reejectmass_fit, sigma_eject )
fits_lost = robust_linefit( alog10(xarr), alog10((relostmassmet[*,nz - 1]-diskgmass_lostmet[*,nz - 1])/sfmass_total/yield), relostmass_fit, sigma_eject )
print,'Expell (total): ',fits_expell
print,'Eject (total): ',fits_eject
print,'Lost (total): ',fits_lost
;oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),color = colore[0],linestyle = 0,psym = -3,thick = thicks[0]
;oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),color = colore[1],linestyle = 0,psym = -3,thick = thicks[0]
;oplot,mlrange,10^fits_lost[0]*mlrange^(fits_lost[1]),color = colore[4],linestyle = 0,psym = -3,thick = thicks[0]
;legend,['Gas Mass Expelled/Stellar Mass Formed','Gas Mass Ejected/Stellar Mass Formed'],psym = symbols,color = colore,/bottom,/left,box=0
;legend,['Metals Expelled','Metals Ejected'],psym = symbols[0:1],color = colore[0:1],/bottom,/left,box=0
legend,['Metals Expelled','Metals Ejected','Metals Removed from Disk'],psym = [symbols[0],symbols[1],symcat(symbols[4])],color = [colore[0],colore[1],colore[4]],/bottom,/left,box=0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Effective Mass Loading*outflow metallicity/ISM metallicity
xarr = vcirc ;xmass
xarrt = vcirct
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_netmass_z.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,(reexpellmass[*,nz - 1] - diskgmass_expell[*,nz - 1])/sfmass_total,ytitle = textoidl('Net mass loss/M_{Stellar, formed}*Z_{wind}/Z_{ISM}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [16,140],/xlog,xstyle = 1,yrange = [0.1,40]
distinct_colors,n_colors = 12
oplot,xarr,10^0.2*(reejectmass[*,nz - 1] - diskgmass[*,nz - 1])/sfmass_total,psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xarr,10^0.2*(reejectmass[*,nz - 1] - diskgmass[*,nz - 1])/sfmass_total,psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xarr,10^0.2*(reexpellmass[*,nz - 1] - diskgmass_expell[*,nz - 1])/sfmass_total,psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xarr,10^0.2*(reexpellmass[*,nz - 1] - diskgmass_expell[*,nz - 1])/sfmass_total,psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xarr,10^0.2*(relostmass[*,nz - 1] - diskgmass_lost[*,nz - 1])/sfmass_total,psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xarr,10^0.2*(relostmass[*,nz - 1] - diskgmass_lost[*,nz - 1])/sfmass_total,psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
fits_expell = robust_linefit( alog10(xarr), alog10(10^0.2*(reexpellmass[*,nz - 1]-diskgmass_expell[*,nz - 1])/sfmass_total/yield), reexpellmass_fit, sigma_expell )
fits_eject = robust_linefit( alog10(xarr), alog10(10^0.2*(reejectmass[*,nz - 1]-diskgmass[*,nz - 1])/sfmass_total/yield), reejectmass_fit, sigma_eject )
fits_lost = robust_linefit( alog10(xarr), alog10(10^0.2*(relostmass[*,nz - 1]-diskgmass_lost[*,nz - 1])/sfmass_total/yield), relostmass_fit, sigma_eject )
print,'Expell (total): ',fits_expell
print,'Eject (total): ',fits_eject
print,'Lost (total): ',fits_lost
;oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),color = colore[0],linestyle = 0,psym = -3,thick = thicks[0]
;oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),color = colore[1],linestyle = 0,psym = -3,thick = thicks[0]
;oplot,mlrange,10^fits_lost[0]*mlrange^(fits_lost[1]),color = colore[4],linestyle = 0,psym = -3,thick = thicks[0]
;legend,['Gas Mass Expelled/Stellar Mass Formed','Gas Mass Ejected/Stellar Mass Formed'],psym = symbols,color = colore,/bottom,/left,box=0
;legend,['Metals Expelled','Metals Ejected'],psym = symbols[0:1],color
;= colore[0:1],/bottom,/left,box=0
oplot,findgen(200),(78/findgen(200))^3.81+0.19 ;Best fit to Tremonti '04 from Peeples+'11
oplot,findgen(200),(111.8/findgen(200))^2.31+1.35 ;PPO4N2 from Peeples+ 2011
oplot,findgen(200),(55.5/findgen(200))^3.04+0.32 ;KK04 from Peeples+2011
legend,['Metals Expelled','Metals Ejected','Metals Removed from Disk'],psym = [symbols[0],symbols[1],symcat(symbols[4])],color = [colore[0],colore[1],colore[4]],/bottom,/left,box=0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Effective Mass Loading*outflow metallicity/ISM
;             metallicity vs. circular velocity
xarr = vcirc ;xmass
xarrt = vcirct
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_netmass_z_vcirc.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],ytitle = textoidl('Net mass loss/M_{Stellar, formed}*Z_{wind}/Z_{ISM}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [16,140],/xlog,xstyle = 1,yrange = [0.01,100]
distinct_colors,n_colors = 12
oplot,xarr,(reejectmassmetr[*,nz - 1] - diskgmassmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xarr,(reejectmassmetr[*,nz - 1] - diskgmassmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xarr,10^0.2*(reejectmass[*,nz - 1] - diskgmass[*,nz - 1])/sfmass_total,psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xarr,(reexpellmassmetr[*,nz - 1] - diskgmass_expellmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xarr,(reexpellmassmetr[*,nz - 1] - diskgmass_expellmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xarr,10^0.2*(reexpellmass[*,nz - 1] - diskgmass_expell[*,nz - 1])/sfmass_total,psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xarr,(relostmassmetr[*,nz - 1] - diskgmass_lostmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xarr,(relostmassmetr[*,nz - 1] - diskgmass_lostmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
oplot,xarr,10^0.2*(relostmass[*,nz - 1] - diskgmass_lost[*,nz - 1])/sfmass_total,psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
fits_expell = robust_linefit( alog10(xarr), alog10((reexpellmassmetr[*,nz - 1]-diskgmass_expellmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1]), reexpellmass_fit, sigma_expell )
fits_eject = robust_linefit( alog10(xarr), alog10((reejectmassmetr[*,nz - 1]-diskgmassmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1]), reejectmass_fit, sigma_eject )
fits_lost = robust_linefit( alog10(xarr), alog10((relostmassmetr[*,nz - 1]-diskgmass_lostmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1]), relostmass_fit, sigma_eject )
print,'Expell (total): ',fits_expell
print,'Eject (total): ',fits_eject
print,'Lost (total): ',fits_lost
legend,['Metals Expelled','Metals Ejected','Metals Removed from Disk'],psym = [symbols[0],symbols[1],symcat(symbols[4])],color = [colore[0],colore[1],colore[4]],/bottom,/left,box=0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Effective Mass Loading*outflow metallicity/ISM
;             metallicity vs. stellar mass
xarr = smass ;xmass
xarrt = smasst
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_netmass_z_smass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,(reexpellmassmetr[*,nz - 1] - diskgmass_expellmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],ytitle = textoidl('Net mass loss/M_{Stellar, formed}*Z_{wind}/Z_{ISM}'),xtitle = textoidl('M_{stellar} [M_{\odot}]'),/nodata,/ylog,xrange = [1e6,2e10],/xlog,xstyle = 1,yrange = [0.1,100]
distinct_colors,n_colors = 12
oplot,xarr,(reejectmassmetr[*,nz - 1] - diskgmassmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xarr,(reejectmassmetr[*,nz - 1] - diskgmassmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xarr,(reexpellmassmetr[*,nz - 1] - diskgmass_expellmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xarr,(reexpellmassmetr[*,nz - 1] - diskgmass_expellmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xarr,(relostmassmetr[*,nz - 1] - diskgmass_lostmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xarr,(relostmassmetr[*,nz - 1] - diskgmass_lostmetr[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
fits_expell = robust_linefit( alog10(xarr), alog10((reexpellmassmet[*,nz - 1]-diskgmass_expellmet[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1]), reexpellmass_fit, sigma_expell )
fits_eject = robust_linefit( alog10(xarr), alog10((reejectmassmet[*,nz - 1]-diskgmassmet[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1]), reejectmass_fit, sigma_eject )
fits_lost = robust_linefit( alog10(xarr), alog10((relostmassmet[*,nz - 1]-diskgmass_lostmet[*,nz - 1])/sfmassr[*,nz - 1]/metallicity_ISMr[*,nz - 1]), relostmass_fit, sigma_eject )
print,'Expell (total): ',fits_expell
print,'Eject (total): ',fits_eject
print,'Lost (total): ',fits_lost
legend,['Metals Expelled','Metals Ejected','Metals Removed from Disk'],psym = [symbols[0],symbols[1],symcat(symbols[4])],color = [colore[0],colore[1],colore[4]],/bottom,/left,box=0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Effective Mass Loading*outflow metallicity/ISM
;             metallicity vs. stellar mass -- Lost at different redshifts
xarr = smass ;xmass
xarrt = smasst
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_netmass_z_smass_lostred.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,(reexpellmassmet - diskgmass_expellmet)/sfmassr/metallicity_ISMr,ytitle = textoidl('\zeta_{net}'),xtitle = textoidl('M_{stellar} [M_{\odot}]'),/nodata,/ylog,xrange = [6e5,2e10],/xlog,xstyle = 1,yrange = [0.3,300]
distinct_colors,n_colors = 12
FOR iz = 0, nz - 1 DO BEGIN
   oplot,xarrt[*,iz],(relostmassmetr[*,iz] - diskgmass_lostmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(symbols[4]),color =  z_colors[iz],symsize = symsize
   oplot,xarrt[*,iz],(relostmassmetr[*,iz] - diskgmass_lostmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
ENDFOR
;FOR iz = 0, nz - 1 DO oplot,xarrt[*,iz],(relostmassmetr[*,iz] - diskgmass_lostmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(symbols[4]),color =  z_colors[iz],symsize = symsize
;FOR iz = 0, nz - 1 DO oplot,xarrt[*,iz],(relostmassmetr[*,iz] - diskgmass_lostmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(sym_outline(symbols[4])),symsize = symsize
legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = symbols[4],color  = z_colors,/right,/top,box = 0
pos_ind = where(((relostmassmetr-diskgmass_lostmetr)/sfmassr/metallicity_ISMr) GT 0)
fits_lost = robust_linefit( alog10(xarrt[pos_ind]), alog10((relostmassmetr[pos_ind]-diskgmass_lostmetr[pos_ind])/sfmassr[pos_ind]/metallicity_ISMr[pos_ind]), relostmass_fit, sigma_eject )
print,'Zwta (all z) : Lost (total): ',fits_lost
oplot,[1e5,1e11],10^fits_lost[0]*[1e5,1e11]^(fits_lost[1]),linestyle = 0,psym = -3,thick = thicks[0]
IF keyword_set(outplot) THEN device, /close ELSE stop

xarr = vcirc ;xmass
xarrt = vcirct
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_netmass_z_vcirc_lostred.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,(reexpellmassmet - diskgmass_expellmet)/sfmassr/metallicity_ISMr,ytitle = textoidl('\zeta_{net}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [10,200],/xlog,xstyle = 1,yrange = [0.3,300]
distinct_colors,n_colors = 12
FOR iz = 0, nz - 1 DO BEGIN
    oplot,xarrt[*,iz],(relostmassmetr[*,iz] - diskgmass_lostmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(symbols[4]),color =  z_colors[iz],symsize = symsize
    oplot,xarrt[*,iz],(relostmassmetr[*,iz] - diskgmass_lostmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
ENDFOR
;FOR iz = 0, nz - 1 DO oplot,xarrt[*,iz],(relostmassmetr[*,iz] - diskgmass_lostmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(symbols[4]),color =  z_colors[iz],symsize = symsize
;FOR iz = 0, nz - 1 DO oplot,xarrt[*,iz],(relostmassmetr[*,iz] - diskgmass_lostmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(sym_outline(symbols[4])),symsize = symsize
legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = symbols[4],color  = z_colors,/left,/bottom,box = 0
pos_ind = where(((relostmassmetr-diskgmass_lostmetr)/sfmassr/metallicity_ISMr) GT 0)
fits_lost = robust_linefit( alog10(xarrt[pos_ind]), alog10((relostmassmetr[pos_ind]-diskgmass_lostmetr[pos_ind])/sfmassr[pos_ind]/metallicity_ISMr[pos_ind]), relostmass_fit, sigma_eject )
print,'Lost (total): ',fits_lost
oplot,[10,1e3],10^fits_lost[0]*[10,1e3]^(fits_lost[1]),linestyle = 0,psym = -3,thick = thicks[0]
oplot,findgen(50),(78/findgen(50))^3.81+0.19,linestyle = 2 ;Best fit to Tremonti '04 from Peeples+'11
oplot,findgen(50),(111.8/findgen(50))^2.31+1.35,linestyle = 1 ;PPO4N2 from Peeples+ 2011
oplot,findgen(50),(55.5/findgen(50))^3.04+0.32,linestyle = 3 ;KK04 from Peeples+2011
vcircarr =  findgen(150) + 50
oplot,vcircarr,(78/vcircarr)^3.81+0.19,linestyle = 2,thick = 10 ;Best fit to Tremonti '04 from Peeples+'11
oplot,vcircarr,(111.8/vcircarr)^2.31+1.35,linestyle = 1,thick = 10 ;PPO4N2 from Peeples+ 2011
oplot,vcircarr,(55.5/vcircarr)^3.04+0.32,linestyle = 3,thick = 10 ;KK04 from Peeples+2011
legend,['Exponential fit to simulations','Fit to Tremonti et al., 2004','Fit to Pettini & Pegal, 2004, NII','Fit to Kobulnicky & Kewley 2004'],linestyle = [0,2,1,3],/top,/right,box = 0

IF keyword_set(outplot) THEN device, /close ELSE stop


;------------ Effective Mass Loading*outflow metallicity/ISM
;             metallicity vs. stellar mass -- Ejected at different redshifts
xarr = smass ;xmass
xarrt = smasst
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_netmass_z_smass_ejectred.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,(reejectmassmet - diskgmassmet)/sfmassr/metallicity_ISMr,ytitle = textoidl('Net mass loss/M_{Stellar, formed}*Z_{wind}/Z_{ISM}'),xtitle = textoidl('M_{stellar} [M_{\odot}]'),/nodata,/ylog,xrange = [6e5,2e10],/xlog,xstyle = 1,yrange = [0.1,2e2]
distinct_colors,n_colors = 12
FOR iz = 0, nz - 1 DO BEGIN
    oplot,xarrt[*,iz],(reejectmassmetr[*,iz] - diskgmassmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(symbols[1]),color =  z_colors[iz],symsize = symsize
    oplot,xarrt[*,iz],(reejectmassmetr[*,iz] - diskgmassmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
ENDFOR
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Effective Mass Loading*outflow metallicity/ISM
;             metallicity vs. stellar mass -- Expelled at different redshifts
xarr = smass ;xmass
xarrt = smasst
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_netmass_z_smass_expellred.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,(reexpellmassmet - diskgmass_expellmet)/sfmassr/metallicity_ISMr,ytitle = textoidl('Net mass loss/M_{Stellar, formed}*Z_{wind}/Z_{ISM}'),xtitle = textoidl('M_{stellar} [M_{\odot}]'),/nodata,/ylog,xrange = [1e6,2e10],/xlog,xstyle = 1,yrange = [0.1,100]
distinct_colors,n_colors = 12
FOR iz = 0, nz - 1 DO BEGIN
    oplot,xarrt[*,iz],(reexpellmassmetr[*,iz] - diskgmass_expellmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(symbols[0]),color =  z_colors[iz],symsize = symsize
    oplot,xarrt[*,iz],(reexpellmassmetr[*,iz] - diskgmass_expellmetr[*,iz])/sfmassr[*,iz]/metallicity_ISMr[*,iz],psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
ENDFOR
IF keyword_set(outplot) THEN device, /close ELSE stop


;------------ Total massloading vs stellar mass
IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_smass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = xsize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = xsize
plot,smass,reexpellmassmet[*,nz - 1]/sfmass_total/yield,yrange = [0.001/yield,0.05/yield],ytitle = textoidl('Effective metal mass loading factor/'+yield_str),xtitle = textoidl('M_*/M')+sunsymbol(),/nodata,/ylog,xrange = [1e7,1e11],/xlog,POSITION=ASPECT(1, Margin=0.10)
distinct_colors,n_colors = 12
oplot,smass,reexpellmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,smass,reexpellmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,smass,reejectmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,smass,reejectmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,smass,relostmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,smass,relostmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
legend,['Metals Expelled','Metals Ejected','Metals Removed from Disk'],psym = [symbols[0],symbols[1],symcat(symbols[4])],color = [colore[0],colore[1],colore[4]],/bottom,/left,box=0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Metal Mass Loading Over Time
xarr = vcirc ;xmass
xarrt = vcirct
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_mass_time.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarr,reejectmassmetr[*,0]/sfmassr[*,0]/yield,psym = symcat(symbols[0]),/xlog,xrange = [10,200],yrange = [0.001/yield,0.1/yield],ytitle = textoidl('\eta_{Z,ejected}/')+yield_str,xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog
    distinct_colors,n_colors = 12
    FOR iz = 0, nz - 1 DO oplot,xarrt[*,iz],reejectmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(ej_psym[iz]),color  = z_colors[iz]
    FOR iz = 0, nz - 1 DO oplot,xarrt[*,iz],reejectmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(sym_outline(ej_psym[iz]),thick = thicks[0])
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;------------ Mass Loading Over Time with current xaxis
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_mass_time_eject.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarrt[*,0],reejectmassmetr[*,0]/sfmassr[*,0]/yield,/xlog,yrange = [0.001/yield,0.1/yield],ytitle = textoidl('\eta_{Z,ejected}/') + yield_str,xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [10,150]
    distinct_colors,n_colors = 12
    vc_arr = findgen(200)
    eta_m15 = vc_arr
    eta_m15[where(vc_arr LT 60)] = 2.9*(vc_arr[where(vc_arr LT 60)]/60)^(-3.2)
    eta_m15[where(vc_arr GE 60)] = 2.9*(vc_arr[where(vc_arr GE 60)]/60)^(-1.0)
;    loadct,0
;    oplot,vc_arr,eta_m15,thick = 2,color = 100,linestyle = 5
;    oplot,vc_arr,4500*vc_arr^(-2.0),thick = 2,color = 100,linestyle = 2 ;Energy driven
;    oplot,vc_arr,75*vc_arr^(-1.0),thick = 2,color = 100,linestyle = 3 ;Momentum driven
;    legend,['Energy-driven','Momentum-driven','Muratov et al. 2015','Fit to data'],linestyle = [2,3,5,0],thick = [2,2,2,2],color = [100,100,100,fgcolor],box = 0,/bottom,/left
;    loadct,39
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xarrt[*,iz],reejectmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(ej_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,xarrt[*,iz],reejectmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(sym_outline(ej_psym[iz]),thick = thicks[0]),symsize = symsize*0.6
;       fits_eject = robust_linefit( alog10(xarrt[*,iz]), alog10(reejectmassr[*,iz]/sfmassr[*,iz]), reejectmass_fit, sigma_eject )
;       print,fits_eject
;       oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),color = z_colors[iz],linestyle = 0,psym = -3,thick = thicks[0]
    ENDFOR
    fits_eject = robust_linefit( alog10(xarrt), alog10(reejectmassmetr/sfmassr/yield), reejectmass_fit, sigma_eject )
    fits_eject_z0 = robust_linefit( alog10(xarrt[*,3]), alog10(reejectmassmetr[*,3]/sfmassr[*,3]/yield), reejectmass_fit, sigma_eject )
    fits_eject_z0_5 = robust_linefit( alog10(xarrt[*,2]), alog10(reejectmassmetr[*,2]/sfmassr[*,2]/yield), reejectmass_fit, sigma_eject )
    fits_eject_z1 = robust_linefit( alog10(xarrt[*,1]), alog10(reejectmassmetr[*,1]/sfmassr[*,1]/yield), reejectmass_fit, sigma_eject )
    fits_eject_z2 = robust_linefit( alog10(xarrt[*,0]), alog10(reejectmassmetr[*,0]/sfmassr[*,0]/yield), reejectmass_fit, sigma_eject )
    print,'Eject: ',fits_eject
    print,'Eject: (z = 0): ',fits_eject_z0
    print,'Eject: (z = 0.5): ',fits_eject_z0_5
    print,'Eject: (z = 1)',fits_eject_z1
    print,'Eject: (z = 2)',fits_eject_z2
    oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),linestyle = 0,psym = -3,thick = thicks[0]
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;----------- Metal Mass Loading Over Time vs stellar mass            --------------------------
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_smass_time_eject.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = xsize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = xsize
    plot,smasst[*,0],reejectmassmetr[*,0]/sfmassr[*,0]/yield,/xlog,yrange = [0.001/yield,0.1/yield],ytitle = textoidl('\eta_{Z,ejected}/') + yield_str,xtitle = textoidl('M_*/M')+sunsymbol(),/nodata,/ylog,xrange = [6e5,1e10],POSITION=ASPECT(1, Margin=0.10)
    distinct_colors,n_colors = 12
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,smasst[*,iz],reejectmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(ej_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,smasst[*,iz],reejectmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(sym_outline(ej_psym[iz]),thick = thicks[0]),symsize = symsize*0.6
    ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/left,/bottom,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF


;------------ Mass Loading Over Time with current xaxis
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_mass_time_expell.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarrt[*,0],reexpellmassmetr[*,0]/sfmassr[*,0]/yield,/xlog,yrange = [0.001/yield,0.1/yield],ytitle = textoidl('\eta_{Z,expelled}/') + yield_str,xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [10,200]
    distinct_colors,n_colors = 12
 ;   oplot,xarr,reexpellmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[0])),symsize = symsize
;    oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),linestyle = 0,psym = -3,thick = 1
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xarrt[*,iz],reexpellmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,xarrt[*,iz],reexpellmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(sym_outline(ex_psym[iz]),thick = thicks[0]),symsize = symsize*0.6
;       fits_expell = robust_linefit( alog10(xarrt[*,iz]), alog10(reexpellmassr[*,iz]/sfmassr[*,iz]), reexpellmass_fit, sigma_expell )
;       print,fits_expell
;       oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),color = z_colors[iz],linestyle = 0,psym = -3,thick = thicks[0]
    ENDFOR
    fits_expell = robust_linefit( alog10(xarrt), alog10(reexpellmassmetr/sfmassr/yield), reexpellmass_fit, sigma_expell )
    print,'Reexpell: ',fits_expell
    oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),linestyle = 0,psym = -3,thick = thicks[0]
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ex_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;----------- Multi-plot mass loading
IF keyword_set(z_cut) THEN BEGIN
   IF keyword_set(outplot) THEN  device,filename = outplot + '_zmassloading_mass_time_multi.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*1.6,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize*1.6
   multiplot,[1,2]

   plot,xarrt[*,0],reejectmassmetr[*,0]/sfmassr[*,0]/yield_str,/xlog,yrange = [0.001/yield,0.1/yield],ytitle = textoidl('\eta_{Z,ejected}/') + yield_str,/nodata,/ylog,xrange = [10,200]
   distinct_colors,n_colors = 12
   FOR iz = 0, nz - 1 DO BEGIN
      oplot,xarrt[*,iz],reejectmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(ej_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
      oplot,xarrt[*,iz],reejectmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(sym_outline(ej_psym[iz]),thick = thicks[0]),symsize = symsize*0.6
   ENDFOR
   fits_eject = robust_linefit( alog10(xarrt), alog10(reejectmassmetr/sfmassr/yield), reejectmass_fit, sigma_eject )
   oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),linestyle = 0,psym = -3,thick = thicks[0]
   legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0

   multiplot
   plot,xarrt[*,0],reexpellmassmetr[*,0]/sfmassr[*,0]/yield,/xlog,yrange = [0.001/yield,0.1/yield],ytitle = textoidl('\eta_{Z,expelled}/') + yield_str,xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [10,200]
   FOR iz = 0, nz - 1 DO BEGIN
      oplot,xarrt[*,iz],reexpellmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
      oplot,xarrt[*,iz],reexpellmassmetr[*,iz]/sfmassr[*,iz]/yield,psym = symcat(sym_outline(ex_psym[iz]),thick = thicks[0]),symsize = symsize*0.6
   ENDFOR
   fits_expell = robust_linefit( alog10(xarrt), alog10(reexpellmassmetr/sfmassr/yield), reexpellmass_fit, sigma_expell )
   oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),linestyle = 0,psym = -3,thick = thicks[0]
   multiplot,/reset
   IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

;------------ Mass expelled/ejected, fraction of disk gas expelled/ejected, mass loading ----------
IF 1 THEN BEGIN
IF keyword_set(outplot) THEN  device,filename = outplot + '_multi_massz.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*2.6,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize*2.6
multiplot,[1,3]
plot,xmass,reexpellmassmet[*,nz - 1],/ylog,/xlog,xrange = xrange,yrange = [1e4,2e9],ytitle = 'Metal Mass [M' + sunsymbol() + ']',/nodata
distinct_colors,n_colors = 12
oplot,xmass,reexpellmassmet[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xmass, reejectmassmet[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, reejectmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xmass, relostmassmet[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize; Metals lost
oplot,xmass, relostmassmet[*,nz - 1],psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),color = fgcolor,symsize = symsize
oplot,xmass,(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(17),color = fgcolor,symsize = symsize
;oplot,xmass,(primehalo_stars_zmass + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[4])),color = fgcolor,symsize = symsize
oplot,xmass, primehalo_stars_zmass, psym = symcat(symbol_s),color = color_s,symsize = symsize ;Metals produced by the main progenitor (consider using currenthalo_stars_zmass if I want metals produced by stars in main halo at z = 0)
oplot,xmass, primehalo_stars_zmass, psym = symcat(sym_outline(symbol_s),thick = thicks[0]),color = fgcolor,symsize = symsize
legend,['Metals expelled','Metals ejected','Metals removed','Metals produced','Metals produced + accreted'],psym = [symbols[0:1],symbols[4],symbol_s,17],color = [colore[0:1],color_s,colore[4],fgcolor],/bottom,/right,box = 0
multiplot

plot,xmass,reexpellmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange = [0.001/yield,0.05/yield],ytitle = textoidl('Effective metal mass loading factor/') + yield_str,/nodata,/ylog
distinct_colors,n_colors = 12
oplot,xmass,reexpellmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
oplot,xmass,reejectmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass,reejectmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xmass,relostmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass,relostmassmet[*,nz - 1]/sfmass_total/yield,psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
legend,['Gas removed from disk'],psym = [symbols[4]],color = [colore[4]],symsize = symsize,/left,/bottom,box = 0
;legend,['Gas Mass Expelled/Stellar Mass Formed','Gas Mass Ejected/Stellar Mass Formed'],psym = symbols,color = colore,/bottom,/left,box=0
multiplot

;XXX Need total metals in the disk
plot,xmass,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange = [0,1],ytitle = textoidl('M_{outflow}/M_{accreted}'),xtitle = xtitle,/nodata,symsize = symsize
distinct_colors,n_colors = 12
oplot,xmass,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,(reexpellmassmet[*,nz - 1] - diskgmass_expellmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[0]),thick = thicks[0]),symsize = symsize
;oplot,xmass,reexpellmass_cool[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[0])),color = colore[0],symsize = symsize
oplot,xmass, (reejectmassmet[*,nz - 1]- diskgmassmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, (reejectmassmet[*,nz - 1]- diskgmassmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[1]),thick = thicks[0]),symsize = symsize
oplot,xmass, (relostmassmet[*,nz - 1]- diskgmassmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass, (relostmassmet[*,nz - 1]- diskgmassmet[*,nz - 1])/(zavail_zmax_total + diskgmassmet_uniq_total),psym = symcat(sym_outline(symbols[4]),thick = thicks[0]),symsize = symsize
;legend,['Gas Mass Expelled','Gas Mass Ejected'],psym = [symbols],color = [colore],box = 0,/bottom,/right;position = [2.2e11,0.4];/bottom,/right
;legend,['Total gas mass expelled','Gas mass expelled -- SN heated','Total gas mass ejected','Gas mass ejected -- SN heated'],psym=[symbols[0],sym_outline(symbols[0]),symbols[1],sym_outline(symbols[1])],color = [colore[0],colore[0],colore[1],colore[1]],/top,/right,box = 0
legend,['Total gas mass ejected','Gas mass ejected -- SN heated'],psym=[symbols[1],sym_outline(symbols[1])],color = [colore[1],colore[1]],/left,/bottom,box = 0
multiplot,/reset
IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

END
