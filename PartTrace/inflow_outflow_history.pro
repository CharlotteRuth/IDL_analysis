;inflow_outflow_history,['/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'],haloid=['1'],/color,/debug

PRO inflow_outflow_history, dirs, haloid = haloid, outplot = outplot, colors = colors, thicks = thicks, linestyle = linestyle, scale = scale,formatthick = formatthick, keys = keys,pmulti = pmulti,allpositive = allpositive,yrange_max = yrange_max,label = label,debug = debug,gmass = gmass

IF keyword_set(allpositive) THEN e_m = 1 ELSE e_m = -1

n = n_elements(dirs)
IF NOT keyword_set(haloid) THEN haloid = strarr(n_elements(dirs)) + '1'

IF keyword_set(outplot) THEN BEGIN 
    formatplot,/outplot,thick = formatthick
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
   formatplot
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF keyword_set(colors) THEN BEGIN
    loadct,39
    if colors[0] eq 1 then  colors = (findgen(n) + 1)*240/n ; else colors = color
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(outplot) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 1
    IF NOT keyword_set(linestyle) THEN linestyle = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(outplot) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 1
    IF NOT keyword_set(linestyle) THEN linestyle = findgen(n)*2   
ENDELSE
IF NOT keyword_set(scale) THEN scale = fltarr(n) + 1.0

ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
z = reverse(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
ticktime_in_z = findgen(7)*2*1e9 ;redshift major plot
tickred_in_z = spline(t,z,ticktime_in_z )

tickred_in_t = reverse(findgen(5))
ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e9

IF keyword_set(pmulti) THEN BEGIN
   IF (keyword_set(outplot)) THEN device,filename = outplot + '_fb.eps',/color,bits_per_pixel= 8,/times,xsize = 36,ysize = 45,xoffset =  2,yoffset =  2 ELSE window,3,xsize = 580,ysize = 800
   multiplot,pmulti,mxTitle = 'Time [Gyr]',myTitle = textoidl('M/M_{vir}'),gap = 0,mxTitSize = 1.5,myTitSize = 1.5,/square;,xgap = 0.01
   xtitle = ''
   ytitle = ''
   xstyle = 1
ENDIF ELSE BEGIN
;   multiplot,/reset
;   multiplot,/default
   xtitle = 'Time [Gyr]'
   ytitle = 'Mass [M'+ sunsymbol() + ']'
   xstyle = 9
ENDELSE
binsize = 5e8
nbins = round(14.0/(binsize/1e9))
FOR i = 0, n - 1 DO BEGIN 
   IF NOT keyword_set(pmulti) THEN $
      IF (keyword_set(outplot)) THEN device,filename = outplot + '_' + strtrim(i,2) + '_fb.eps',/color,bits_per_pixel= 8,/times,ysize = 16,xsize = 18,xoffset =  2,yoffset =  2 ELSE window,3,xsize = 712,ysize = 600
;   print,dirs[i],haloid[i]
   filebase = 'grp' + haloid[i] + '.mass_metal_track.dat'
   cd,dirs[i]
   spawn,'ls ' + dirs[i] + '*.grp' + haloid[i] + '.haloid.dat',files_haloid
   readcol,files_haloid,files,halos,format='a,l',/silent
   readcol,dirs[i] + filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,ox,fe,mHI,mH2,mcoldg,/silent ;,mH2,H2frac,Pressure,r25,SFR
   scale[i] = 1.0/mtot[n_elements(mtot) - 1]
   nhalostat = n_elements(halo)
   spawn,'ls ' + dirs[i] + '/*_2merge.starlog.fits',files_sl
   IF file_test(files_sl[0]) THEN BEGIN
      slall = mrdfits(files_sl[0],1,/silent)
      print,files_sl[0]
   ENDIF ELSE BEGIN
      spawn,'ls ' + dirs[i] + '/*.starlog',files_sl
      slall = rstarlog(files_sl[0],/molecularH)
      print,files_sl[0]
   ENDELSE
   spawn,'ls ' + dirs[i] + 'h*param',pfile
   units = tipsyunits(pfile[0],/silent)
   IF max(time) LT 0.5 THEN time = time*units.timeunit/1e9
   spawn,'ls ' + dirs[i] + '*512/*cosmo**512',file
   IF keyword_set(debug) THEN rtipsy,file,h,g,d,s ELSE rtipsy,file,h,g,d,s,/justhead
   IF keyword_set(debug) THEN BEGIN
      spawn,'ls ' + dirs[i] + '*512/*512.amiga.grp',file_grp
      readarr,file_grp,h,grp,/ascii
      grp_star = grp[h.ngas + h.ndark:h.n - 1]
      grp_gas = grp[0:h.ngas - 1]
      spawn,'ls ' + dirs[i] + '*512/*512.iord',file_iord
      readarr,file_iord,  h,iord,/ascii,type = 'long'
      iord_star = iord[h.ngas + h.ndark:h.n - 1]
      iord_gas = iord[0:h.ngas - 1]
;      spawn,'ls ' + dirs[i] + '*512/*512.massform',file_massform
;      readarr,file_massform,  h,massform,/ascii
;      massform_star = massform[h.ngas + h.ndark:h.n - 1]
;   spawn,'ls ' + dirs[i] + '*512/*512.coolontime',file_coolon
;   readarr,file_coolon,h,coolon,part = 'gas',/ascii
      spawn,'ls ' + dirs[i] + '*512/*512.amiga.stat',file_stat
      stat = read_stat_struc_amiga(file_stat)
      main = where(stat.group EQ haloid[i])
      satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
      inhalo = fltarr(n_elements(grp))
      FOR j = 0, n_elements(satellites) - 1 DO BEGIN
         test = where(grp EQ stat[satellites[j]].group, ntest)
         IF ntest NE 0 THEN inhalo[test] = 1
      ENDFOR
      inhalo_star = inhalo[h.ngas + h.ndark:h.n - 1]
      inhalo_gas = inhalo[0:h.ngas - 1]
   ENDIF

   IF keyword_set(debug) THEN BEGIN
;   iords = mrdfits(dirs[i] + '/grp' + haloid[i] + '.allgas.iord.fits')
; ------------------- Inflow on to Halo -----------------------------
   print,'Read accretion files'
   inflowz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccr_z.fits') ;'accr_rvir_z';'.accrz.fits'
   inflowm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reaccr.fits') ;'mass_at_accr_rvir';'.mass_at_accr.fits'
   inflowi = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccr_iord.fits')
   inflowi_sf = inflowi[where(inflowm LE 0)]
   match,inflowi_sf,slall.iordergas,ind1,ind2
   inflowm[(where(inflowm LE 0))[ind1]] = slall[ind2].massform*units.massunit
   inflowt = z_to_t(inflowz)
   inflowmcum  = weighted_histogram(inflowt, weight =  inflowm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)  
   nbins = n_elements(inflowmcum)

; Need unique accretion onto halo here
   ind_inflow_uniq = uniq(inflowi,sort(inflowi))
   inflowz_uniq = inflowz[ind_inflow_uniq]
   inflowm_uniq = inflowm[ind_inflow_uniq]
   inflowi_uniq = inflowi[ind_inflow_uniq]
   inflowt_uniq = z_to_t(inflowz_uniq)    
   inflowmcum_uniq = weighted_histogram(inflowt_uniq, weight =  inflowm_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
ENDIF

; ------------------- Inflow Disk -----------------------------
   inflowdiskz_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdiskall_z.fits')
   inflowdiskm_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reaccrdiskall.fits')
;   inflowdiskm[where(inflowdiskm LT 0)] = 0
   inflowdiski_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdiskall_iord.fits')  
   inflowdiski_all_sf = inflowdiski_all[where(inflowdiskm_all LE 0)]
   match,inflowdiski_all_sf,slall.iordergas,ind1,ind2
   inflowdiskm_all[(where(inflowdiskm_all LE 0))[ind1]] = slall[ind2].massform*units.massunit
   inflowdiskt_all = z_to_t(inflowdiskz_all)
   inflowmdiskcum_all  = weighted_histogram(inflowdiskt_all, weight =  inflowdiskm_all,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)

IF keyword_set(debug) THEN BEGIN
;Unique accretion onto disk
   ind_inflowdisk_all_uniq = uniq(inflowdiski_all[sort(inflowdiski_all)])
   inflowdiskz_all_uniq = inflowdiskz_all[ind_inflowdisk_all_uniq]
   inflowdiskm_all_uniq = inflowdiskm_all[ind_inflowdisk_all_uniq]
   inflowdiski_all_uniq = inflowdiski_all[ind_inflowdisk_all_uniq]
   inflowdiskt_all_uniq = z_to_t(inflowdiskz_all_uniq)
   inflowmdiskcum_all_uniq = weighted_histogram(inflowdiskt_all_uniq, weight =  inflowdiskm_all_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
;-------------------- Inflow Disk (reaccretion of all gas that leaves
;                     disk by being heated -------
   inflowdiskz_heat = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdiskheat_z.fits')
   inflowdiski_heat = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdiskheat_iord.fits')
   inflowdiskm_heat = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reaccrdiskheat.fits')
;   inflowdiskm_heat[where(inflowdiskm_heat LT 0)] = 0
   inflowdiski_heat_sf = inflowdiski_heat[where(inflowdiskm_heat LE 0)]
   match,inflowdiski_heat_sf,slall.iordergas,ind1,ind2
   inflowdiskm_heat[(where(inflowdiskm_heat LE 0))[ind1]] = slall[ind2].massform*units.massunit
   inflowdiskt_heat = z_to_t(inflowdiskz_heat)
   inflowmdiskcum_heat  = weighted_histogram(inflowdiskt_heat, weight =  inflowdiskm_heat,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)

;Unique cooling onto the disk
   ind_inflowdisk_heat_uniq = uniq(inflowdiski_heat,sort(inflowdiski_heat))
   inflowdiskz_heat_uniq = inflowdiskz_heat[ind_inflowdisk_heat_uniq]
   inflowdiskm_heat_uniq = inflowdiskm_heat[ind_inflowdisk_heat_uniq]
   inflowdiski_heat_uniq = inflowdiski_heat[ind_inflowdisk_heat_uniq]
   inflowdiskt_heat_uniq = z_to_t(inflowdiskz_heat_uniq)
   inflowmdiskcum_heat_uniq = weighted_histogram(inflowdiskt_heat_uniq, weight =  inflowdiskm_heat_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   ENDIF
;-------------------- Inflow Disk (reaccretion of all gas that is ejected from disk by supernova
   inflowdiskz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdisk_z.fits',/silent)
   inflowdiski = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdisk_iord.fits',/silent)
   inflowdiskm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reaccrdisk.fits',/silent)
;   inflowdiskm_heat[where(inflowdiskm_heat LT 0)] = 0
   inflowdiski_sf = inflowdiski[where(inflowdiskm LE 0)]
   match,inflowdiski_sf,slall.iordergas,ind1,ind2
   inflowdiskm[(where(inflowdiskm LE 0))[ind1]] = slall[ind2].massform*units.massunit
   inflowdiskt = z_to_t(inflowdiskz)
   inflowmdiskcum  = weighted_histogram(inflowdiskt, weight =  inflowdiskm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)

;Unique cooling onto the disk
   ind_inflowdisk_uniq = uniq(inflowdiski,sort(inflowdiski))
   inflowdiskz_uniq = inflowdiskz[ind_inflowdisk_uniq]
   inflowdiskm_uniq = inflowdiskm[ind_inflowdisk_uniq]
   inflowdiski_uniq = inflowdiski[ind_inflowdisk_uniq]
   inflowdiskt_uniq = z_to_t(inflowdiskz_uniq)
   inflowmdiskcum_uniq = weighted_histogram(inflowdiskt_uniq, weight =  inflowdiskm_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
;------------------ Expelled Material -------------------------
;Expelled is expelled because of supernova feedback
;    print,'Read files of gas expullsion'
;    expellm_all  = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_expell.fits',0)
;    expellz_all  = mrdfits(dirs[i] + '/grp' + haloid[i] + '.expell_z.fits',0)
;    expellm      = expellm_all[    where(expellz_all  NE 99.0)]
;    expellz      = expellz_all[    where(expellz_all  NE 99.0)] 
;    expelli      = iords      [    where(expellz_all  NE 99.0)]
;    expellt      = z_to_t(expellz) 
;    expellmcum  = weighted_histogram(expellt, weight =  expellm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)

;------------------ Expelled, Ejected and Heated Material -----------------
   reejecti = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reeject_iord.fits',0,/silent)
   reejectm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reeject.fits',0,/silent)
;   reejectm[where(reejectm LT 0)] = 0
   reejectz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reeject_z.fits',0,/silent)
   reejectt = z_to_t(reejectz)
;    reejecti_uniq = reejecti[uniq(reejecti,sort(reejecti))]
   reejectmcum  = weighted_histogram(reejectt, weight =  reejectm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
    
   reexpelli = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reexpell_iord.fits',0,/silent)
   reexpellm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reexpell.fits',0,/silent) 
;   reexpellm[where(reexpellm LT 0)] = 0
   reexpellz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reexpell_z.fits',0,/silent)
   reexpellt = z_to_t(reexpellz)
;    reexpelli_uniq = reexpelli[uniq(reexpelli,sort(reexpelli))]
   reexpellmcum  = weighted_histogram(reexpellt, weight =  reexpellm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)

   IF keyword_set(debug) THEN BEGIN
      reheati = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reheat_iord.fits',0)
      reheatm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reheat.fits',0) 
;   reheatm[where(reheatm LT 0)] = 0
      reheatz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reheat_z.fits',0)
      reheatt = z_to_t(reheatz)
;    reheati_uniq = reheati[uniq(reheati,sort(reheati))]
      reheatmcum  = weighted_histogram(reheatt, weight =  reheatm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)

      ind_reheat_uniq = uniq(reheati,sort(reheati))
      reheatz_uniq = reheatz[ind_reheat_uniq]
      reheatm_uniq = reheatm[ind_reheat_uniq]
      reheati_uniq = reheati[ind_reheat_uniq]
      reheatt_uniq = z_to_t(reheatz_uniq)
      reheatmcum_uniq = weighted_histogram(reheatt_uniq, weight =  reheatm_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   ENDIF

      relosti = mrdfits(dirs[i] + '/grp' + haloid[i] + '.relost_iord.fits',0)
      relostm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_relost.fits',0) 
;   relostm[where(relostm LT 0)] = 0
      relostz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.relost_z.fits',0)
      relostt = z_to_t(relostz)
;    relosti_uniq = relosti[uniq(relosti,sort(relosti))]
      relostmcum  = weighted_histogram(relostt, weight =  relostm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)

IF keyword_set(debug) THEN BEGIN
      
;------------------ Outflow Material --------------------------
;Outflow gas is gas that is lost because of stripping or wierd Amiga
;halo definition
      print,'Read files of outflows'
      outflowm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reoutflow.fits')
;   outflowm[where(outflowm LT 0)] = 0
      outflowi = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reoutflow_iord.fits')  
      outflowz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reoutflow_z.fits')
      outflowt = z_to_t(outflowz) 
      outflowmcum = weighted_histogram(outflowt,weight =  outflowm, min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
      
; Need unique outflow halo here
      ind_outflow_uniq = uniq(outflowi,sort(outflowi))
      outflowz_uniq = outflowz[ind_outflow_uniq]
      outflowm_uniq = outflowm[ind_outflow_uniq]
      outflowi_uniq = outflowi[ind_outflow_uniq]
      outflowt_uniq = z_to_t(outflowz_uniq)    
      outflowmcum_uniq = weighted_histogram(outflowt_uniq, weight =  outflowm_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
   ENDIF

IF keyword_set(debug) THEN BEGIN
;------------------- Stars at z = 0 ----------------------
   startfile = files[n_elements(files) - 1]
   rtipsy,startfile,h,g,d,s,/justhead
   readarr,startfile + '.iord',h,startiord,/ascii,type = 'long'
   readarr,startfile + '.amiga.grp',h,startgrp,/ascii,type = 'long'
   stat = read_stat_struc_amiga(startfile+'.amiga.stat')
   main = where(stat.group EQ halos[n_elements(files) - 1])
   satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
   inhalo = fltarr(n_elements(startgrp))
   FOR j = 0, n_elements(satellites) - 1 DO BEGIN
      test = where(startgrp EQ stat[satellites[j]].group, ntest)
      IF ntest NE 0 THEN inhalo[test] = 1
   ENDFOR
   inhalo_star_start = inhalo[h.ngas + h.ndark:h.n - 1]
   startiord_star = startiord[h.ngas + h.ndark:h.n - 1]
   startiord_star = startiord_star[where(inhalo_star_start EQ 1)]
   match2,slall.iorderstar,iord_star[where(inhalo_star EQ 1)],ind1,ind2
   slgal = slall[where(ind1 NE -1)]
;   match2,slgal.iordergas,inflowi_uniq,ind3,ind4
;   accr_star_iord = slgal[where(ind3 EQ -1)].iorderstar
;   match2,accr_star_iord,startiord_star,ind5,ind6
;   accr_star_iord = accr_star_iord[where(ind5 EQ -1)]
;   insitu_star_iord = [slgal[where(ind3 NE -1)].iorderstar,startiord_star]
;   insitu_star_iord = insitu_star_iord[uniq(insitu_star_iord,sort(insitu_star_iord))]

;   print,total(slgal[where(ind3 EQ -1)].massform)*units.massunit,total(slgal[where(ind3 NE -1)].massform)*units.massunit,total(slgal.massform)*units.massunit,mstar[nhalostat - 1],format='(E,E,E,E)'
;   accrstarcum = weighted_histogram(slgal[where(ind3 EQ -1)].timeform*units.timeunit/1e9, weight = slgal[where(ind3 EQ -1)].massform*units.massunit,min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
;   insitustarcum = weighted_histogram(slgal[where(ind3 NE -1)].timeform*units.timeunit/1e9, weight = slgal[where(ind3 NE -1)].massform*units.massunit,min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
;   sgal_igasord = slgal[where(ind3 NE -1)].iordergas
  
; ------------------- Formed Stars -----------------------------
   sl = mrdfits(dirs[i] + 'starlog.cut.'+haloid[i] + '.fits',1)
   formstarm = sl.massform*units.massunit
   formstart = sl.timeform
   formstarcum = weighted_histogram(sl.timeform, weight = formstarm,min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
ENDIF

; ------------------- Inflow Stars ----------------------------- 
 ;  IF file_test(dirs[i] + '/grp' + haloid[i] + '.accrstars_z.fits') THEN BEGIN
      inflowstarz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accrstars_z.fits',/silent)
      inflowstarm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accrstars_m.fits',/silent)
      inflowstari = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accrstars_i.fits',/silent)
      inflowstari = inflowstari[where(inflowstarz NE 99)]
      inflowstarm = inflowstarm[where(inflowstarz NE 99)]
      inflowstarz = inflowstarz[where(inflowstarz NE 99)]
      inflowstart = z_to_t(inflowstarz)
      inflowmstarcum  = weighted_histogram(inflowstart, weight =  inflowstarm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)   
      track_stars = 1

IF keyword_set(debug) THEN BEGIN
      match2,slgal.iorderstar,inflowstari,ind12,ind13
      insitu_star_iord = slgal[where(ind12 EQ -1)].iorderstar
      insitustarcum = weighted_histogram(slgal[where(ind12 NE -1)].timeform*units.timeunit/1e9, weight = slgal[where(ind12 NE -1)].massform*units.massunit,min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
      ENDIF
 ;  ENDIF ELSE BEGIN
 ;     track_stars = 0
;       mstarform = spline(timearr,formstarcum,time)
 ;     inflowmstarcum = fltarr(n_elements(timearr))
 ;  ENDELSE


   IF keyword_set(debug) THEN BEGIN
;      plot,time,mstar,yrange = [0,2.3e8]
;      oplot,timearr,formstarcum,color = 100
;      oplot,timearr,accrstarcum,color = 254
;      oplot,timearr,insitustarcum,color = 120
;      oplot,timearr,(insitustarcum + accrstarcum),color = 190
;      print,((insitustarcum[nbins - 1] + accrstarcum[nbins - 1]) - mstar[nhalostat - 1])/mstar[nhalostat - 1]

;Check that all gas currently in the halo was accreted
      match2,inflowi_uniq,iord_gas[where(inhalo_gas EQ 1)],ind5,ind6
      help,where(ind6 EQ -1) ;Should be the empty array

      ind_outflow_uniq = uniq(outflowi,sort(outflowi))
      outflowz_uniq = outflowz[ind_outflow_uniq]
      outflowm_uniq = outflowm[ind_outflow_uniq]
      outflowi_uniq = outflowi[ind_outflow_uniq]
      outflowt_uniq = z_to_t(outflowz_uniq)
      match,outflowi_uniq,inflowi_uniq,ind7,ind8
      timediff = outflowt_uniq[ind7] - inflowt_uniq[ind8]
      print,minmax(timediff)
      
      outhalo_z0i = outflowi_uniq[ind7(where(timediff GT 0))]
      match,outhalo_z0i,iord_gas,ind9,ind10
      help,where(inhalo_gas[ind10] EQ 1) ;Should be the empty array

      print,total(inflowm_uniq[ind8(where(timediff GT 0))] - outflowm_uniq[ind7(where(timediff GT 0))])

      match2,outflowi_uniq,inflowi_uniq,ind11,ind12
      inhalo_z0i = [inflowi_uniq[ind8(where(timediff LT 0))],inflowi_uniq[where(ind12 EQ -1)]]
      inhalo_z0m = [inflowm_uniq[ind8(where(timediff LT 0))],inflowm_uniq[where(ind12 EQ -1)]]
   ENDIF

;    fbase = strmid(files[i],0,strlen(files[i]) - 7)
;    rtipsy,fbase,h,g,d,s,/justhead
;    readarr,fbase + '.iord',h,iord,/ascii,part = 'gas'
;    readarr,fbase + '.amiga.grp',h,grp,/ascii,part = 'gas'
;    match,iord,outflowiord,iordind,ofiordind
;    iord0 = iord[iordind];Particles once accreted that still exist
;    grp0 = grp[iordind]
;    outflowiord0 = outflowiord[ofiordind]
;    outflowm0 =    outflowm[ofiordind]
;    outflowt0 =    outflowt[ofiordind]
;    outflowm00 =   outflowm0[where(grp0 NE 1)]
;    outflowt00 =   outflowt0[where(grp0 NE 1)]

   IF 0 THEN BEGIN
;-------------------- Bin Outflow and Inflows according to the same
;                     steps in the data file
      a_expellm  = fltarr(n_elements(mtot))
      a_expellm_cum  = fltarr(n_elements(mtot))
      a_outflowm = fltarr(n_elements(mtot))
      a_outflowm_cum = fltarr(n_elements(mtot))
      a_inflowm  = fltarr(n_elements(mtot))
      a_inflowm_cum  = fltarr(n_elements(mtot))
      a_time = [0,time]
      
      FOR j = 0, n_elements(mtot) - 1 DO BEGIN 
         t0 = a_time[j] 
         t1 = a_time[j + 1]  
         IF (where(expellt  GT t0 and expellt  LE t1))[0] NE -1 THEN a_expellm[j] = total(expellm[where(expellt GT t0 and expellt LE t1)]) 
         IF (where(expellt  LE t1))[0] NE -1 THEN a_expellm_cum[j] = total(expellm[where(expellt LE t1)]) 
         IF (where(outflowt GT t0 and outflowt LE t1))[0] NE -1 THEN a_outflowm[j] = total(outflowm[where(outflowt GT t0 and outflowt LE t1)]) 
         IF (where(outflowt LE t1))[0] NE -1 THEN a_outflowm_cum[j] = total(outflowm[where(outflowt LE t1)]) 
         IF (where(inflowt  GT t0 and inflowt  LE t1))[0] NE -1 THEN a_inflowm[j] = total(inflowm[where(inflowt GT t0 and inflowt LE t1)]) 
         IF (where(inflowt  LE t1))[0] NE -1 THEN a_inflowm_cum[j] = total(inflowm[where(inflowt LE t1)]) 
      ENDFOR
      IF i eq 0 THEN $
         a_array0 = [[time],[a_inflowm],[a_inflowm_cum],[a_outflowm],[a_outflowm_cum],[a_expellm],[a_expellm_cum]] ELSE $
            a_array1 = [[time],[a_inflowm],[a_inflowm_cum],[a_outflowm],[a_outflowm_cum],[a_expellm],[a_expellm_cum]]
   ENDIF
   
   IF keyword_set(yrange_max) THEN yrange = [0,yrange_max[i]] ELSE $
      IF keyword_set(allpositive) THEN yrange = [0,[max([mgas+mstar,(inflowmdiskcum + inflowmstarcum),inflowmdiskcum_all])]]*scale[i] ELSE  yrange = [-1*max([outflowmcum,relostmcum]),max([inflowmcum,inflowmdiskcum_all])]*scale[i]
   print,i,mtot[n_elements(mtot) - 1],yrange[1],yrange[1]/mtot[n_elements(mtot) - 1]
;   IF i EQ 0 THEN BEGIN
   IF keyword_set(debug) THEN stop
   plot,time,(mgas+mstar)*scale[i],ytitle = ytitle,xstyle =xstyle,xrange=[0.1,ageUniverse/1e9],xtitle = xtitle,/nodata,xmargin = [17,17],ymargin = [6,6],yrange = yrange,ystyle = 1,xcharsize = 1,ycharsize = 1
;        axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
   IF NOT keyword_set(pmulti) THEN axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = n_elements(ticktime_in_t) - 1
   IF NOT keyword_set(allpositive) THEN oplot,[0,14],[0,0],thick = 2
;   ENDIF

;----------------- Mass within halo --------------
   oplot,time,(mgas+mstar)*scale[i],                            linestyle = 5,thick = thicks[i];,color = colors[i] ;total baryons within rvir
   IF keyword_set(debug) THEN BEGIN
      oplot,timearr,inflowmcum*scale[i],                linestyle = 0,thick = thicks[i],color = 185;colors[i] ;,psym = -4 ;gas mass accreted into the virial mass (reaccr)
      oplot,timearr,inflowmcum_uniq*scale[i],                linestyle = 1,thick = thicks[i],color = 185;colors[i] ;gas mass accreted into the virial mass -- once (reaccr)
      oplot,timearr,e_m*outflowmcum*scale[i],          linestyle = 0,thick = thicks[i],color = 185;colors[i] ;gas mass out (reoutflow_z.fits, mass_at_reoutflow.fits)
      oplot,timearr,(inflowmcum - outflowmcum + inflowmstarcum)*scale[i],linestyle = 0,thick = thicks[i],color = 185;colors[i] ;gas accreted minus gas out
   ENDIF

;==================================
;stop
;   IF i EQ 0 THEN BEGIN
;   yrange = [-1*max(relostmcum*scale[i]),max(inflowmdiskcum_all*scale[i])]
;      plot,time,mgas+mstar,ytitle = 'Baryonic Mass [M' + sunsymbol() + ']',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',/nodata,xmargin = [17,17],ymargin = [6,6],yrange = yrange
;        axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
;      axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = n_elements(ticktime_in_t) - 1        
;      oplot,[0,14],[0,0],thick = 2
;   ENDIF

;--------------------- Disk --------------
   oplot,timearr,(inflowmdiskcum_all + inflowmstarcum)*scale[i],      linestyle = 0,thick = thicks[i],color = fgcolor
   oplot,timearr,e_m*relostmcum*scale[i],           linestyle = 0,thick = thicks[i],color = 110 ;colors[i] + 20
   IF keyword_set(debug) THEN BEGIN
      ghalo = g[where(inhalo_gas NE 0)]
      IF keyword_set(gmass) THEN mdisk = gmass[i] ELSE mdisk = total(g[where(ghalo.dens*units.rhounit GE 0.1 AND ghalo.tempg LT 1.2e4)].mass)*units.massunit
 ;     oplot,time,(mgas+mstar)*scale[i],               linestyle = 5,thick = thicks[i] ;,color = colors[i] ;total baryons within rvir
      oplot,time,(mcoldg + mstar)*scale[i],           linestyle = 5,thick = thicks[i] ;,color = colors[i] + 20 ;baryonic mass
      oplot,[max(time),max(time)],[mdisk + mstar[n_elements(mstar) -1],mdisk + mstar[n_elements(mstar) -1]]*scale[i],psym = 7,symsize = 2
      oplot,timearr,inflowmdiskcum_all*scale[i],      linestyle = 0,thick = thicks[i],color = 225;colors[i] + 20
;      oplot,timearr,inflowmdiskcum_all_uniq*scale[i], linestyle = 1,thick = thicks[i],color = 225;colors[i] + 20
      oplot,timearr,(inflowmdiskcum_all  - relostmcum + inflowmstarcum)*scale[i],          linestyle = 2,thick = thicks[i],color = 225;colors[i] + 20
      oplot,timearr,inflowmdiskcum_heat*scale[i],     linestyle = 0,thick = thicks[i],color = 210;colors[i] + 50
;      oplot,timearr,inflowmdiskcum_heat_uniq*scale[i],linestyle = 1,thick = thicks[i],color = 210;colors[i] + 50
      oplot,timearr,-1*reheatmcum*scale[i],           linestyle = 0,thick = thicks[i],color = 210;colors[i] + 50
      oplot,timearr,(inflowmdiskcum_heat - reheatmcum + inflowmstarcum)*scale[i],linestyle = 2,thick = thicks[i],color = 210;colors[i] + 50
    ENDIF

;--------------------- Ejecta and Expelled --------------
 ;   oplot,timearr,(inflowmdiskcum + inflowmstarcum)*scale[i],            linestyle = 0,thick = thicks[i],color = 60;254 ;gas mass into disk 'accrdisk_z.fits'
;    oplot,timearr,(inflowmdiskcum_uniq + inflowmstarcum)*scale[i],       linestyle = 1,thick = thicks[i],color = 60;254 ;gas mass into disk -- once 'accrdisk_z.fits'
    oplot,timearr,e_m*reejectmcum*scale[i],          linestyle = 0,thick = thicks[i],color = 254 ;gas mass ejected form disk (multiple times) 'reeject_z'
    oplot,timearr,e_m*reexpellmcum*scale[i],         linestyle = 0,thick = thicks[i],color = 50 ;gas mass expelled from disk by supernova 'reexpell_z'
    IF keyword_set(debug) THEN $
       oplot,timearr,(inflowmdiskcum - reejectmcum + inflowmstarcum)*scale[i],linestyle = 5,thick = thicks[i],color = 60;colors[i] + 100 ;gas + stars in - gas out

;------------------ Stars ----------------------
;Note that there will be a difference between the current mass of stars and the mass of stars formed
    IF keyword_set(debug) THEN BEGIN
       oplot,time,mstar*scale[i],                               linestyle = 0,thick = thicks[i],color = 100
       oplot,timearr,(insitustarcum  + inflowmstarcum)*scale[i],linestyle = 2,thick = thicks[i],color = 100;colors[i]+150
;    IF track_stars THEN oplot,timearr,(inflowmdiskcum + inflowmstarcum - reejectmcum)*scale[i],linestyle = 0,thick = thicks[i],color = colors[i]+50 ;gas + stars in to disk - gas out of disk (need to add something that will show multiple accretions onto disk)
    ENDIF
       
;Baryonic mass will be different from inflow - outflow because of gas that is reaccreted later
;    oplot,time,(mcoldg + mstar)*scale[i],                linestyle = 1,thick = thicks[i],color = colors[i]+100 ;baryonic disk mass. Possibly different because of stripping, lack of rediskaccretion, and different definitions of cold gas?
;    oplot,time,mgas*scale[i],linestyle = 2,thicks = 2,color = colors[i]
;    oplot,time,mstar*scale[i],linestyle = 2,thicks = 2,color = colors[i]
;    oplot,time,mstar*scale[i] + mgas*scale[i],linestyle = 0,thicks = 2,color = colors[i]
;    oplot,time,cuminflow*scale[i] - (mstar*scale[i] + mgas*scale[i]),linestyle = 1,thicks = 2,color = colors[i]
;    IF i EQ 0 THEN outflowratio0 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]] ELSE outflowratio1 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]]
    IF keyword_set(label) THEN legend,[label[i]],box = 0,/left,/top
    IF i EQ 0 THEN $
       IF keyword_set(keys) THEN legend,keys,color = [fgcolor,110,254,50,fgcolor],linestyle = [0,0,0,0,2],/top,/right,box = 0,position = [14,0.26];[0.293480,0.95];,position = [385,90]
;    print,'Position: ',position
    IF keyword_set(pmulti) THEN multiplot
    IF keyword_set(debug) AND NOT keyword_set(outplot) THEN stop
ENDFOR
;legend,['Gas Accreted','Gas Accreted to Disk','Stars Accerted','Gas Lost','Gas Expelled','Gas Accreted - Gas Lost','Baryonic Mass'],linestyle = [2,5,4,2,3,0,1]
;legend,['Gas Accreted','Gas Accreted to Disk','Gas Lost','Gas Expelled','Gas Accreted - Gas Lost','Disk Gas Accreted - Disk Gas Lost'],linestyle = [2,5,2,3,0,0]
;legend,['Virial Mass','Gas Accreted to Disk','Gas Accreted Once to Disk','Ejected Gas','Expelled Gas'],linestyle = [5,0,1,0,0],color = [fgcolor,254,254,60,20]
multiplot,/reset
IF keyword_set(outplot) THEN device,/close ELSE stop
END
