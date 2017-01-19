PRO ejectz_v_mass,dirs,files,halo = halo,finalstep = finalstep,colors = colors,outplot = outplot,stellarmass = stellarmass,z_cut = z_cut,symbols = symbols,formatthick = formatthick,z_colors = z_colors,gmass = gmass,accr_star = accr_star,rewrite = rewrite

formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
    IF keyword_set(stellarmass) THEN outplot = outplot + '_sm_' ELSE outplot = outplot + '_vm'
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

diskgmass_total  = fltarr(n)
halogmass_total  = fltarr(n)
allgmass_total   = fltarr(n)
sfmass_total     = fltarr(n)
smass_accr_total = fltarr(n)

relostmass = fltarr(n,nz)
;reheatmass = fltarr(n,nz)
;reheatmassr = fltarr(n,nz)
;reejectmass_rv2  = fltarr(n,nz)
;reejectmassr_rv2  = fltarr(n,nz)
;reejectmass_rv5  = fltarr(n,nz)
;reejectmassr_rv5  = fltarr(n,nz)
reejectmass  = fltarr(n,nz)
reejectmassr  = fltarr(n,nz)
;reejectmass_cool  = fltarr(n,nz)
;reejectmass_coolr  = fltarr(n,nz)
reexpellmass = fltarr(n,nz)
reexpellmassr = fltarr(n,nz)
;reexpellmass_cool = fltarr(n,nz)
;reexpellmass_coolr = fltarr(n,nz)
diskgmass = fltarr(n,nz)
;diskgmass_h = fltarr(n,nz)
;diskgmass_rv2 = fltarr(n,nz)
;diskgmass_rv5 = fltarr(n,nz)
diskgmass_lost = fltarr(n,nz)
halogmass = fltarr(n,nz)
sfmassr = fltarr(n,nz)
reoutflowmass = fltarr(n,nz)
reoutflowmassr = fltarr(n,nz)
reoutflowmass_uniq = fltarr(n,nz)

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
   FOR iz = 0, n_elements(z_bins) - 1 DO BEGIN
      temp = min(abs(time_t - t_bins[iz]),ind_t)
      vmasst[i,iz] = mtot_t[ind_t]
      smasst[i,iz] = mstar_t[ind_t] 
      temp = min(abs(align.z - z_bins[iz]),ind_t)
      stat = read_stat_struc_amiga(dirs[i] + '/' + align[ind_t].file + '.amiga.stat')
      ind = (where(stat.group EQ align[ind_t].haloid))[0]
;      vcirct[i,iz] = stat[ind].vc
      vcirct[i,iz] = sqrt(6.67e-11*stat[ind].m_tot*1.9891e30/(3.08567758e19*stat[ind].rvir))/1000.
   ENDFOR
   IF keyword_set(rewrite) OR NOT file_test(dirs[i] + '/grp' + halo[i] + '.ejectz_quant.txt') THEN $
      ejectz_v_mass_data,dirs[i],files[i],halo[i],z_bins,zmax,mlbin
   openr,lun,dirs[i] + '/grp' + halo[i] + '.ejectz_quant.txt',/get_lun
   readf,lun,allgmass_total_temp,halogmass_total_temp,diskgmass_total_temp,relost_total_temp,reeject_total,reexpell_total_temp,sfmass_total_temp,smass_accr_total_temp
   close,lun
   free_lun,lun
   allgmass_total[i] = allgmass_total_temp
   halogmass_total[i] = halogmass_total_temp
   diskgmass_total[i] = diskgmass_total_temp
   sfmass_total[i] = sfmass_total_temp
   smass_accr_total[i] = smass_accr_total_temp
   readcol,dirs[i] + '/grp' + halo[i] + '.ejectz_quant.txt',z_bins_temp,relostmass_temp,reejectmass_temp,reejectmassr_temp,reexpellmass_temp,reexpellmassr_temp,diskgmass_lost_temp,diskgmass_temp,halogmass_temp,sfmassr_temp

   relostmass[i,*] = relostmass_temp
;   reheatmass[i,*] = reheatmass_temp
;   reheatmassr[i,*] = reheatmassr_temp
;   reejectmass_rv2[i,*] = reejectmass_rv2_temp
;   reejectmassr_rv2[i,*] = reejectmassr_rv2_temp
;   reejectmass_rv5[i,*] = reejectmass_rv5_temp
;   reejectmassr_rv5[i,*] = reejectmassr_rv5_temp
   reejectmass[i,*] = reejectmass_temp
   reejectmassr[i,*] = reejectmassr_temp
;   reejectmass_cool[i,*] = reejectmass_cool_temp
;   reejectmass_coolr[i,*] = reejectmass_coolr_temp
   reexpellmass[i,*] = reexpellmass_temp
   reexpellmassr[i,*] = reexpellmassr_temp
;   reexpellmass_cool[i,*] = reexpellmass_cool_temp
;   reexpellmass_coolr[i,*] = reexpellmass_coolr_temp
   diskgmass_lost[i,*] = diskgmass_lost_temp
   diskgmass[i,*] = diskgmass_temp
;   diskgmass_h[i,*] = diskgmass_h_temp
;   diskgmass_rv2[i,*] = diskgmass_rv2_temp
;   diskgmass_rv5[i,*] = diskgmass_rv5_temp
   halogmass[i,*] = halogmass_temp
   sfmassr[i,*] = sfmassr_temp

;For Rachel
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
;

;plot,xmass,((halogmass_total + smass_accr_total) - (smass + gmass_hot))/reexpellmass[*,3],psym = 4,/xlog
;plot,xmass,(halogmass_total + smass_accr_total) - (smass + gmass_hot),psym = 4,/xlog,/ylog
;oplot,xmass,reexpellmass[*,3],psym = 5

   plot_alllost = 1
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
;stop
ENDFOR

IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*254/n else colors = colors
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(symbols) THEN symbols = [14,15,17,18,16,34];[4,5]
    colore = [50,254,180,200,110,30]
    ;expelled, ejected, rvir/2,rvir/5,heated,lost
    IF keyword_set(z_cut) THEN BEGIN
        IF NOT keyword_set(z_colors) THEN z_colors = (findgen(nz) + 1)/n_elements(z_bins)*254
        IF NOT keyword_set(z_psym)   THEN z_psym   =  fltarr( nz) + 16;4
    ENDIF
    color_s = 150;fgcolor
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
    colore =  [50,254,210,30,120,30]
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

;------------ Expelled and Ejected Metal Mass as a Function of Halo Mass
IF keyword_set(outplot) THEN  device,filename = outplot + '_eject_expell_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,reexpellmass[*,nz - 1],psym = symcat(symbols[0]),/ylog,/xlog,xrange = xrange,yrange = [1e6,1e12],xtitle = xtitle,ytitle = 'Mass [M' + sunsymbol() + ']',symsize = symsize,/nodata
oplot,xmass,reexpellmass[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmass[*,nz - 1],psym = symcat(sym_outline(symbols[0])),color = fgcolor,symsize = symsize
;oplot,xmass,reexpellmass_cool[*,nz - 1],psym = symcat(sym_outline(symbols[0])),color = colore[0],symsize = symsize
oplot,xmass, reejectmass[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, reejectmass[*,nz - 1],psym = symcat(sym_outline(symbols[1])),color = fgcolor,symsize = symsize
;oplot,xmass, reejectmass_cool[*,nz - 1],psym = symcat(sym_outline(symbols[1])),color = colore[1],symsize = symsize
;oplot,xmass, reejectmass_rv2[*,nz - 1],psym = symcat(symbols[2]),color = colore[2],symsize = symsize
;oplot,xmass, reejectmass_rv2[*,nz - 1],psym = symcat(sym_outline(symbols[2])),color = fgcolor,symsize = symsize
oplot,xmass,smass,psym = symcat(symbol_s),color = color_s,symsize = symsize ;smass
;oplot,xmass,diskgmass[*,nz - 1],psym = symcat(symbol_a),color = color_a,symsize = symsize
oplot,xmass,smass,psym = symcat(sym_outline(symbol_s)),color = fgcolor,symsize = symsize
oplot,xmass,halogmass[*,nz - 1],psym = symcat(9),color = color_a,symsize = symsize
;,/ylog,/xlog,xrange = [1e8,1e11],yrange = [1e8,1e11],xtitle = xtitle,ytitle = 'Mass of Gas Ejected'
;oplot,xmass,reheatmass[*,nz - 1],psym = symcat(symbols(4)),color  = colore[4],symsize = symsize  
;oplot,xmass,reheatmass[*,nz - 1],psym = symcat(sym_outline(symbols(4))),symsize = symsize
;oplot,xmass,reoutflowmass[*,nz - 1],psym = symcat(14),color = fgcolor,symsize = symsize  
legend,['Gas mass expelled','Gas mass ejected','Stellar mass','Gas accreted to halo'],psym = [symbols[0:1],symbol_s,9],color = [colore[0:1],color_s,color_a],/bottom,/right,box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Expelled and Ejected Metal Frac as a Function of Halo Mass
readcol,'~/Datafiles/HIcubes/moster.stars.z0',mhalo,mstar_mhalo,mstar,logmstar

IF keyword_set(outplot) THEN  device,filename = outplot + '_eject_expell_fbar.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(halogmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(symbols[0]),/ylog,/xlog,xrange = xrange,xtitle = xtitle,yrange = [4e-3,1.3],ytitle = textoidl('M/(f_B M_{vir})'),symsize = symsize,/nodata
oplot,xrange,[1,1],linestyle = 1
oplot,xmass,(smass + gmass_hot)/(f_bar*vmass),psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,(smass + gmass_hot)/(f_bar*vmass),psym = symcat(sym_outline(symbols[0])),color = fgcolor,symsize = symsize
oplot,xmass,(halogmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(sym_outline(symbols[0])),color = colore[0],symsize = symsize
oplot,xmass,(smass + gmass)/(f_bar*vmass),psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass,(smass + gmass)/(f_bar*vmass),psym = symcat(sym_outline(symbols[1])),color = fgcolor,symsize = symsize
oplot,xmass,(diskgmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(sym_outline(symbols[1])),color = colore[1],symsize = symsize
oplot,xmass,smass/(f_bar*vmass),psym = symcat(symbol_s),color = color_s,symsize = symsize ;smass
oplot,xmass,smass/(f_bar*vmass),psym = symcat(sym_outline(symbol_s)),color = fgcolor,symsize = symsize
;oplot,xmass,diskgmass[*,nz - 1],psym = symcat(symbol_a),color = color_a,symsize = symsize
oplot,xmass,halogmass[*,nz - 1],psym = symcat(9),color = color_a,symsize = symsize
;,/ylog,/xlog,xrange = [1e8,1e11],yrange = [1e8,1e11],xtitle = xtitle,ytitle = 'Mass of Gas Ejected'
;oplot,10^mhalo,10^logmstar/(f_bar*10^mhalo) ;Plot Moster
legend,['Gas accreted to halo','Baryons in halo','Gas accreted to disk','Disk baryons','Stars'],psym = [sym_outline(symbols[0]),symbols[0],sym_outline(symbols[1]),symbols[1],symbol_s],color = [colore[0],colore[0],colore[1],colore[1],color_s],/bottom,/right,box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Supression as a Function of Halo Mass
IF keyword_set(outplot) THEN  device,filename = outplot + '_eject_expell_fsuppress.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(halogmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(symbols[0]),/xlog,xrange = xrange,xtitle = xtitle,yrange = [0,1.1],ytitle = textoidl('Unsuppressed Mass Fraction'),symsize = symsize,/nodata
oplot,xrange,[1,1],linestyle = 1
oplot,xmass,smass/(smass + gmass),psym = symcat(14),color = color_s,symsize = symsize
oplot,xmass,smass/(smass + gmass),psym = symcat(sym_outline(14)),color = fgcolor,symsize = symsize
oplot,xmass,(smass + gmass)/(diskgmass_total + smass_accr_total),psym = symcat(17),color = colore[1],symsize = symsize
oplot,xmass,(smass + gmass)/(diskgmass_total + smass_accr_total),psym = symcat(sym_outline(17)),color = fgcolor,symsize = symsize
;oplot,xmass,(smass + gmass)/(diskgmass_lost[*,nz - 1] + smass_accr_total),psym = symcat(sym_outline(17)),color = colore[1],symsize = symsize
oplot,xmass,(halogmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(sym_outline(symbols[0])),color = colore[0],symsize = symsize ;Fraction of avaliable gas accreted to halo
oplot,xmass,(diskgmass_total + smass_accr_total)/(halogmass_total + smass_accr_total),psym = 2,color = colore[0],symsize = symsize
legend,['Accreted to halo/'+textoidl('f_B M_{vir}'),'Accreted to disk/accreted to halo','Baryons in disk/unique accretions to disk','Stellar mass/disk mass'],psym = [sym_outline(symbols[0]),2,17,14],color = [colore[0],colore[0],colore[1],color_s],/top,/left,box = 0
;oplot,xmass,(smass + gmass_hot)/(f_bar*vmass),psym = symcat(symbols[0]),color = colore[0],symsize = symsize
;oplot,xmass,(smass + gmass_hot)/(f_bar*vmass),psym = symcat(sym_outline(symbols[0])),color = fgcolor,symsize = symsize
;oplot,xmass,(smass + gmass)/(f_bar*vmass),psym = symcat(symbols[1]),color = colore[1],symsize = symsize
;oplot,xmass,(smass + gmass)/(f_bar*vmass),psym = symcat(sym_outline(symbols[1])),color = fgcolor,symsize = symsize
;oplot,xmass,(diskgmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(sym_outline(symbols[1])),color = colore[1],symsize = symsize
;oplot,xmass,smass/(f_bar*vmass),psym = symcat(symbol_s),color = color_s,symsize = symsize ;smass
;oplot,xmass,smass/(f_bar*vmass),psym = symcat(sym_outline(symbol_s)),color = fgcolor,symsize = symsize
;oplot,xmass,diskgmass[*,nz - 1],psym = symcat(symbol_a),color = color_a,symsize = symsize
;oplot,xmass,halogmass[*,nz - 1],psym = symcat(9),color = color_a,symsize = symsize
;,/ylog,/xlog,xrange = [1e8,1e11],yrange = [1e8,1e11],xtitle = xtitle,ytitle = 'Mass of Gas Ejected'
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Multiplot with suppresion ------------
IF keyword_set(outplot) THEN  device,filename = outplot + '_eject_expell_fbar_multi.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*1.6,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = 1.6*ysize
multiplot,[1,2]

plot,xmass,(halogmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(symbols[0]),/ylog,/xlog,xrange = xrange,yrange = [4e-3,1.3],ytitle = textoidl('M/(f_B M_{vir})'),symsize = symsize,/nodata
oplot,xrange,[1,1],linestyle = 1
oplot,xmass,(smass + gmass_hot)/(f_bar*vmass),psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,(smass + gmass_hot)/(f_bar*vmass),psym = symcat(sym_outline(symbols[0])),color = fgcolor,symsize = symsize
oplot,xmass,(halogmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(sym_outline(symbols[0])),color = colore[0],symsize = symsize
oplot,xmass,(smass + gmass)/(f_bar*vmass),psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass,(smass + gmass)/(f_bar*vmass),psym = symcat(sym_outline(symbols[1])),color = fgcolor,symsize = symsize
oplot,xmass,(diskgmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(sym_outline(symbols[1])),color = colore[1],symsize = symsize
oplot,xmass,smass/(f_bar*vmass),psym = symcat(symbol_s),color = color_s,symsize = symsize ;smass
oplot,xmass,smass/(f_bar*vmass),psym = symcat(sym_outline(symbol_s)),color = fgcolor,symsize = symsize
oplot,xmass,halogmass[*,nz - 1],psym = symcat(9),color = color_a,symsize = symsize
legend,['Gas Accreted to Halo','Baryons in Halo','Gas Accreted to Disk','Disk Baryons','Stars'],psym = [sym_outline(symbols[0]),symbols[0],sym_outline(symbols[1]),symbols[1],symbol_s],color = [colore[0],colore[0],colore[1],colore[1],color_s],/bottom,/right,box = 0
multiplot

plot,xmass,(halogmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(symbols[0]),/xlog,xrange = xrange,xtitle = xtitle,yrange = [0,1.2],ytitle = textoidl('Unsuppressed Mass Fraction'),symsize = symsize,/nodata
oplot,xrange,[1,1],linestyle = 1
oplot,xmass,smass/(smass + gmass),psym = symcat(14),color = color_s,symsize = symsize
oplot,xmass,smass/(smass + gmass),psym = symcat(sym_outline(14)),color = fgcolor,symsize = symsize
oplot,xmass,(smass + gmass)/(diskgmass_total + smass_accr_total),psym = symcat(17),color = colore[1],symsize = symsize
oplot,xmass,(smass + gmass)/(diskgmass_total + smass_accr_total),psym = symcat(sym_outline(17)),color = fgcolor,symsize = symsize
;oplot,xmass,(smass + gmass)/(diskgmass_lost[*,nz - 1] + smass_accr_total),psym = symcat(sym_outline(17)),color = colore[1],symsize = symsize
oplot,xmass,(halogmass_total + smass_accr_total)/(f_bar*vmass),psym = symcat(sym_outline(symbols[0])),color = colore[0],symsize = symsize ;Fraction of avaliable gas accreted to halo
oplot,xmass,(diskgmass_total + smass_accr_total)/(halogmass_total + smass_accr_total),psym = 2,color = colore[0],symsize = symsize
legend,['Accreted to halo/'+textoidl('f_B M_{vir}'),'Accreted to disk/accreted to halo','Baryons in disk/unique accretions to disk','Stellar mass/disk mass'],psym = [sym_outline(symbols[0]),2,17,14],color = [colore[0],colore[0],colore[1],color_s],/top,/left,box = 0
multiplot,/reset
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------------------- Plots of halo outflow for Rachel -----------------
IF keyword_set(outplot) THEN  device,filename = outplot + '_eject_expell_halooutflow.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(halogmass_total + smass_accr_total) - (smass + gmass_hot),psym = 4,/xlog,xrange = xrange,yrange = [2e7,6e10],xtitle = xtitle,ytitle = textoidl('Outflow mass from halos'),symsize = symsize,/ylog,/nodata
oplot,xmass,reexpellmass[*,3],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmass[*,3],psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass,reoutflowmass_uniq[*,3],psym = symcat(14),color = fgcolor,symsize = symsize
oplot,xmass,reoutflowmass_uniq[*,3],psym = symcat(sym_outline(14)),symsize = symsize
oplot,xmass,(halogmass_total + smass_accr_total) - (smass + gmass_hot),psym = symcat(sym_outline(14)),color = 254,symsize = symsize
legend,['Difference between current halo and accreted','Expelled from disk','Uniq. part. exiting halo'],psym = [4,14,14],/left,/top,box = 0,color = [254,colore[0],fgcolor]
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_eject_expell_halooutflowfrac2.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,((halogmass_total + smass_accr_total) - (smass + gmass_hot))/reoutflowmass_uniq[*,3],psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange  = [0.4,1.2],xtitle = xtitle,ytitle = textoidl('Diff in halo mass/halo outflow mass'),symsize = symsize
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_eject_expell_halooutflowfrac.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,((halogmass_total + smass_accr_total) - (smass + gmass_hot))/reexpellmass[*,3],psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange  = [1,3.5],xtitle = xtitle,ytitle = textoidl('Diff in halo mass/expelled mass'),symsize = symsize
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Expelled Mass as a Function of Halo Mass and Time
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_expell_mass_time.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xmass,reexpellmass[*,0],/ylog,/xlog,xrange = xrange,yrange = [5e5,2e11],xtitle = xtitle,ytitle = 'Mass of Gas Expelled [M' + sunsymbol() + ']',/nodata ;yrange = [1e3,1e11]
;    FOR iz = 0, nz - 1 DO oplot,xmass,reexpellmass[*,iz],psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize  
;    FOR iz = 0, nz - 1 DO oplot,xmass,reexpellmass[*,iz],psym = symcat(sym_outline(ex_psym[iz])),symsize = symsize        
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xmass,reexpellmass[*,iz],psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize  
       oplot,xmass,reexpellmass[*,iz],psym = symcat(sym_outline(ex_psym[iz])),symsize = symsize
    ENDFOR
    ;for Rachel
    oplot,xmass,reoutflowmass[*,3],psym = symcat(14),symsize = symsize
    ;
    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = ex_psym,color  = z_colors,/right,/bottom,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
ENDIF

;------------ Ejected Mass as a Function of Halo Mass and Time
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_eject_mass_time.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xmass,reejectmass[*,0],/ylog,/xlog,xrange = xrange,yrange = [5e5,1e11],xtitle = xtitle,ytitle = 'Mass of Gas Ejected [M' + sunsymbol() + ']',/nodata;yrange = [1e4,1e11]
;    FOR iz = 0, nz - 1 DO oplot,xmass,reejectmass[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz] ,symsize = symsize   
;    FOR iz = 0, nz - 1 DO oplot,xmass,reejectmass[*,iz],psym = symcat(sym_outline(ej_psym[iz])),symsize = symsize         
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xmass,reejectmass[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz] ,symsize = symsize 
       oplot,xmass,reejectmass[*,iz],psym = symcat(sym_outline(ej_psym[iz])),symsize = symsize
    ENDFOR
    ;for Rachel
    oplot,xmass,reoutflowmass[*,3],psym = symcat(14),symsize = symsize
    ;
    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = ej_psym,color  = z_colors,/right,/bottom,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
ENDIF

;------------ Fraction of Gas Mass Lost that is Expelled & Ejected
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_lostgas_eject_expell_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,reexpellmass[*,nz - 1]/relostmass[*,nz - 1],/xlog,xrange = xrange,yrange = [0,0.6],xtitle = xtitle,ytitle = 'Outflow Mass/All Gas Mass Lost',/nodata
oplot,xmass,reexpellmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass, reejectmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, reejectmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[1])),symsize = symsize
;oplot,xmass, reejectmass_rv2[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(symbols[2]),color = colore[2],symsize = symsize
;oplot,xmass, reejectmass_rv2[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[2])),symsize = symsize
;oplot,xmass, reejectmass_rv5[*,nz - 1]/diskgmass_rv5[*,nz - 1],psym = symcat(symbols[3]),color = colore[3],symsize = symsize
;oplot,xmass, reejectmass_rv5[*,nz - 1]/diskgmass_rv5[*,nz - 1],psym = symcat(sym_outline(symbols[3])),symsize = symsize
;oplot,xmass, reheatmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize
;oplot,xmass, reheatmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[4])),symsize = symsize
;oplot,xmass,sfmass/diskgmass[*,nz - 1],psym = symcat(symbol_s),color = color_s,symsize = symsize
;,/ylog,/xlog,xrange = [1e8,1e11],yrange = [1e8,1e11],xtitle = xtitle,ytitle = 'Mass of Gas Ejected'
;legend,['Gas Mass Expelled','Gas Mass Ejected','Stellar Mass Formed'],psym = [symbols,symbol_s],color = [colore,color_s],box = 0,/bottom,/right;position = [2.2e11,0.4];/bottom,/right
legend,['Gas Mass Expelled','Gas Mass Ejected'],psym = [symbols[0:1]],color = [colore[0:1]],box = 0,/top,/right;position = [2.2e11,0.4];/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop



;------------ Fraction of Gas Mass Expelled & Ejected
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_eject_expell_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,reexpellmass[*,nz - 1]/diskgmass_lost[*,nz - 1],/xlog,xrange = xrange,yrange = [0,1.0],xtitle = xtitle,ytitle = 'Outflow Mass/Disk Gas Mass',/nodata
oplot,xmass,reexpellmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass, reejectmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, reejectmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(sym_outline(symbols[1])),symsize = symsize

;oplot,xmass,sfmass/diskgmass[*,nz - 1],psym = symcat(symbol_s),color = color_s,symsize = symsize
;,/ylog,/xlog,xrange = [1e8,1e11],yrange = [1e8,1e11],xtitle = xtitle,ytitle = 'Mass of Gas Ejected'
;legend,['Gas Mass Expelled','Gas Mass Ejected','Stellar Mass Formed'],psym = [symbols,symbol_s],color = [colore,color_s],box = 0,/bottom,/right;position = [2.2e11,0.4];/bottom,/right
;IF plot_alllost THEN BEGIN
;   oplot,xmass,relostmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[5]),color = colore[5],symsize = symsize
;   oplot,xmass,relostmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(sym_outline(symbols[5])),symsize = symsize
;ENDIF
legend,['Gas Mass Expelled','Gas Mass Ejected'],psym = [symbols[0:1]],color = [colore[0:1]],box = 0,/top,/right;position = [2.2e11,0.4];/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop


;------------ Fraction of Gas Mass Expelled
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_expelled_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,reexpellmass[*,nz - 1]/diskgmass_lost[*,nz - 1],/xlog,xrange = xrange,yrange = [0,0.4],xtitle = xtitle,ytitle = 'Fraction of Gas Mass Expelled',/nodata;yrange = [0,1]
oplot,xmass,reexpellmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(sym_outline(symbols[0])),symsize = symsize
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Fraction of Gas Mass Expelled Over Time
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_expelled_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xmass,reexpellmass/diskgmass_lost,psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange = [0,0.51],xtitle = xtitle,ytitle = 'Fraction of Gas Mass Expelled',/nodata,symsize = symsize; yrange = [0,1]
;    FOR iz = 0, nz - 1 DO oplot,xmass,reexpellmass[*,iz]/diskgmass[*,iz],psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize
;    FOR iz = 0, nz - 1 DO oplot,xmass,reexpellmass[*,iz]/diskgmass[*,iz],psym = symcat(sym_outline(ex_psym[iz])),symsize = symsize    
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xmass,reexpellmass[*,iz]/diskgmass_lost[*,iz],psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize
       oplot,xmass,reexpellmass[*,iz]/diskgmass_lost[*,iz],psym = symcat(sym_outline(ex_psym[iz])),symsize = symsize
    ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = ex_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

;------------ Fraction of Gas Mass Ejected
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_ejected_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,reejectmass[*,nz - 1]/diskgmass_lost[*,nz - 1],/xlog,xrange = xrange,yrange = [0,0.5],xtitle = xtitle,ytitle = 'Fraction of Gas Mass Ejected',/nodata;yrange = [0,1]
oplot,xmass, reejectmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, reejectmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(sym_outline(symbols[1])),symsize = symsize
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Fraction of Gas Mass Ejected Over Time
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_ejected_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xmass,reejectmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange = [0,0.65],xtitle = xtitle,symsize = symsize,ytitle = 'Fraction of Gas Mass Ejected',/nodata; yrange = [0,1]
;    FOR iz = 0, nz - 1 DO oplot,xmass,reejectmass[*,iz]/diskgmass[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz],symsize = symsize
;    FOR iz = 0, nz - 1 DO oplot,xmass,reejectmass[*,iz]/diskgmass[*,iz],psym = symcat(sym_outline(ej_psym[iz])),symsize = symsize
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xmass,reejectmass[*,iz]/diskgmass_lost[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz],symsize = symsize
       oplot,xmass,reejectmass[*,iz]/diskgmass_lost[*,iz],psym = symcat(sym_outline(ej_psym[iz])),symsize = symsize
    ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

;------------ Fraction of Ejected Gas Expelled
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_expelled_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,reejectmass[*,nz - 1]/relostmass[*,nz - 1],ytitle = textoidl('M/M_{lost}'),xtitle = xtitle,yrange = [0,0.6],xrange = xrange,psym = symcat(18),/nodata,/xlog
loadct,0
IF keyword_set(outplot) THEN bar_thick = 14 ELSE bar_thick = 4
FOR i = 0,n_elements(xmass) - 1 DO oplot,[xmass[i],xmass[i]],[0,1.0],color = 100,thick = bar_thick
IF keyword_set(colors) THEN loadct,39
FOR i = 0,n_elements(xmass) - 1 DO oplot,[xmass[i],xmass[i]],[0,reejectmass[i,nz - 1]/relostmass[i,nz - 1]],color = colore[1],thick = bar_thick
FOR i = 0,n_elements(xmass) - 1 DO oplot,[xmass[i],xmass[i]],[0,reexpellmass[i,nz - 1]/relostmass[i,nz - 1]],color = colore[0],thick = bar_thick
oplot,xmass,reejectmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(symbols[1]),symsize = symsize,color = colore[1]
oplot,xmass,reejectmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,xmass,reexpellmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(symbols[0]),symsize = symsize,color = colore[0]
oplot,xmass,reexpellmass[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[0])),symsize = symsize
;oplot,xmass,reejectmass_rv2[*,nz - 1]/reheatmass[*,nz - 1],psym = symcat(symbols[2]),symsize = symsize,color = colore[2]
;oplot,xmass,reejectmass_rv2[*,nz - 1]/reheatmass[*,nz - 1],psym = symcat(sym_outline(symbols[2])),symsize = symsize
;legend,['Ejected from disk','Ejected beyond 0.2*R' + textoidl('_{vir}')],psym = [16,9],box = 0,/right
legend,['Gas Expelled','Gas Ejected'],psym = symbols[0:1],color = colore[0:1],/top,/right,box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Fraction of Ejected Gas Expelled (2)
IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_expelled_mass2.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,reexpellmass[*,nz - 1]/reejectmass[*,nz - 1],psym = symcat(9),/xlog,xrange = xrange,yrange = [0,1.0],xtitle = xtitle,ytitle = 'Ratio of Ejected to Expelled Mass',symsize = symsize,/nodata
oplot,xmass,reexpellmass[*,nz - 1]/reejectmass[*,nz - 1],symsize = symsize,psym = symcat(9),color = colore[1]
oplot,xmass,reexpellmass[*,nz - 1]/relostmass[*,nz - 1],symsize = symsize,psym = symcat(symbol_a),color = colore[4]
;oplot,xmass,reexpellmass[*,nz - 1]/relostmass[*,nz - 1],symsize = symsize,psym = symcat(symbol_a),color = colore[5]
legend,['Lost','Unbound from disk','Beyond R' + textoidl('_{vir}')],psym = [symbol_a,9,36],color = [colore[4],colore[1],colore[2]],box = 0
;oplot,xmass,reexpellmass[*,nz - 1]/reejectmass_rv5[*,nz - 1],symsize = symsize,psym = symcat(35)
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Fraction of Ejected Gas Expelled Over Time
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_frac_gas_expelled_mass_time.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xmass,reexpellmass[*,nz - 1]/reejectmass[*,nz - 1],/xlog,xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = 'Fraction of Gas Mass Ejected that is Expelled',/nodata;yrange = [0,1]
;    FOR iz = 0, nz - 1 DO oplot,xmass,reexpellmass[*,iz]/reejectmass[*,iz],psym = symcat(z_psym[iz]),color  = z_colors[iz],symsize = symsize
;    FOR iz = 0, nz - 1 DO oplot,xmass,reexpellmass[*,iz]/reejectmass[*,iz],psym = symcat(sym_outline(z_psym[iz])),symsize = symsize
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xmass,reexpellmass[*,iz]/reejectmass[*,iz],psym = symcat(z_psym[iz]),color  = z_colors[iz],symsize = symsize
       oplot,xmass,reexpellmass[*,iz]/reejectmass[*,iz],psym = symcat(sym_outline(z_psym[iz])),symsize = symsize
    ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)') + ' ',psym = z_psym,color  = z_colors,/right,/bottom,box=0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
ENDIF

mlrange = [1e9,1e12]
mlrange = [10,1000]
;------------ Mass Loading
xarr = vcirc ;xmass
xarrt = vcirct
IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xarr,reexpellmass[*,nz - 1]/sfmass_total,yrange = [0.1,100],ytitle = textoidl('Effective mass loading factor'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [25,200],/xlog
oplot,xarr,reexpellmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xarr,reexpellmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xarr,reejectmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xarr,reejectmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[1])),symsize = symsize
;oplot,xarr,reejectmass_rv5[*,nz - 1]/sfmass_total,psym = symcat(symbols[3]),color = colore[3],symsize = symsize
;oplot,xarr,reejectmass_rv5[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[3])),symsize = symsize
;oplot,xarr,reheatmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[4]),color = colore[4],symsize = symsize
;oplot,xarr,reheatmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[4])),symsize = symsize
oplot,xarr,relostmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[5]),color = colore[5],symsize = symsize
oplot,xarr,relostmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[5])),symsize = symsize
fits_expell = robust_linefit( alog10(xarr), alog10(reexpellmass[*,nz - 1]/sfmass_total), reexpellmass_fit, sigma_expell )
fits_eject = robust_linefit( alog10(xarr), alog10(reejectmass[*,nz - 1]/sfmass_total), reejectmass_fit, sigma_eject )
print,'Expell (total): ',fits_expell
print,'Eject (total): ',fits_eject
;oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),color = colore[0],linestyle = 0,psym = -3,thick = thicks[0]
;oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),color = colore[1],linestyle = 0,psym = -3,thick = thicks[0]
;legend,['Gas Mass Expelled/Stellar Mass Formed','Gas Mass Ejected/Stellar Mass Formed'],psym = symbols,color = colore,/bottom,/left,box=0
legend,['Gas Expelled','Gas Ejected'],psym = symbols[0:1],color = colore[0:1],/bottom,/left,box=0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Total massloading vs stellar mass
IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_smass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = xsize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = xsize
plot,smass,reexpellmass[*,nz - 1]/sfmass_total,yrange = [0.1,100],ytitle = textoidl('Effective mass loading factor'),xtitle = textoidl('M_*/M')+sunsymbol(),/nodata,/ylog,xrange = [1e7,1e11],/xlog,POSITION=ASPECT(1, Margin=0.10)
oplot,smass,reexpellmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,smass,reexpellmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,smass,reejectmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,smass,reejectmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[1])),symsize = symsize
;oplot,smass,reejectmass_rv2[*,nz - 1]/sfmass_total,psym = symcat(symbols[2]),color = colore[2],symsize = symsize
;oplot,smass,reejectmass_rv2[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[2])),symsize = symsize
;oplot,smass,reejectmass_rv5[*,nz - 1]/sfmass_total,psym = symcat(symbols[3]),color = colore[3],symsize = symsize
;oplot,smass,reejectmass_rv5[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[3])),symsize = symsize
;oplot,smass,reheatmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[4]),color = colore[4],symsize = symsize
;oplot,smass,reheatmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[4])),symsize = symsize
oplot,smass,relostmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[5]),color = colore[5],symsize = symsize
oplot,smass,relostmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[5])),symsize = symsize
legend,['Gas Expelled','Gas Ejected','Gas Lost from Disk'],psym = [symbols[0],symbols[1],symcat(symbols[4])],color = [colore[0],colore[1],colore[4]],/bottom,/left,box=0
IF keyword_set(outplot) THEN device, /close ELSE stop

;------------ Mass Loading Over Time
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_mass_time.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarr,reejectmassr[*,0]/sfmassr[*,0],psym = symcat(symbols[0]),/xlog,xrange = [20,200],yrange = [0.09,50],ytitle = textoidl('\eta_{ejected}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog
    FOR iz = 0, nz - 1 DO oplot,xarr,reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz]
    FOR iz = 0, nz - 1 DO oplot,xarr,reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(ej_psym[iz]))
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xarr,reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz]
       oplot,xarr,reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(ej_psym[iz]))
       ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;------------ Mass Loading Over Time with current xaxis
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_mass_time_eject.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarrt[*,0],reejectmassr[*,0]/sfmassr[*,0],/xlog,yrange = [0.09,200],ytitle = textoidl('\eta_{ejected}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [10,150]
    vc_arr = findgen(200)
    eta_m15 = vc_arr
    eta_m15[where(vc_arr LT 60)] = 2.9*(vc_arr[where(vc_arr LT 60)]/60)^(-3.2)
    eta_m15[where(vc_arr GE 60)] = 2.9*(vc_arr[where(vc_arr GE 60)]/60)^(-1.0)
    loadct,0
    oplot,vc_arr,eta_m15,thick = 2,color = 100,linestyle = 5
    oplot,vc_arr,4500*vc_arr^(-2.0),thick = 2,color = 100,linestyle = 2 ;Energy driven
    oplot,vc_arr,75*vc_arr^(-1.0),thick = 2,color = 100,linestyle = 3 ;Momentum driven
    legend,['Energy-driven','Momentum-driven','Muratov et al. 2015','Fit to data'],linestyle = [2,3,5,0],thick = [2,2,2,2],color = [100,100,100,fgcolor],box = 0,/bottom,/left
    loadct,39
;    oplot,xarr,reejectmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[1])),symsize = symsize
;    oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),linestyle = 0,psym = -3,thick = thicks[0]
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xarrt[*,iz],reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,xarrt[*,iz],reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(ej_psym[iz])),symsize = symsize*0.6
;       fits_eject = robust_linefit( alog10(xarrt[*,iz]), alog10(reejectmassr[*,iz]/sfmassr[*,iz]), reejectmass_fit, sigma_eject )
;       print,fits_eject
;       oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),color = z_colors[iz],linestyle = 0,psym = -3,thick = thicks[0]
    ENDFOR
    fits_eject = robust_linefit( alog10(xarrt), alog10(reejectmassr/sfmassr), reejectmass_fit, sigma_eject )
    print,'Eject: ',fits_eject
    oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),linestyle = 0,psym = -3,thick = thicks[0]
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;----------- Mass Loading Over Time vs stellar mass            --------------------------
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_smass_time_eject.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = xsize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = xsize
    plot,smasst[*,0],reejectmassr[*,0]/sfmassr[*,0],/xlog,yrange = [0.1,100],ytitle = textoidl('\eta_{ejected}'),xtitle = textoidl('M_*/M')+sunsymbol(),/nodata,/ylog,xrange = [1e6,1e10],POSITION=ASPECT(1, Margin=0.10)
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,smasst[*,iz],reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,smasst[*,iz],reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(ej_psym[iz])),symsize = symsize*0.6
;       oplot,smasst[*,iz],reheatmassr[*,iz]/sfmassr[*,iz],psym = symcat(symbols[4]),color  = z_colors[iz],symsize = symsize*0.6
;       oplot,smasst[*,iz],reheatmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(symbols[4])),symsize = symsize*0.6,color  = z_colors[iz]
    ENDFOR
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/left,/bottom,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;------------ Mass Loading Over Time with current xaxis
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_mass_time_eject_rv2.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarrt[*,0],reejectmassr_rv2[*,0]/sfmassr[*,0],/xlog,yrange = [0.09,50],ytitle = textoidl('\eta_{ejected}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [20,200]
 ;   oplot,xarr,reejectmass_rv2[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[2])),symsize = symsize
    FOR iz = 0, nz - 1 DO BEGIN ;h986.2
       oplot,xarrt[*,iz],reejectmassr_rv2[*,iz]/sfmassr[*,iz],psym = symcat(symbols[2]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,xarrt[*,iz],reejectmassr_rv2[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(symbols[2])),symsize = symsize*0.6
    ENDFOR
    fits_eject = robust_linefit( alog10(xarrt), alog10(reejectmassr_rv2/sfmassr), reejectmass_fit, sigma_eject )
    print,'Eject (rvir/2): ',fits_eject
    oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),linestyle = 0,psym = -3,thick = thicks[0]
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;------------ Mass Loading Over Time with current xaxis
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_mass_time_eject_rv5.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarrt[*,0],reejectmassr_rv5[*,0]/sfmassr[*,0],/xlog,yrange = [0.09,50],ytitle = textoidl('\eta_{ejected}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [20,200]
;    oplot,xarr,reejectmass_rv5[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[3])),symsize = symsize
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xarrt[*,iz],reejectmassr_rv5[*,iz]/sfmassr[*,iz],psym = symcat(symbols[3]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,xarrt[*,iz],reejectmassr_rv5[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(symbols[3])),symsize = symsize*0.6
    ENDFOR
    fits_eject = robust_linefit( alog10(xarrt), alog10(reejectmassr_rv5/sfmassr), reejectmass_fit, sigma_eject )
    print,'Eject (rvir/5): ',fits_eject
    oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),linestyle = 0,psym = -3,thick = thicks[0]
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;------------ Mass Loading Over Time with current xaxis
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_mass_time_lost.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarrt[*,0],relostmassr[*,0]/sfmassr[*,0],/xlog,yrange = [0.9,150],ytitle = textoidl('\eta_{lost}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [20,200]
;    oplot,xarr,relostmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[4])),symsize = symsize
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xarrt[*,iz],relostmassr[*,iz]/sfmassr[*,iz],psym = symcat(symbols[4]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,xarrt[*,iz],relostmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(symbols[4])),symsize = symsize*0.6
    ENDFOR
    fits_eject = robust_linefit( alog10(xarrt), alog10(relostmassr/sfmassr), reejectmass_fit, sigma_eject )
    print,'Lost: ',fits_eject
    oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),linestyle = 0,psym = -3,thick = thicks[0]
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;------------ Mass Loading Over Time for material lost from halo with current xaxis
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_mass_time_halooutflow.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarrt[*,0],reoutflowmassr[*,0]/sfmassr[*,0],/xlog,yrange = [0.9,200],ytitle = textoidl('\eta_{HHT}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [10,200]
;    oplot,xarr,relostmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[4])),symsize = symsize
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xarrt[*,iz],reoutflowmassr[*,iz]/sfmassr[*,iz],psym = symcat(symbols[4]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,xarrt[*,iz],reoutflowmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(symbols[4])),symsize = symsize*0.6
    ENDFOR
    fits_eject = robust_linefit( alog10(xarrt), alog10(reoutflowmassr/sfmassr), reejectmass_fit, sigma_eject )
    print,'Halo Outflow: ',fits_eject
    oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),linestyle = 0,psym = -3,thick = thicks[0]
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;------------ Mass Loading Over Time with current xaxis
IF keyword_set(z_cut) THEN BEGIN
    IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_mass_time_expell.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
    plot,xarrt[*,0],reexpellmassr[*,0]/sfmassr[*,0],/xlog,yrange = [0.09,50],ytitle = textoidl('\eta_{expelled}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [20,200]
 ;   oplot,xarr,reexpellmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[0])),symsize = symsize
;    oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),linestyle = 0,psym = -3,thick = 1
    FOR iz = 0, nz - 1 DO BEGIN
       oplot,xarrt[*,iz],reexpellmassr[*,iz]/sfmassr[*,iz],psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
       oplot,xarrt[*,iz],reexpellmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(ex_psym[iz])),symsize = symsize*0.6
;       fits_expell = robust_linefit( alog10(xarrt[*,iz]), alog10(reexpellmassr[*,iz]/sfmassr[*,iz]), reexpellmass_fit, sigma_expell )
;       print,fits_expell
;       oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),color = z_colors[iz],linestyle = 0,psym = -3,thick = thicks[0]
    ENDFOR
    fits_expell = robust_linefit( alog10(xarrt), alog10(reexpellmassr/sfmassr), reexpellmass_fit, sigma_expell )
    print,'Reexpell: ',fits_expell
    oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),linestyle = 0,psym = -3,thick = thicks[0]
    legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ex_psym,color  = z_colors,/right,/top,box = 0
    IF keyword_set(outplot) THEN device, /close ELSE stop    
 ENDIF

;----------- Multi-plot mass loading
IF keyword_set(z_cut) THEN BEGIN
   IF keyword_set(outplot) THEN  device,filename = outplot + '_massloading_mass_time_multi.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*1.6,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize*1.6
   multiplot,[1,2]

   plot,xarrt[*,0],reejectmassr[*,0]/sfmassr[*,0],/xlog,yrange = [0.09,50],ytitle = textoidl('\eta_{ejected}'),/nodata,/ylog,xrange = [20,200]
   FOR iz = 0, nz - 1 DO BEGIN
      oplot,xarrt[*,iz],reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(ej_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
      oplot,xarrt[*,iz],reejectmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(ej_psym[iz])),symsize = symsize*0.6
   ENDFOR
   fits_eject = robust_linefit( alog10(xarrt), alog10(reejectmassr/sfmassr), reejectmass_fit, sigma_eject )
   oplot,mlrange,10^fits_eject[0]*mlrange^(fits_eject[1]),linestyle = 0,psym = -3,thick = thicks[0]
   legend,'z = ' + string(z_bins_legend,format = '(A3)'),psym = ej_psym,color  = z_colors,/right,/top,box = 0

   multiplot
   plot,xarrt[*,0],reexpellmassr[*,0]/sfmassr[*,0],/xlog,yrange = [0.09,50],ytitle = textoidl('\eta_{expelled}'),xtitle = textoidl('V_{circ} [km s^{-1}]'),/nodata,/ylog,xrange = [20,200]
   FOR iz = 0, nz - 1 DO BEGIN
      oplot,xarrt[*,iz],reexpellmassr[*,iz]/sfmassr[*,iz],psym = symcat(ex_psym[iz]),color  = z_colors[iz],symsize = symsize*0.6
      oplot,xarrt[*,iz],reexpellmassr[*,iz]/sfmassr[*,iz],psym = symcat(sym_outline(ex_psym[iz])),symsize = symsize*0.6
   ENDFOR
   fits_expell = robust_linefit( alog10(xarrt), alog10(reexpellmassr/sfmassr), reexpellmass_fit, sigma_expell )
   oplot,mlrange,10^fits_expell[0]*mlrange^(fits_expell[1]),linestyle = 0,psym = -3,thick = thicks[0]
   multiplot,/reset
   IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF

;------------ Mass expelled/ejected, fraction of disk gas expelled/ejected, mass loading ----------
IF keyword_set(outplot) THEN  device,filename = outplot + '_multi_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize*2.6,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize*2.6
multiplot,[1,3]
plot,xmass,reexpellmass[*,nz - 1],/ylog,/xlog,xrange = xrange,yrange = [1e6,1e11],ytitle = 'Mass [M' + sunsymbol() + ']',/nodata
oplot,xmass,reexpellmass[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmass[*,nz - 1],psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass, reejectmass[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, reejectmass[*,nz - 1],psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,xmass,smass,psym = symcat(symbol_s),color = color_s,symsize = symsize ;smass
;,/ylog,/xlog,xrange = [1e8,1e11],yrange = [1e8,1e11],xtitle = xtitle,ytitle = 'Mass of Gas Ejected'
oplot,xmass,smass,psym = symcat(sym_outline(symbol_s)),color = fgcolor,symsize = symsize
legend,['Gas ejected','Gas expelled beyond '+textoidl('R_{vir}'),'Stars'],psym = [symbols[1],symbols[0],symbol_s],color = [colore[1],colore[0],color_s],/bottom,/right,box = 0
multiplot

plot,xmass,reexpellmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange = [0.1,100],ytitle = textoidl('Effective mass loading factor'),/nodata,/ylog
oplot,xmass,reexpellmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[0]),color = colore[0],symsize = symsize
oplot,xmass,reexpellmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass,reejectmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass,reejectmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,xmass,relostmass[*,nz - 1]/sfmass_total,psym = symcat(symbols[4]),color = colore[4],symsize = symsize
oplot,xmass,relostmass[*,nz - 1]/sfmass_total,psym = symcat(sym_outline(symbols[4])),symsize = symsize
legend,['Gas lost from disk'],psym = [symbols[4]],color = [colore[4]],symsize = symsize,/left,/bottom,box = 0
;legend,['Gas Mass Expelled/Stellar Mass Formed','Gas Mass Ejected/Stellar Mass Formed'],psym = symbols,color = colore,/bottom,/left,box=0
multiplot

plot,xmass,reexpellmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[0]),/xlog,xrange = xrange,yrange = [0.06,0.5],ytitle = textoidl('M_{outflow}/M_{accreted}'),xtitle = xtitle,/nodata,symsize = symsize,/ylog
;oplot,xmass,reexpellmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[0]),color = colore[0],symsize = symsize
;oplot,xmass,reexpellmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(sym_outline(symbols[0])),symsize = symsize
;oplot,xmass,reexpellmass_cool[*,nz - 1]/relostmass[*,nz - 1],psym = symcat(sym_outline(symbols[0])),color = colore[0],symsize = symsize
oplot,xmass, reejectmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[1]),color = colore[1],symsize = symsize
oplot,xmass, reejectmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(sym_outline(symbols[1])),symsize = symsize
;oplot,xmass, reheatmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(symbols[4]),color = colore[4],symsize = symsize
;oplot,xmass, reheatmass[*,nz - 1]/diskgmass_lost[*,nz - 1],psym = symcat(sym_outline(symbols[4])),symsize = symsize
;legend,['Gas Mass Expelled','Gas Mass Ejected'],psym = [symbols],color = [colore],box = 0,/bottom,/right;position = [2.2e11,0.4];/bottom,/right
;legend,['Total gas mass expelled','Gas mass expelled -- SN heated','Total gas mass ejected','Gas mass ejected -- SN heated'],psym=[symbols[0],sym_outline(symbols[0]),symbols[1],sym_outline(symbols[1])],color = [colore[0],colore[0],colore[1],colore[1]],/top,/right,box = 0
legend,['Total gas mass ejected','Gas mass ejected -- SN heated'],psym=[symbols[1],sym_outline(symbols[1])],color = [colore[1],colore[1]],/left,/bottom,box = 0
multiplot,/reset
IF keyword_set(outplot) THEN device, /close ELSE stop

END
