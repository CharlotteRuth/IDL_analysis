;inflow_outflow_historyz,['/nobackupp8/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'],haloid=['1'],/color,/debug

;relost -- leaves the disk
;reheat -- leaves the disk AND was affected by supernovae
;inflowmdiskcum -- all gas that enters the disk



PRO inflow_outflow_historyz, dirs, haloid = haloid, outplot = outplot, colors = colors, thicks = thicks, linestyle = linestyle, scale = scale,formatthick = formatthick, keys = keys,pmulti = pmulti,allpositive = allpositive,yrange_max = yrange_max,label = label,debug = debug,gmass = gmass,plot_halo = plot_halo
!EXCEPT = 0 ;Supress math errors
zsolar = 0.0130215
;----------------------- Set parameters ------------------
IF keyword_set(allpositive) THEN e_m = 1 ELSE e_m = -1

n = n_elements(dirs)
IF NOT keyword_set(haloid) THEN haloid = strarr(n_elements(dirs)) + '1'

; *************************** Plot Formatting ********************
IF keyword_set(outplot) THEN formatplot,/outplot,thick = formatthick ELSE formatplot
IF keyword_set(colors) THEN BEGIN
    distinct_colors,n_colors = 12                            ;loadct,39
    white = 13 ;white = 255
    black = 0 ;black = 0
    IF keyword_set(outplot) THEN BEGIN
        fgcolor = black
        bgcolor = white
    ENDIF ELSE BEGIN
        fgcolor = white
        bgcolor = black
    ENDELSE
;    cheat = 230
;    cexpell = 60
;    cpristine = 30
;    cstar =130
;    ceject = 85
    cheat = 6                   ;150 ;yellow/pea green
    cexpell = 2                 ;60 ;blue
    ceject = 11                 ;254 ;red
    cpristine = 1              ;20 ;purple/blue
    cstar = 5                   ;210 ;green
    cheat_halo = 7              ; gold/pea green
;    cexpell_halo = 3            ; cyan
    csim = 4
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(outplot) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 1
    IF NOT keyword_set(linestyle) THEN linestyle = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    IF keyword_set(outplot) THEN BEGIN
        fgcolor = 0
        bgcolor = 255
    ENDIF ELSE BEGIN
        fgcolor = 255
        bgcolor = 0
    ENDELSE
    colors = findgen(5)*fgcolor/5
    cheat = colors[0]
    cexpell = colors[1]
    ceject = colors[2]
    cpristine = colors[3]
    cstar = colors[4]
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(outplot) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 1
    IF NOT keyword_set(linestyle) THEN linestyle = findgen(n)*2   
ENDELSE
IF NOT keyword_set(scale) THEN scale = fltarr(n) + 1.0

;--------------------------- Set (multi)plotting parameters --------------
IF keyword_set(pmulti) THEN BEGIN
   IF (keyword_set(outplot)) THEN device,filename = outplot + '_fbz.eps',/color,bits_per_pixel= 8,/times,xsize = 36,ysize = 45,xoffset =  2,yoffset =  2 ELSE window,3,xsize = 580,ysize = 800
   multiplot,pmulti,mxTitle = 'Time [Gyr]',myTitle = textoidl('Cumulative M_{Z}/M_{Z, total produced}'),gap = 0,mxTitSize = 1.5,myTitSize = 1.5,/square;,xgap = 0.01;myTitle = 'Cumulative Metal Mass [M' + sunsymbol() + ']'
   xtitle = ''
   ytitle = ''
   xstyle = 1
ENDIF ELSE BEGIN
;   multiplot,/resetx
;   multiplot,/default
   xtitle = 'Time [Gyr]'
   ytitle = 'Mass [M'+ sunsymbol() + ']'
   xstyle = 9
   !x.charsize=1.5    
   !y.charsize=1.5
ENDELSE
binsize = 5e8
nbins = round(14.0/(binsize/1e9))


;-------------------------- Conversion between redshift and time ----------------
ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
z = reverse(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
ticktime_in_z = findgen(7)*2*1e9 ;redshift major plot
tickred_in_z = spline(t,z,ticktime_in_z )
tickred_in_t = reverse(findgen(5))
ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e

range = [[-5.5e5,3e5],[-2.5e6,1.5e6],[-5e6,3e6],[-1.75e8,1.5e+8],[ -1.4e+09,1.5e+09]]
range = [[-3.3,3],[-3,2],[-2.5,2],[-3.3,3],[-3.3,3]]
;******************** Iterate through files ***************
;At i = 6
;(/nobackupp8/crchrist/MolecH/h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/9)
;7/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/2
;8/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/1
;11/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/1
;14/nobackupp8/crchrist/MolecH/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/2
;Got to 16
;print,i,dirs[i],haloid[i]
FOR i = 0, n - 1 DO BEGIN ;Iterate through the files
;iarray = [1,12,15,17]
;FOR i_iarr = 0, n_elements(iarray) - 1 DO BEGIN
;    i = iarray[i_iarr]
   print,dirs[i],haloid[i]
   cd,dirs[i]
   halodat = mrdfits('grp'+haloid[i]+'.alignment.fits',1,/silent)
   filebase = 'grp' + haloid[i] + '.mass_metal_track.dat'
   spawn,'ls ' + dirs[i] + '*.grp' + haloid[i] + '.haloid.dat',files_haloid
   spawn,'ls ' + dirs[i] + 'h*param',pfile
   units = tipsyunits(pfile[0],/silent)
   IF max(halodat.time) LT 0.5 THEN halodat.time = halodat.time*units.timeunit/1e9

;--------------- Read Halo Information at each step -------------
   IF keyword_set(debug) THEN readcol,files_haloid[0],files,halos,format='a,l',/silent
   IF keyword_set(debug) THEN BEGIN
      readcol,dirs[i] + filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,ox,fe,mHI,mH2,mcoldg,/silent ;,mH2,H2frac,Pressure,r25,SFR
      IF n_elements(time) NE n_elements(halodat.time) THEN print,'Number of times not equal';stop
      time = halodat.time
   ENDIF ELSE time = halodat.time
   readcol,'grp' + haloid[i] + '.sfmetals.txt',zstars,oxstars,festars,zstars_in,oxstars_in,festars_in,/silent
   readcol,dirs[i] + '/grp' + haloid[i] + '.metals.txt',zmetals,ox,fe,coldg,zmetals_H,ox_H,fe_H,mHI_m,mH2_m,Hgas,/silent 
   IF keyword_set(plot_halo) THEN readcol,dirs[i] + '/grp' + haloid[i] + '.halometals.txt',zmetals_ha,ox_ha,fe_ha,/silent
   IF n_elements(zstars) NE n_elements(zmetals) THEN stop
   zstars_cum = fltarr(n_elements(zstars))
   FOR j = 0, n_elements(zstars) - 1 DO zstars_cum[j] = total(zstars[0:j])
;   IF keyword_set(debug) THEN scale[i] = 1.0/mtot[n_elements(mtot) - 1] ;Paramter to scale by virial mass
   nhalostat = n_elements(halo)

;-------------- Read in information about final step ---------
   spawn,'ls ' + dirs[i] + '*512/*cosmo**512',file
   rtipsy,file,h,g,d,s
   spawn,'ls ' + dirs[i] + '*512/*512.iord',file_iord
   readarr,file_iord,  h,iord,/ascii,type = 'long'
   iord_star = iord[h.ngas + h.ndark:h.n - 1]
   iord_gas = iord[0:h.ngas - 1]
   spawn,'ls ' + dirs[i] + '*512/*512.OxMassFrac',file_ox
   readarr,file_ox,  h,ox,/ascii
   ox_star = ox[h.ngas + h.ndark:h.n - 1]
   spawn,'ls ' + dirs[i] + '*512/*512.FeMassFrac',file_fe
   readarr,file_fe,  h,fe,/ascii
   fe_star = fe[h.ngas + h.ndark:h.n - 1]
   zmetal_star = 2.09*ox_star + 1.06*fe_star
   s.metals = zmetal_star
   g.zmetal = 2.09*ox[0:h.ngas - 1] + 1.06*fe[0:h.ngas - 1]
   IF keyword_set(debug) THEN BEGIN
;      rtipsy,file,h,g,d,s
       spawn,'ls ' + dirs[i] + '*512/*512.amiga.grp',file_grp
       readarr,file_grp,h,grp,/ascii,type = 'long'
       grp_star = grp[h.ngas + h.ndark:h.n - 1]
       grp_gas = grp[0:h.ngas - 1]
;      spawn,'ls ' + dirs[i] + '*512/*512.massform',file_massform
;      readarr,file_massform,  h,massform,/ascii
;      massform_star = massform[h.ngas + h.ndark:h.n - 1]
;   spawn,'ls ' + dirs[i] + '*512/*512.coolontime',file_coolon
;   readarr,file_coolon,h,coolon,part = 'gas',/ascii
       spawn,'ls ' + dirs[i] + '*512/*512.amiga.stat',file_stat
       stat = read_stat_struc_amiga(file_stat[0])
       main = where(stat.group EQ haloid[i])
       satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
       inhalo = fltarr(n_elements(grp))
       FOR j = 0, n_elements(satellites) - 1 DO BEGIN
           test = where(grp EQ stat[satellites[j]].group, ntest)
           IF ntest NE 0 THEN inhalo[test] = 1
       ENDFOR
;1 if in halo or satellite at z = 0, 0 otherwise
       inhalo_star = inhalo[h.ngas + h.ndark:h.n - 1] 
       inhalo_gas = inhalo[0:h.ngas - 1]
       
;------------ Read in information about stars 
; ------------------- Formed Stars -----------------------------
;Read in stars formed within halo
       sl = mrdfits(dirs[i] + 'starlog.cut.'+haloid[i] + '.fits',1,/silent)
;Remove duplicates from starlog file
       sl = sl[uniq(sl.iorderstar,sort(sl.iorderstar))]
       match2,sl.iorderstar,iord_star,ind1,ind2
       formstarcum = weighted_histogram(sl.timeform, weight = sl.massform*units.massunit,min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
       formstarzcum = weighted_histogram(sl.timeform, weight = sl.massform*units.massunit*s[ind1].metals,min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
       
;Determine what gas particles stars that formed in the main progenitor
;formed from
       spawn,'ls ' + dirs[i] + '*512/*512.igasorder',file_igasorder
       readarr,file_igasorder,h,igasorder,/ascii,type = 'long',part = 'star'
       match,sl.iorderstar,iord_star,ind1,ind2
       sl[ind1].iordergas = igasorder[ind2]
   ENDIF ELSE rtipsy,file[0],h,g,d,s,/justhead

; ------------------- Inflow on to Halo -----------------------------
   print,'Read accretion files'
   inflowz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccr_z.fits',/silent)          ;'accr_rvir_z';'.accrz.fits'
   inflowm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reaccr.fits',/silent)    ;'mass_at_accr_rvir';'.mass_at_accr.fits'
   inflowi = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccr_iord.fits',/silent)
   inflow_data = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccr_history.fits',1,/silent)
   inflowi_sf = inflowi[where(inflowm LE 0)]
;      match,inflowi_sf,slall.iordergas,ind1,ind2
;      inflowm[(where(inflowm LE 0))[ind1]] = slall[ind2].massform*units.massunit
   inflowm[where(inflowm LE 0)] = 0 ;Set mass to zero when it comes in as a star
   inflowt = z_to_t(inflowz)
   inflowmcum  = weighted_histogram(inflowt, weight =  inflowm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
   inflowzcum  = weighted_histogram(inflowt, weight =  inflowm*inflow_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   nbins = n_elements(inflowmcum)
   
; Need unique accretion onto halo here
   ind_inflow_uniq = rem_dup(inflowi) ;uniq(inflowi,sort(inflowi))
   inflowz_uniq = inflowz[ind_inflow_uniq]
   inflowm_uniq = inflowm[ind_inflow_uniq]
   inflowi_uniq = inflowi[ind_inflow_uniq]
   inflow_data_uniq = inflow_data[ind_inflow_uniq]
   inflowt_uniq = z_to_t(inflowz_uniq)    
   inflowmcum_uniq = weighted_histogram(inflowt_uniq, weight =  inflowm_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
   inflowzcum_uniq = weighted_histogram(inflowt_uniq, weight =  inflowm_uniq*inflow_data_uniq.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)


; ------------------- Inflow Disk -----------------------------
   inflowdiskz_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdiskall_z.fits',/silent)
   inflowdiskm_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reaccrdiskall.fits',/silent)
   inflowdiski_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdiskall_iord.fits',/silent)
   inflowdisk_all_data = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdiskall_history.fits',1,/silent) 
;   inflowdiski_all_sf = inflowdiski_all[where(inflowdiskm_all LE 0)]
;   match,inflowdiski_all_sf,slall.iordergas,ind1,ind2
;   inflowdiskm_all[(where(inflowdiskm_all LE 0))[ind1]] = slall[ind2].massform*units.massunit
   inflowdiskm_all[where(inflowdiskm_all LE 0)] = 0
   inflowdiskt_all = z_to_t(inflowdiskz_all)
   inflowmdiskcum_all  = weighted_histogram(inflowdiskt_all, weight =  inflowdiskm_all,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   inflowzdiskcum_all  = weighted_histogram(inflowdiskt_all, weight =  inflowdiskm_all*inflowdisk_all_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   inflowmdiskhist_all  = weighted_histogram(inflowdiskt_all, weight =  inflowdiskm_all,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05)
   inflowzdiskhist_all  = weighted_histogram(inflowdiskt_all, weight =  inflowdiskm_all*inflowdisk_all_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05)

;Unique accretion onto disk
;Select for those particles the first time they are accreted onto the disk
   ind_inflowdisk_all_uniq = rem_dup(inflowdiski_all);uniq(inflowdiski_all,sort(inflowdiski_all))
   inflowdiski_all_uniq = inflowdiski_all[ind_inflowdisk_all_uniq]
   inflowdiskz_all_uniq = inflowdiskz_all[ind_inflowdisk_all_uniq]
   inflowdiskm_all_uniq = inflowdiskm_all[ind_inflowdisk_all_uniq]
   inflowdiskmet_all_uniq_disk = inflowdisk_all_data[ind_inflowdisk_all_uniq].metallicity ;Metallicity when accreted to the disk
; Metallicity when accreted to the halo
   match2,inflowi_uniq,inflowdiski_all_uniq,ind_halo,ind_disk  
   IF n_elements(inflowdiski_all_uniq) NE n_elements(where(ind_halo NE -1)) THEN BREAK ;Sanity check that all particle accreted to disk are also accreted to the halo
   inflowdiskmet_all_uniq = inflow_data_uniq[ind_halo[where(ind_halo NE -1)]].metallicity
   inflowdiskt_all_uniq = z_to_t(inflowdiskz_all_uniq)

   inflowmdiskcum_all_uniq = weighted_histogram(inflowdiskt_all_uniq, weight =  inflowdiskm_all_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
   inflowzdiskcum_all_uniq =      weighted_histogram(inflowdiskt_all_uniq, weight =  inflowdiskm_all_uniq*inflowdiskmet_all_uniq, min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
   inflowzdiskcum_all_uniq_disk = weighted_histogram(inflowdiskt_all_uniq, weight =  inflowdiskm_all_uniq*inflowdiskmet_all_uniq_disk, min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
   inflowmdiskhist_all_uniq = weighted_histogram(inflowdiskt_all_uniq, weight =  inflowdiskm_all_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05) 
   inflowzdiskhist_all_uniq = weighted_histogram(inflowdiskt_all_uniq, weight =  inflowdiskm_all_uniq*inflowdiskmet_all_uniq, min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05)
   inflowzdiskhist_all_uniq_disk = weighted_histogram(inflowdiskt_all_uniq, weight =  inflowdiskm_all_uniq*inflowdiskmet_all_uniq_disk, min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05) 
   IF keyword_set(debug) THEN BEGIN
;-------------------- Inflow Disk (reaccretion of all gas that leaves
;                     disk by being heated) -------
      inflowdiskz_heat = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdiskheat_z.fits',/silent)
      inflowdiski_heat = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdiskheat_iord.fits',/silent)
      inflowdiskm_heat = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reaccrdiskheat.fits',/silent)
;   inflowdiskm_heat[where(inflowdiskm_heat LT 0)] = 0
      inflowdiski_heat_sf = inflowdiski_heat[where(inflowdiskm_heat LE 0)]
;      match,inflowdiski_heat_sf,slall.iordergas,ind1,ind2
;      inflowdiskm_heat[(where(inflowdiskm_heat LE 0))[ind1]] = slall[ind2].massform*units.massunit
      inflowdiskm_heat[where(inflowdiskm_heat LE 0)]= 0
      inflowdiskt_heat = z_to_t(inflowdiskz_heat)
      inflowmdiskcum_heat  = weighted_histogram(inflowdiskt_heat, weight =  inflowdiskm_heat,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   ENDIF

;-------------------- Inflow Disk (reaccretion of all gas that is ejected from disk by supernova
;   inflowdiskz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdisk_z.fits',/silent)
;   inflowdiski = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdisk_iord.fits',/silent)
;   inflowdiskm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reaccrdisk.fits',/silent)
;Read in the following file which only includes information about the
;particles on their reaccretion after being ejected.
   inflowdisk_data = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccrdisk_history.fits',1,/silent)
   inflowdiski_sf = inflowdisk_data[where(inflowdisk_data.mass LE 0)].iord
;   match,inflowdiski_sf,slall.iordergas,ind1,ind2
;   inflowdisk_data[(where(inflowdisk_data.mass LE 0))[ind1]].mass = slall[ind2].massform*units.massunit
   inflowdisk_data[(where(inflowdisk_data.mass LE 0))].mass = 0
   inflowdiskt = z_to_t(inflowdisk_data.red)
   inflowmdiskcum  = weighted_histogram(inflowdiskt, weight = inflowdisk_data.mass,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   inflowzdiskcum  = weighted_histogram(inflowdiskt, weight = inflowdisk_data.mass*inflowdisk_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   inflowmdiskhist  = weighted_histogram(inflowdiskt, weight = inflowdisk_data.mass,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05)
   inflowzdiskhist  = weighted_histogram(inflowdiskt, weight = inflowdisk_data.mass*inflowdisk_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05)

;Unique cooling onto the disk
;   ind_inflowdisk_uniq = uniq(inflowdisk_data.iord,sort(inflowdisk_data.iord))
;   inflowdiski_uniq = inflowdisk_data[ind_inflowdisk_uniq].iord
;   inflowdiskz_uniq = fltarr(n_elements(inflowdiski_uniq))
;   inflowdiskm_uniq = fltarr(n_elements(inflowdiski_uniq)) 
;   inflowdiskmet_uniq = fltarr(n_elements(inflowdiski_uniq))
;   FOR ind = 0, n_elements(inflowdiski_all_uniq) - 1 DO BEGIN
 ;Iterate through, saving the values that correspond to the highest redshift
;      indiord = where(inflowdiski_all EQ inflowdiski_all_uniq[ind])
;      temp = max(inflowdiskz_all[indiord],maxz_ind)
;      inflowdiskz_uniq[ind] = inflowdiskz_all[indiord[maxz_ind]]
;      inflowdiskm_uniq[ind] = inflowdiskm_all[indiord[maxz_ind]]
;      inflowdiskmet_uniq[ind] = inflowdisk_all_data[indiord[maxz_ind]].metallicity
;   ENDFOR
;   inflowdiskt_uniq = z_to_t(inflowdiskz_uniq)
;   inflowmdiskcum_uniq = weighted_histogram(inflowdiskt_uniq, weight =  inflowdiskm_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
;   inflowzdiskcum_uniq = weighted_histogram(inflowdiskt_uniq, weight =  inflowdiskm_uniq*inflowdiskmet_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)

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
  IF (where(reejectm LT 0))[0] NE -1 THEN stop ; reejectm[where(reejectm LT 0)] = 0
   reejectz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reeject_z.fits',0,/silent)
   reejectt = z_to_t(reejectz)
   reeject_data = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reeject_halo.fits',1,/silent)
   reejectmcum = weighted_histogram(reejectt, weight =  reejectm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   reejectzcum = weighted_histogram(reejectt, weight =  reejectm*reeject_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
   reejectmhist = weighted_histogram(reejectt, weight =  reejectm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05)
   reejectzhist = weighted_histogram(reejectt, weight =  reejectm*reeject_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05) 

;   IF keyword_set(debug) THEN BEGIN
;      reheati = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reheat_iord.fits',0)
;      reheatm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reheat.fits',0) 
;      reheatm[where(reheatm LT 0)] = 0
;      reheatz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reheat_z.fits',0)
;      reheatt = z_to_t(reheatz)
;      reheatmcum  = weighted_histogram(reheatt, weight =  reheatm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
;      ind_reheat_uniq = uniq(reheati,sort(reheati))
;      reheatz_uniq = reheatz[ind_reheat_uniq]
;      reheatm_uniq = reheatm[ind_reheat_uniq]
;      reheati_uniq = reheati[ind_reheat_uniq]
;      reheatt_uniq = z_to_t(reheatz_uniq)
;      reheatmcum_uniq = weighted_histogram(reheatt_uniq, weight =  reheatm_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
;   ENDIF

;---------------------- All gas lost from the disk ---------------
   relosti = mrdfits(dirs[i] + '/grp' + haloid[i] + '.relost_iord.fits',0,/silent)
   relostm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_relost.fits',0,/silent) 
   IF (where(relostm LT 0))[0] NE -1 THEN stop ;relostm[where(relostm LT 0)] = 0
   relostz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.relost_z.fits',0,/silent)
   relostt = z_to_t(relostz)
   relost_data = mrdfits(dirs[i] + '/grp' + haloid[i] + '.relost_history.fits',1,/silent)
   relostmcum  = weighted_histogram(relostt, weight =  relostm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   relostzcum = weighted_histogram(relostt, weight =  relostm*relost_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
   relostmhist  = weighted_histogram(relostt, weight =  relostm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05)
   relostzhist = weighted_histogram(relostt, weight =  relostm*relost_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,binsize = 0.05)

; ------------------ Outflow from halo -------------
      reexpelli = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reexpell_iord.fits',0,/silent)
      reexpellm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reexpell.fits',0,/silent) 
;   reexpellm[where(reexpellm LT 0)] = 0
      reexpellz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reexpell_z.fits',0,/silent)
      eject_expell_iord = mrdfits(dirs[i] + '/grp' + haloid[i] + '.eject_expell_iord.fits',/silent)
      reexpell_data = reeject_data[eject_expell_iord]
      reexpellt = z_to_t(reexpellz)
;    reexpelli_uniq = reexpelli[uniq(reexpelli,sort(reexpelli))]
      reexpellmcum  = weighted_histogram(reexpellt, weight =  reexpellm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
      reexpellzcum = weighted_histogram(reexpellt, weight =  reexpellm*reexpell_data.metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)

;Compare the metallicity at the time gas was lost to when it was expelled
      IF keyword_set(debug) THEN BEGIN
         i_reexpell_data_r = reverse(indgen(n_elements(reexpell_data),/long))
         i_relost_data_r = reverse(indgen(n_elements(relost_data),/long))
         reexpell_data_r=reexpell_data[i_reexpell_data_r]
         relost_data_r = relost_data[i_relost_data_r]
         match2,reexpell_data_r.iord,relost_data_r.iord,inda,indb
;         plot,reexpellt[i_reexpell_data_r] - relostt[i_relost_data_r[inda[where(inda ne -1)]]],reexpell_data_r.metallicity/relost_data_r[inda[where(inda ne -1)]].metallicity,psym = 3
;         contour_plus_orig,reexpellt[i_reexpell_data_r] - relostt[i_relost_data_r[inda[where(inda ne -1)]]],alog10(reexpell_data_r.metallicity/relost_data_r[inda[where(inda ne -1)]].metallicity),xbinsize = 0.2,ybinsize = 0.1,threshold = 20,nlevels = 254
         oplot,[0,14],[0,0]
         print,'Total metals expelled: ',total(reexpell_data_r.metallicity*reexpell_data_r.mass)
         print,'Total metals lost and later expelled: ',total(relost_data_r[inda[where(inda ne -1)]].metallicity*relost_data_r[inda[where(inda ne -1)]].mass)
      ENDIF

   IF keyword_set(plot_halo) THEN BEGIN
;---------------- Find all those gas particles that were reaccreted
;                 to the halo after being expelled ------
      ind_inflow = indgen(n_elements(inflowi))
      match2,ind_inflow,ind_inflow_uniq,ind1,ind2
      reinflowi = inflowi[ind_inflow(where(ind1 eq -1))] ;Find all iords that are being reaccreted for a second or more time
      reinflowz = inflowz[ind_inflow(where(ind1 eq -1))]
      reinflowt = inflowt[ind_inflow(where(ind1 eq -1))]
      reinflowm = inflowm[ind_inflow(where(ind1 eq -1))]
      reinflow_data = inflow_data[ind_inflow(where(ind1 eq -1))]
      match,reinflowi,reexpelli,ind3,ind4 ;Find iords of all gas particles that were expelled and accreted a second time (i.e., were potentiall reaccreted to the halo after being expelled
      ind_reinflow_expell=[0]
      FOR j = 0,n_elements(reexpelli[ind4]) -1 DO BEGIN
         IF (where(inflowz LT reexpellz[ind4[j]] AND inflowi EQ reexpelli[ind4[j]]))[0] NE -1 THEN $ ;If there is an instance of reaccretion after expulsion
            ind_reinflow_expell = [ind_reinflow_expell,(where(inflowz LT reexpellz[ind4[j]] AND inflowi EQ reexpelli[ind4[j]]))[0]] ;Save the index to the first one
;         IF (where(inflowz LT reexpellz[ind4[j]] AND inflowi EQ reexpelli[ind4[j]]))[0] NE -1 THEN $
;            print,reexpelli[ind4[j]],reexpellz[ind4[j]],inflowz[where(inflowi EQ reexpelli[ind4[j]])],inflowz[(where(inflowz LT reexpellz[ind4[j]] AND inflowi EQ reexpelli[ind4[j]]))[0]]
 ;        stop
      ENDFOR
;      FOR j = 0,n_elements(reexpelli[ind4]) -1 DO BEGIN
;         IF (where(reinflowz LT reexpellz[ind4[j]] AND reinflowi EQ reexpelli[ind4[j]]))[0] NE -1 THEN $ ;If there is an instance of reaccretion after expulsion
 ;           ind_reinflow_expell = [ind_reinflow_expell,(where(reinflowz LT reexpellz[ind4[j]] AND reinflowi EQ reexpelli[ind4[j]]))[0]] ;Save the index to the first one
;         print,reexpelli[ind4[j]],reexpellz[ind4[j]],reinflowz[where(reinflowi EQ reexpelli[ind4[j]])]
;         stop
;      ENDFOR

      IF n_elements(ind_reinflow_expell) NE 1 THEN BEGIN
         ind_reinflow_expell = ind_reinflow_expell[1:n_elements(ind_reinflow_expell)-1]
         reinflow_expellmcum  = weighted_histogram(inflowt[ind_reinflow_expell], weight =  inflowm[ind_reinflow_expell],  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
         reinflow_expellzcum  = weighted_histogram(inflowt[ind_reinflow_expell], weight =  inflowm[ind_reinflow_expell]*inflow_data[ind_reinflow_expell].metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
;         reinflow_expellmcum  = weighted_histogram(reinflowt[ind_reinflow_expell], weight =  reinflowm[ind_reinflow_expell],  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
;         reinflow_expellzcum  = weighted_histogram(reinflowt[ind_reinflow_expell], weight =  reinflowm[ind_reinflow_expell]*reinflow_data[ind_reinflow_expell].metallicity,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
      ENDIF ELSE BEGIN
         reinflow_expellmcum = fltarr(n_elements(timearr))
         reinflow_expellzcum = fltarr(n_elements(timearr))
      ENDELSE
   ENDIF


   IF keyword_set(plot_halo) OR keyword_set(debug) THEN BEGIN
;------------------ Outflow Material --------------------------
;Outflow gas is gas that is lost because of stripping or wierd Amiga halo definition
      print,'Read files of outflows'
      outflowm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reoutflow.fits',/silent)
;   outflowm[where(outflowm LT 0)] = 0
      outflowi = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reoutflow_iord.fits',/silent)  
      outflowz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reoutflow_z.fits',/silent)
      outflowt = z_to_t(outflowz) 
      outflow_data = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reoutflow_history.fits',1,/silent)
      outflowmcum = weighted_histogram(outflowt,weight =  outflowm, min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
      outflowzcum = weighted_histogram(outflowt,weight =  outflowm*outflow_data.metallicity, min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
; Need unique outflow halo here
      ind_outflow_uniq = rem_dup(outflowi);uniq(outflowi,sort(outflowi))
      outflowz_uniq = outflowz[ind_outflow_uniq]
      outflowm_uniq = outflowm[ind_outflow_uniq]
      outflowi_uniq = outflowi[ind_outflow_uniq]
      outflowt_uniq = z_to_t(outflowz_uniq)    
      outflowmcum_uniq = weighted_histogram(outflowt_uniq, weight =  outflowm_uniq,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05) 
   ENDIF

;Stars accreted onto the halo or any of the satellites
; ------------------- Inflow Stars ----------------------------- 
   inflowstarz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accrstars_z.fits',/silent)
   inflowstarm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accrstars_m.fits',/silent)
   inflowstari = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accrstars_i.fits',/silent)
;Since we don't current track stars leaving the halo, only look
;at the first instance of accretion (i.e., don't allow multiple accretions)
   inflowstarz = inflowstarz[uniq(inflowstari,sort(inflowstari))]
   inflowstarm = inflowstarm[uniq(inflowstari,sort(inflowstari))]
   inflowstari = inflowstari[uniq(inflowstari,sort(inflowstari))]
   match2,inflowstari,iord_star,ind1,ind2
   IF n_elements(inflowstari) NE n_elements(where(ind2 NE -1)) THEN BREAK ;Sanity check that all stars accreted still exist
   inflowstarmetal = zmetal_star[ind1]
   inflowstari = inflowstari[where(inflowstarz NE 99)]
   inflowstarm = inflowstarm[where(inflowstarz NE 99)]
   inflowstarz = inflowstarz[where(inflowstarz NE 99)]
   inflowstarmetal = inflowstarmetal[where(inflowstarz NE 99)]
   inflowstart = z_to_t(inflowstarz)
   inflowmstarcum  = weighted_histogram(inflowstart, weight =  inflowstarm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)   
   inflowzstarcum  = weighted_histogram(inflowstart, weight =  inflowstarm*inflowstarmetal,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)   
   track_stars = 1
     
;      IF keyword_set(debug) THEN BEGIN        
;   spawn,'ls ' + dirs[i] + '/*_2merge.starlog.fits',files_sl
;   IF file_test(files_sl[0]) THEN BEGIN
;      slall = mrdfits(files_sl[0],1,/silent)
;      print,files_sl[0]
;   ENDIF ELSE BEGIN
;      spawn,'ls ' + dirs[i] + '/*.starlog',files_sl
;      slall = rstarlog(files_sl[0],/molecularH)
;      print,files_sl[0]
;      IF slall EQ -1 THEN 
;   ENDELSE
;   slall = mrdfits(dirs[i] + 'starlog.' + haloid[i] + '.fits')
;   match2,slall.iorderstar,iord_star[where(inhalo_star EQ 1)],ind1,ind2
;   slgal = slall[where(ind1 NE -1)]
;   match2,slgal.iordergas,inflowi_uniq,ind3,ind4
;   accr_star_iord = slgal[where(ind3 EQ -1)].iorderstar
;8/31/17 -- I deleted where startiord_star was calcluated because it
;           was buggy. Maybe see what inhalo_star can do?
;   match2,accr_star_iord,startiord_star,ind5,ind6
;   accr_star_iord = accr_star_iord[where(ind5 EQ -1)]
;   insitu_star_iord = [slgal[where(ind3 NE -1)].iorderstar,startiord_star]
;   insitu_star_iord = insitu_star_iord[uniq(insitu_star_iord,sort(insitu_star_iord))]

;   print,total(slgal[where(ind3 EQ -1)].massform)*units.massunit,total(slgal[where(ind3 NE -1)].massform)*units.massunit,total(slgal.massform)*units.massunit,mstar[nhalostat - 1],format='(E,E,E,E)'
;   accrstarcum = weighted_histogram(slgal[where(ind3 EQ -1)].timeform*units.timeunit/1e9, weight = slgal[where(ind3 EQ -1)].massform*units.massunit,min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
;   insitustarcum = weighted_histogram(slgal[where(ind3 NE -1)].timeform*units.timeunit/1e9, weight = slgal[where(ind3 NE -1)].massform*units.massunit,min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
;   sgal_igasord = slgal[where(ind3 NE -1)].iordergas

;         match2,slgal.iorderstar,inflowstari,ind12,ind13
;         insitu_star_iord = slgal[where(ind12 EQ -1)].iorderstar
;         insitustarcum = weighted_histogram(slgal[where(ind12 NE -1)].timeform*units.timeunit/1e9, weight = slgal[where(ind12 NE -1)].massform*units.massunit,min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum,binsize = 0.05)
 
 ;  ENDIF ELSE BEGIN
 ;     track_stars = 0
;       mstarform = spline(timearr,formstarcum,time)
 ;     inflowmstarcum = fltarr(n_elements(timearr))
 ;  ENDELSE


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

   IF i EQ 7 OR i EQ 11 THEN BEGIN ;h516
      mint = z_to_t(3.130)      ;Tracing only is robust back to z= 3
      temp = min(abs(timearr -(mint)[0]),ind_mint_timearr)
      temp = min(abs(time -(mint)[0]),ind_mint_time)
      miss_disk_out = (zstars_cum[ind_mint_time] + inflowzdiskcum_all[ind_mint_timearr] - relostzcum[ind_mint_timearr]) - zmetals[ind_mint_time]
      IF keyword_set(plot_halo) THEN miss_halo_out = (zstars_cum[ind_mint_time] + inflowzcum[ind_mint_timearr] - outflowzcum[ind_mint_time]) - zmetals_ha[ind_mint_time]
   ENDIF ELSE BEGIN
      IF i EQ 9 OR i EQ 14 OR i EQ 16 THEN BEGIN ;603
         mint = z_to_t(3.436)   ;Tracing only is robust back to z= 3
         temp = min(abs(timearr -(mint)[0]),ind_mint_timearr)
         temp = min(abs(time -(mint)[0]),ind_mint_time)
         miss_disk_out = (zstars_cum[ind_mint_time] + inflowzdiskcum_all[ind_mint_timearr] - relostzcum[ind_mint_timearr]) - zmetals[ind_mint_time]
         IF keyword_set(plot_halo) THEN miss_halo_out = (zstars_cum[ind_mint_time] + inflowzcum[ind_mint_timearr] - outflowzcum[ind_mint_time]) - zmetals_ha[ind_mint_time]
      ENDIF ELSE BEGIN
         miss_disk_out = 0  ;zstars_cum[ind_mint_time] - zmetals[ind_mint_time]
         miss_halo_out = 0  ;zstars_cum[ind_mint_time] - zmetals_ha[ind_mint_time]
         ind_mint_timearr = 0
         ind_mint_time = 0
      ENDELSE
   ENDELSE
   indtarr = indgen(n_elements(timearr) - ind_mint_timearr) + ind_mint_timearr; ind_mint_timearr :n_elements(timearr)-1
;***************************** Start Plotting ********************
;***************************** History of metals in disk ***********
   IF NOT keyword_set(pmulti) THEN $
      IF (keyword_set(outplot)) THEN device,filename = outplot + '_' + strtrim(i,2) + '_fbzmass.eps',/color,bits_per_pixel= 8,/times,ysize = 16,xsize = 18,xoffset =  2,yoffset =  2 ELSE window,3,xsize = 712,ysize = 600  
   plot,timearr,inflowzdiskcum_all,xrange=[0.1,ageUniverse/1e9],xtitle = xtitle,/nodata,xmargin = [17,17],ymargin = [6,6],ystyle = 1,xcharsize = 1,ycharsize = 1,xstyle = xstyle,yrange = range[*,i/4];yrange = [-1*max([inflowzdiskcum_all,relostzcum]) - miss_disk_out,max([inflowzdiskcum_all,relostzcum])]
   IF NOT keyword_set(pmulti) THEN axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = n_elements(ticktime_in_t) - 1,ytitle = 'Metal Mass [M' + sunsymbol() + ']'
   oplot,[0,14],[0,0]
   yscale = max(zstars_cum)
; Total accretion onto the disk (pristine + reaccretion) (blue solid)
;   oplot,timearr,inflowzdiskcum_all,color = cheat,linestyle = 0,thick = thicks[i]
; "Pristine" accretion onto the disk (purple dot-dashed)
;   oplot,timearr,inflowzdiskcum_all_uniq,color = cpristine,linestyle = 3,thick = thicks[i]
;   oplot,timearr,inflowzstarcum,color = cpristine,linestyle = 4,thick = thicks[i]  
   oplot,timearr,(inflowzstarcum+inflowzdiskcum_all_uniq)/yscale,color = cpristine,linestyle = 3,thick = thicks[i]
; Reaccretion of gas onto disk (blue dashed)
   oplot,timearr,(inflowzdiskcum_all - inflowzdiskcum_all_uniq)/yscale,color = cheat,thick = thicks[i] ,linestyle = 0
; Total metals accreted from ejected gas (red solid)
   oplot,timearr,inflowzdiskcum/yscale,color = ceject,linestyle = 0,thick = thicks[i]
; Metals reaccreted from ejected gas (red dashed)
;   oplot,timearr,inflowzdiskcum,color = 254,thick = thicks[i],linestyle = 2
; Total metals (re)lost from disk (blue solid)
   oplot,timearr[indtarr],(-1*relostzcum[indtarr] - miss_disk_out)/yscale,color = cheat,thick = thicks[i]
; Metals (re)ejected from disk (red solid)
   oplot,timearr[indtarr],(-1*reejectzcum[indtarr] - miss_disk_out)/yscale,color = ceject,thick = thicks[i]
; Net metals (re)ejected from disk (red dotted)
;   oplot,timearr[indtarr],inflowzdiskcum[indtarr]-reejectzcum[indtarr] - miss_disk_out,color = ceject,linestyle = 2,thick = thicks[i]
; Net metals (re)lost from disk (blue dotted): inflow - pristine inflow
   oplot,timearr[indtarr],((inflowzdiskcum_all[indtarr]-inflowzstarcum[indtarr] - inflowzdiskcum_all_uniq[indtarr])-relostzcum[indtarr]- miss_disk_out)/yscale,color = cheat,linestyle = 2,thick = thicks[i]

;   oplot,time,zmetals,thick = thicks[i],linestyle = 2 ;disk metals (gas)
;   oplot,timearr,formstarzcum,thick = thicks[i],linestyle = 2
;   oplot,time,zstars_in,thick = thicks[i],linestyle = 1 ;disk metals (stars)
   oplot,time,(zmetals + zstars_in)/yscale,thick = thicks[i] ;disk metals total
   oplot,time,zstars_cum/yscale,color = cstar,thick = thicks[i],linestyle = 5 ;Metals produced by stars
   zstars_cum_spline = spline(time[uniq(time)],zstars_cum[uniq(time)],timearr)
   zstars_cum_spline[where(timearr lt min(time))] = 0
   ;Metals produced + metals (re)accreted + metals lost
   ;Metals produced (zstar_cum_spline) + metals accreted as stars (inflowzstarcum)+ accretion to disk (inflowzdiskcum_all)- lost from the disk (relostzcum)
   IF keyword_set(debug) THEN oplot,timearr[indtarr],(zstars_cum_spline[indtarr] + inflowzstarcum[indtarr]+ inflowzdiskcum_all[indtarr] - relostzcum[indtarr] - miss_disk_out)/yscale,linestyle = 3,thick = thicks[i],color = csim
   ;Metals produced (zstars_cum_spline) + accretion to the disk (inflowzdiskcum_all)
;   oplot,timearr[indtarr],zstar_cum_spline[indtarr] + inflowzdiskcum_all[indtarr] - inflowzdiskcum_all_uniq[indtarr] - relostzcum[indtarr] - miss_disk_out,linestyle = 3,thick = thicks[i],color = 100
; Metals expelled from the halo (blue solid)
   oplot,timearr,-1*reexpellzcum/yscale,color = cexpell,thick = thicks[i]
; Net metals expelled from halo, i.e. expelled - reaccreted to halo (blue dotted)
;      oplot,timearr,-1*reexpellzcum + reinflow_expellzcum,color = 60,linestyle = 2,thick = thicks[i]

   IF keyword_set(label) THEN legend,[label[i]],box = 0,/right,/top
;   IF i EQ 0 OR NOT keyword_set(pmulti) THEN legend,['Total','Initial','Reaccreted','Remain out'],color = [fgcolor,fgcolor,fgcolor,fgcolor],linestyle = [0,3,2,1],box = 0,/left,/top,thick = [thicks[i],thicks[i],thicks[i],thicks[i]]
;   IF i EQ 0 OR NOT keyword_set(pmulti) THEN legend,['Total','Initial','Remain out'],color = [fgcolor,fgcolor,fgcolor],linestyle = [0,3,2],box = 0,/left,/top,thick = [thicks[i],thicks[i],thicks[i]]
;   IF i EQ 0 OR NOT keyword_set(pmulti) THEN legend,['ISM + Stars','First Accretion','Produced','Ejected/Reaccreted','Heated/Reaccreted'],color = [fgcolor,20,210,254,60],linestyle = [0,0,0,0,0],box = 0,/left,/bottom,thick = [thicks[i],thicks[i],thicks[i],thicks[i],thicks[i]]
;   IF i EQ 0 OR NOT keyword_set(pmulti) THEN legend,['Disk (solid),
;   ISM (dashed), Stars (dotted)','First
;   Accretion','Produced','Expelled Beyond
;   Rvir','Heated/Reaccreted'],color =
;   [fgcolor,cpristine,cstar,cexpell,cheat],linestyle =
;   [0,0,0,0,0],box = 0,/left,/bottom,thick =
;   [thicks[i],thicks[i],thicks[i],thicks[i],thicks[i]]
   IF i EQ 0 OR NOT keyword_set(pmulti) THEN legend,['Disk','First accretion','Produced'],color = [fgcolor,cpristine,cstar],linestyle = [0,3,5],box = 0,/left,/top,thick = [thicks[i],thicks[i],thicks[i]]
   IF i EQ 0 OR NOT keyword_set(pmulti) THEN legend,[textoidl('Expelled beyond R_{vir}'),'Outflow/Reaccreted','Removed from disk/Reaccreted','Net removed from disk'],color = [cexpell,ceject,cheat,cheat],linestyle = [0,0,0,2],box = 0,/left,/bottom,thick = [thicks[i],thicks[i],thicks[i],thicks[i]]

   IF i MOD 4 EQ 3 THEN print,'Max/min: ',max(inflowzdiskcum_all - inflowzdiskcum_all_uniq),min(-1*relostzcum[indtarr] - miss_disk_out),max(inflowzdiskcum_all - inflowzdiskcum_all_uniq),min(-1*relostzcum[indtarr] - miss_disk_out)/yscale
 
  IF keyword_set(pmulti) THEN multiplot ELSE $
      IF keyword_set(outplot) THEN device,/close ;ELSE stop

;************** Plot History of Halo Metal Content *************************
   IF NOT keyword_set(pmulti) AND keyword_set(plot_halo) THEN BEGIN
      IF (keyword_set(outplot)) THEN device,filename = outplot + '_' + strtrim(i,2) + '_fbzmass_halo.eps',/color,bits_per_pixel= 8,/times,ysize = 16,xsize = 18,xoffset =  2,yoffset =  2 ELSE window,7,xsize = 712,ysize = 600  
      plot,timearr,inflowzcum,ytitle = 'Metal Mass [M' + sunsymbol() + ']',yrange = [-1*max([reexpellzcum,outflowzcum])-miss_halo_out,max([zstars_cum,zstars_cum_spline + (inflowzcum - inflowzcum_uniq) - reexpellzcum,zmetals_ha])],xrange=[0.1,ageUniverse/1e9],xtitle = xtitle,/nodata,xmargin = [17,17],ymargin = [6,6],ystyle = 1,xcharsize = 2,ycharsize = 2,xstyle = xstyle
      axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = n_elements(ticktime_in_t) - 1
      oplot,[0,14],[0,0]

      oplot,time,zmetals_ha + zstars_in,thick = thicks[i] ;Total Halo metals
      oplot,time,zstars_in,thick = thicks[i],linestyle = 1 ;halo metals (stars)
      oplot,time,zmetals_ha,thick = thicks[i],linestyle = 2 ;halo metals (gas)

      oplot,time,zstars_cum,color = cstar,thick = thicks[i] ; Metals produced by stars
; Total Metals (re)accreted onto the halo from both gas and stars (pea green solid)
      oplot,timearr,inflowzcum+inflowzstarcum,color = cheat_halo,thick = thicks[i]
; "Pristine" and stellar accretion onto the halo (purple/blue dot-dashed)
      oplot,timearr,inflowzcum_uniq+inflowzstarcum,color = cpristine,thick = thicks[i],linestyle = 3
; Reaccretion of metals onto the halo (pea green dashed)
      oplot,timearr,inflowzcum - inflowzcum_uniq,color = cheat_halo,thick = thicks[i],linestyle = 2
; Reaccretion of metals after they were expelled from the virial radius (blue solid)
      oplot,timearr,reinflow_expellzcum,color = cexpell,thick = thicks[i],linestyle = linestyle[2]

; Metals outflow from the halo (pea green solid)
      oplot,timearr[indtarr],-1*outflowzcum[indtarr] - miss_halo_out,color = cheat_halo,thick = thicks[i]
; Net Metals outflowing from the halo (pea green dashed)
      oplot,timearr[indtarr],-1*outflowzcum[indtarr] - miss_halo_out + inflowzcum[indtarr] - inflowzcum_uniq[indtarr],color = cheat_halo,linestyle = 2,thick = thicks[i]
; Metals expelled from the halo (blue solid)
      oplot,timearr,-1*reexpellzcum,color = cexpell,thick = thicks[i]
; Net metals expelled (blue dashed)
      oplot,timearr,-1*reexpellzcum + reinflow_expellzcum,color = cexpell,linestyle = 2,thick = thicks[i]
; Metals produced and accreted minus net metals expelled (cyan dot-dashed) -- does not include gas that leaves the halo that was never in the disk) 
;Metals produced (zstars_cum_spline) + Metals accreted through pristine accretion (inflowzcum_uniq) + Metals reaccreted as stars (inflowzstarcum) - Net metals expelled (reexpellzcum - reinflow_expellzcum)
; Expelled material must once have been in the disk so net metals
; expelled is != net metals lost from the halo
;      oplot,timearr[indtarr],$
;            zstars_cum_spline[indtarr] + inflowzcum_uniq[indtarr] + inflowzstarcum[indtarr] + $
;            reinflow_expellzcum[indtarr] - reexpellzcum[indtarr] - miss_halo_out,thick = thicks[i],linestyle = 3,color = 100

;Metals produced and accreted minus net metals outflowing (cyan dot-dashed)
;Metals produced + inflow - outflow - start
      oplot,timearr[indtarr],$
            zstars_cum_spline[indtarr] + inflowzcum[indtarr]     + inflowzstarcum[indtarr] - $
            outflowzcum[indtarr] - miss_halo_out,thick = thicks[i],linestyle = 3,color = csim

      IF keyword_set(label) THEN legend,[label[i]],box = 0,/right,/top
      IF i EQ 0 OR NOT keyword_set(pmulti) THEN legend,['Total','Initial','Reaccreted','Remain out'],color = [fgcolor,fgcolor,fgcolor,fgcolor],linestyle = [0,3,2,1],box = 0,/left,/top,thick = [thicks[i],thicks[i],thicks[i],thicks[i]]
      IF i EQ 0 OR NOT keyword_set(pmulti) THEN legend,['Metals in halo','[Sim] Metals in halo','First accretion','Metals produced','Lost/reaccreted to halo','Expelled/reaccreted to halo'],thick = [thicks[i],thicks[i],thicks[i],thicks[i],thicks[i],thicks[i]],color = [fgcolor,csim,cpristine,cstar,cheat_halo,cexpell],linestyle = [0,0,3,0,0,0],box = 0,/left,/bottom
      IF keyword_set(outplot) THEN device,/close ;ELSE stop
   ENDIF

   IF NOT keyword_set(pmulti) THEN BEGIN
;***** Plot average metallicity of gas being accreted or ejected as a function of time
      IF (keyword_set(outplot)) THEN device,filename = outplot + '_' + strtrim(i,2) + '_fbz.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 ELSE window,2,xsize = 600,ysize = 400
      maxplot = max(reejectzhist[where(reejectmhist NE 0)]/reejectmhist[where(reejectmhist NE 0)])
      plot,timearr[where(reejectmhist NE 0)],reejectzhist[where(reejectmhist NE 0)]/reejectmhist[where(reejectmhist NE 0)]/zsolar,ytitle = 'Metallicity [Z/Z' + sunsymbol() + ']',xrange=[0.1,ageUniverse/1e9],yrange = [0,maxplot/zsolar],xtitle = xtitle,/nodata;,xmargin = [17,17],ymargin = [6,6],ystyle = 1,xcharsize = 1,ycharsize = 1;,xstyle = 9
 ;     axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = n_elements(ticktime_in_t) - 1
;Metallicity of gas during initial external accretion
;    oplot,timearr[where(inflowmdiskhist_all_uniq NE 0)],inflowzdiskhist_all_uniq[where(inflowmdiskhist_all_uniq NE 0)]/inflowmdiskhist_all_uniq[where(inflowmdiskhist_all_uniq NE 0)]/zsolar,linestyle = 3,color = 20,thick = thicks[i]
;Metallicity of all gas accreted to the disk
;    oplot,timearr[where(inflowmdiskhist_all NE 0)],inflowzdiskhist_all[where(inflowmdiskhist_all NE 0)]/inflowmdiskhist_all[where(inflowmdiskhist_all NE 0)],color = 50,linestyle = 2
;Metallicity of all gas that returns to the disk after being lost (i.e., I'm ignoring the metallicity of the gas particles accreted for the first time.)
      inflowmdiskhist_all_return = inflowmdiskhist_all - inflowmdiskhist_all_uniq
      inflowzdiskhist_all_return = inflowzdiskhist_all - inflowzdiskhist_all_uniq_disk
      oplot,timearr[where(inflowmdiskhist_all NE 0)],inflowzdiskhist_all[where(inflowmdiskhist_all NE 0)]/inflowmdiskhist_all[where(inflowmdiskhist_all NE 0)]/zsolar,linestyle = 5,color = cpristine,thick = thicks[i]
      oplot,timearr[where(inflowmdiskhist_all_return NE 0)],inflowzdiskhist_all_return[where(inflowmdiskhist_all_return NE 0)]/inflowmdiskhist_all_return[where(inflowmdiskhist_all_return NE 0)]/zsolar,linestyle = 2,color = cheat,thick = thicks[i]
;Metallicity of all gas removed from the disk
      oplot,timearr[where(relostmhist NE 0)],relostzhist[where(relostmhist NE 0)]/relostmhist[where(relostmhist NE 0)]/zsolar,color = cheat,thick = thicks[i]
;    oplot,timearr[where(inflowmdiskhist NE 0)],inflowzdiskhist[where(inflowmdiskhist NE 0)]/inflowmdiskhist[where(inflowmdiskhist NE 0)],color = 254,linestyle = 2
;Metallicity of all gas reaccreted after being ejected
      oplot,timearr[where(inflowmdiskhist NE 0)],inflowzdiskhist[where(inflowmdiskhist NE 0)]/inflowmdiskhist[where(inflowmdiskhist NE 0)]/zsolar,linestyle = 2,color = ceject,thick = thicks[i]
;Metallicity of the ejected gas
      oplot,timearr[where(reejectmhist NE 0)],reejectzhist[where(reejectmhist NE 0)]/reejectmhist[where(reejectmhist NE 0)]/zsolar,color = ceject,thick = thicks[i]
      oplot,time[where(coldg NE 0)],zmetals[where(coldg NE 0)]/coldg[where(coldg NE 0)]/zsolar,thick = thicks[i]
      IF i EQ 1 THEN legend,['Disk Gas','Accreted','Ejected/Reaccreted','Lost/Reaccreted'],linestyle = [0,5,0,0],color = [fgcolor,cpristine,cheat,ceject],/left,/top,box=0,thick = [thicks[i],thicks[i],thicks[i],thicks[i]]
;      legend,['Disk
;      Gas','Ejected/Reaccreted','Lost/Reaccreted','Initial
;      Accretion'],linestyle = [0,0,0,3],color =
;      [fgcolor,85,254,20],/left,/top,box=0,thick =
;      [thicks[i],thicks[i],thicks[i],thicks[i]]
      IF i EQ 1 THEN $
        IF keyword_set(label) THEN legend,[label[i]],box = 0,/right,/top
      IF i NE 1 THEN $
        IF keyword_set(label) THEN legend,[label[i]],box = 0,/right,/bottom

;Calculate the average metallicity enhancement of the outflow
      ;Only complete analysis since z = 3.5
      reejectmhist_tcut = reejectmhist[where(timearr GE 1.89322)]
      reejectzhist_tcut = reejectzhist[where(timearr GE 1.89322)]
      timearr_tcut = timearr[where(timearr GE 1.89322)]
      diskmetallicity = zmetals[where(coldg NE 0)]/coldg[where(coldg NE 0)]
      diskmetallicity_spline = spline(time[where(coldg NE 0)],diskmetallicity,timearr_tcut[where(reejectmhist_tcut NE 0)])
      metallicityenhance = (reejectzhist_tcut[where(reejectmhist_tcut NE 0)]/reejectmhist_tcut[where(reejectmhist_tcut NE 0)])/diskmetallicity_spline
      avemetallicityenhance = total(metallicityenhance*reejectmhist_tcut[where(reejectmhist_tcut NE 0)])/total(reejectmhist_tcut[where(reejectmhist_tcut NE 0)])
      print,'Average metallicity enhancement of ejecta: ',avemetallicityenhance
;      oplot,timearr_tcut[where(reejectmhist_tcut NE 0)],diskmetallicity_spline/zsolar*avemetallicityenhance,color = csim,thick = thicks[i]
      IF keyword_set(outplot) THEN device,/close ELSE stop
   ENDIF
      
   IF NOT keyword_set(pmulti) AND keyword_set(debug) THEN BEGIN
;*********** Plot total mass accreted/expelled from disk. Used here for
;                debugging. See inflow_outflow_history.pro for a
;                dedicated program -----------------------------------
      IF (keyword_set(outplot)) THEN device,filename = outplot + '_' + strtrim(i,2) + '_fbmass.eps',/color,bits_per_pixel= 8,/times,ysize = 16,xsize = 18,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 712,ysize = 600
      IF keyword_set(yrange_max) THEN $
         IF keyword_set(allpositive) THEN yrange = [0,yrange_max[i]] ELSE yrange = [-1*yrange_max[i],yrange_max[i]] ELSE $
            IF keyword_set(allpositive) THEN yrange = [0,[max([mgas+mstar,(inflowmdiskcum + inflowmstarcum),inflowmdiskcum_all])]]*scale[i] ELSE  yrange = [-1*max([outflowmcum,relostmcum]),max([inflowmcum,inflowmdiskcum_all])]*scale[i]
      print,i,mtot[n_elements(mtot) - 1],yrange[1],yrange[1]/mtot[n_elements(mtot) - 1]

      plot,time,(mgas+mstar)*scale[i],ytitle = ytitle,xstyle = 9,xrange=[0.1,ageUniverse/1e9],xtitle = 'Mass [M'+ sunsymbol() + ']',/nodata,xmargin = [17,17],ymargin = [6,6],yrange = yrange,ystyle = 1,xcharsize = 1,ycharsize = 1
      axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = n_elements(ticktime_in_t) - 1
      IF NOT keyword_set(allpositive) THEN oplot,[0,14],[0,0],thick = 2
      
;----------------- Mass within halo --------------
      IF 1 THEN BEGIN
         oplot,time,(mgas+mstar)*scale[i],                            linestyle = 0,thick = thicks[i] ;,color = colors[i] ;total baryons within rvir
         oplot,timearr,(inflowmcum+inflowzstarcum)*scale[i],                linestyle = 0,thick = thicks[i],color = cheat_halo ;colors[i] ;,psym = -4 ;gas mass accreted into the virial mass (reaccr)
         oplot,timearr,(inflowmcum - inflowmcum_uniq)*scale[i],linestyle = 2,thick = thicks[i],color = cheat_halo ;Just gas reaccreted onto the halo
         oplot,timearr,inflowmcum_uniq*scale[i],                linestyle = 1,thick = thicks[i],color = cpristine ;colors[i] ;gas mass accreted into the virial mass -- once (reaccr)
         oplot,timearr,e_m*outflowmcum*scale[i],          linestyle = 0,thick = thicks[i],color = cheat_halo ;colors[i] ;gas mass out (reoutflow_z.fits, mass_at_reoutflow.fits)
         oplot,timearr,(inflowmcum - outflowmcum + inflowmstarcum)*scale[i],linestyle = 3,thick = thicks[i],color = csim ;colors[i] ;gas accreted minus gas out
      ENDIF
      
;==================================
;   IF i EQ 0 THEN BEGIN
;   yrange = [-1*max(relostmcum*scale[i]),max(inflowmdiskcum_all*scale[i])]
;      plot,time,mgas+mstar,ytitle = 'Baryonic Mass [M' + sunsymbol() + ']',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',/nodata,xmargin = [17,17],ymargin = [6,6],yrange = yrange
;        axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
;      axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = n_elements(ticktime_in_t) - 1        
;      oplot,[0,14],[0,0],thick = 2
;   ENDIF

;--------------------- Disk --------------
      oplot,timearr,e_m*relostmcum*scale[i],           linestyle = 0,thick = thicks[i],color = cheat ;colors[i] + 20
;      ghalo = g[where(inhalo_gas NE 0)]
;      IF keyword_set(gmass) THEN mdisk = gmass[i] ELSE mdisk = total(g[where(ghalo.dens*units.rhounit GE 0.1 AND ghalo.tempg LT 1.2e4)].mass)*units.massunit
;     oplot,time,(mgas+mstar)*scale[i],               linestyle = 5,thick = thicks[i] ;,color = colors[i] ;total baryons within rvir
     oplot,time,(coldg + mstar)*scale[i],           linestyle = 0,thick = thicks[i] ;,color = colors[i] + 20 ;baryonic mass
;     oplot,[max(time),max(time)],[mdisk + mstar[n_elements(mstar) -1],mdisk + mstar[n_elements(mstar) -1]]*scale[i],psym = 7,symsize = 2
      oplot,timearr,inflowmdiskcum_all*scale[i],      linestyle = 0,thick = thicks[i],color = cheat ;colors[i] + 20
      oplot,timearr,inflowmdiskcum_all_uniq*scale[i], linestyle = 3,thick = thicks[i],color = cpristine;colors[i] + 20
;      oplot,timearr,(inflowmdiskcum_all  - relostmcum + inflowmstarcum)*scale[i],          linestyle = 2,thick = thicks[i],color = 225 ;colors[i] + 20
;      oplot,timearr,inflowmdiskcum_heat*scale[i],     linestyle = 0,thick = thicks[i],color = 210 ;colors[i] + 50
;      oplot,timearr,inflowmdiskcum_heat_uniq*scale[i],linestyle = 1,thick = thicks[i],color = 210;colors[i] + 50
;      oplot,timearr,-1*reheatmcum*scale[i],           linestyle = 0,thick = thicks[i],color = 210 ;colors[i] + 50
;      oplot,timearr,(inflowmdiskcum_heat - reheatmcum + inflowmstarcum)*scale[i],linestyle = 2,thick = thicks[i],color = 210 ;colors[i] + 50
      oplot,timearr,(inflowmstarcum + inflowmdiskcum_all - relostmcum)*scale[i],color = csim,linestyle = 3,thick = thicks[i]
;--------------------- Ejecta and Expelled --------------
                                ;   oplot,timearr,(inflowmdiskcum + inflowmstarcum)*scale[i],            linestyle = 0,thick = thicks[i],color = 60;254 ;gas mass into disk 'accrdisk_z.fits'
;    oplot,timearr,(inflowmdiskcum_uniq + inflowmstarcum)*scale[i],       linestyle = 1,thick = thicks[i],color = 60;254 ;gas mass into disk -- once 'accrdisk_z.fits'
      oplot,timearr,inflowmdiskcum*scale[i],linestyle = 0,thick = thicks[i],color = ceject
      oplot,timearr,e_m*reejectmcum*scale[i],          linestyle = 0,thick = thicks[i],color = ceject ;gas mass ejected form disk (multiple times) 'reeject_z'
      oplot,timearr,e_m*reexpellmcum*scale[i],         linestyle = 0,thick = thicks[i],color = cexpell ;gas mass expelled from disk by supernova 'reexpell_z'
;      oplot,timearr,(inflowmdiskcum - reejectmcum + inflowmstarcum)*scale[i],linestyle = 5,thick = thicks[i],color = 60 ;colors[i] + 100 ;gas + stars in - gas out
      
;------------------ Stars ----------------------
;Note that there will be a difference between the current mass of stars and the mass of stars formed
;     oplot,time,mstar*scale[i],                               linestyle = 0,thick = thicks[i],color = 100
;     oplot,timearr,(insitustarcum  + inflowmstarcum)*scale[i],linestyle = 2,thick = thicks[i],color = 100 ;colors[i]+150
;    IF track_stars THEN oplot,timearr,(inflowmdiskcum + inflowmstarcum - reejectmcum)*scale[i],linestyle = 0,thick = thicks[i],color = colors[i]+50 ;gas + stars in to disk - gas out of disk (need to add something that will show multiple accretions onto disk)
       
;Baryonic mass will be different from inflow - outflow because of gas that is reaccreted later
;    oplot,time,(coldg + mstar)*scale[i],                linestyle = 1,thick = thicks[i],color = colors[i]+100 ;baryonic disk mass. Possibly different because of stripping, lack of rediskaccretion, and different definitions of cold gas?
;    oplot,time,mgas*scale[i],linestyle = 2,thicks = 2,color = colors[i]
;    oplot,time,mstar*scale[i],linestyle = 2,thicks = 2,color = colors[i]
;    oplot,time,mstar*scale[i] + mgas*scale[i],linestyle = 0,thicks = 2,color = colors[i]
;    oplot,time,cuminflow*scale[i] - (mstar*scale[i] + mgas*scale[i]),linestyle = 1,thicks = 2,color = colors[i]
;    IF i EQ 0 THEN outflowratio0 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]] ELSE outflowratio1 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]]
      IF keyword_set(label) THEN legend,[label[i]],box = 0,/left,/top
      legend,['Gas lost/reaccreted to disk','Gas ejected/reaccreted to disk','Gas expelled','First accretion to disk','Mass in disk'],color = [cheat,ceject,cexpell,cpristine,fgcolor],linestyle = [0,0,0,3,0],/bottom,/left,box = 0 ,thick = [thicks[i],thicks[i],thicks[i],thicks[i],thicks[i]]
;*     legend,['Gas lost/reaccreted to disk','Gas ejected/reaccreted to disk','Pristine Accretion','Mass in disk'],color = [60,254,20,fgcolor],linestyle = [0,0,2,0],/bottom,/left,box = 0 ,thick = [thicks[i],thicks[i],thicks[i],thicks[i]];[0.293480,0.95];,position = [385,90]
;    print,'Position: ',position 
;      IF debug THEN legend,['mgas+mstar','inflow','inflow, uniq','in-out','in+in','relost','mdisk','in disk','in-out disk','heat','reheat','in-heat','rereject','reexpell','in-reeject','star','instar-outstar'],color = [fgcolor,185,185,185,185,fgcolor,110,fgcolor,225,210,210,210,254,50,60,100,100],linestyle = [5,0,1,0,0,0,0,5,0,0,0,2,0,0,5,0,2],/left,/bottom,box = 0
      IF keyword_set(outplot) THEN device,/close ELSE stop
  ENDIF

;****************** Debugging by comparing particles accreted/ejected
;                  to those in final halo *******************
   IF keyword_set(debug) AND NOT keyword_set(pmulti) THEN BEGIN
;Confirm that there is no overlap between stars accreted to halo and
;stars formed in it (I think this could happen if a star is lost from halo,
;then reaccreted)
       match,inflowstari,sl.iorderstar,ind1,ind2
       IF ind1[0] NE -1 OR ind2[0] NE -1 THEN BEGIN
           print,'ERROR: ',n_elements(ind1),' stars are both accreted and formed within the halo'
           stop
       ENDIF
;Check that all stars that are in the halo at z = 0 either formed
;within it or were accreted onto it
       star_form_inflow_iord = [sl.iorderstar,inflowstari]
       star_form_inflow_iord = star_form_inflow_iord[uniq(star_form_inflow_iord,sort(star_form_inflow_iord))]
       match2,iord_star[where(inhalo_star EQ 1)],star_form_inflow_iord,ind1,ind2
       IF (where(ind1 EQ -1))[0] NE -1 THEN print,'Potential problem: ',strtrim(n_elements(where(ind1 EQ -1)),2),' stars in z = 0 halo neither accreted nor formed in halo. Equates to ',strtrim(total(s[where(ind1 EQ -1)].mass*units.massunit*s[where(ind1 EQ -1)].metals),2),' Msol metal mass'
       IF (where(ind2 EQ -1))[0] NE -1 THEN print,'Potential problem: ',strtrim(n_elements(where(ind2 EQ -1)),2),' accreted/formed stars not in halo at z = 0. Equates to ',strtrim(total(s[where(ind1 EQ -1)].mass*units.massunit*s[where(ind1 EQ -1)].metals),2),' Msol metal mass'

;Check that all gas currently in the halo was accreted and listed in *.reaccr_iord.fits
       match2,inflowi_uniq,iord_gas[where(inhalo_gas EQ 1)],ind5,ind6
;ind6 EQ -1 should be the empty array if all gas in the halo at z = 0
;appears in *.reaccr_iord.fits
       IF (where(ind6 EQ -1))[0] NE -1 THEN BEGIN
           print,'ERROR: Z = 0 halo has particles that do not appear in *.reaccr_iord.fits'
           stop
       ENDIF

;Calculate uniq particles that are accreted onto halo, finding the instance when
;they last accrete the halo (contrast wtih inflowi_uniq, which find
;first accretion)
       ind_inflow_uniq_r = rem_dup(inflowi,inflowt) ;uniq(inflowi,sort(inflowi))
       inflowz_uniq_r = inflowz[ind_inflow_uniq_r]
       inflowm_uniq_r = inflowm[ind_inflow_uniq_r]
       inflowi_uniq_r = inflowi[ind_inflow_uniq_r]
       inflowt_uniq_r = inflowt[ind_inflow_uniq_r] 
       inflow_data_uniq_r = inflow_data[ind_inflow_uniq_r]
       
;Calculate uniq particles that leave halo, finding the instance when
;they last leave the halo
       ind_outflow_uniq_r = rem_dup(outflowi,outflowt) ;uniq(outflowi,sort(outflowi))
       outflowz_uniq_r = outflowz[ind_outflow_uniq_r]
       outflowm_uniq_r = outflowm[ind_outflow_uniq_r]
       outflowi_uniq_r = outflowi[ind_outflow_uniq_r]
       outflowt_uniq_r = outflowt[ind_outflow_uniq_r]
       outflow_data_uniq_r = outflow_data[ind_outflow_uniq_r]
       
       match,outflowi_uniq_r,inflowi_uniq_r,ind7,ind8
       IF n_elements(ind7) NE n_elements(outflowi_uniq_r) THEN BEGIN
           print,'ERROR: some outflowing from halo particles are not inflowing to halo'
           stop
       ENDIF
       timediff = outflowt_uniq_r[ind7] - inflowt_uniq_r[ind8]
       print,'Min/Max time between last halo outflow and last halo accretion (+ = out, - = in): ',strtrim(minmax(timediff),2)
       window,4
;       plot,inflowt_uniq_r[ind8],outflowt_uniq_r[ind7],psym = 3,xtitle = 'Time at last inflow',ytitle = 'Time at last outflow'

;Iords of particles that were once accreted to the halo but removed by z = 0
       outhalo_z0i = outflowi_uniq_r[ind7(where(timediff GT 0))]
       match,outhalo_z0i,iord_gas,ind9,ind10
       IF (where(inhalo_gas[ind10] EQ 1))[0] NE -1 THEN BEGIN ;Should be the empty array
           print,'ERROR: gas particles thought to have left halo are still in it'
           stop
       ENDIF

;Find an array of all gas particles in halo at z = 0 along with all
;gas particles that spawned stars in the halo and were deleted
       match2,iord_gas,sl.iordergas,ind9a,ind10a
       IF (where(ind10a EQ -1))[0] NE -1 THEN BEGIN
           sl_deleted = sl[where(ind10a EQ -1)]
;Removed duplicates because the same gas particle can form multiple stars 
           sl_deleted = sl_deleted[uniq(sl_deleted.iordergas,sort(sl_deleted.iordergas))]
           iord_gas_inhalo = [iord_gas[where(inhalo_gas EQ 1)],sl_deleted.iordergas] 
       ENDIF ELSE iord_gas_inhalo = iord_gas[where(inhalo_gas EQ 1)]
;Compare with array of all gas particles accreted and never removed or
;accreted after being removed
       match2,outflowi_uniq_r,inflowi_uniq_r,ind11,ind12
       inhalo_z0i = [outflowi_uniq_r[ind7(where(timediff LT 0))],inflowi_uniq_r[where(ind12 EQ -1)]]
       match2,iord_gas_inhalo,inhalo_z0i,ind13,ind14
       IF (where(ind14 EQ -1))[0] NE -1 THEN BEGIN ;Particles that were accreted and not ejected are not found in the halo at z = 0
           print,'Potential problem: ',strtrim(n_elements(where(ind14 EQ -1)),2),' gas particles that were accreted and not listed as lost are not in it at z = 0'
           match,iord_gas,inhalo_z0i[where(ind14 EQ -1)],ind_z0,ind_parttrace
           IF ind_z0[0] EQ -1 THEN print,'None of the gas particles still exist' ELSE print,strtrim(n_elements(ind_z0),2),' of the gas particles still exist'
           print,'Perhaps these particles formed stars that were then lost from the halo?'
           stop
           ;print,inhalo_z0i[where(ind14 EQ -1)]
           ;print,grp_star[where(igasorder EQ (inhalo_z0i[where(ind14 EQ -1)])[0])]
       ENDIF
       IF (where(ind13 EQ -1))[0] NE -1 THEN BEGIN ;Particles found in halo at z = 0 that were not accreted
           print,'Potential problem: ',strtrim(n_elements(where(ind13 EQ -1)),2),' gas particles / gas particles that spawned stars in the halo at z = 0 are not listed as accreted/reaccreted post expulsion'
           match2,iord_gas[where(inhalo_gas EQ 1)],inhalo_z0i,ind15,ind16
           IF (where(ind15 EQ -1))[0] NE -1 THEN print,'This problem includes ',n_elements(where(ind15 EQ -1)),' z = 0 gas particles'
           match2,sl_deleted.iordergas,inhalo_z0i,ind17,ind18
           IF (where(ind17 EQ -1))[0] NE -1 THEN BEGIN
               print,'This problem includes ',strtrim(n_elements(where(ind17 EQ -1)),2),' z = 0 star-spawning gas particles or ',strtrim(total(sl_deleted[where(ind17 EQ -1)].massform)*units.massunit,2),' solar masses'
               slmissing = sl_deleted[where(ind17 EQ -1)]
               match,slmissing.iorderstar,iord_star,ind19,ind20
               plot,s[where(inhalo_star EQ 1)].x,s[where(inhalo_star EQ 1)].y,psym = 3
               IF n_elements(ind20) EQ 1 THEN oplot,[s[ind20].x,s[ind20].x],[s[ind20].y,s[ind20].y],psym = 5,color = 7 $ 
               ELSE oplot,s[ind20].x,s[ind20].y,psym = 5,color = 7
               print,'Metal mass: ',strtrim(total(s[ind20].metals*s[ind20].mass*units.massunit),2)
               stop
           ENDIF
       ENDIF
;Match current gas particles in halo to those accreted or accreted and either
;not lost or reaccreted after being lost
       iord_gas_inhalo = iord_gas[where(inhalo_gas EQ 1)]
       g_inhalo = g[where(inhalo_gas EQ 1)]
       s_inhalo = s[where(inhalo_star EQ 1)]
       match,iord_gas_inhalo,inhalo_z0i,ind21,ind22
       IF n_elements(ind21) NE n_elements(iord_gas_inhalo) OR n_elements(ind21) NE n_elements(g_inhalo) THEN BEGIN
           print,'ERROR: Some gas particles in halo at z = 0 are not accounted for in particle tracking census'
           stop
       ENDIF
       match,inflowi_uniq_r,iord_gas,ind23,ind24
       match,sl.iorderstar,iord_star,ind25,ind26
       print,'Metal mass contained within all gas once in the halo or stars formed from that gas: ',strtrim(total(g[ind24].mass*units.massunit*g[ind24].zmetal)+total(s[where(ind1 EQ -1)].mass*units.massunit*s[where(ind1 EQ -1)].metals),2)

       inhalo_z0data = [outflow_data_uniq_r[ind7(where(timediff LT 0))],inflow_data_uniq_r[where(ind12 EQ -1)]]
       inhalo_z0timeaccr = [(inflowt_uniq_r[ind8])[where(timediff LT 0)],inflowt_uniq_r[where(ind12 EQ -1)]]
       outhalo_z0timeaccr = (outflowt_uniq_r[ind7])[where(timediff GT 0)]
       outhalo_z0i = (outflowi_uniq_r[ind7])[where(timediff GT 0)]
       outhalo_z0data = (outflow_data_uniq_r[ind7])[where(timediff GT 0)]
       match,outhalo_z0i,iord_gas,ind27,ind28
;       window,4
;Plot current metallicity vs metallicity at time of accretion. Should
;not be a 45 degree line because gas should have become more enriched
;       plot,inhalo_z0data[ind22].metallicity,g_inhalo[ind21].zmetal,psym = 3,xtitle = 'Metallicity at (re)accretion',ytitle = 'z = 0 metallicity',/xlog,/ylog
;       plot,inhalo_z0timeaccr[ind22],g_inhalo[ind21].zmetal - inhalo_z0data[ind22].metallicity,psym = 3,xtitle = 'Time of last (re)accretion',ytitle = 'Current - accreted metallicity',xcharsize = 2,ycharsize = 2
       plot,outhalo_z0timeaccr[ind27],g[ind28].zmetal - outhalo_z0data[ind27].metallicity,psym = 3,xtitle = 'Time of last (re)accretion',ytitle = 'Current - accreted metallicity',xcharsize = 2,ycharsize = 2
       oplot,[0,14],[0,0]
       print,'Metal diffused into (+)/from (-) halo: ',strtrim(total(g[ind28].zmetal*g[ind28].mass*units.massunit- outhalo_z0data[ind27].metallicity*outhalo_z0data[ind27].mass),2),', or ',strtrim(abs(100*total(g[ind28].zmetal*g[ind28].mass*units.massunit - outhalo_z0data[ind27].metallicity*outhalo_z0data[ind27].mass)/total(outhalo_z0data[ind27].metallicity*outhalo_z0data[ind27].mass)),2),'% of total metals in gas'

       print,'z = 0 Metal mass in gas: ',strtrim(total(g_inhalo[ind21].zmetal*g_inhalo[ind21].mass*units.massunit),2),', z = 0 Metal mass in stars: ',strtrim(total(s_inhalo.metals*s_inhalo.mass*units.massunit),2),', z = 0 Metal mass in halo: ',strtrim(total(g_inhalo[ind21].zmetal*g_inhalo[ind21].mass*units.massunit)+total(s_inhalo.metals*s_inhalo.mass*units.massunit),2)
       print,'z = 0 mass in gas: ',strtrim(total(g_inhalo[ind21].mass*units.massunit),2),', z = 0 mass in stars: ',strtrim(total(s_inhalo.mass*units.massunit),2),', z = 0 mass in halo: ',strtrim(total(g_inhalo[ind21].mass*units.massunit)+total(s_inhalo.mass*units.massunit),2)

;       print,'Mass gain (+)/loss (-) before reaccretion: ',strtrim(total(inflowm_uniq[ind8(where(timediff GT 0))] - outflowm_uniq[ind7(where(timediff GT 0))]))
       match2,outflowi_uniq,inflowi_uniq,ind11,ind12
       inhalo_z0i = [inflowi_uniq[ind8(where(timediff LT 0))],inflowi_uniq[where(ind12 EQ -1)]]
       inhalo_z0m = [inflowm_uniq[ind8(where(timediff LT 0))],inflowm_uniq[where(ind12 EQ -1)]]

;*************************** Debug by tallying disk stars
;Confirm that there is no overlap between stars accreted to halo and
;stars formed in it (I think this could happen if a star is lost from halo,
;then reaccreted)
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
       scale_factor = 1  ;Because at redshift zero
       iz0 = n_elements(halodat) - 1
       az = reform(a[iz0,*])
       az0 = az[0]
       az1 = az[1]
       az2 = az[2]
       ax = [az2/sqrt(az0*az0 + az2*az2),0,-1.0*az0/sqrt(az0*az0 + az2*az2)]
       ay = crossp(az,ax)       ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
       basis = [[ax],[ay],[az]]
       gpos = [[g.x*units.lengthunit - center[iz0,2]*scale_factor],[g.y*units.lengthunit - center[iz0,3]*scale_factor],[g.z*units.lengthunit - center[iz0,4]*scale_factor]] ;*scale_factor
       gpos = transpose(transpose(basis)#transpose(gpos))
       g.x = gpos[*,0]
       g.y = gpos[*,1]
       g.z = gpos[*,2]
       indisk = fltarr(h.ngas)
       indisk[where(g.dens*units.rhounit/scale_factor^3 GE 0.1 AND g.tempg LE 1.2e4 AND abs(g.z) LE 3.0)] = 1
       indisk = indisk*inhalo_gas ;ensure that disk gas is within the halo
       iord_gas_indisk = iord_gas[where(indisk EQ 1)]
;Calculate uniq particles that are accreted onto disk, finding the instance when
;they last accrete the disk (contrast wtih inflowdiski_uniq, which find
;first accretion)
       ind_inflowdisk_all_uniq_r = rem_dup(inflowdiski_all,inflowdiskt_all) ;uniq(inflowi,sort(inflowi))
       inflowdiskz_all_uniq_r = inflowdiskz_all[ind_inflowdisk_all_uniq_r]
       inflowdiskm_all_uniq_r = inflowdiskm_all[ind_inflowdisk_all_uniq_r]
       inflowdiski_all_uniq_r = inflowdiski_all[ind_inflowdisk_all_uniq_r]
       inflowdiskt_all_uniq_r = inflowdiskt_all[ind_inflowdisk_all_uniq_r] 
       inflowdisk_all_data_uniq_r = inflowdisk_all_data[ind_inflowdisk_all_uniq_r]
       
;Calculate uniq particles that leave disk, finding the instance when
;they last leave the disk
       ind_relost_uniq_r = rem_dup(relosti,relostt) ;uniq(relosti,sort(relosti))
       relostz_uniq_r = relostz[ind_relost_uniq_r]
       relostm_uniq_r = relostm[ind_relost_uniq_r]
       relosti_uniq_r = relosti[ind_relost_uniq_r]
       relostt_uniq_r = relostt[ind_relost_uniq_r]
       relost_data_uniq_r = relost_data[ind_relost_uniq_r]
       
       match,relosti_uniq_r,inflowdiski_all_uniq_r,ind7,ind8
       IF n_elements(ind7) NE n_elements(relosti_uniq_r) THEN BEGIN
           print,'ERROR: some particles lost from disk were not accreted to the disk'
           stop
       ENDIF
       timediff = relostt_uniq_r[ind7] - inflowdiskt_all_uniq_r[ind8]
       print,'Min/Max time between last time lost from disk and last time reaccreted (+ = out, - = in)',strtrim(minmax(timediff),2)
;       window,1
;       plot,inflowdiskt_all_uniq_r[ind8],relostt_uniq_r[ind7],psym = 3,xtitle = 'Time at last inflow',ytitle = 'Time at last relost'

;Iords of particles that were once accreted to the disk but removed by z = 0
       outdisk_z0i = relosti_uniq_r[ind7(where(timediff GT 0))]
       match,outdisk_z0i,iord_gas_indisk,ind9,ind10
       IF (where(inhalo_gas[ind10] EQ 1))[0] NE -1 THEN BEGIN ;Should be the empty array
           print,'ERROR: gas particles thought to have left disk are still in it'
           stop
       ENDIF

;Find an array of all gas particles in disk at z = 0 along with all
;gas particles that spawned stars and were deleted
       match2,iord_gas,sl.iordergas,ind9a,ind10a
       IF (where(ind10a EQ -1))[0] NE -1 THEN BEGIN
           sl_deleted = sl[where(ind10a EQ -1)]
;Removed duplicates because the same gas particle can form multiple stars 
           sl_deleted = sl_deleted[uniq(sl_deleted.iordergas,sort(sl_deleted.iordergas))]
           iord_gas_indisk = [iord_gas_indisk,sl_deleted.iordergas] 
       ENDIF
;Find array of all gas particles accreted and never removed or
;accreted after being removed
       match2,relosti_uniq_r,inflowdiski_all_uniq_r,ind11,ind12
       indisk_z0i = [relosti_uniq_r[ind7(where(timediff LT 0))],inflowdiski_all_uniq_r[where(ind12 EQ -1)]]
;Plot all gas particles 
;thought because of cycling to be in the disk at z = 0
       match,iord_gas,indisk_z0i,ind_z0,ind_parttrace
       plot,g[ind_z0].x,g[ind_z0].y,psym = 3 
;known to be in the disk at z = 0
       oplot,g[where(indisk EQ 1)].x,g[where(indisk EQ 1)].y,psym = 3,color = 7
;Match iord_gas_indisk (all particles in disk at z= 0 + all that
;accreted and never removed or accreted after being removed
;formed stars) against indisk_z0i (all particles that were 
       match2,iord_gas_indisk,indisk_z0i,ind13,ind14
       IF (where(ind14 EQ -1))[0] NE -1 THEN BEGIN ;Particles that were accreted and not ejected are not found in the disk at z = 0
           print,'Potential problem: ',strtrim(n_elements(where(ind14 EQ -1)),2),' gas particles that were accreted and not listed as lost are not in it at z = 0'
           match,iord_gas,indisk_z0i[where(ind14 EQ -1)],ind_z0,ind_parttrace
           IF ind_z0[0] EQ -1 THEN print,'None of the gas particles still exist' ELSE BEGIN
               print,strtrim(n_elements(ind_z0),2),' of the gas particles still exist'
               IF n_elements(ind_z0) EQ 1 THEN oplot,[g[ind_z0].x,g[ind_z0].x],[g[ind_z0].y,g[ind_z0].y],psym = 3,color = 14$ 
               ELSE oplot,g[ind_z0].x,g[ind_z0].y,psym = 3,color = 14
           ENDELSE
           print,'Perhaps some of these particles formed stars that were then lost from the halo?'
           stop
           ;print,inhalo_z0i[where(ind14 EQ -1)]
           ;print,grp_star[where(igasorder EQ (inhalo_z0i[where(ind14 EQ -1)])[0])]
       ENDIF
       IF (where(ind13 EQ -1))[0] NE -1 THEN BEGIN ;Particles found in disk at z = 0 that were not included in accreted gas
           print,'Potential problem: ',strtrim(n_elements(where(ind13 EQ -1)),2),' gas particles/ gas particles that spawned stars, which are in the halo at z = 0 are not listed as accreted/reaccreted post expulsion'
           match2,iord_gas[where(indisk EQ 1)],indisk_z0i,ind15,ind16
           IF (where(ind15 EQ -1))[0] NE -1 THEN print,'This problem includes ',strtrim(n_elements(where(ind15 EQ -1)),2),' z = 0 gas particles'
           match2,sl_deleted.iordergas,indisk_z0i,ind17,ind18
           IF (where(ind17 EQ -1))[0] NE -1 THEN BEGIN
               print,'This problem includes ',strtrim(n_elements(where(ind17 EQ -1))),' z = 0 star-spawning gas particles or ',strtrim(total(sl_deleted[where(ind17 EQ -1)].massform)*units.massunit,2),' solar masses'
               slmissing = sl_deleted[where(ind17 EQ -1)]
               match,slmissing.iorderstar,iord_star,ind19,ind20
               plot,s[where(inhalo_star EQ 1)].x,s[where(inhalo_star EQ 1)].y,psym = 3
               IF n_elements(ind20) EQ 1 THEN oplot,[s[ind20].x,s[ind20].x],[s[ind20].y,s[ind20].y],psym = 5,color = 7 $ 
               ELSE oplot,s[ind20].x,s[ind20].y,psym = 5,color = 7
               print,'Metal mass: ',strtrim(total(s[ind20].metals*s[ind20].mass*units.massunit),2)
               stop
           ENDIF
       ENDIF
       match,iord_gas,indisk_z0i,ind13a,ind14a 
       mdisk_z0i = [relostm_uniq_r[ind7(where(timediff LT 0))],inflowdiskm_all_uniq_r[where(ind12 EQ -1)]]
       zdisk_z0i = [relostz_uniq_r[ind7(where(timediff LT 0))],inflowdiskz_all_uniq_r[where(ind12 EQ -1)]]
       plot,g[ind13a].mass*units.massunit,mdisk_z0i[ind14a],psym = 3

;Match current gas particles in halo to those accreted or and either
;not lost or reaccreted after being lost
       iord_gas_inhalo = iord_gas[where(inhalo_gas EQ 1)]
       g_inhalo = g[where(inhalo_gas EQ 1)]
       s_inhalo = s[where(inhalo_star EQ 1)]
       match,iord_gas_inhalo,inhalo_z0i,ind21,ind22
       IF n_elements(ind21) NE n_elements(iord_gas_inhalo) OR n_elements(ind21) NE n_elements(g_inhalo) THEN BEGIN
           print,'ERROR: Some gas particles in halo at z = 0 are not accounted for in particle tracking census'
           stop
       ENDIF
       inhalo_z0data = [relost_data_uniq_r[ind7(where(timediff LT 0))],inflow_data_uniq_r[where(ind12 EQ -1)]]
;       window,4
;Plot current metallicity vs metallicity at time of accretion. Should
;not be a 45 degree line because gas should have become more enriched
 ;      plot,inhalo_z0data[ind22].metallicity,g_inhalo[ind21].zmetal,psym = 3,xtitle = 'Metallicity at (re)accretion',ytitle = 'z = 0 metallicity',/xlog,/ylog
       print,'z = 0 Metal mass in gas: ',strtrim(total(g_inhalo[ind21].zmetal*g_inhalo[ind21].mass*units.massunit),2),', z = 0 Metal mass in stars: ',strtrim(total(s_inhalo.metals*s_inhalo.mass*units.massunit),2),', z = 0 Metal mass in halo',strtrim(total(g_inhalo[ind21].zmetal*g_inhalo[ind21].mass*units.massunit)+total(s_inhalo.metals*s_inhalo.mass*units.massunit),2)
       print,'z = 0 mass in gas: ',strtrim(total(g_inhalo[ind21].mass*units.massunit),2),', z = 0 mass in stars: ',strtrim(total(s_inhalo.mass*units.massunit),2),', z = 0 mass in halo',strtrim(total(g_inhalo[ind21].mass*units.massunit)+total(s_inhalo.mass*units.massunit),2)
   ENDIF
;   stop
ENDFOR
;legend,['Gas Accreted','Gas Accreted to Disk','Stars Accerted','Gas Lost','Gas Expelled','Gas Accreted - Gas Lost','Baryonic Mass'],linestyle = [2,5,4,2,3,0,1]
;legend,['Gas Accreted','Gas Accreted to Disk','Gas Lost','Gas Expelled','Gas Accreted - Gas Lost','Disk Gas Accreted - Disk Gas Lost'],linestyle = [2,5,2,3,0,0]
;legend,['Virial Mass','Gas Accreted to Disk','Gas Accreted Once to Disk','Ejected Gas','Expelled Gas'],linestyle = [5,0,1,0,0],color = [fgcolor,254,254,60,20]
multiplot,/reset
IF keyword_set(outplot) THEN device,/close ELSE stop
END
