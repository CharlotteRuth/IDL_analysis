PRO plot_inflow_outflow, dirs, haloid = haloid, outplot = outplot, colors = colors, thicks = thicks, linestyle = linestyle, scale = scale,formatthick = formatthick, keys = keys
; psym = psym, symsize = symsize

filebase = 'steps/masstrack.dat'
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
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2;1
    IF NOT keyword_set(linestyle) THEN linestyle = fltarr(n) 
;    IF NOT keyword_set(psym) THEN psym =-1*( fltarr(n) + 4)
;    IF NOT keyword_set(symsize) THEN symsize = fltarr(n) + 2
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2;1 ;findgen(n)*2
    IF NOT keyword_set(linestyle) THEN linestyle = findgen(n)*2   
;    IF NOT keyword_set(psym) THEN psym = -1*(fltarr(n) + 4)
;    IF NOT keyword_set(symsize) THEN symsize = fltarr(n) + 2
ENDELSE
IF NOT keyword_set(scale) THEN scale = fltarr(n) + 1.0

ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
z = reverse(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
ticktime_in_z = findgen(7)*2*1e9 ;redshift major plot
tickred_in_z = spline(t,z,ticktime_in_z )

tickred_in_t = reverse(findgen(5))
ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e9

IF (keyword_set(outplot)) THEN device,filename = outplot + '_fb.eps',/color,bits_per_pixel= 8,/times,ysize = 16,xsize = 18,xoffset =  2,yoffset =  2 ELSE window,9,xsize = 712,ysize = 600
binsize = 5e8
nbins = round(14.0/(binsize/1e9))

FOR i = 0, n - 1 DO BEGIN
    cd,dirs[i]
    readcol,dirs[i] + filebase,halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR
    spawn,'ls ' + dirs[i] + '*512/*cosmo**512',file
    rtipsy,file,h,g,d,s,/justhead
    spawn,'ls ' + dirs[i] + '*512/*512.coolontime',file_coolon
    readarr,file_coolon,h,coolon,part = 'gas',/ascii
    spawn,'ls ' + dirs[i] + '*512/*512.iord',file_iord
    readarr,file_iord,  h,iord,  part = 'gas',/ascii
    spawn,'ls ' + dirs[i] + 'h*param',pfile
    units = tipsyunits(pfile[0])


    iords = mrdfits(dirs[i] + '/grp1.allgas.iord.fits')
; ------------------- Inflow -----------------------------
    print,'Read accretion files'
    inflowz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accr_rvir_z.fits');'.accrz.fits'
    inflowm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_accr_rvir.fits');'.mass_at_accr.fits'
    inflowi = iords
;    inflowt = (ageUniverse - wmap3_lookback(inflowz))/1e9  

    uniqz = inflowz[uniq(inflowz,sort(inflowz))]
    uniqt = (ageUniverse - wmap3_lookback(uniqz))/1e9

    inflowt = z_to_t(uniqz,uniqt,inflowz)
    inflowmcum  = weighted_histogram(inflowt, weight =  inflowm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)   

; ------------------- Inflow Disk -----------------------------
    inflowdiskz_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accrdisk_z.fits')
    inflowdiskm_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_accrdisk.fits')
    inflowdiskz     = inflowdiskz_all[ where(inflowdiskz_all NE 99.0)]
    inflowdiskm     = inflowdiskm_all[ where(inflowdiskz_all NE 99.0)]
    inflowdiski     = iords          [ where(inflowdiskz_all NE 99.0)] 
;    inflowdiskt = (ageUniverse - wmap3_lookback(inflowdiskz))/1e9  
    inflowdiskt = z_to_t(uniqz,uniqt,inflowdiskz)
    inflowmdiskcum  = weighted_histogram(inflowdiskt, weight =  inflowdiskm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)

; ------------------- Inflow Stars -----------------------------
    IF file_test(dirs[i] + '/grp' + haloid[i] + '.accrstars_z.fits') THEN BEGIN
       inflowstarz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accrstars_z.fits')
       inflowstarm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.accrstars_m.fits')
;    inflowstart = (ageUniverse - wmap3_lookback(inflowstarz))/1e9  
       inflowstart = z_to_t(uniqz,uniqt,inflowstarz)
       inflowmstarcum  = weighted_histogram(inflowstart, weight =  inflowstarm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)   
       track_stars = 1
    ENDIF ELSE track_stars = 0

;------------------ Expelled Material -------------------------
;Expelled is expelled because of supernova feedback
    print,'Read files of gas expullsion'
    expellm_all  = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_expell.fits',0)
    expellz_all  = mrdfits(dirs[i] + '/grp' + haloid[i] + '.expell_z.fits',0)
    expellm      = expellm_all[    where(expellz_all  NE 99.0)]
    expellz      = expellz_all[    where(expellz_all  NE 99.0)] 
    expelli      = iords      [    where(expellz_all  NE 99.0)]
;    expellt      = (ageUniverse - wmap3_lookback(expellz))/1e9
    expellt      = z_to_t(uniqz,uniqt,expellz) 

;Outflow gas is gas that is lost because of stripping or wierd Amiga
;halo definition
    print,'Read files of outflows'
;    outflowiord_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.outflow_iord.fits')  
    outflowz_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.outflow_z.fits')
    outflowm_all = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_outflow.fits')
;    outflowiord  = outflowiord_all[where(outflowz_all NE 99.0)]
    outflowz     = outflowz_all[   where(outflowz_all NE 99.0)]
    outflowm     = outflowm_all[   where(outflowz_all NE 99.0)]
    outflowi     = iords       [   where(outflowz_all NE 99.0)]
;    outflowt     = (ageUniverse - wmap3_lookback(outflowz))/1e9
    outflowt = z_to_t(uniqz,uniqt,outflowz) 

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
    
    expellmcum  = weighted_histogram(expellt, weight =  expellm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)
    outflowmcum = weighted_histogram(outflowt,weight =  outflowm, min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)

;------------------ (Re)Accreted, Expelled, Ejected, and Outflowing Material -----------------
    reejecti = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reeject_iord.fits',0)
    reejectm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reeject.fits',0)
    reejectz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reeject_z.fits',0)
    reejectt = z_to_t(uniqz,uniqt,reejectz)
    reejecti_uniq = reejecti[uniq(reejecti,sort(reejecti))]
    reejectmcum  = weighted_histogram(reejectt, weight =  reejectm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)
    
    reexpelli = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reexpell_iord.fits',0)
    reexpellm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reexpell.fits',0) 
    reexpellz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reexpell_z.fits',0)
    reexpellt = z_to_t(uniqz,uniqt,reexpellz)
    reexpelli_uniq = reexpelli[uniq(reexpelli,sort(reexpelli))]
    reexpellmcum  = weighted_histogram(reexpellt, weight =  reexpellm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)

    IF file_test(dirs[i] + '/grp' + haloid[i] + '.reaccr_iord.fits') THEN BEGIN
       reaccri = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccr_iord.fits',0)
       reaccrm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reaccr.fits',0)
       reaccrz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reaccr_z.fits',0)

       uniqz = reaccrz[uniq(reaccrz,sort(reaccrz))]
       uniqt = (ageUniverse - wmap3_lookback(uniqz))/1e9

       reaccrt = z_to_t(uniqz,uniqt,reaccrz)
       reaccri_uniq = reaccri[uniq(reaccri,sort(reaccri))]
       reaccrmcum  = weighted_histogram(reaccrt, weight =  reaccrm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)  
   
       reoutflowi = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reoutflow_iord.fits',0)
       reoutflowm = mrdfits(dirs[i] + '/grp' + haloid[i] + '.mass_at_reoutflow.fits',0)
       reoutflowz = mrdfits(dirs[i] + '/grp' + haloid[i] + '.reoutflow_z.fits',0)
       reoutflowt = z_to_t(uniqz,uniqt,reoutflowz)
       reoutflowi_uniq = reoutflowi[uniq(reoutflowi,sort(reoutflowi))]
       reoutflowmcum  = weighted_histogram(reoutflowt, weight =  reoutflowm,  min = 1e-6, max = round(ageUniverse/1e9), locations = timearr,/cum)
       
       track_reout = 1
    ENDIF ELSE track_reout = 0

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

    yrange = [-2.5e10,5e10]
    yrange = [-2e9,3e9]
    yrange = [-3e10,3.5e10] ;h986
    yrange = [-3.5e10,5.5e10]
;    yrange = [-2e7,4e7]
    IF i EQ 0 THEN BEGIN
        plot,timearr,expellmcum*scale[i],ytitle = 'Gas  Mass [M' + sunsymbol() + ']',xstyle =9,xrange=[0,ageUniverse/1e9],xtitle = 'Time [Gyr]',yrange = yrange,/nodata,ystyle = 9,xmargin = [17,17],ymargin = [6,6]
        axis,ystyle = 1,yaxis = 1,ytitle = textoidl('8 X M_{star}') + '[M' + sunsymbol() + ']'
        axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtitle = 'z',xtickv = ticktime_in_t,xticks = n_elements(ticktime_in_t) - 1        
        oplot,[0,14],[0,0],thick = 2
    ENDIF
    oplot,timearr,inflowmcum*scale[i],                linestyle = 2,thick = thicks[i],color = colors[i] ;gas mass in
    IF track_reout THEN oplot,timearr,reaccrmcum*scale[i],                linestyle = 2,thick = thicks[i],color = colors[i] + 50 ;gas mass in (multiple times)
    oplot,timearr,inflowmdiskcum*scale[i],            linestyle = 5,thick = thicks[i],color = colors[i] ;gas mass into disk
;    oplot,timearr,inflowmstarcum*scale[i],            linestyle = 4,thick = thicks[i],color = colors[i] ;star mass in
    oplot,timearr,-1.0*expellmcum*scale[i],           linestyle = 3,thick = thicks[i],color = colors[i] ;gas mass expelled
    oplot,timearr,-1.0*reejectmcum*scale[i],           linestyle = 3,thick = thicks[i],color = colors[i]; + 50 ;gas mass ejected form disk (multiple times)
    oplot,timearr,-1.0*outflowmcum*scale[i],          linestyle = 2,thick = thicks[i],color = colors[i] ;gas mass out
    IF track_reout THEN oplot,timearr,-1.0*reoutflowmcum*scale[i],          linestyle = 2,thick = thicks[i],color = colors[i] + 50 ;gas mass out
    IF track_reout AND track_stars THEN BEGIN
       oplot,timearr,(reaccrmcum + inflowmstarcum - reoutflowmcum)*scale[i],linestyle = 0,thick = thicks[i],color = colors[i] + 50 ;gas + stars in - gas out
       oplot,time,(mgas + mstar)*scale[i],               linestyle = 1,thick = thicks[i],color = colors[i] ;baryonic mass
    ENDIF
;Baryonic mass will be different from inflow - outflow because of gas that is reaccreted later
;    IF track_stars THEN oplot,timearr,(inflowmdiskcum + inflowmstarcum - reejectmcum)*scale[i],linestyle = 0,thick = thicks[i],color = colors[i]+50 ;gas + stars in to disk - gas out of disk (need to add something that will show multiple accretions onto disk)
;    oplot,time,(mcoldg + mstar)*scale[i],                linestyle = 1,thick = thicks[i],color = colors[i]+100 ;baryonic disk mass. Possibly different because of stripping, lack of rediskaccretion, and different definitions of cold gas?
;    oplot,time,mgas*scale[i],linestyle = 2,thicks = 2,color = colors[i]
;    oplot,time,mstar*scale[i],linestyle = 2,thicks = 2,color = colors[i]
;    oplot,time,mstar*scale[i] + mgas*scale[i],linestyle = 0,thicks = 2,color = colors[i]
;    oplot,time,cuminflow*scale[i] - (mstar*scale[i] + mgas*scale[i]),linestyle = 1,thicks = 2,color = colors[i]
;    IF i EQ 0 THEN outflowratio0 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]] ELSE outflowratio1 = [[outflowtbin],[cumoutflow/cuminflow],[outflowmbin]]

stop
ENDFOR
IF keyword_set(keys) THEN legend,keys,color = colors,/bottom
;legend,['Gas Accreted','Gas Accreted to Disk','Stars Accreted','Gas Lost','Gas Expelled','Gas Accreted - Gas Lost','Baryonic Mass'],linestyle = [2,5,4,2,3,0,1]
legend,['Gas Accreted','Gas Accreted to Disk','Gas Lost','Gas Expelled','Gas Accreted - Gas Lost','Disk Gas Accreted - Disk Gas Lost'],linestyle = [2,5,2,3,0,0]
IF keyword_set(outplot) THEN device,/close ELSE stop
END
