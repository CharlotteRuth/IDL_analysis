;Altered 1/23 to determine metallicity of reaccretion by using
;inflowzdiskcum_all_uniq rather than inflowzdiskcum_all_uniq_disk. I have only
;updated the data in *.inflow_outflow_z.dat for the four galaxies in
;the paper. Will have to rerun inflow_outflow_historyz.pro to get data
;for more galaxies if needed. Check the header of
;*.inflow_outflow_z.dat to know for sure if the file has been
;updated. Updated files will list inflowzdiskcum_all_uniq rather than inflowzdiskcum_all_uniq_disk

PRO inflow_outflow_metallicity_plot,dirs,haloid = haloid,outplot = outplot,colors = colors,label = label,pmulti = pmulti

!EXCEPT = 0 ;Supress math errors
zsolar = 0.0130215

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

IF keyword_set(pmulti) THEN BEGIN
    mxtitle = 'Time [Gyr]'
    mytitle = 'Metallicity [Z/Z' + sunsymbol() + ']'    
    IF (keyword_set(outplot)) THEN device,filename = outplot + '_inflow_outflow_z.eps',/color,bits_per_pixel= 8,/times,ysize = 32,xsize = 32,xoffset =  2,yoffset =  2 ELSE window,2,xsize = 800,ysize = 800
    multiplot,pmulti,gap = 0.04,mxTitle = mxtitle,myTitle = mytitle,doyaxis = 1,doxaxis = 1,mxTitSize = !x.charsize,myTitSize = !y.charsize,mxTitOffset = !x.charsize*1.5,myTitOffset = !y.charsize*1.2;,YTICKFORMAT='(F6.3)'
ENDIF ELSE BEGIN
    xtitle = 'Time [Gyr]'
    ytitle = 'Metallicity [Z/Z' + sunsymbol() + ']'
ENDELSE

ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
z = reverse(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
ticktime_in_z = findgen(7)*2*1e9 ;redshift major plot
tickred_in_z = spline(t,z,ticktime_in_z )
tickred_in_t = reverse(findgen(5))
ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e

binsize = 0.5
timearr_bin = findgen(14/binsize + 1)*binsize
reejectmhist_bin = fltarr(14/binsize)
reejectzhist_bin = fltarr(14/binsize)
inflowmdiskhist_all_bin = fltarr(14/binsize)
inflowzdiskhist_all_bin = fltarr(14/binsize)
inflowmdiskhist_all_uniq_bin = fltarr(14/binsize)
inflowzdiskhist_all_uniq_bin = fltarr(14/binsize)
relostmhist_bin = fltarr(14/binsize)
relostzhist_bin = fltarr(14/binsize)
inflowmdiskhist_bin = fltarr(14/binsize)
inflowzdiskhist_bin = fltarr(14/binsize)
reejectmhist_bin = fltarr(14/binsize)
reejectzhist_bin = fltarr(14/binsize)

ageUniverse = 13.7346*1e9 ;wmap3_lookback(100)
z = reverse(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)
ticktime_in_z = findgen(7)*2*1e9 ;redshift major plot
tickred_in_z = spline(t,z,ticktime_in_z )
tickred_in_t = reverse(findgen(5))
ticktime_in_t = (ageUniverse - wmap3_lookback(tickred_in_t))/1e9

FOR i = 0, n - 1 DO BEGIN ;Iterate through the files
    print,i,dirs[i],haloid[i]
   spawn,'ls ' + dirs[i] + 'h*param',pfile
   units = tipsyunits(pfile[0],/silent)
   outbase = strsplit(pfile,'.param',/extract,/regex)
;    openr,1,outbase[0] + '.grp' + haloid[i] + '.inflow_outflow_z.dat'
   readcol,dirs[i] + '/grp' + haloid[i] + '.metals.txt',zmetals,ox,fe,coldg,zmetals_H,ox_H,fe_H,mHI_m,mH2_m,Hgas,/silent 
   IF file_test(dirs[i] + '/grp' + haloid[i] + '.mass_metal_track.dat') THEN BEGIN
       readcol,dirs[i] + '/grp' + haloid[i] + '.mass_metal_track.dat',halo,time,z,mtot,mgas,mstar,mdark,metal,ox,fe,mHI,mH2,mcoldg,/silent
       print,max(mtot)
       IF max(time) LT 1 THEN time = (time*units.timeunit)/1e9
   ENDIF
   IF outbase EQ '/nobackupp8/crchrist/MolecH/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK' THEN BEGIN ; AND haloid[i] eq '2' THEN BEGIN
       coldg[49] = 0 ;Problems with these steps so ignore
       coldg[100] = 0
   ENDIF
   IF NOT file_test(outbase[0] + '.grp' + haloid[i] + '.inflow_outflow_z.dat') THEN CONTINUE
   readcol,outbase[0] + '.grp' + haloid[i] + '.inflow_outflow_z.dat',timearr,reejectmhist,reejectzhist,inflowmdiskhist_all,inflowzdiskhist_all,inflowmdiskhist_all_uniq,inflowzdiskhist_all_uniq,relostmhist,relostzhist,inflowmdiskhist,inflowzdiskhist,reejectmhist,reejectzhist ;,format = '(e,e,e,e,e,e,e,e,e,e,e,e,e)'
    FOR ibin = 0, n_elements(timearr_bin) - 2 DO BEGIN
        ind = where(timearr GE timearr_bin[ibin] AND timearr LT timearr_bin[ibin + 1])
        reejectmhist_bin[ibin] = total(reejectmhist[ind]) ;mass of ejected material
        reejectzhist_bin[ibin] = total(reejectzhist[ind]) ;metal mass of ejected material
        inflowmdiskhist_all_bin[ibin] = total(inflowmdiskhist_all[ind]) ;mass of all accretion onto disk
        inflowzdiskhist_all_bin[ibin] = total(inflowzdiskhist_all[ind]) ;metal mass of all accretion onto disk
        inflowmdiskhist_all_uniq_bin[ibin] = total(inflowmdiskhist_all_uniq[ind]) ;external mass accretion
        inflowzdiskhist_all_uniq_bin[ibin] = total(inflowzdiskhist_all_uniq[ind]) ;metal mass accreted externally
        relostmhist_bin[ibin] = total(relostmhist[ind])
        relostzhist_bin[ibin] = total(relostzhist[ind])
        inflowmdiskhist_bin[ibin] = total(inflowmdiskhist[ind]) ;reaccretion after ejection
        inflowzdiskhist_bin[ibin] = total(inflowzdiskhist[ind])
        reejectmhist_bin[ibin] = total(reejectmhist[ind])
        reejectzhist_bin[ibin] = total(reejectzhist[ind])
    ENDFOR
    
    IF NOT keyword_set(pmulti) THEN $
      IF (keyword_set(outplot)) THEN device,filename = outplot + '_' + strtrim(i,2) + '_inflow_outflow_z.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 ELSE window,2,xsize = 600,ysize = 400
    minplot = min(inflowzdiskhist_all_bin[where(inflowmdiskhist_all_bin NE 0)]/inflowmdiskhist_all_bin[where(inflowmdiskhist_all_bin NE 0)]/zsolar)
    maxplot = max(reejectzhist_bin[where(reejectmhist_bin NE 0)]/reejectmhist_bin[where(reejectmhist_bin NE 0)]/zsolar)*1.1
    IF maxplot LT 0.15 THEN ytickname = ['0.00','0.02','0.04','0.06','0.08','0.10','0.12','0.14'] $
    ELSE IF maxplot LT 0.27 THEN ytickname = ['0.00','0.05','0.10','0.15','0.20','0.25'] ELSE $
      IF maxplot LT 1 THEN ytickname = ['0.0','0.2','0.4','0.6','0.8','1.0']
    plot,timearr_bin[where(reejectmhist_bin NE 0)]+ binsize/2,reejectzhist_bin[where(reejectmhist_bin NE 0)]/reejectmhist_bin[where(reejectmhist_bin NE 0)]/zsolar,xrange=[0.1,ageUniverse/1e9],ytitle = ytitle,xtitle = xtitle,/nodata,xstyle = 9,yrange = [0,maxplot],xtickname = ['2','4','6','8','10','12'],ytickname = ytickname,ystyle = 1;,ytext_color = 'red';,ystyle = 1
    print,minplot,maxplot
    axis,yaxis = 0;,;,color = 'red'
    axis,xstyle = 1,xaxis = 1,XTICKFORMAT = 'z_from_time_axes',xtickv = ticktime_in_t,xticks = n_elements(ticktime_in_t) - 1
;All accretion
    oplot,timearr_bin[where(inflowmdiskhist_all_bin NE 0)]+ binsize/2,inflowzdiskhist_all_bin[where(inflowmdiskhist_all_bin NE 0)]/inflowmdiskhist_all_bin[where(inflowmdiskhist_all_bin NE 0)]/zsolar,linestyle = 5,color = cpristine,thick = thicks[i]
;Pristine accretion
;    oplot,timearr_bin[where(inflowmdiskhist_all_uniq_bin NE 0)] + binsize/2,inflowzdiskhist_all_uniq_bin[where(inflowmdiskhist_all_uniq_bin NE 0)]/inflowmdiskhist_all_uniq_bin[where(inflowmdiskhist_all_uniq_bin NE 0)],linestyle = 5,color = cpristine,thick = thicks[i]
;Reaccreted material    
    inflowmdiskhist_all_return_bin = inflowmdiskhist_all_bin - inflowmdiskhist_all_uniq_bin
    inflowzdiskhist_all_return_bin = inflowzdiskhist_all_bin - inflowzdiskhist_all_uniq_bin
    oplot,timearr_bin[where(inflowmdiskhist_all_return_bin NE 0)]+ binsize/2,inflowzdiskhist_all_return_bin[where(inflowmdiskhist_all_return_bin NE 0)]/inflowmdiskhist_all_return_bin[where(inflowmdiskhist_all_return_bin NE 0)]/zsolar,linestyle = 2,color = cheat,thick = thicks[i]
;Material lost from the disk
    oplot,timearr_bin[where(relostmhist_bin NE 0)]+ binsize/2,relostzhist_bin[where(relostmhist_bin NE 0)]/relostmhist_bin[where(relostmhist_bin NE 0)]/zsolar,color = cheat,thick = thicks[i]
;Reaccretion of ejected material
;    oplot,timearr_bin[where(inflowmdiskhist_bin NE 0)]+ binsize/2,inflowzdiskhist_bin[where(inflowmdiskhist_bin NE 0)]/inflowmdiskhist_bin[where(inflowmdiskhist_bin NE 0)]/zsolar,linestyle = 2,color = ceject,thick = thicks[i]
;Metallicity of the ejected gas
    oplot,timearr_bin[where(reejectmhist_bin NE 0)]+ binsize/2,reejectzhist_bin[where(reejectmhist_bin NE 0)]/reejectmhist_bin[where(reejectmhist_bin NE 0)]/zsolar,color = ceject,thick = thicks[i]
    IF file_test(dirs[i] + '/grp' + haloid[i] + '.mass_metal_track.dat') THEN oplot,time[where(coldg NE 0)],zmetals[where(coldg NE 0)]/coldg[where(coldg NE 0)]/zsolar,thick = thicks[i]
    IF i EQ 1 THEN legend,['Disk Gas','All Accreted','Reaccreted','Removed from disk','Ejected from disk'],linestyle = [0,5,5,0,0],color = [fgcolor,cpristine,cheat,cheat,ceject],/right,/bottom,box=0,thick = [thicks[i],thicks[i],thicks[i],thicks[i],thicks[i]]
    IF keyword_set(label) THEN legend,[label[i]],box = 0,/left,/top
    IF keyword_set(pmulti) THEN multiplot ELSE $
      IF keyword_set(outplot) THEN device,/close ELSE stop
ENDFOR
IF keyword_set(pmulti) THEN BEGIN
    multiplot,/reset
    !P.multi = 0
    IF keyword_set(outplot) THEN device,/close
ENDIF

END
