PRO bulge_sfh,dir,file,haloid = haloid,key = key,step0 = step0,step1 = step1, step2 = step2,binsize = binsize,outplot = outplot,massloc = massloc,label = label
print,dir+file+'.'+haloid
legend = 1

f_bar = 0.16510
xrange = [0,13.75]
massloc_lr = massloc

formatplot,outplot = outplot,thick = formatthick
;loadct,39
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
    thick = 2
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 600
    ysize = 400
    thick = 1
ENDELSE
device,decomposed = 0

IF NOT keyword_set(haloid) THEN haloid = '1'
IF NOT keyword_set(key) THEN key = file
IF NOT keyword_set(step0) THEN step0 = '00512'
IF NOT keyword_set(binsize) THEN binsize = 0.05
binsize2 = 0.2

readcol,dir + '/grp' + haloid + '.mass_metal_track.dat',halo_t,time_t,z_t,mtot_t,mgas_t,mstar_t,mdark_t,metal_t,ox_t,fe_t,HI_t,H2_t,coldg_t,/silent
IF keyword_set(step1) THEN BEGIN
   spawn,'ls ' + dir + '*' + step1 + '/*cosmo*' + step1,tfile1
   rtipsy,tfile1,h1,g,d,s,/justhead
   time1 = z_to_t((1 - h1.time)/h1.time )
ENDIF

IF keyword_set(step2) THEN BEGIN
   spawn,'ls ' + dir + '*' + step2 + '/*cosmo*' + step2,tfile2
   rtipsy,tfile2,h2,g,d,s,/justhead
   time2 = z_to_t((1 - h2.time)/h2.time )
ENDIF

spawn,'ls ' + dir + 'h*param',pfile
units = tipsyunits(pfile[0],/silent)
spawn,'ls ' + dir + '*' + step0 + '/*cosmo*' + step0 + '.halo.' + haloid + '.std',tfile
rtipsy,tfile,h,g,d,s
;spawn,'ls ' + dir + '*' + step0 + '/*' + step0 + '.amiga.grp',file_grp
;readarr,file_grp,h,grp,/ascii
;grp_star = grp[h.ngas + h.ndark:h.n - 1]
;grp_gas = grp[0:h.ngas - 1]
spawn,'ls ' + dir + '*' + step0 + '/*' + step0 + '.halo.' + haloid + '.iord',file_iord
readarr,file_iord,  h,iord,/ascii,type = 'long'
iord_star = iord[h.ngas + h.ndark:h.n - 1]
iord_gas = iord[0:h.ngas - 1]

spawn,'ls ' + dir + '/' + file + '.starlog.fits',files_sl
IF file_test(files_sl[0]) THEN BEGIN
   sl = mrdfits(files_sl[0],1,/silent)
   print,files_sl[0]
ENDIF ELSE BEGIN
   spawn,'ls ' + dir + '/*_2merge.starlog.fits',files_sl
   IF file_test(files_sl[0]) THEN BEGIN
      sl = mrdfits(files_sl[0],1,/silent)
      print,files_sl[0]
   ENDIF ELSE BEGIN
      spawn,'ls ' + dir + '/*.starlog',files_sl
      sl = rstarlog(files_sl[0],/molecularH)
      print,files_sl[0]
   ENDELSE
ENDELSE
match,sl.iorderstar,iord_star,ind0,ind1 ;stars in the halo
sl = sl[ind0]

inflowstarz = mrdfits(dir + '/grp' + haloid + '.accrstars_z.fits',/silent)
inflowstari = mrdfits(dir + '/grp' + haloid + '.accrstars_i.fits',/silent)
;inflowstarm = mrdfits(dir + '/grp' + haloid + '.accrstars_m.fits',/silent)
inflowstari = inflowstari[where(inflowstarz NE 99)]
;inflowstarm = inflowstarm[where(inflowstarz NE 99)]
inflowstarz = inflowstarz[where(inflowstarz NE 99)]
;inflowstart = z_to_t(inflowstarz)

match2,iord_star,inflowstari,ind0,ind1
ind_insitu = where(ind0 EQ -1)
ind_exsitu = where(ind1 EQ -1)

relosti = mrdfits(dir + '/grp' + haloid + '.relost_iord.fits',0,/silent)
relostm = mrdfits(dir + '/grp' + haloid + '.mass_at_relost.fits',0,/silent) 
relostz = mrdfits(dir + '/grp' + haloid + '.relost_z.fits',0,/silent)
relostt = z_to_t(relostz)

reejecti = mrdfits(dir + '/grp' + haloid + '.reeject_rvir0.2_iord.fits',0,/silent)
reejectm = mrdfits(dir + '/grp' + haloid + '.mass_at_reeject_rvir0.2.fits',0,/silent)
reejectz = mrdfits(dir + '/grp' + haloid + '.reeject_rvir0.2_z.fits',0,/silent)
reejectt = z_to_t(reejectz)

sfh = weighted_histogram(sl.timeform*units.timeunit/1e9,weight = sl.massform*units.massunit/(binsize*1e9),min = 1e-6,max = 13.75,binsize = binsize,locations = timearr)
sfh_insitu = weighted_histogram(sl[ind_insitu].timeform*units.timeunit/1e9,weight = sl[ind_insitu].massform*units.massunit/(binsize*1e9),min = 1e-6,max = 13.75,binsize = binsize,locations = timearr)
sfh2_insitu = weighted_histogram(sl[ind_insitu].timeform*units.timeunit/1e9,weight = sl[ind_insitu].massform*units.massunit/(binsize2*1e9),min = 1e-6,max = 13.75,binsize = binsize2,locations = timearr2)
relosthm  = weighted_histogram(relostt, weight =  relostm/(binsize2*1e9),  min = 1e-6, max = 13.75, locations = timearr2,binsize = binsize2)
reejecthm  = weighted_histogram(reejectt, weight =  reejectm/(binsize2*1e9),  min = 1e-6, max = 13.75, locations = timearr2,binsize = binsize2)

!y.style = 0
!x.style = 1
;---------------------------- SFH ----------------------
IF keyword_set(outplot) THEN  device,filename = outplot + key + haloid + '_sfh.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
;yrange = [0,max([sfh,relosthm])]
yrange = [0,max([sfh,reejecthm])]
plot,timearr,sfh,xtitle = 'Time [Gyr]',ytitle = 'Rate [M' + sunsymbol() + textoidl(' yr^{-1}]'),psym = 10,thick = 1,yrange = yrange,title = file + '.' + haloid,xrange = xrange
oplot,timearr,sfh_insitu,color = 254,psym = 10,thick = thick
;oplot,timearr2,relosthm,color = 60,psym = 10,thick = 1
oplot,timearr2,reejecthm,color = 60,psym = 10,thick = thick
IF keyword_set(time1) THEN oplot,[time1,time1],[0,1e5],linestyle = 2
IF keyword_set(time2) THEN oplot,[time2,time2],[0,1e5],linestyle = 2
IF keyword_set(legend) THEN $
   legend,['Star Formation','In Situ Star Formation','Outflow'],linestyle = [0,0,0],color = [fgcolor,254,60],/top,/right,box = 0
IF keyword_set(outplot) THEN device, /close ;ELSE stop

;--------------------------- M_500 -------------------------
readcol,dir + '/grp' + haloid + '.mass500.txt',haloids,time,z,gmass,smass,dmass,slope,gmass1,smass1,dmass1,slope1
IF keyword_set(outplot) THEN  device,filename = outplot + key + haloid + '_m500.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
yrange = [0,max([dmass1,gmass1,smass1])]
plot,time,[gmass1 + smass1],xtitle = 'Time [Gyr]',ytitle = textoidl('M_{baryon} [M') + sunsymbol() + ']',yrange = yrange,psym = -4,thick = thick,/nodata,xrange = xrange
oplot,time,[gmass1],color = 60,psym = -4,thick = thick
;oplot,time,[gmass + smass],color = 60,psym = -4,thick = thick
oplot,time,smass1,color = 254,psym = -4,thick = thick
oplot,time,[dmass1],psym = -4,thick = thick
;oplot,time,[gmass+smass+dmass],psym = -4,thick = thick
IF keyword_set(time1) THEN oplot,[time1,time1],[0,1e11],linestyle = 2
IF keyword_set(time2) THEN oplot,[time2,time2],[0,1e11],linestyle = 2
;IF keyword_set(legend) THEN $
;   legend,['Total Mass','Baryons','Stars'],linestyle = [0,0,0],color = [fgcolor,60,254],/bottom,/right,box = 0
IF keyword_set(legend) THEN $
   legend,['Dark Matter','Cold Gas','Stars'],linestyle = [0,0,0],color = [fgcolor,60,254],right = massloc,box = 0;,bottom = massloc
IF keyword_set(outplot) THEN device, /close ;ELSE stop

;------------------------ SFH + m_500 **********************
IF keyword_set(outplot) THEN  device,filename = outplot + key + haloid + '_sfh_m500.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = 1.5*ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = 1.5*ysize
multiplot,[1,2];,mxtitle = 'Time [Gyr]'
yrange = [0,max([sfh,reejecthm])]
plot,timearr,sfh,ytitle = 'Rate [M' + sunsymbol() + textoidl(' yr^{-1}]'),psym = 10,thick = 1,yrange = yrange,xrange = xrange;,title = file + '.' + haloid
oplot,timearr,sfh_insitu,color = 254,psym = 10,thick = thick
;oplot,timearr2,relosthm,color = 60,psym = 10,thick = 1
oplot,timearr2,reejecthm,color = 60,psym = 10,thick = thick
IF keyword_set(time1) THEN oplot,[time1,time1],[0,1e5],linestyle = 2
IF keyword_set(time2) THEN oplot,[time2,time2],[0,1e5],linestyle = 2
IF keyword_set(legend) THEN $
   legend,['Star Formation','In Situ Star Formation','Outflow'],linestyle = [0,0,0],color = [fgcolor,254,60],/top,/right,box = 0
IF keyword_set(label) THEN $
   legend,label,/top,/left,box = 0
multiplot

yrange = [0,max([dmass1,gmass1,smass1])]
plot,time,[gmass1 + smass1],ytitle = textoidl('M_{center} [M') + sunsymbol() + ']',yrange = yrange,psym = -4,thick = thick,/nodata,xrange = xrange,xtitle = 'Time [Gyr]'
oplot,time,[gmass1],color = 60,psym = -4,thick = thick
;oplot,time,[gmass + smass],color = 60,psym = -4,thick = thick
oplot,time,smass1,color = 254,psym = -4,thick = thick
;axis,/yaxis,ytitle = textoidl('M_{total} [M') + sunsymbol() + ']',yrange = yrange
oplot,time,[dmass1],psym = -4,thick = thick
;oplot,time,[gmass+smass+dmass],psym = -4,thick = thick
IF keyword_set(time1) THEN oplot,[time1,time1],[0,1e11],linestyle = 2
IF keyword_set(time2) THEN oplot,[time2,time2],[0,1e11],linestyle = 2
;IF keyword_set(legend) THEN $
;   legend,['Total Mass','Baryons','Stars'],linestyle = [0,0,0],color = [fgcolor,60,254],/bottom,/right,box = 0
IF keyword_set(legend) THEN $
   legend,['Dark Matter','Cold Gas','Stars'],linestyle = [0,0,0],color = [fgcolor,60,254],right = massloc,box = 0;bottom = massloc

IF keyword_set(outplot) THEN  device,/close ELSE stop
multiplot,/reset

IF 0 THEN BEGIN
;------------------------ SFH2 + m_500
!X.MARGIN = [12,3]
IF keyword_set(outplot) THEN  device,filename = outplot + key + haloid + '_sfh_m5002.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = 2*ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = 2*ysize
multiplot,[1,3];,mxtitle = 'Time [Gyr]'
yrange = [0,max([sfh,reejecthm])]
plot,timearr,sfh,ytitle = 'Rate [M' + sunsymbol() + textoidl(' yr^{-1}]'),psym = 10,thick = 1,yrange = yrange,xrange = xrange;,title = file + '.' + haloid
oplot,timearr,sfh_insitu,color = 254,psym = 10,thick = thick
;oplot,timearr2,relosthm,color = 60,psym = 10,thick = 1
oplot,timearr2,reejecthm,color = 60,psym = 10,thick = thick
IF keyword_set(time1) THEN oplot,[time1,time1],[0,1e5],linestyle = 2
IF keyword_set(time2) THEN oplot,[time2,time2],[0,1e5],linestyle = 2
multiplot

yrange = [0,max(gmass1+smass1)]
plot,time,[gmass1 + smass1],ytitle = textoidl('M_{center} [M') + sunsymbol() + ']',yrange = yrange,psym = -4,thick = thick,/nodata,xrange = xrange
oplot,time,[gmass1 + smass1],color = 60,psym = -4,thick = thick
oplot,time,smass1,color = 254,psym = -4,thick = thick
IF keyword_set(time1) THEN oplot,[time1,time1],[0,1e10],linestyle = 2
IF keyword_set(time2) THEN oplot,[time2,time2],[0,1e10],linestyle = 2
multiplot

yrange = [0,max(dmass1+gmass1+smass1)]
plot,time,[dmass1 + gmass1 + smass1],ytitle = textoidl('M_{total} [M') + sunsymbol() + ']',yrange = yrange,psym = -4,thick = thick,/nodata,xrange = xrange,xtitle = 'Time [Gyr]'
oplot,time,[dmass1 + gmass1 + smass1],psym = -4,thick = thick
IF keyword_set(time1) THEN oplot,[time1,time1],[0,1e12],linestyle = 2
IF keyword_set(time2) THEN oplot,[time2,time2],[0,1e12],linestyle = 2
IF keyword_set(outplot) THEN  device,/reset ;ELSE stop
multiplot,/reset
ENDIF

;----------------------------------------------
multiplot,[1,1]
IF keyword_set(outplot) THEN  device,filename = outplot +  key + haloid + '_ml_time.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
ml_lost = relosthm/sfh2_insitu
ml_eject = reejecthm/sfh2_insitu
;yrange = [min([ml_lost[where(ml_lost NE 0)],ml_eject[where(ml_eject NE 0)]]),max([ml_lost,ml_eject])]
yrange = [min(ml_eject[where(ml_eject NE 0)]),max(ml_eject)]
plot,timearr2,ml_eject,xtitle = 'Time [Gyr]',ytitle = textoidl('\eta'),psym = 4,/ylog,yrange = yrange,/nodata,title = file + '.' + haloid
oplot,timearr2,ml_eject,psym = 4,color = 254
;oplot,timearr2,ml_lost,psym = 4,color = 60
IF keyword_set(outplot) THEN device, /close ;ELSE stop

;------------------------------------------------
marr2 = spline(time_t,mtot_t,timearr2)
IF keyword_set(outplot) THEN  device,filename = outplot +  key + haloid + '_ml_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
ml_lost = relosthm/sfh2_insitu
ml_eject = reejecthm/sfh2_insitu
yrange = [min([ml_lost[where(ml_lost NE 0)],ml_eject[where(ml_lost NE 0)]]),max([ml_lost,ml_eject])]
mxrange = [min(marr2[where(marr2 GT 0 AND ml_eject GT 0)]),max(marr2)]
plot,marr2,ml_eject,xtitle = textoidl('M_{vir} [M') + sunsymbol() + ']',ytitle = textoidl('\eta'),psym = 4,/ylog,yrange = yrange,/xlog,/nodata,title = file + '.' + haloid,xrange = mxrange
oplot,marr2,ml_eject,psym = 4,color = 254
;oplot,marr2,ml_lost,psym = 4,color = 60
IF keyword_set(outplot) THEN device, /close ;ELSE stop

openw,lun,dir + '/grp' + haloid + '.massloading.txt',/get_lun
printf,lun,'Time [Gyr]','M_{vir} [M_{sun}]','ML_{eject}','ML_{lost}',format='(A18,A18,A18,A18)'
printf,lun,transpose([[timearr2],[marr2],[ml_eject],[ml_lost]])
close,lun

END
