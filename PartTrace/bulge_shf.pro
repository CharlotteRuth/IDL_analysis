PRO bulge_sfh,dir,file,haloid = haloid,step = step,step1 = step1, step2 = step2,binsize = binsize
loadct,39
formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 600
    ysize = 400
ENDELSE
device,decomposed = 0

IF NOT keyword_set(haloid) THEN haloid = '1'
IF NOT keyword_set(step) THEN step = '00512'
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
spawn,'ls ' + dir + '*' + step + '/*cosmo*' + step + '.halo.' + haloid + '.std',tfile
rtipsy,tfile,h,g,d,s
;spawn,'ls ' + dir + '*' + step + '/*' + step + '.amiga.grp',file_grp
;readarr,file_grp,h,grp,/ascii
;grp_star = grp[h.ngas + h.ndark:h.n - 1]
;grp_gas = grp[0:h.ngas - 1]
spawn,'ls ' + dir + '*' + step + '/*' + step + '.halo.' + haloid + '.iord',file_iord
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

IF keyword_set(outplot) THEN  device,filename = outplot + '_sfh.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
yrange = [0,max([sfh,relosthm])]
plot,timearr,sfh,xtitle = 'Time [Gyr]',ytitle = 'M' + sunsymbol() + textoidl(' yr^{-1}'),psym = 10,thick = 1,yrange = yrange
oplot,timearr,sfh_insitu,color = 254,psym = 10,thick = 1
;oplot,timearr2,relosthm,color = 60,psym = 10,thick = 1
oplot,timearr2,reejecthm,color = 60,psym = 10,thick = 1
IF keyword_set(time1) THEN oplot,[time1,time1],[0,1e5],linestyle = 2
IF keyword_set(time2) THEN oplot,[time2,time2],[0,1e5],linestyle = 2
IF keyword_set(outplot) THEN device, /close ELSE stop


IF keyword_set(outplot) THEN  device,filename = outplot + '_ml_time.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
ml_lost = relosthm/sfh2_insitu
ml_eject = reejecthm/sfh2_insitu
;yrange = [min([ml_lost[where(ml_lost NE 0)],ml_eject[where(ml_lost NE 0)]]),max([ml_lost,ml_eject])]
yrange = [min(ml_eject[where(ml_lost NE 0)]),max(ml_eject)]
plot,timearr2,ml_eject,xtitle = 'Time [Gyr]',ytitle = textoidl('\eta'),psym = 4,/ylog,yrange = yrange,/nodata
oplot,timearr2,ml_eject,psym = 4,color = 254
;oplot,timearr2,ml_lost,psym = 4,color = 60
IF keyword_set(outplot) THEN device, /close ELSE stop

marr2 = spline(time_t,mtot_t,timearr2)
IF keyword_set(outplot) THEN  device,filename = outplot + '_ml_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
ml_lost = relosthm/sfh2_insitu
ml_eject = reejecthm/sfh2_insitu
yrange = [min([ml_lost[where(ml_lost NE 0)],ml_eject[where(ml_lost NE 0)]]),max([ml_lost,ml_eject])]
plot,marr2,ml_eject,xtitle = textoidl('M_{vir} [M') + sunsymbol() + ']',ytitle = textoidl('\eta'),psym = 4,/ylog,yrange = yrange,/xlog,/nodata
oplot,marr2,ml_eject,psym = 4,color = 254
;oplot,marr2,ml_lost,psym = 4,color = 60
IF keyword_set(outplot) THEN device, /close ELSE stop

openw,lun,dir + '/grp' + haloid + '.massloading.txt',/get_lun
printf,lun,'Time [Gyr]','M_{vir} [M_{sun}]','ML_{eject}','ML_{lost}',format='(A18,A18,A18,A18)'
printf,lun,transpose([[timearr2],[marr2],[ml_eject],[ml_lost]])
close,lun

END
