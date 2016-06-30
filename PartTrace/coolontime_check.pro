PRO coolontime_check,dirs,files,halo,finalstep = finalstep,colors = colors
zmax = 3.130 ;Instituted because I am frequently missing steps below 80.
ageUniverse = wmap3_lookback(1000)
IF NOT keyword_set(finalstep) THEN finalstep = '00512'

n = n_elements(dirs)

formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 800
    ysize = 500
 ENDELSE
IF NOT keyword_set(colors) THEN colors = (findgen(n) + 1)*255/n
loadct,39
device,decomposed = 0

FOR i = 6, n_elements(dirs) - 1 DO BEGIN
    units = tipsyunits(dirs[i] + files[i] + '.param')
    rtipsy,dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep,h,g,d,s
    readarr,dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.iord',h,iord,/ascii,part = 'gas',type = 'long'
    readarr,dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.coolontime',h,coolontime,/ascii,part = 'gas',type = 'float'
    coolontime = coolontime*units.timeunit

    reejecti = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_iord.fits',0)
    reejectm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reeject.fits',0)
    reejectz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_z.fits',0)
    reejectt = z_to_t(reejectz)
    reejecti_uniq = reejecti[uniq(reejecti,sort(reejecti))]
    reejectt_uniq = reejectt[uniq(reejecti,sort(reejecti))]

    reexpelli = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_iord.fits',0)
    reexpellm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reexpell.fits',0) 
    reexpellz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_z.fits',0)
    reexpellt = z_to_t(reexpellz)
    reexpelli_uniq = reexpelli[uniq(reexpelli,sort(reexpelli))]
    reexpellt_uniq = reexpellt[uniq(reexpelli,sort(reexpelli))]

    match,iord,reejecti_uniq,indt,indej
    cooloff_ej = coolontime[indt] - reejectt_uniq[indej]*1e9
    match,iord,reexpelli_uniq,indt,index
    cooloff_ex = coolontime[indt] - reexpellt_uniq[index]*1e9
;    IF i EQ 0 THEN 
    histogramp,alog10(cooloff_ej),nbins  = 200,min = 5,max = 10.5,xtitle = 'Log( Cool Off Time [Years])',/normalize
    histogramp,alog10(cooloff_ej),nbins  = 200,min = 5,max = 10.5,color = colors[i],/overplot,/normalize
    histogramp,alog10(cooloff_ex),nbins  = 200,min = 5,max = 10.5,color = colors[i],linestyle = 2,/overplot,/normalize

    diskaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdisk.fits',0)
    diskaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0)
    diskaccri = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_iord.fits',0)
    diskaccrt = z_to_t(diskaccrz)

    IF 0 THEN BEGIN
    ind_longcool = where(cooloff_ex GT 1e9)
;    iord_longcool = iord[indt[ind_longcool]]
    reexpelli_longcool = reexpelli_uniq[index[ind_longcool]]
    reexpellt_longcool = reexpellt_uniq[index[ind_longcool]]
    nreacc = 0
    FOR j = 0, n_elements(iord_longcool) - 1 DO $
       IF max(diskaccrt[where(diskaccri EQ reexpelli_longcool[j])]) GT reexpellt_longcool[j] THEN nreacc = nreacc + 1

    ind_longcool = where(cooloff_ej GT 1e9)
;    iord_longcool = iord[indt[ind_longcool]]
    reejecti_longcool = reejecti_uniq[indej[ind_longcool]]
    reejectt_longcool = reejectt_uniq[indej[ind_longcool]]
    nreacc = 0
    FOR j = 0, n_elements(iord_longcool) - 1 DO $
       IF max(diskaccrt[where(diskaccri EQ reejecti_longcool[j])]) GT reejectt_longcool[j] THEN nreacc = nreacc + 1

    ENDIF


    stop

 ENDFOR

END
