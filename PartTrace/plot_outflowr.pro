;Plot the radius of the particles at the step after being ejected from
;the disk

PRO plot_outflowr,dirs,filenames,finalid = finalid,outplot = outplot
IF NOT keyword_set(finalid) THEN finalid = '1' 
n = n_elements(dirs)
finalstep = '00512'
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

loadct,39
colors = (findgen(n) + 1)*254/n
medianr = findgen(n)
halomass = findgen(n)

FOR i = 0, n - 1 DO BEGIN
   cd,dirs[i]
   align = mrdfits('grp' + finalid[i] + '.alignment.fits',1)
;   vcstep = fltarr(n_elements(align))
   readcol,'grp' + finalid[i] + '.vc.dat',halo,time,redshift,vcstep,format='(d4)'
   ejectdata = mrdfits('grp' + finalid[i] + '.reeject_halo.fits',1)
   ejectz = mrdfits(dirs[i] + 'grp' + finalid[i] + '.reeject_z.fits')
   r = sqrt(ejectdata.x^2 + ejectdata.y^2 + ejectdata.x^2)
   normalr = fltarr(n_elements(r))
   FOR j = 0, n_elements(align) - 1 DO $
         IF (where(ejectz EQ align[j].z))[0] NE -1 THEN $
            normalr[where(ejectz EQ align[j].z)] = float(r[where(ejectz EQ align[j].z)]/align[j].rvir)
   IF i EQ 0 THEN normalr_all = normalr ELSE normalr_all = [normalr_all,normalr]
   IF i EQ 0 THEN histogramp,normalr,nbins = 100,min = 1e-6,max = 1,/normalize
   histogramp,normalr,nbins = 100,min = 1e-6,max = 1,color = colors[i],/overplot,/normalize
   print,'Median r [Rvir]',median(normalr)
   medianr[i]= median(normalr)
;   halomass[i] = align[n_elements(align) - 1].mvir
  stat = read_stat_struc_amiga(dirs[i] + filenames[i] + '.' + finalstep + '/' + filenames[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group eq finalid[i]))[0]
   halomass[i] = stat[ind].m_tot
ENDFOR
stop
END
