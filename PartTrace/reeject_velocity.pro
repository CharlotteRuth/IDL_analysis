;------------------------- Metallicity of feedback comes from ---------
PRO reeject_velocity,dirs,filenames,finalid = finalid,scale = scale,normalize = normalize,outplot = outplot,keys = keys,colors = colors,linestyles = linestyles,thicks = thicks,label = label,ctables = ctables,formatthick = formatthick,expell = expell,vscale = vscale,vs_step = vs_step,medv_ej = medv_ej,stdv_ej = stdv_ej,medv_exp = medv_exp,stdv_exp = stdv_exp,vfinal = vfinal
formatplot,outplot = outplot,thick = formatthick
finalstep = '00512'
n = n_elements(dirs)
IF NOT keyword_set(finalid) THEN finalid = '1'
IF NOT keyword_set(scale) THEN scale = fltarr(n) + 1
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
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF colors[0] EQ 1 THEN  colors = (findgen(n) + 1)*240/n ELSE colors = colors
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = findgen(n)*2  
ENDELSE
vc = fltarr(n_elements(dirs))
mvir = fltarr(n_elements(dirs))
medv_ej = fltarr(n_elements(dirs))
stdv_ej = fltarr(n_elements(dirs))
medv_exp = fltarr(n_elements(dirs))
stdv_exp = fltarr(n_elements(dirs))

IF keyword_set(outplot) THEN BEGIN
    IF keyword_set(normalize)  THEN $
      IF keyword_set(vscale) OR keyword_set(vs_step) THEN device,filename = outplot + '_fbvscale_norm.eps',/encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8 ELSE device,filename = outplot + '_fbv_norm.eps',/encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8 $
    ELSE $
      IF keyword_set(vscale) OR keyword_set(vs_step) THEN device,filename = outplot + '_fbvscale.eps',     /encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8 ELSE device,filename = outplot + '_fbv.eps',     /encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8
 ENDIF ELSE window,0

FOR i = 0, n_elements(dirs)-1 DO BEGIN
   loadct,ctables[i]
   cd,dirs[i]
   align = mrdfits('grp' + finalid[i] + '.alignment.fits',1)
   splitbase = strsplit(filenames[i],'.')
   base = strmid(filenames[i],0,splitbase[n_elements(splitbase) - 1] - 1)
   stat = read_stat_struc_amiga(dirs[i] + filenames[i] + '.' + finalstep + '/' + filenames[i] + '.' +  finalstep + '.amiga.stat')
   ind = (where(stat.group eq finalid[i]))[0]
;   vc[i] = stat[ind].vc
   vc[i] = sqrt(6.67e-11*stat[ind].m_tot*1.9891e30/(3.08567758e19*stat[ind].rvir))/1000.
   IF keyword_set(vfinal) THEN BEGIN
;   plot,vc,vfinal/vc,psym = 4,xtitle = 'V_circ from virial mass and radius',ytitle = 'v_btf/v_circ',yrange = [1,2]
      vc = vfinal
   ENDIF
   IF keyword_set(vs_step) THEN BEGIN
      vcstep = fltarr(n_elements(align))
;      readcol,'grp' + finalid[i] + '.vc_amiga.dat',halo,time,redshift,vcstep,format='(d4)'
      readcol,'grp' + finalid[i] + '.vc.dat',halo,time,redshift,vcstep,format='(d4)'
 ;     stop
   ENDIF  
   mvir[i] = stat[ind].m_tot
   totalmass = 1
;  eject_gas_character,dirhead[i],ejecthistory,expellhistory,totalmass = totalmass
   ejecthistory =  mrdfits(dirs[i] + 'grp' + finalid[i] + '.reeject_halo.fits',1) ;properties of the ejected gas in the step immediately after it is in the disk
;   ejecthistory =  mrdfits(dirs[i] + 'grp' + finalid[i] + '.reeject_disk.fits',1)
   ejectz = mrdfits(dirs[i] + 'grp' + finalid[i] + '.reeject_z.fits')
;   expellhistory = mrdfits(dirs[i] + 'grp' + finalid[i] + '.expell_disk.fits',1)
   velocity = sqrt(ejecthistory.vx^2 + ejecthistory.vy^2 + ejecthistory.vz^2)
   r = sqrt(ejecthistory.x^2 + ejecthistory.y^2 + ejecthistory.z^2)
   energy = velocity
;   stop
   IF keyword_set(vs_step) THEN BEGIN
       normalvelocity = fltarr(n_elements(velocity))
       FOR j = 0, n_elements(vcstep) - 1 DO $
         IF (where(ejectz EQ align[j].z))[0] NE -1 THEN $
         normalvelocity[where(ejectz EQ align[j].z)] = float(velocity[where(ejectz EQ align[j].z)]/vcstep[j]) ELSE IF j GT n_elements(align.z)  - 100 THEN print,'problem '
   ENDIF
   IF file_test(dirs[i]+'grp'+finalid[i]+'.eject_expell_iord.fits') THEN BEGIN
       eject_expell_iord = mrdfits(dirs[i]+'grp'+finalid[i]+'.eject_expell_iord.fits')
       expellhistory = ejecthistory[eject_expell_iord]
       velocity_expell = sqrt(expellhistory.vx^2 + expellhistory.vy^2 + expellhistory.vz^2)
   ENDIF
   IF keyword_set(vs_step) THEN BEGIN
       normalvelocity_expell = fltarr(n_elements(velocity_expell))
       FOR j = 0, n_elements(vcstep) - 1 DO $
         IF (where(ejectz[eject_expell_iord] EQ align[j].z))[0] NE -1 THEN $
           normalvelocity_expell[where(ejectz[eject_expell_iord] EQ align[j].z)] = velocity_expell[where(ejectz[eject_expell_iord] EQ align[j].z)]/vcstep[j]
    ENDIF
;   stop
   mcolor = (alog10(mvir[i]) - 9.5)/2.5*256
   print,mcolor

   IF keyword_set(vscale) OR keyword_set(vs_step) THEN BEGIN
       IF keyword_set(vscale) THEN BEGIN 
           velocity = velocity/vc[i]
           velocity_expell = velocity_expell/vc[i]
       ENDIF
       vmax = 4;5;3
       ymax = 0.035;5;0.028;0.045;0.028
       xtitle = textoidl('V/V_{circ}')
   ENDIF ELSE BEGIN
       vmax = 300
       ymax = 0.17
       xtitle = 'Velocity [km/s]'
    ENDELSE
;   stop

;   IF 0 THEN BEGIN
   IF keyword_set(vs_step) THEN BEGIN
       IF i EQ 0 AND     keyword_set(normalize) THEN histogramp,normalvelocity, nbins = 100,max = vmax,weight = ejecthistory.mass*scale[i]/total(ejecthistory.mass*scale[i]),xrange = [0,vmax],xtitle = xtitle,ytitle = '1/M dM/dr',title = label, /nodata,yrange = [0,ymax]
       IF i EQ 0 AND NOT keyword_set(normalize) THEN histogramp,normalvelocity, nbins = 100,max = vmax,weight = ejecthistory.mass*scale[i],            xrange = [0,vmax],xtitle = xtitle,ytitle = 'dM/dr'    ,title = label, /nodata ;yrange=[0,2e9] ;,xmargin = [16,3]
       IF                keyword_set(normalize) THEN histogramp,normalvelocity, nbins = 100,max = vmax,weight = ejecthistory.mass*scale[i]/total(ejecthistory.mass*scale[i]),/overplot,color = mcolor,thick = thicks[i],linestyle = linestyles[i]
       IF NOT            keyword_set(normalize) THEN histogramp,normalvelocity, nbins = 100,max = vmax,weight = ejecthistory.mass*scale[i],            /overplot,color = mcolor,thick = thicks[i],linestyle = linestyles[i]
;   histogramp,                                              alog10(expellhistory.metallicity/zsolar),nbins = 100,max = xrange_fb[1],weight = expellhistory.mass*scale[i],           /overplot,color = colors[i],thick = thicks[i],linestyle = 2
       IF                keyword_set(normalize) AND keyword_set(expell) THEN histogramp,normalvelocity_expell, nbins = 100,max = vmax,weight = expellhistory.mass*scale[i]/total(ejecthistory.mass*scale[i]),/overplot,color = mcolor,thick = thicks[i],linestyle = 2
       IF NOT            keyword_set(normalize) AND keyword_set(expell) THEN histogramp,normalvelocity_expell, nbins = 100,max = vmax,weight = expellhistory.mass*scale[i],            /overplot,color = mcolor,thick = thicks[i],linestyle = 2
   ENDIF ELSE BEGIN       
       IF i EQ 0 AND     keyword_set(normalize) THEN histogramp,velocity, nbins = 100,max = vmax,weight = ejecthistory.mass*scale[i]/total(ejecthistory.mass*scale[i]),xrange = [0,vmax],xtitle = xtitle,ytitle = '1/M dM/dr',title = label, /nodata,yrange = [0,ymax]
       IF i EQ 0 AND NOT keyword_set(normalize) THEN histogramp,velocity, nbins = 100,max = vmax,weight = ejecthistory.mass*scale[i],            xrange = [0,vmax],xtitle = xtitle,ytitle = 'dM/dr'    ,title = label, /nodata ;yrange=[0,2e9] ;,xmargin = [16,3]
       IF                keyword_set(normalize) THEN histogramp,velocity, nbins = 100,max = vmax,weight = ejecthistory.mass*scale[i]/total(ejecthistory.mass*scale[i]),/overplot,color = mcolor,thick = thicks[i],linestyle = linestyles[i]
       IF NOT            keyword_set(normalize) THEN histogramp,velocity, nbins = 100,max = vmax,weight = ejecthistory.mass*scale[i],            /overplot,color = mcolor,thick = thicks[i],linestyle = linestyles[i]
;   histogramp,                                              alog10(expellhistory.metallicity/zsolar),nbins = 100,max = xrange_fb[1],weight = expellhistory.mass*scale[i],           /overplot,color = colors[i],thick = thicks[i],linestyle = 2
       IF                keyword_set(normalize) AND keyword_set(expell) THEN histogramp,velocity_expell, nbins = 100,max = vmax,weight = expellhistory.mass*scale[i]/total(ejecthistory.mass*scale[i]),/overplot,color = mcolor,thick = thicks[i],linestyle = 1
       IF NOT            keyword_set(normalize) AND keyword_set(expell) THEN histogramp,velocity_expell, nbins = 100,max = vmax,weight = expellhistory.mass*scale[i],            /overplot,color = mcolor,thick = thicks[i],linestyle = 2
   ENDELSE
;   ENDIF

   IF 0 THEN BEGIN
      y_eject =  weighted_histogram(normalvelocity,        nbins = 100,max = vmax,min = 1e-6,weight =  ejecthistory.mass*scale[i]/total(ejecthistory.mass*scale[i]),location = xarray)
      y_expell = weighted_histogram(normalvelocity_expell, nbins = 100,max = vmax,min = 1e-6,weight = expellhistory.mass*scale[i]/total(ejecthistory.mass*scale[i]),location = xarray)   
      IF i EQ 0 THEN plot,xarray,y_expell/y_eject,/nodata
      oplot,xarray,y_expell/y_eject,color = mcolor,thick = thicks[i],linestyle = linestyles[i]
      oplot,[0,vmax],[0.9,0.9]
   ENDIF

   IF keyword_set(vscale)  OR keyword_set(vs_step) THEN oplot,[1,1],[0,1] ELSE oplot,[stat[ind].vc,stat[ind].vc],[0,1],color = mcolor,thick = thicks[i],linestyle = 1
   IF keyword_set(keys) THEN legend,keys,color = colors,linestyle = linestyles,thick = thicks,/top,/right,box = 0
 ;  stop

   IF keyword_set(vs_step) THEN BEGIN
       medv_ej[i] = median(normalvelocity) 
       stdv_ej[i] = stdev(normalvelocity)
       medv_exp[i] = median(normalvelocity_expell) 
       stdv_exp[i] = stdev(normalvelocity_expell)
   ENDIF ELSE BEGIN
       medv_ej[i] = median(velocity)
       stdv_ej[i] = stdev(velocity)
       medv_exp[i] = median(velocity_expell)
       stdv_exp[i] = stdev(velocity_expell)
   ENDELSE
;   oplot,[medv_ej[i],medv_ej[i]],[0,1],color = mcolor,linestyle = 3
   IF keyword_set(expell) THEN oplot,[medv_exp[i],medv_exp[i]],[0,1],color = mcolor,linestyle = 2
   print,max(weighted_histogram(velocity,weight = totalmass/1e8,max = vmax,nbins = 100))
   undefine,totalmass
   undefine,ejecthistory
;   undefine,expellhistory
ENDFOR    
cgColorbar,range=[9.5,12],/vertical,position=[0.90, 0.34, 0.92, 0.90],title = 'Log Virial Mass',divisions = 5
IF keyword_set(outplot) THEN device,/close ELSE stop
END
