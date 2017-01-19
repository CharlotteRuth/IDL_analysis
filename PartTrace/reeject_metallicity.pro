;------------------------- Metallicity of feedback comes from ---------
PRO reeject_metallicity,dir,filenames,finalid = finalid,scale = scale,normalize = normalize,outplot = outplot,keys = keys,colors = colors,linestyles = linestyles,thicks = thicks,label = label,ctables = ctables,formatthick = formatthick, absolute = absolute,modez = modez
formatplot,outplot = outplot,thick = formatthick
n = n_elements(dir)
zsolar  =  0.0130215
;metalscale = 2 ;The code was run the metallicity bug where zmetal
;just is ox+fe.  The actual metallicity is about twice as high (z =
;1.06*fe + 2.09*ox). Updated 8/13/16 except this is already accounted
;for in find_history.pro
IF keyword_set(absolute) THEN xrange_fb = [-2.5,0.5] ELSE xrange_fb = [-0.5,1.0]
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
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
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
modez = fltarr(n_elements(dir))
IF keyword_set(absolute) THEN outext = '_fbz_abs' ELSE outext = '_fbz'
IF keyword_set(normalize) THEN outext = outext + '_norm.eps' ELSE outext = outext + '.eps'

IF keyword_set(outplot) THEN device,filename = outplot + outext,/encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8 ELSE window,0,ysize=ysize,xsize=xsize

FOR i = 0, n_elements(dir)-1 DO BEGIN
   loadct,ctables[i]
   cd,dir[i]
   splitbase = strsplit(filenames[i],'.')
   base = strmid(filenames[i],0,splitbase[n_elements(splitbase) - 1] - 1)

   totalmass = 1
;  eject_gas_character,dirhead[i],ejecthistory,expellhistory,totalmass = totalmass
   ejecthistory =  mrdfits(dir[i] + 'grp' + finalid[i] + '.reeject_halo.fits',1)
   ejectz = mrdfits(dir[i] + 'grp' + finalid[i] + '.reeject_z.fits')
   print,dir[i] + 'grp' + finalid[i] + '.reeject_z.fits'

   align = mrdfits('grp' + finalid[i] + '.alignment.fits',1)
   readcol,'grp' + finalid[i] + '.metals.txt',z,ox,fe,coldgas,z_H,ox_H,fe_H,H2,HI,H
   stepmetallicity = z/coldgas

   normalmetallicity = fltarr(n_elements(ejecthistory))
   FOR j = 0, n_elements(z) - 1 DO $
       normalmetallicity[where(ejectz EQ align[j].z)] = ejecthistory[where(ejectz EQ align[j].z)].metallicity/stepmetallicity[j]

;   expellhistory = mrdfits(dir[i] + 'grp' + finalid[i] + '.expell_disk.fits',1)

   IF keyword_set(absolute) THEN BEGIN
       IF i EQ 0 AND     keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,xrange = xrange_fb,xtitle = 'Log(Z/Z'+sunsymbol()+')',ytitle = '1/M dM/dr',title = label, /nodata,yrange = [0,0.08]
       IF i EQ 0 AND NOT keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            xrange = xrange_fb,xtitle = 'Log(Z/Z'+sunsymbol()+')',ytitle = 'dM/dr'    ,title = label, /nodata ;yrange=[0,2e9] ;,xmargin = [16,3]
       IF                keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,/overplot,color = colors[i],thick = thicks[i],linestyle = linestyles[i]
       IF NOT            keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            /overplot,color = colors[i],thick = thicks[i],linestyle = linestyles[i]
;   histogramp,                                              alog10(expellhistory.metallicity/zsolar),nbins = 100,max = xrange_fb[1],weight = expellhistory.mass*scale[i],           /overplot,color = colors[i],thick = thicks[i],linestyle = 2
       print,max(weighted_histogram(alog10(ejecthistory.metallicity/zsolar), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],locations = yaxis),maxindex)
;       oplot,[yaxis[maxindex],yaxis[maxindex]],[0,1],color = colors[i],thick = thicks[i],linestyle = 1
;       modez[i] = 10^(yaxis[maxindex])
;       medianz = total(ejecthistory.metallicity/zsolar*ejecthistory.mass*scale)/total(ejecthistory.mass*scale)
       medianz = median(ejecthistory.metallicity/zsolar)
       oplot,[alog10(medianz),alog10(medianz)],[0,1],color = colors[i],thick = thicks[i],linestyle = 1
       modez[i] = medianz
   ENDIF ELSE BEGIN
       IF i EQ 0 AND     keyword_set(normalize) THEN histogramp,alog10(normalmetallicity), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,xrange = xrange_fb,xtitle = textoidl('Log(Z_{eject}/Z_{gas})'),ytitle = '1/M dM/dr',title = label, /nodata,yrange = [0,0.08]
       IF i EQ 0 AND NOT keyword_set(normalize) THEN histogramp,alog10(normalmetallicity), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            xrange = xrange_fb,xtitle = textoidl('Log(Z_{eject}/Z_{gas})'),ytitle = 'dM/dr'    ,title = label, /nodata ;yrange=[0,2e9] ;,xmargin = [16,3]
       oplot,[0,0],[0,1],linestyle = 2
       IF                keyword_set(normalize) THEN histogramp,alog10(normalmetallicity), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,/overplot,color = colors[i],thick = thicks[i],linestyle = linestyles[i]
       IF NOT            keyword_set(normalize) THEN histogramp,alog10(normalmetallicity), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            /overplot,color = colors[i],thick = thicks[i],linestyle = linestyles[i]
;   histogramp,                                              alog10(expellhistory.metallicity/zsolar),nbins = 100,max = xrange_fb[1],weight = expellhistory.mass*scale[i],           /overplot,color = colors[i],thick = thicks[i],linestyle = 2
       print,max(weighted_histogram(alog10(normalmetallicity),nbins = 100, min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],locations = yaxis),maxindex)
;       oplot,[alog10(mean(normalmetallicity)),alog10(mean(normalmetallicity))],[0,1],color = colors[i],thick = thicks[i],linestyle = 2
;       oplot,[yaxis[maxindex],yaxis[maxindex]],[0,1],color = colors[i],thick = thicks[i],linestyle = 1
;       modez[i] = 10^(yaxis[maxindex])
       medianz = median(normalmetallicity)
       oplot,[alog10(medianz),alog10(medianz)],[0,1],color = colors[i],thick = thicks[i],linestyle = 1
       modez[i] = medianz
   ENDELSE
   undefine,totalmass
   undefine,ejecthistory
   undefine,expellhistory
ENDFOR    
IF keyword_set(keys) THEN legend,keys,color = colors,linestyle = linestyles,thick = thicks,/top,/left,box = 0
IF keyword_set(outplot) THEN device,/close ELSE stop
END
