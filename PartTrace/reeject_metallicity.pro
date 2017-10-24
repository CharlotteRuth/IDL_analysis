;------------------------- Metallicity of feedback comes from ---------
PRO reeject_metallicity,dir,filenames,finalid = finalid,scale = scale,normalize = normalize,outplot = outplot,keys = keys,colors = colors,linestyles = linestyles,thicks = thicks,label = label,ctables = ctables,formatthick = formatthick, absolute = absolute,avez = avez,red_avez = red_avez
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
;    loadct,11
;    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
;    rainbow_colors
;    loadcv,80,/reverse;,/noqual
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(outplot) THEN thicks = fltarr(n) + 5 ELSE thicks = fltarr(n) + 2
;       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 4
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) 
ENDIF ELSE BEGIN
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(outplot) THEN thicks = fltarr(n) + 5 ELSE thicks = fltarr(n) + 2
;       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = findgen(n)*2  
ENDELSE
avez_struct = {meanz: 0.0, medianz: 0.0, modez: 0.0, stdz: 0.0}
avez = replicate(avez_struct,n_elements(dir))
IF keyword_set(absolute) THEN outext = '_fbz_abs' ELSE outext = '_fbz'
IF keyword_set(normalize) THEN outext = outext + '_norm.eps' ELSE outext = outext + '.eps'

IF keyword_set(outplot) THEN device,filename = outplot + outext,/encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8 ELSE window,0,ysize=ysize,xsize=xsize

z_bins = [2,1,0.5,-1e-10]
ageUniverse = wmap3_lookback(1000)
t_bins = (ageUniverse - wmap3_lookback(z_bins))/1e9
mlbin = 1.0 ;in billion
red_avez_struct = {mvir:0.0, mstar:0.0, ZISM:0.0, Zeject:0.0} 
red_avez = replicate(red_avez_struct,n_elements(dir),n_elements(t_bins)) ;galaxy, redshift bins, (virial mass, stellar mass, average ISM metallicity, average ejecta metallicity)

FOR i = 0, n_elements(dir)-1 DO BEGIN
;   loadct,ctables[i]
   loadct,0
   cd,dir[i]
   splitbase = strsplit(filenames[i],'.')
   base = strmid(filenames[i],0,splitbase[n_elements(splitbase) - 1] - 1)

   totalmass = 1
;  eject_gas_character,dirhead[i],ejecthistory,expellhistory,totalmass = totalmass
   ejecthistory =  mrdfits(dir[i] + 'grp' + finalid[i] + '.reeject_halo.fits',1)
   ejectz = mrdfits(dir[i] + 'grp' + finalid[i] + '.reeject_z.fits')
   ejectt = z_to_t(ejectz)
   print,dir[i] + 'grp' + finalid[i] + '.reeject_z.fits'

   align = mrdfits('grp' + finalid[i] + '.alignment.fits',1)
   readcol,'grp' + finalid[i] + '.mass_metal_track.dat',halo_t,time_t,z_t,mtot_t,mgas_t,mstar_t,mdark_t,metal_t,ox_t,fe_t,HI_t,H2_t,coldg_t,/silent
   readcol,'grp' + finalid[i] + '.metals.txt',z,ox,fe,coldgas,z_H,ox_H,fe_H,H2,HI,H
   stepmetallicity = z/coldgas
   IF keyword_set(colors) THEN linecolor = (alog10(mtot_t[n_elements(mtot_t) - 1]) - 9.5)/2.5*(256-13) + 13  ELSE linecolor = 0
   normalmetallicity = fltarr(n_elements(ejecthistory))
   FOR j = 0, n_elements(z) - 1 DO $
       normalmetallicity[where(ejectz EQ align[j].z)] = ejecthistory[where(ejectz EQ align[j].z)].metallicity/stepmetallicity[j]

;   expellhistory = mrdfits(dir[i] + 'grp' + finalid[i] + '.expell_disk.fits',1)
  FOR iz = 0, n_elements(z_bins) - 1 DO BEGIN
     IF t_bins[iz] EQ n_elements(t_bins) - 1 THEN mlbin2 = mlbin ELSE mlbin2 = mlbin/2
     indtime = where(time_t GE t_bins[iz] - mlbin2 AND time_t LE t_bins[iz] + mlbin2)
     ;Average virial mass
     red_avez[i,iz].mvir = mean(mtot_t[indtime])
     ;Average stellar mass
     red_avez[i,iz].mstar = mean(mstar_t[indtime])
     ;Average ISM metallicity
     red_avez[i,iz].ZISM = total(z[indtime])/total(coldgas[indtime])/zsolar
     indeject = where(ejectt GE t_bins[iz] - mlbin2 AND ejectt LE t_bins[iz] + mlbin2)
     red_avez[i,iz].Zeject = total(ejecthistory[indeject].metallicity*ejecthistory[indeject].mass)/total(ejecthistory[indeject].mass)/zsolar
  ENDFOR
  loadct,0
   IF keyword_set(absolute) THEN BEGIN
       IF i EQ 0 AND     keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,xrange = xrange_fb,xtitle = 'Log(Z/Z'+sunsymbol()+')',ytitle = '1/M dM/dr',title = label, /nodata,yrange = [0,0.08]
       IF i EQ 0 AND NOT keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            xrange = xrange_fb,xtitle = 'Log(Z/Z'+sunsymbol()+')',ytitle = 'dM/dr'    ,title = label, /nodata ;yrange=[0,2e9] ;,xmargin = [16,3]
       loadcv,80,/reverse;,/noqual
       IF                keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,/overplot,color = linecolor,thick = thicks[i],linestyle = linestyles[i]
       IF NOT            keyword_set(normalize) THEN histogramp,alog10(ejecthistory.metallicity/zsolar), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            /overplot,color = linecolor,thick = thicks[i],linestyle = linestyles[i]
       medianz = median(ejecthistory.metallicity)/zsolar
       meanz = total(ejecthistory.metallicity*ejecthistory.mass)/total(ejecthistory.mass)/zsolar
       temp = max(weighted_histogram(alog10(ejecthistory.metallicity/zsolar),weight = ejecthistory.mass*scale[i],min = xrange_fb[0],max = xrange_fb[1],locations = locations,nbins = 100),maxindex)
       modez = locations[maxindex]
;       oplot,[alog10(medianz),alog10(medianz)],[0,1],color = linecolor,thick = thicks[i],linestyle = 1
       avez[i].medianz = medianz
       avez[i].meanz = meanz
       avez[i].modez = 10^modez
       avez[i].stdz = stdev(ejecthistory.metallicity/zsolar)
   ENDIF ELSE BEGIN
       IF i EQ 0 AND     keyword_set(normalize) THEN histogramp,alog10(normalmetallicity), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,xrange = xrange_fb,xtitle = textoidl('Log(Z_{eject}/Z_{gas})'),ytitle = '1/M dM/dr',title = label, /nodata,yrange = [0,0.08]
       IF i EQ 0 AND NOT keyword_set(normalize) THEN histogramp,alog10(normalmetallicity), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            xrange = xrange_fb,xtitle = textoidl('Log(Z_{eject}/Z_{gas})'),ytitle = 'dM/dr'    ,title = label, /nodata ;yrange=[0,2e9] ;,xmargin = [16,3]
       oplot,[0,0],[0,1],linestyle = 2
       loadcv,80,/reverse;,/noqual
       IF                keyword_set(normalize) THEN histogramp,alog10(normalmetallicity), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i], /normalize,/overplot,color = linecolor,thick = thicks[i],linestyle = linestyles[i]
       IF NOT            keyword_set(normalize) THEN histogramp,alog10(normalmetallicity), nbins = 100,min = xrange_fb[0],max = xrange_fb[1],weight = ejecthistory.mass*scale[i],            /overplot,color = linecolor,thick = thicks[i],linestyle = linestyles[i]
       medianz = median(normalmetallicity)
       meanz = total(normalmetallicity*ejecthistory.mass)/total(ejecthistory.mass)
       temp = max(weighted_histogram(alog10(normalmetallicity),weight = ejecthistory.mass*scale[i],min = xrange_fb[0],max = xrange_fb[1],locations = locations,nbins = 100),maxindex)
       modez = 10^locations[maxindex]
;       oplot,[alog10(medianz),alog10(medianz)],[0,1],color = linecolor,thick = thicks[i],linestyle = 1
       avez[i].medianz = medianz
       avez[i].meanz = meanz
       avez[i].modez = modez
       avez[i].stdz = stdev(normalmetallicity)
    ENDELSE
   print,'Metallicity enhancement (mean, median, mode, stdev): ',avez[i].meanz,avez[i].medianz,avez[i].modez,avez[i].stdz
   undefine,totalmass
   undefine,ejecthistory
   undefine,expellhistory
ENDFOR    
cgColorbar,range=[9.5,12],/vertical,position=[0.90, 0.22, 0.92, 0.9],title = 'Log Virial Mass [M'+sunsymbol()+']',divisions = 5,tcharsize = 1.5,color = "black",ncolors = 256-13,bottom = 13 ;256 ;"Black"
IF keyword_set(keys) THEN legend,keys,color = colors,linestyle = linestyles,thick = thicks,/top,/left,box = 0
IF keyword_set(outplot) THEN device,/close ELSE stop
END
