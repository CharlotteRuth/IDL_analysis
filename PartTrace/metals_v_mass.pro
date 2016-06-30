;Plots the mass of metals in different states

PRO metals_v_mass,dirs,files,halo,finalstep = finalstep,colors = colors,outplot = outplot,stellarmass = stellarmass,symbols = symbols,formatthick = formatthick

n = n_elements(dirs)

IF NOT keyword_set(finalstep) THEN finalstep = '00512'
IF keyword_set(stellarmass) THEN BEGIN
   xtitle = 'Stellar Mass [M' + sunsymbol() + ']'
   xrange = [1e6,1e11]
ENDIF ELSE BEGIN
   xtitle = 'Virial Mass [M' + sunsymbol() + ']'
   xrange = [3e9,1e12]
ENDELSE

formatplot,outplot = outplot,thick = formatthick

IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
    IF keyword_set(stellarmass) THEN outplot = outplot + '_sm_' ELSE outplot = outplot + '_vm'
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 800
    ysize = 500
ENDELSE
device,decomposed = 0
IF keyword_set(colors) THEN BEGIN
    loadct,39
;    IF n_elements(colors) EQ 0 AND colors[0] EQ 0  THEN colors = [50,254,fgcolor]
    colors = [50,254,120,fgcolor]
ENDIF ELSE BEGIN
    loadct,0    
    colors = [0,0,0,0]
ENDELSE
IF NOT keyword_set(symbols) THEN symbols = [15,14,46,16] ;Halo/Ever, ISM, Stars, ISM+stars
IF keyword_set(formatthick) THEN BEGIN
   thick = 4
   IF NOT keyword_set(symsize) THEN symsize = 2
ENDIF ELSE BEGIN
   thick = 2
   IF NOT keyword_set(symsize) THEN symsize = 1.5
ENDELSE

vmass = fltarr(n)

s_mass = fltarr(n)
s_oxmass = fltarr(n)
s_femass = fltarr(n)

shalo_mass = fltarr(n)
shalo_oxmass = fltarr(n)
shalo_femass = fltarr(n)

srvir_mass = fltarr(n)
srvir_oxmass = fltarr(n)
srvir_femass = fltarr(n)

ism_mass = fltarr(n)
ism_oxmass = fltarr(n)
ism_femass = fltarr(n)

accr_mass = fltarr(n)
accr_oxmass = fltarr(n)
accr_femass = fltarr(n)

accrh_mass = fltarr(n)
accrh_oxmass = fltarr(n)
accrh_femass = fltarr(n)

accr_halo_mass = fltarr(n)
accr_halo_oxmass = fltarr(n)
accr_halo_femass = fltarr(n)

accrh_halo_mass = fltarr(n)
accrh_halo_oxmass = fltarr(n)
accrh_halo_femass = fltarr(n)

grvir_mass = fltarr(n)
grvir_oxmass = fltarr(n)
grvir_femass = fltarr(n)

prevtfile = ''
close,/all
FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   tfile = dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep
   print,tfile + '.' + halo[i]
   IF (tfile NE prevtfile) THEN BEGIN  
      units = tipsyunits(dirs[i] + files[i] + '.param')
      rtipsy,tfile,h,g0,d0,s0
      readarr,tfile + '.iord',h,iord,/ascii,type = 'long'
      iordgas = iord[0:h.ngas - 1]
      iordstar = iord[h.ngas + h.ndark: h.n - 1]
;    readarr,tfile + '.OxMassFrac',h,Ox,/ascii,type = 'float'
      read_tipsy_arr, tfile + '.OxMassFrac', h, Ox,type = 'float'
      ox_g = ox[0:h.ngas - 1]
      ox_s = ox[h.ngas + h.ndark: h.n - 1]
;    readarr,tfile + '.FeMassFrac',h,Fe,/ascii,type = 'float'
      read_tipsy_arr,tfile + '.FeMassFrac',h,Fe,type = 'float'
      fe_g = fe[0:h.ngas - 1]
      fe_s = fe[h.ngas + h.ndark: h.n - 1]
      IF file_test(dirs[i] + files[i] + '.starlog') THEN sl = rstarlog(dirs[i] + files[i] + '.starlog',/molecularH) ELSE BEGIN
         readarr,tfile + '.igasorder',h,igasorder,/ascii,type = 'long',part = 'star'
         sl = {iordergas:0L,iorderstar:0L}
         sl = replicate(sl,h.nstar)
         sl.iordergas = igasorder
         sl.iorderstar = iordstar 
      ENDELSE
   ENDIF
   prevtfile = tfile
   
   stat = read_stat_struc_amiga(tfile + '.amiga.stat')
   ind = (where(stat.group eq halo[i]))[0]
   vmass[i] = stat[ind].m_tot
   rvir = stat[ind].rvir
 
   halodat = mrdfits(dirs[i] + '/grp' + halo[i] + '.alignment.fits',1)
   stepi = where(halodat.file eq files[i]+'.'+finalstep+'/'+files[i]+'.'+finalstep)
;    scale = 1.0/(1.0 + halodat[stepi].z)
   center = fltarr(3)
   center[0] = halodat[stepi].xc ;*1000.0 - units.lengthunit/2.0 ;halodat is in Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
   center[1] = halodat[stepi].yc ;*1000.0 - units.lengthunit/2.0
   center[2] = halodat[stepi].zc ;*1000.0 - units.lengthunit/2.0
   az = reform([[halodat[stepi].xa],[halodat[stepi].ya],[halodat[stepi].za]])
   ax = [az[2]/sqrt(az[0]*az[0] + az[2]*az[2]),0,-1.0*az[0]/sqrt(az[0]*az[0] + az[2]*az[2])]
   ay = crossp(az,ax)           ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
   basis = [[ax],[ay],[az]]
   gpos = [[g0.x*units.lengthunit - center[0]*units.lengthunit],$
           [g0.y*units.lengthunit - center[1]*units.lengthunit],$
           [g0.z*units.lengthunit - center[2]*units.lengthunit]]
   gpos = transpose(transpose(basis)#transpose(gpos))
   g = g0
   g.x = gpos[*,0]
   g.y = gpos[*,1]
   g.z = gpos[*,2]
   spos = [[s0.x*units.lengthunit - center[0]*units.lengthunit],$
           [s0.y*units.lengthunit - center[1]*units.lengthunit],$
           [s0.z*units.lengthunit - center[2]*units.lengthunit]]
   spos = transpose(transpose(basis)#transpose(spos))
   s = s0
   s.x = spos[*,0]
   s.y = spos[*,1]
   s.z = spos[*,2]
   
;;    diskaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdisk.fits',0)
;;    diskaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0)
   diskaccri = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_iord.fits',0)
;;    diskaccrt = z_to_t(diskaccrz)
;;    diskearlym = mrdfits(dirs[i] + '/grp' + halo[i] + '.earlydisk_mass.fits',0)
   diskearlyi = mrdfits(dirs[i] + '/grp' + halo[i] + '.earlydisk_iord.fits',0)
;;    diskaccrm1 = [diskearlym,diskaccrm[uniq(diskaccri,sort(diskaccri))]]
   diskaccri1 = [diskearlyi,diskaccri]
   diskaccri1 = diskaccri[uniq(diskaccri1,sort(diskaccri1))]

;   haloaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccr.fits',0)
;   haloaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccr_z.fits',0)
   haloaccri = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccr_iord.fits',0)
;   haloaccrt = z_to_t(diskaccrz)
;   haloearlym = mrdfits(dirs[i] + '/grp' + halo[i] + '.earlyhalo_mass.fits',0)
   haloearlyi = mrdfits(dirs[i] + '/grp' + halo[i] + '.earlyhalo_iord.fits',0) 
;   halogmass1[i] = total(haloaccrm[uniq(haloaccri,sort(haloaccri))]) + total(haloearlym) 
   haloaccri1 = [haloearlyi,haloaccri]
   haloaccri1 = haloaccri1[uniq(haloaccri1,sort(haloaccri1))]

;------------------------------------------------------------
   match,diskaccri1,iordgas,ind0,diskacc_ind
   match,haloaccri1,iordgas,ind1,haloacc_ind
   g_gal = g[diskacc_ind] ;Gas ever accreted onto disk
   g_halo = g[haloacc_ind] ;Gas ever accreted onto halo

   match2,diskaccri1,sl.iordergas,ind2,sl_ind
   match,sl[where(sl_ind NE -1)].iorderstar,iordstar,ind4,saccr_ind ;Stars formed from material ever accreted onto disk
   accrrvir_ind = where(sqrt(s[saccr_ind].x*s[saccr_ind].x + s[saccr_ind].y*s[saccr_ind].y + s[saccr_ind].z*s[saccr_ind].z) LE rvir)
   s_gal = s[saccr_ind[accrrvir_ind]] ;Stars formed from material ever accreted onto disk within a virial radius
   s_gal_ind = saccr_ind[accrrvir_ind]

   match2,haloaccri1,sl.iordergas,ind3,sl_ind2
   match,sl[where(sl_ind2 NE -1)].iorderstar,iordstar,ind5,shalo_ind ;Stars formed from material ever accreted onto disk
   halorvir_ind = where(sqrt(s[shalo_ind].x*s[shalo_ind].x + s[shalo_ind].y*s[shalo_ind].y + s[shalo_ind].z*s[shalo_ind].z) LE rvir)
   s_halo = s[shalo_ind[halorvir_ind]] ;Stars formed from material ever accreted onto halo within a virial radius
   s_halo_ind = shalo_ind[halorvir_ind]
    
   srvir_ind = where(sqrt(s.x*s.x + s.y*s.y + s.z*s.z) LE rvir) ;Stars within a virial radius
   srvir = s[srvir_ind]

   grvir_ind = where(sqrt(g.x*g.x + g.y*g.y + g.z*g.z) LE rvir) ;Gas within a virial radius
   grvir = g[grvir_ind]

   match,iordgas[grvir_ind],diskaccri1,ind5,ind6 ;Gas ever accreted onto the disk within rvir
   grvir_acc_ind = grvir_ind[ind5]
   grvir_acc = g[grvir_ind[ind5]]

   match,iordgas[grvir_ind],haloaccri1,ind7,ind8 ;Gas ever accreted on the halo within rvir
   grvir_halo_ind = grvir_ind[ind7]
   grvir_halo = g[grvir_ind[ind7]]

   gism_ind = where(g.dens*units.rhounit/(h.time)^3 GE 0.1 AND g.tempg LE 1.2e4 AND abs(g.z) LE 10.0)
   match,iordgas[gism_ind],iordgas[grvir_acc_ind],ind9,ind10 ;Gas ever accreted onto disk, now within rvir as ISM
   gism = g[gism_ind[ind9]]
   gism_ind = gism_ind[ind9]

;--------------------------------
;Stars formed from material ever accreted onto disk within a virial radius  
;   stop
   weight_gal = integrate_ssp((max(s0.tform) - s_gal.tform)*units.timeunit,s_gal.metals)
   s_mass[i] = total(s_gal.mass*weight_gal)*units.massunit
   s_oxmass[i] = total(s_gal.mass*weight_gal*ox_s[s_gal_ind])*units.massunit
   s_femass[i] = total(s_gal.mass*weight_gal*fe_s[s_gal_ind])*units.massunit
;   stop

;Stars formed from material ever accreted onto halo within a virial radius  
   weight_halo = integrate_ssp((max(s0.tform) - s_halo.tform)*units.timeunit,s_halo.metals)
   shalo_mass[i] = total(s_halo.mass*weight_halo)*units.massunit
   shalo_oxmass[i] = total(s_halo.mass*weight_halo*ox_s[s_halo_ind])*units.massunit
   shalo_femass[i] = total(s_halo.mass*weight_halo*fe_s[s_halo_ind])*units.massunit
  
;Stars within a virial radius
   weight_rvir = integrate_ssp((max(s0.tform) - srvir.tform)*units.timeunit,srvir.metals)
   srvir_mass[i] = total(srvir.mass*weight_rvir)*units.massunit
   srvir_oxmass[i] = total(srvir.mass*weight_rvir*ox_s[srvir_ind])*units.massunit
   srvir_femass[i] = total(srvir.mass*weight_rvir*fe_s[srvir_ind])*units.massunit

;Gas ever accreted onto disk, now within rvir as ISM
   ism_mass[i] = total(gism.mass)*units.massunit
   ism_oxmass[i] = total(gism.mass*ox_g[gism_ind])*units.massunit
   ism_femass[i] = total(gism.mass*fe_g[gism_ind])*units.massunit

;Material that was ever accreted onto the disk as gas    
   accr_mass[i] = total(g_gal.mass)*units.massunit;total(s_gal.mass)*units.massunit + 
   accr_oxmass[i] =  total(g_gal.mass*ox_g[diskacc_ind])*units.massunit ;total(s_gal.mass*ox_s[s_gal_ind])*units.massunit +
   accr_femass[i] = total(g_gal.mass*fe_g[diskacc_ind])*units.massunit ;total(s_gal.mass*fe_s[s_gal_ind])*units.massunit + 

;Material that was ever accreted onto the halo as gas    
   accrh_mass[i] = total(g_halo.mass)*units.massunit;total(s_gal.mass)*units.massunit + 
   accrh_oxmass[i] =  total(g_halo.mass*ox_g[haloacc_ind])*units.massunit ;total(s_gal.mass*ox_s[s_gal_ind])*units.massunit +
   accrh_femass[i] = total(g_halo.mass*fe_g[haloacc_ind])*units.massunit ;total(s_gal.mass*fe_s[s_gal_ind])*units.massunit + 

;Gas that was ever accreted onto the disk and now is within rvir
   accr_halo_mass[i] = total(grvir_acc.mass)*units.massunit ;total(s_gal.mass)*units.massunit
   accr_halo_oxmass[i] = total(grvir_acc.mass*ox_g[grvir_acc_ind])*units.massunit ;total(s_gal.mass*ox_s[ind5])*units.massunit
   accr_halo_femass[i] = total(grvir_acc.mass*fe_g[grvir_acc_ind])*units.massunit ;total(s_gal.mass*fe_s[ind5])*units.massunit

;Gas that was ever accreted onto the halo and now is within rvir
   accrh_halo_mass[i] = total(grvir_halo.mass)*units.massunit ;total(s_gal.mass)*units.massunit
   accrh_halo_oxmass[i] = total(grvir_halo.mass*ox_g[grvir_halo_ind])*units.massunit ;total(s_gal.mass*ox_s[ind5])*units.massunit
   accrh_halo_femass[i] = total(grvir_acc.mass*fe_g[grvir_acc_ind])*units.massunit ;total(s_gal.mass*fe_s[ind5])*units.massunit

;Gas that is now within rvir
   grvir_mass[i] = total(grvir.mass)*units.massunit ;total(srvir.mass)*units.massunit
   grvir_oxmass[i] = total(grvir.mass*ox_g[grvir_ind])*units.massunit ;total(srvir.mass*ox_s[srvir_ind])*units.massunit
   grvir_femass[i] = total(grvir.mass*fe_g[grvir_ind])*units.massunit ;total(srvir.mass*fe_s[srvir_ind])*units.massunit


   IF 0 THEN BEGIN
       plot,s_halo.x,s_halo.y,psym = 3
       oplot,srvir.x,srvir.y,psym = 3,color = 50
       oplot,s_halo.x,s_halo.y,psym = 3,color = 255
       oplot,s_gal.x,s_gal.y,psym = 3,color = 254
       stop

       histogramp,srvir.x,nbins = 100
       histogramp,srvir.x,nbins = 100,color = 50,/overplot
       histogramp,s_halo.x,nbins = 100,/overplot
       histogramp,s_gal.x,nbins = 100,color = 254,/overplot
    ENDIF
ENDFOR
;Material that was ever accreted onto the disk as gas
accr_zmass = 2.09*accr_oxmass + 1.06*accr_femass
;Material that was ever accreted onto the halo as gas
accrh_zmass = 2.09*accrh_oxmass + 1.06*accrh_femass
;Gas that was ever accreted onto the disk and now is within rvir
accr_halo_zmass = 2.09*accr_halo_oxmass + 1.06*accr_halo_femass
;Gas that was ever accreted onto the disk and now is within rvir
accrh_halo_zmass = 2.09*accrh_halo_oxmass + 1.06*accrh_halo_femass
;Gas that is now within rvir
grvir_zmass = 2.09*grvir_oxmass + 1.06*grvir_femass
;Gas ever accreted onto disk, now within rvir as ISM
ism_zmass = 2.09*ism_oxmass + 1.06*ism_femass
;Stars formed from material ever accreted onto disk within a virial radius  
s_zmass = 2.09*s_oxmass + 1.06*s_femass
;Stars formed from material ever accreted onto halo within a virial radius
shalo_zmass = 2.09*shalo_oxmass + 1.06*shalo_femass
;Stars within a virial radius
srvir_zmass = 2.09*srvir_oxmass + 1.06*srvir_femass

IF keyword_set(stellarmass) THEN xmass = srvir_mass ELSE xmass = vmass

;Fraction of stellar metals that are in stars formed from accreted material
IF keyword_set(outplot) THEN  device,filename = outplot + '_stellar_insitu_accr.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,s_zmass/srvir_zmass,/xlog,psym = symcat(sym_outline(symbols[0])),xrange = xrange,yrange= [0,1],xtitle = xtitle,ytitle = "In Situ/(Accreted + In Situ)"
oplot,xmass,s_mass/srvir_mass,psym = symcat(symbols[0])
oplot,xmass,shalo_mass/srvir_mass,psym = symcat(symbols[0]),color = 100
oplot,xmass,shalo_zmass/srvir_zmass,psym = symcat(sym_outline(symbols[0])),color = 100
legend,["Stellar Mass","Stellar Metals"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of metals in the halo that were added to once-disk-accreted material
IF keyword_set(outplot) THEN  device,filename = outplot + '_halo_accr.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accr_halo_zmass + s_zmass)/(accr_zmass + s_zmass),/xlog,psym = symcat(sym_outline(symbols[0])),xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = "Gas and Stars in Halo from Disk/Gas Ever in Disk"
oplot,xmass,(accr_halo_mass + s_mass)/(accr_mass + s_mass),psym = symcat(symbols[0])
legend,["Gas Mass","Gas Metals"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of metals in the halo that were added to once-halo-accreted material
IF keyword_set(outplot) THEN  device,filename = outplot + '_halo_accrh.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accrh_halo_zmass + shalo_zmass)/(accrh_zmass + shalo_zmass),/xlog,psym = symcat(sym_outline(symbols[0])),xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = "Gas and Stars in Halo from Halo/Gas Ever in Halo"
oplot,xmass,(accrh_halo_mass + shalo_mass)/(accrh_mass + shalo_mass),psym = symcat(symbols[0])
legend,["Gas Mass","Gas Metals"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of mass within rvir that is bound/was accreted to disk
IF keyword_set(outplot) THEN  device,filename = outplot + '_halo_rvir.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,accrh_halo_mass/grvir_mass,xtitle = xtitle,xrange = xrange,/xlog,psym = symcat(sym_outline(symbols[0])),yrange = [0,1],ytitle = "Accreted Gas within Rvir/Gas within Rvir"
oplot,xmass,accr_halo_mass/grvir_mass,psym = symcat(symbols[0])
legend,["Accreted to Disk","Accreted to Halo"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop


;Fraction of metals in the halo that were added to once-accreted material
IF keyword_set(outplot) THEN  device,filename = outplot + '_halo_accrh.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,accrh_halo_zmass/accrh_zmass,/xlog,psym = symcat(sym_outline(symbols[0])),xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = "Gas in Halo/Gas Ever in Halo"
oplot,xmass,accrh_halo_mass/accrh_mass,psym = symcat(symbols[0])
legend,["Gas Mass","Gas Metals"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop


;Fraction of metals produced in accreted material that are still within rvir
IF keyword_set(outplot) THEN  device,filename = outplot + '_zfrac_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accr_zmass + s_zmass)/(accr_zmass + s_zmass),/xlog,xrange = xrange,yrange = [0,0.8],/nodata,xtitle = xtitle,ytitle = textoidl('M_z/M_z avaliable')
oplot,xmass,(accr_halo_zmass + s_zmass)/(accrh_zmass + s_zmass),psym = symcat(symbols[0]),color = colors[0],symsize = symsize
oplot,xmass,(accr_halo_zmass + s_zmass)/(accrh_zmass + s_zmass),psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass,(s_zmass + ism_zmass)/(accrh_zmass + s_zmass),psym = symcat(symbols[3]),color = 190,symsize = symsize
oplot,xmass,(s_zmass + ism_zmass)/(accrh_zmass + s_zmass),psym = symcat(sym_outline(symbols[3])),symsize = symsize
oplot,xmass,ism_zmass/(accrh_zmass + s_zmass),psym = symcat(symbols[1]),color = colors[1],symsize = symsize
oplot,xmass,ism_zmass/(accrh_zmass + s_zmass),psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,xmass,s_zmass/(accrh_zmass + s_zmass),psym = symcat(symbols[2]),color = colors[2],symsize = symsize
oplot,xmass,s_zmass/(accrh_zmass + s_zmass),psym = symcat(sym_outline(symbols[2])),symsize = symsize
legend,["Within Rvir","Stars + ISM","Stars","ISM"],color = [colors[0],190,colors[2],colors[1]],psym = [symbols[0],symbols[3],symbols[2],symbols[1]],box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of metals produced in accreted material that are still within rvir
IF keyword_set(outplot) THEN  device,filename = outplot + '_zfrac_smass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,shalo_mass,(accr_zmass + s_zmass)/(accr_zmass + s_zmass),/xlog,xrange = xrange,yrange = [0,1],/nodata,xtitle = 'Stellar Mass [M' + sunsymbo() + ']',ytitle = textoidl('M_z/M_z avaliable')
oplot,shalo_mass,(accr_halo_zmass + s_zmass)/(accr_zmass + s_zmass),psym = symcat(symbols[0]),color = colors[0],symsize = symsize
oplot,shalo_mass,(accr_halo_zmass + s_zmass)/(accr_zmass + s_zmass),psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,shalo_mass,(s_zmass + ism_zmass)/(accr_zmass + s_zmass),psym = symcat(symbols[3]),color = colors[3],symsize = symsize
oplot,shalo_mass,(s_zmass + ism_zmass)/(accr_zmass + s_zmass),psym = symcat(sym_outline(symbols[3])),symsize = symsize
oplot,shalo_mass,ism_zmass/(accr_zmass + s_zmass),psym = symcat(symbols[1]),color = colors[1],symsize = symsize
oplot,shalo_mass,ism_zmass/(accr_zmass + s_zmass),psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,shalo_mass,s_zmass/(accr_zmass + s_zmass),psym = symcat(symbols[2]),color = colors[2],symsize = symsize
oplot,shalo_mass,s_zmass/(accr_zmass + s_zmass),psym = symcat(sym_outline(symbols[2])),symsize = symsize
legend,["Within Rvir","Stars + ISM","Stars","ISM"],color = [colors[0],colors[3],colors[2],colors[1]],psym = [symbols[0],symbols[3],symbols[2],symbols[1]],box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop


;Fraction of metals produced in accreted material that are still within rvir
IF keyword_set(outplot) THEN  device,filename = outplot + '_zfrach_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accrh_zmass + shalo_zmass)/(accrh_zmass + shalo_zmass),/xlog,xrange = xrange,yrange = [0,1],/nodata,xtitle = xtitle,ytitle = textoidl('M_z/M_z Halo Accretion')
oplot,xmass,(accrh_halo_zmass + shalo_zmass)/(accrh_zmass + shalo_zmass),psym = symcat(symbols[0]),color = colors[0],symsize = symsize
oplot,xmass,(accrh_halo_zmass + shalo_zmass)/(accrh_zmass + shalo_zmass),psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass,(shalo_zmass + ism_zmass)/(accrh_zmass + shalo_zmass),psym = symcat(symbols[3]),color = colors[3],symsize = symsize
oplot,xmass,(shalo_zmass + ism_zmass)/(accrh_zmass + shalo_zmass),psym = symcat(sym_outline(symbols[3])),symsize = symsize
oplot,xmass,ism_zmass/(accrh_zmass + shalo_zmass),psym = symcat(symbols[1]),color = colors[1],symsize = symsize
oplot,xmass,ism_zmass/(accrh_zmass + shalo_zmass),psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,xmass,shalo_zmass/(accrh_zmass + shalo_zmass),psym = symcat(symbols[2]),color = colors[2],symsize = symsize
oplot,xmass,shalo_zmass/(accrh_zmass + shalo_zmass),psym = symcat(sym_outline(symbols[2])),symsize = symsize
legend,["Within Rvir","Stars + ISM","Stars","ISM"],color = [colors[0],colors[3],colors[2],colors[1]],psym = [symbols[0],symbols[3],symbols[2],symbols[1]],box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_massfrac_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accr_mass + s_mass)/(accr_mass + s_mass),/xlog,xrange = xrange,yrange = [0,1],/nodata,xtitle = xtitle,ytitle = textoidl('Mass/Mass Accreted to Disk As Gas')
oplot,xmass,(accr_halo_mass + s_mass)/(accr_mass + s_mass),psym = symcat(symbols[0]),color = colors[0],symsize = symsize
oplot,xmass,(accr_halo_mass + s_mass)/(accr_mass + s_mass),psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass,(s_mass + ism_mass)/(accr_mass + s_mass),psym = symcat(symbols[3]),color = colors[3],symsize = symsize
oplot,xmass,(s_mass + ism_mass)/(accr_mass + s_mass),psym = symcat(sym_outline(symbols[3])),symsize = symsize
oplot,xmass,ism_mass/(accr_mass + s_mass),psym = symcat(symbols[1]),color = colors[1],symsize = symsize
oplot,xmass,ism_mass/(accr_mass + s_mass),psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,xmass,s_mass/(accr_mass + s_mass),psym = symcat(symbols[2]),color = colors[2],symsize = symsize
oplot,xmass,s_mass/(accr_mass + s_mass),psym = symcat(sym_outline(symbols[2])),symsize = symsize
legend,["Within Rvir","Stars + ISM","Stars","ISM"],color = [colors[0],colors[3],colors[2],colors[1]],psym = [symbols[0],symbols[3],symbols[2],symbols[1]],box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

IF keyword_set(outplot) THEN  device,filename = outplot + '_ztotal_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,shalo_mass,(accrh_zmass + shalo_zmass),xtitle = 'Stellar Mass [M' + sunsymbol() + ']',ytitle = textoidl("Log M_{Z}/M") + sunsymbol(),/nodata,/xlog,/ylog,xrange = [1e6,5e10],yrange = [1e2,1e9];[min(ism_zmass),max(accr_zmass)]
oplot,shalo_mass,(accrh_zmass + shalo_zmass),psym = symcat(sym_outline(symbols[3])),color = colors[3],symsize = symsize
;oplot,shalo_mass,(shalo_zmass + ism_zmass),psym = symcat(symbols[3]),color = colors[3],symsize = symsize
;oplot,shalo_mass,(shalo_zmass + ism_zmass),psym = symcat(sym_outline(symbols[3])),color = fgcolor,symsize = symsize
oplot,shalo_mass,shalo_zmass,psym = symcat(symbols[2]),color = colors[2],symsize = symsize
oplot,shalo_mass,shalo_zmass,psym = symcat(sym_outline(symbols[2])),color = fgcolor,symsize = symsize
oplot,shalo_mass,ism_zmass,psym = symcat(symbols[1]),color = colors[1],symsize = symsize
oplot,shalo_mass,ism_zmass,psym = symcat(sym_outline(symbols[1])),color = fgcolor,symsize = symsize
legend,["Ever in Galaxy","Disk","Stars","ISM"],color = [colors[3],colors[3],colors[2],colors[1]],psym = [sym_outline(symbols[3]),symbols[3],symbols[2],symbols[1]]
IF keyword_set(outplot) THEN device, /close ELSE stop

END
