;Plots the mass of metals in different states
;.r metals_produced.pro
;.r /home1/crchrist/IDL/IDL_analysis/procedures/legend.pro
PRO metals_v_mass,dirs,files,halo,finalstep = finalstep,colors = colors,outplot = outplot,stellarmass = stellarmass,symbols = symbols,formatthick = formatthick

n = n_elements(dirs)
zsolar  =  0.0130215
IF NOT keyword_set(finalstep) THEN finalstep = '00512'
IF finalstep EQ '00512' THEN redshift  = '0'
IF finalstep EQ '00128' THEN redshift  = '2'
IF keyword_set(stellarmass) THEN BEGIN
   xtitle = 'Stellar Mass [M' + sunsymbol() + ']'
   IF finalstep EQ '00128' THEN xrange = [5e5,1e10] ELSE xrange = [1e6,1e11]
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
    IF keyword_set(stellarmass) THEN outplot = outplot + '_sm' ELSE outplot = outplot + '_vm'
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 800
    ysize = 500
ENDELSE
device,decomposed = 0
IF keyword_set(colors) THEN BEGIN
;    loadct,39
    distinct_colors,n_colors = 12
;    IF n_elements(colors) EQ 0 AND colors[0] EQ 0  THEN colors = [50,254,fgcolor]
    gascolor = 10
    starcolor = 5
    diskcolor = 7
    cgmcolor = 2
    colors = [cgmcolor,gascolor,starcolor,diskcolor]
ENDIF ELSE BEGIN
    loadct,0    
    colors = [0,0,0,0]
ENDELSE
IF NOT keyword_set(symbols) THEN symbols = [14,15,46,16] ;Halo/Ever, ISM, Stars, ISM+stars
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
srvir_remmass = fltarr(n)
srvir_remz = fltarr(n)

ism_mass = fltarr(n)
ism_oxmass = fltarr(n)
ism_femass = fltarr(n)

HIH2_mass = fltarr(n)
HIH2_oxmass = fltarr(n)
HIH2_femass = fltarr(n)

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

gr150_mass = fltarr(n)
gr150_oxmass = fltarr(n)
gr150_femass = fltarr(n)

primehalo_stars_oxmass = fltarr(n)
primehalo_stars_femass = fltarr(n)
currenthalo_stars_oxmass = fltarr(n)
currenthalo_stars_femass = fltarr(n)

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
      IF file_test(tfile + '.massform') THEN read_tipsy_arr,tfile + '.massform',h,massform,part = 'star',type = 'float' ELSE BEGIN
;          rtipsy,tfile,h,g,d,s
;          massform = s.mass
          match,iordstar,sl.iorderstar,ind0,ind1
          massform = sl[ind1].massform
      ENDELSE
   ENDIF
   prevtfile = tfile
   
   readcol,dirs[i] + files[i] + '.grp' + halo[i] + '.haloid.dat', file_steps, halo_steps, format='a,l'
   ind_step = where(file_steps EQ files[i] + '.' + finalstep + '/' + files[i] + '.' + finalstep)
   IF ind_step EQ -1 THEN stop
   halo_step = halo_steps[ind_step]
   stat = read_stat_struc_amiga(tfile + '.amiga.stat')
   ind = (where(stat.group eq halo_step[0]))[0]
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

 ;  inhalo = [0]
 ;  IF keyword_set(centralhalo) THEN BEGIN
 ;     test = where(gpart.grp[i] EQ stat[main].group, ntest)
 ;     IF ntest NE 0 THEN inhalo = [inhalo,test]
 ;  ENDIF ELSE BEGIN
 ;     FOR j = 0, n_elements(satellites) - 1 DO BEGIN
 ;        test = where(gpart.grp[i] EQ stat[satellites[j]].group, ntest)
 ;        IF ntest NE 0 THEN inhalo = [inhalo,test]
 ;     ENDFOR
 ;  ENDELSE
   
;;    diskaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdisk.fits',0)
   diskaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0)
   diskaccri = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_iord.fits',0)
   IF keyword_set(redshift) THEN diskaccri = diskaccri[where(diskaccrz GE (1/h.time-1)-0.0022)] ;allow for slight inaccuracies in redshift
;;    diskaccrt = z_to_t(diskaccrz)
;;    diskearlym = mrdfits(dirs[i] + '/grp' + halo[i] + '.earlydisk_mass.fits',0)
   diskearlyi = mrdfits(dirs[i] + '/grp' + halo[i] + '.earlydisk_iord.fits',0)
;;    diskaccrm1 = [diskearlym,diskaccrm[uniq(diskaccri,sort(diskaccri))]]
   diskaccri1 = [long(diskearlyi),long(diskaccri)]
   diskaccri1 = diskaccri1[uniq(diskaccri1,sort(diskaccri1))]

;   haloaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccr.fits',0)
   haloaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccr_z.fits',0)
   haloaccri = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccr_iord.fits',0)
   IF keyword_set(redshift) THEN haloaccri = haloaccri[where(haloaccrz GE (1/h.time-1)-0.0022)]
;   haloaccrt = z_to_t(diskaccrz)
;   haloearlym = mrdfits(dirs[i] + '/grp' + halo[i] + '.earlyhalo_mass.fits',0)
   haloearlyi = mrdfits(dirs[i] + '/grp' + halo[i] + '.earlyhalo_iord.fits',0) 
;   halogmass1[i] = total(haloaccrm[uniq(haloaccri,sort(haloaccri))]) + total(haloearlym) 
   haloaccri1 = [long(haloearlyi),long(haloaccri)]
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
   massform_srvir = massform[srvir_ind]

; Calculate the metals produced by the stars currently in the halo
   dDelta = 0.000691208*units.timeunit ;516, 799 0.000691088, 603, 986 (0.000691208), 239 (0.000691208)
   tmax = dDelta*float(finalstep)
   srvir.tform = srvir.tform*units.timeunit
   which_imf = 0
   z_snii_srvir = zsnovaii(0, tmax + dDelta, srvir, massform_srvir*units.massunit, which_imf)
   srvir_remmass[i] = total(z_snii_srvir.stellarremmass)
   srvir_remz[i] = total(z_snii_srvir.stellarremz)
   z_snia_srvir = zsnovaia(0, tmax + dDelta, srvir, massform_srvir*units.massunit)
   currenthalo_stars_oxmass[i] = z_snii_srvir.oxmassloss + z_snia_srvir.oxmassloss
   currenthalo_stars_femass[i] = z_snii_srvir.femassloss + z_snia_srvir.femassloss
   readcol,dirs[i] + '/grp' + halo[i] + '.sfmetals.txt',zstars,oxstars,festars ;This doesn't work at higher redshifts
   primehalo_stars_oxmass[i] = total(oxstars)
   primehalo_stars_femass[i] = total(festars)

   grvir_ind = where(sqrt(g.x*g.x + g.y*g.y + g.z*g.z) LE rvir) ;Gas within a virial radius
   grvir = g[grvir_ind]

   gr150_ind = where(sqrt(g.x*g.x + g.y*g.y + g.z*g.z) LE 150) ;Gas within 150 kpc from the center
   gr150 = g[gr150_ind]

   match,iordgas[grvir_ind],diskaccri1,ind5,ind6 ;Gas ever accreted onto the disk within rvir
   grvir_acc_ind = grvir_ind[ind5]
   grvir_acc = g[grvir_ind[ind5]]

   match,iordgas[grvir_ind],haloaccri1,ind7,ind8 ;Gas ever accreted on the halo within rvir
   grvir_halo_ind = grvir_ind[ind7]
   grvir_halo = g[grvir_ind[ind7]]
   IF n_elements(ind7) NE n_elements(iordgas[grvir_ind]) THEN BEGIN ;Why are particles not being included in accretion onto the halo?
       print,'i: ',strtrim(i,2),' ',tfile,'halo: ',halo[i]
       print,'# particles in halo: ',strtrim(n_elements(iordgas[grvir_ind]),2),', # accreted particles in halo: ',strtrim(n_elements(ind7),2),', difference: ',strtrim(n_elements(iordgas[grvir_ind]) - n_elements(ind7),2),' (',strtrim(100*(n_elements(iordgas[grvir_ind]) - n_elements(ind7))/float(n_elements(iordgas[grvir_ind])),2),'%)'
       match2,iordgas[grvir_ind],haloaccri1,ind7_2,ind8_2
       missing_iord = (iordgas[grvir_ind])[where(ind7_2 EQ -1)] ;iords of Particles in the halo that were not included in the accretion
                                ;Some of these particles are probably being
                                ;missed because accretion only counts
                                ;if it is in the main halo or a
                                ;satellite of the main halo (just
                                ;within rvir doesn't count)
                                ;The rest are likely missed because
                                ;they are hot enough to be unbound
                                ;from the halo
       read_tipsy_arr,tfile + '.amiga.grp',h,grp,type = 'int',part = 'gas'
       grp_missing = (grp[grvir_ind])[where(ind7_2 EQ -1)] 
       print,grp_missing[uniq(grp_missing,sort(grp_missing))]
;print,missing_iord[where(grp_missing eq halo[i])]
       missing_g = (g[grvir_ind])[where(ind7_2 EQ -1)] ;Gas particles in the halo that were not included in accretion
       haloaccri_all =  mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccr_iord.fits',0)
       haloaccri1_all = [haloearlyi,haloaccri_all]
       haloaccrz1_all = [fltarr(n_elements(haloearlyi)) + 9,haloaccrz]
       match2,iordgas[grvir_ind],haloaccri1_all,ind7_all_2,ind8_all_2
       missing_iord_all = (iordgas[grvir_ind])[where(ind7_all_2 EQ -1)] ;iords of Particles in the halo that were not included in the accretion
       missing_g_all = (g[grvir_ind])[where(ind7_all_2 EQ -1)] ;Gas particles in the halo that were not included in accretion
       plot,grvir.x,grvir.y,psym =3 ;Particles in halo
       oplot,g_halo.x,g_halo.y,psym = 3,color = 3 ;Particles accreted to halo
       oplot,grvir_halo.x,grvir_halo.y,psym = 3,color = 6 ;Particles accreted to the halo, currently in the halo
       oplot,missing_g.x,missing_g.y,psym = 3,color = 12 ;Particles that are missing from the accreted particles
       oplot,missing_g_all.x,missing_g_all.y,psym = 3,color = 8 ;Particles that are missing from the particles accreted at any time
;       stop
   ENDIF

   gdisk_ind = where(abs(g.z) LE 3.0 AND sqrt(g.x*g.x + g.y*g.y + g.z*g.z) LE rvir) ;A physical cut of the disk to be used when scaling by HI and H2 (e.g. HIH2mass)
   gdisk = g[gdisk_ind]
   gism_ind = where(g.dens*units.rhounit/(h.time)^3 GE 0.1 AND g.tempg LE 1.2e4 AND abs(g.z) LE 3.0 AND sqrt(g.x*g.x + g.y*g.y + g.z*g.z) LE rvir)
;Matching the gas that meets the definition of ISM is unnecessary and
;potentially problematic. Better to just do away with the requirment
;that we have it being accreted onto the ISM. -- CRC 7/20/17
   ;match,iordgas[gism_ind],iordgas[grvir_acc_ind],ind9,ind10 ;Gas ever accreted onto disk, now within rvir as ISM
   gism = g[gism_ind] ;[ind9]]-- CRC 7/20/17
   ;gism_ind = gism_ind[ind9]-- CRC 7/20/17

   readarr,tfile+'.H2',h,H2,/ascii,type = 'float',part = 'gas'
   readarr,tfile+'.HI',h,HI,/ascii,type = 'float',part = 'gas'
   ;readarr,tfile+'.HeI',h,HeI,/ascii,type = 'float',part = 'gas'
   print,'Stellar mass: ',alog10(total(srvir.mass)*units.massunit)
   print,'Disk gass mass: ',total(gism.mass*units.massunit)
   ;print,'HI + H2 + HeI: ',total((H2[gism_ind]*2 + HI[gism_ind] + HeI[gism_ind])*gism.mass*units.massunit)
   print,'Metal mass in ISM: ',total(gism.mass*(2.09*ox_g[gism_ind] + 1.06*fe_g[gism_ind]))*units.massunit
   metalmassinrvir = total(grvir_acc.mass*(2.09*ox_g[grvir_acc_ind] + 1.06*fe_g[grvir_acc_ind]))*units.massunit 
   print,'Metal mass within gas in Rvir: ',metalmassinrvir
   print,'Fraction of metals in Rvir to metals in ISM: ',(total(gism.mass*(2.09*ox_g[gism_ind] + 1.06*fe_g[gism_ind]))*units.massunit)/metalmassinrvir
;   stop
;--------------------------------
;Stars formed from material ever accreted onto disk within a virial radius  
;   stop
;   weight_gal = integrate_ssp((max(s0.tform) - s_gal.tform)*units.timeunit,s_gal.metals) ;in ms_mass.pro
   weight_gal = s_gal.mass*0 + 1
   s_mass[i] = total(s_gal.mass*weight_gal)*units.massunit
   s_oxmass[i] = total(s_gal.mass*weight_gal*ox_s[s_gal_ind])*units.massunit
   s_femass[i] = total(s_gal.mass*weight_gal*fe_s[s_gal_ind])*units.massunit
;   stop

;Stars formed from material ever accreted onto halo within a virial radius  
;   weight_halo = integrate_ssp((max(s0.tform) - s_halo.tform)*units.timeunit,s_halo.metals)
   weight_halo = s_halo.mass*0 + 1
   shalo_mass[i] = total(s_halo.mass*weight_halo)*units.massunit
   shalo_oxmass[i] = total(s_halo.mass*weight_halo*ox_s[s_halo_ind])*units.massunit
   shalo_femass[i] = total(s_halo.mass*weight_halo*fe_s[s_halo_ind])*units.massunit
  
;Stars within a virial radius
;   weight_rvir = integrate_ssp((max(s0.tform) - srvir.tform)*units.timeunit,srvir.metals)
   weight_rvir = srvir.mass*0 + 1
   srvir_mass[i] = total(srvir.mass*weight_rvir)*units.massunit
   srvir_oxmass[i] = total(srvir.mass*weight_rvir*ox_s[srvir_ind])*units.massunit
   srvir_femass[i] = total(srvir.mass*weight_rvir*fe_s[srvir_ind])*units.massunit

;Gas now within rvir as ISM
   ism_mass[i] = total(gism.mass)*units.massunit
   ism_oxmass[i] = total(gism.mass*ox_g[gism_ind])*units.massunit
   ism_femass[i] = total(gism.mass*fe_g[gism_ind])*units.massunit

;Gas scaled by HI + H2 fraction with abs(z) < 3
   maxHfrac = max(HI + H2)
   HIH2_mass[i] = total(gdisk.mass*(HI[gdisk_ind] + H2[gdisk_ind])/maxHfrac)*units.massunit
   HIH2_oxmass[i] = total(gdisk.mass*(HI[gdisk_ind] + H2[gdisk_ind])/maxHfrac*ox_g[gdisk_ind])*units.massunit
   HIH2_femass[i] = total(gdisk.mass*(HI[gdisk_ind] + H2[gdisk_ind])/maxHfrac*fe_g[gdisk_ind])*units.massunit

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
   accrh_halo_femass[i] = total(grvir_halo.mass*fe_g[grvir_halo_ind])*units.massunit ;total(s_gal.mass*fe_s[ind5])*units.massunit

;Gas that is now within rvir
   grvir_mass[i] = total(grvir.mass)*units.massunit ;total(srvir.mass)*units.massunit
   grvir_oxmass[i] = total(grvir.mass*ox_g[grvir_ind])*units.massunit ;total(srvir.mass*ox_s[srvir_ind])*units.massunit
   grvir_femass[i] = total(grvir.mass*fe_g[grvir_ind])*units.massunit ;total(srvir.mass*fe_s[srvir_ind])*units.massunit

;Gas that is now within 150 kpc
   gr150_mass[i] = total(gr150.mass)*units.massunit ;total(sr150.mass)*units.massunit
   gr150_oxmass[i] = total(gr150.mass*ox_g[gr150_ind])*units.massunit ;total(sr150.mass*ox_s[sr150_ind])*units.massunit
   gr150_femass[i] = total(gr150.mass*fe_g[gr150_ind])*units.massunit ;total(sr150.mass*fe_s[sr150_ind])*units.massunit


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
;Gas that was ever accreted onto the halo and now is within rvir
accrh_halo_zmass = 2.09*accrh_halo_oxmass + 1.06*accrh_halo_femass
;Gas that is now within rvir
grvir_zmass = 2.09*grvir_oxmass + 1.06*grvir_femass
;Gas that is now within 150 kpc of the center
gr150_zmass = 2.09*gr150_oxmass + 1.06*gr150_femass
;Gas now within rvir as ISM
ism_zmass = 2.09*ism_oxmass + 1.06*ism_femass
;Gas now within rvir and 3 kpc of disk plane, scaled by HI+H2 abundance
HIH2_zmass = 2.09*HIH2_oxmass + 1.06*HIH2_femass
;Stars formed from material ever accreted onto disk within a virial radius  
s_zmass = 2.09*s_oxmass + 1.06*s_femass
;Stars formed from material ever accreted onto halo within a virial radius
shalo_zmass = 2.09*shalo_oxmass + 1.06*shalo_femass
;Stars within a virial radius
srvir_zmass = 2.09*srvir_oxmass + 1.06*srvir_femass
;Metals produced by all stars in the virial radius at z = 0
currenthalo_stars_zmass = 2.09*currenthalo_stars_oxmass + 1.06*currenthalo_stars_femass
;Metals produced by stars in the main halo
primehalo_stars_zmass = 2.09*primehalo_stars_oxmass + 1.06*primehalo_stars_femass

writecol,'~/Zproduced_z0.txt',files,halo,currenthalo_stars_zmass,primehalo_stars_zmass,(accrh_zmass + srvir_zmass - srvir_remz),(srvir_zmass - srvir_remz),ism_zmass,(accr_halo_zmass + s_zmass),format = '(A,I,E,E,E,E,E,E)'
;readcol,'~/Zproduced_z0.txt',files,halo,currenthalo_stars_zmass,primehalo_stars_zmass,accrh_shalo_zmass,s_zmass,ism_zmass,accr_halo_s_zmass,format = '(A,I,F,F,F,F,F,F)'

;;;;;;;;;;;;; Stars Plots ;;;;;;;;;;;;;;;;
IF keyword_set(stellarmass) THEN xmass = srvir_mass ELSE xmass = vmass
;Verify stellar metallicity. Match with Gallazzi, Panel 8
IF keyword_set(outplot) THEN  device,filename = outplot + '_stellar_metallicity.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
logstellarmetallicity = alog10((srvir_zmass - srvir_remz)/(srvir_mass - srvir_remmass)/zsolar)
logstellarmetallicity_obs = (logstellarmetallicity+0.16)/1.08 ;Peeples 2014, Eq 7
plot,xmass,logstellarmetallicity,/xlog,psym = symcat(symbols[0]),xrange = xrange,yrange= [-1.6,0],xtitle = xtitle,ytitle = 'Stellar Metallicity'
oplot,xmass,logstellarmetallicity_obs,psym = symcat(sym_outline(symbols[0])),color = 100
legend,['Stellar Metallicity','Stellar Metallicity from Peeples+ 14, Eq 7'],psym = [symbols[0],sym_outline(symbols[0])],color = [fgcolor,100]
IF keyword_set(outplot) THEN device, /close ELSE stop

zmetals_predict_SNII = 10^(1.0146*alog10(srvir_mass - srvir_remmass) + alog10(0.030) + 0.1091)
zmetals_predict_SNIa = 10^(1.043*alog10(srvir_mass - srvir_remmass) - 2.683)
;zmetals_predict_SNII = (srvir_massform)*0.030
zmetals_predict = (srvir_mass - srvir_remmass)^1.004166*0.059

IF keyword_set(outplot) THEN  device,filename = outplot + '_totalmetals_comp.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,(srvir_mass - srvir_remmass),currenthalo_stars_zmass,psym = symcat(sym_outline(symbols[0])),/xlog,xtitle = 'Current Stellar Mass [M' + sunsymbol() + ']',ytitle = 'Metal Mass',/ylog,yrange = [1e2,1e11],xrange = xrange
oplot,(srvir_mass - srvir_remmass),zmetals_predict,psym = 5
oplot,(srvir_mass - srvir_remmass),(srvir_zmass - srvir_remz),psym = symcat(sym_outline(symbols[0])),color = 254
oplot,(srvir_mass - srvir_remmass),ism_zmass,psym = symcat(sym_outline(symbols[0])),color = 60
legend,['Metals produced by stars','Metals predicted','Metals in stars','Metals in ISM'],psym = [sym_outline(symbols[0]),5,sym_outline(symbols[0]),sym_outline(symbols[0])],color = [fgcolor,fgcolor,254,60]
IF keyword_set(outplot) THEN device, /close ELSE stop


;Fraction of stellar metals that are in stars formed from accreted material
IF keyword_set(outplot) THEN  device,filename = outplot + '_stellar_insitu_accr.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,s_zmass/srvir_zmass,/xlog,psym = symcat(sym_outline(symbols[0])),xrange = xrange,yrange= [0,1],xtitle = xtitle,ytitle = "In Situ/(Accreted + In Situ)"
oplot,xmass,s_mass/srvir_mass,psym = symcat(symbols[0])
oplot,xmass,shalo_mass/srvir_mass,psym = symcat(symbols[0]),color = 100
oplot,xmass,shalo_zmass/srvir_zmass,psym = symcat(sym_outline(symbols[0])),color = 100
legend,["Stellar Mass","Stellar Metals"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right,box = 0
legend,["In Situ (Accreted to Disk)","In Situ (Accreted to Halo)"],psym = [symbols[0],symbols[0]],color = [fgcolor,100],/bottom,/left,box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of metals in the halo that were added to once-disk-accreted material
;Fraction of metals accreted onto disk that were accreted onto the halo
IF keyword_set(outplot) THEN  device,filename = outplot + '_halo_accr.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accr_halo_zmass + s_zmass)/(accr_zmass + s_zmass),/xlog,psym = symcat(sym_outline(symbols[0])),xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = "Gas and Stars in Halo from Disk/Gas Ever in Disk"
oplot,xmass,(accr_halo_mass + s_mass)/(accr_mass + s_mass),psym = symcat(symbols[0])
legend,["Gas Mass","Gas Metals"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of metals in the halo that were added to once-halo-accreted material
;Fraction of metals that were accreted to halo that are currently in
;the halo
IF keyword_set(outplot) THEN  device,filename = outplot + '_halo_accrh.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accrh_halo_zmass + shalo_zmass)/(accrh_zmass + shalo_zmass),/xlog,psym = symcat(sym_outline(symbols[0])),xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = "Gas and Stars in Halo/Gas Ever in Halo"
oplot,xmass,(accrh_halo_mass + shalo_mass)/(accrh_mass + shalo_mass),psym = symcat(symbols[0])
legend,["Gas Mass","Gas Metals"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of mass within rvir that is bound/was accreted to disk
;Why isn't accrh_halo_mass/grvir_mass always 1? (Check, in
;particular, for low mass halos)
IF keyword_set(outplot) THEN  device,filename = outplot + '_halo_rvir.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,accrh_halo_mass/grvir_mass,xtitle = xtitle,xrange = xrange,/xlog,psym = symcat(sym_outline(symbols[0])),yrange = [0,1],ytitle = "Accreted Gas within Rvir/Gas within Rvir"
oplot,xmass,accr_halo_mass/grvir_mass,psym = symcat(symbols[0])
legend,["Accreted to Disk","Accreted to Halo"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop


;Fraction of metals in the halo that were added to once-accreted material
IF keyword_set(outplot) THEN  device,filename = outplot + '_halo_accrh_frac.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,accrh_halo_zmass/accrh_zmass,/xlog,psym = symcat(sym_outline(symbols[0])),xrange = xrange,yrange = [0,1],xtitle = xtitle,ytitle = "Gas in Halo/Gas Ever in Halo"
oplot,xmass,accrh_halo_mass/accrh_mass,psym = symcat(symbols[0])
legend,["Gas Mass","Gas Metals"],psym = [symbols[0],sym_outline(symbols[0])],/bottom,/right
IF keyword_set(outplot) THEN device, /close ELSE stop

;Comparing different ways to total the metals in the halo
IF keyword_set(outplot) THEN  device,filename = outplot + '_totalmetals.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,primehalo_stars_zmass,/xlog,psym = symcat(15),xrange = xrange,xtitle = xtitle,ytitle = "Total Metal Mass",/ylog,yrange = [1e4,1e9],/nodata;[min([primehalo_stars_zmass,currenthalo_stars_zmass,(accrh_zmass + s_zmass)]),max([primehalo_stars_zmass,currenthalo_stars_zmass,(accrh_zmass + s_zmass)])]
if redshift EQ 0 THEN oplot,xmass,primehalo_stars_zmass,psym = symcat(15) ;only valid at z = 0
oplot,xmass,currenthalo_stars_zmass,psym = symcat(sym_outline(15))
oplot,xmass,accrh_zmass + s_zmass,psym = symcat(4)
legend,["Metals Produced by Stars in Main Progenitor","Metals Produced by Stars in Final Halo","Metals in All Material Ever in Halo"],psym = [15,sym_outline(15),4],/bottom,/right,box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop
print,'Current/Ever: ',currenthalo_stars_zmass/primehalo_stars_zmass
print,'All Material/Ever: ',(accrh_zmass + s_zmass)/primehalo_stars_zmass

;Fraction of metals produced in accreted material that are still within rvir
;Add in metals accreted
IF keyword_set(outplot) THEN  device,filename = outplot + '_zfrac_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accr_zmass + s_zmass)/(accr_zmass + s_zmass),/xlog,xrange = xrange,yrange = [0,1],/nodata,xtitle = xtitle,ytitle = textoidl('M_z/M_z avaliable')
oplot,xmass,(accrh_halo_zmass + s_zmass)/(accrh_zmass + s_zmass),psym = symcat(symbols[0]),color = colors[0],symsize = symsize
oplot,xmass,(accrh_halo_zmass + s_zmass)/(accrh_zmass + s_zmass),psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass,ism_zmass/(accrh_zmass + s_zmass),psym = symcat(symbols[1]),color = colors[1],symsize = symsize
oplot,xmass,ism_zmass/(accrh_zmass + s_zmass),psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,xmass,s_zmass/(accrh_zmass + s_zmass),psym = symcat(symbols[2]),color = colors[2],symsize = symsize
oplot,xmass,s_zmass/(accrh_zmass + s_zmass),psym = symcat(sym_outline(symbols[2])),symsize = symsize
oplot,xmass,(s_zmass + ism_zmass)/(accrh_zmass + s_zmass),psym = symcat(symbols[3]),color = colors[3],symsize = symsize
oplot,xmass,(s_zmass + ism_zmass)/(accrh_zmass + s_zmass),psym = symcat(sym_outline(symbols[3])),symsize = symsize
legend,["Within Rvir","Stars + ISM","Stars","ISM"],color = [colors[0],colors[3],colors[2],colors[1]],psym = [symbols[0],symbols[3],symbols[2],symbols[1]],box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of metals produced in accreted material that are still within rvir
;Add in metals accreted
IF keyword_set(outplot) THEN  device,filename = outplot + '_zfrac_mass_hist.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,alog10(xmass),(accr_zmass + s_zmass)/(accr_zmass + s_zmass),xrange = alog10(xrange),yrange = [0,1],/nodata,xtitle = textoidl('Log(M_*/M')+sunsymbol()+')',ytitle = textoidl('M_z/M_z avaliable')
binsize = 0.1
xcoords = [[alog10(xmass) - binsize/2],[alog10(xmass) - binsize/2],[alog10(xmass) + binsize/2],[alog10(xmass) + binsize/2]]
frac_inhalo = (accrh_halo_zmass + s_zmass)/(accrh_zmass + s_zmass)
frac_ism = ism_zmass/(accrh_zmass + s_zmass)
frac_star = s_zmass/(accrh_zmass + s_zmass)
;FOR i = 0, n_elements(xmass) - 1 DO polyfill,xcoords[i,*],[0,frac_star[i],frac_star[i],0],color = colors[2],/line_fill,orientation = 45,spacing = 0.1
;FOR i = 0, n_elements(xmass) - 1 DO oplot,xcoords[i,*],[0,frac_star[i],frac_star[i],0],color = colors[2],thick = 5
;FOR i = 0, n_elements(xmass) - 1 DO polyfill,xcoords[i,*],[frac_star[i],frac_star[i]+frac_ism[i],frac_star[i]+frac_ism[i],frac_star[i]],color = colors[1],/line_fill,orientation = 60,spacing = 0.1
;FOR i = 0, n_elements(xmass) - 1 DO oplot,xcoords[i,*],[frac_star[i],frac_star[i]+frac_ism[i],frac_star[i]+frac_ism[i],frac_star[i]],color = colors[1],thick = 5
;FOR i = 0, n_elements(xmass) - 1 DO polyfill,xcoords[i,*],[frac_star[i]+frac_ism[i],frac_inhalo[i],frac_inhalo[i],frac_star[i]+frac_ism[i]],color = colors[0],/line_fill,orientation = 30,spacing = 0.1
;FOR i = 0, n_elements(xmass) - 1 DO oplot,xcoords[i,*],[frac_star[i]+frac_ism[i],frac_inhalo[i],frac_inhalo[i],frac_star[i]+frac_ism[i]],color = colors[0],thick = 5

FOR i = 0, n_elements(xmass) - 1 DO BEGIN &  $
polyfill,xcoords[i,*],[0,frac_star[i],frac_star[i],0],color = colors[2],/line_fill,orientation = 45,spacing = 0.1 &  $
oplot,xcoords[i,*],[0,frac_star[i],frac_star[i],0],color = colors[2],thick = 5 &  $
polyfill,xcoords[i,*],[frac_star[i],frac_star[i]+frac_ism[i],frac_star[i]+frac_ism[i],frac_star[i]],color = colors[1],/line_fill,orientation = 60,spacing = 0.1 &  $
oplot,xcoords[i,*],[frac_star[i],frac_star[i]+frac_ism[i],frac_star[i]+frac_ism[i],frac_star[i]],color = colors[1],thick = 5 &  $
polyfill,xcoords[i,*],[frac_star[i]+frac_ism[i],frac_inhalo[i],frac_inhalo[i],frac_star[i]+frac_ism[i]],color = colors[0],/line_fill,orientation = 30,spacing = 0.1 &  $
oplot,xcoords[i,*],[frac_star[i]+frac_ism[i],frac_inhalo[i],frac_inhalo[i],frac_star[i]+frac_ism[i]],color = colors[0],thick = 5 &  $
ENDFOR


legend,["Within Rvir","Stars","Disk Gas"],color = [colors[0],colors[2],colors[1]],linestyl = [0,0,0],box = 0
IF keyword_set(redshift) THEN legend,['z = ' + strtrim(redshift,2)],box = 0,/right,/top
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of metals produced by stars in current main halo in that are still within rvir
IF keyword_set(outplot) THEN  device,filename = outplot + '_zcurrfrac_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accrh_zmass + s_zmass)/(currenthalo_stars_zmass),/xlog,xrange = xrange,yrange = [0,1],/nodata,xtitle = xtitle,ytitle = textoidl('M_z/M_z avaliable')
oplot,xmass,(accrh_halo_zmass + s_zmass)/(currenthalo_stars_zmass),psym = symcat(symbols[0]),color = colors[0],symsize = symsize
oplot,xmass,(accrh_halo_zmass + s_zmass)/(currenthalo_stars_zmass),psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass,ism_zmass/(currenthalo_stars_zmass),psym = symcat(symbols[1]),color = colors[1],symsize = symsize
oplot,xmass,ism_zmass/(currenthalo_stars_zmass),psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,xmass,(srvir_zmass - srvir_remz)/(currenthalo_stars_zmass),psym = symcat(symbols[2]),color = colors[2],symsize = symsize
oplot,xmass,(srvir_zmass - srvir_remz)/(currenthalo_stars_zmass),psym = symcat(sym_outline(symbols[2])),symsize = symsize
oplot,xmass,(srvir_zmass - srvir_remz + ism_zmass)/(currenthalo_stars_zmass),psym = symcat(symbols[3]),color = colors[3],symsize = symsize
oplot,xmass,(srvir_zmass - srvir_remz + ism_zmass)/(currenthalo_stars_zmass),psym = symcat(sym_outline(symbols[3])),symsize = symsize
legend,["Within Rvir","Stars + Disk Gas","Stars","Disk Gas"],color = [colors[0],colors[3],colors[2],colors[1]],psym = [symbols[0],symbols[3],symbols[2],symbols[1]],box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop

;Fraction of metals produced by all stars in the virial radius at z = 0 are still within rvir
;Add in metals accreted
IF keyword_set(outplot) THEN  device,filename = outplot + '_zcurrfrac_mass_hist.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,alog10(xmass),(accr_zmass + s_zmass)/(currenthalo_stars_zmass),xrange = alog10(xrange),yrange = [0,1],/nodata,xtitle = textoidl('Log(M_*/M')+sunsymbol()+')',ytitle = textoidl('M_z/M_z avaliable')
binsize = 0.1
xcoords = [[alog10(xmass) - binsize/2],[alog10(xmass) - binsize/2],[alog10(xmass) + binsize/2],[alog10(xmass) + binsize/2]]
frac_inhalo = (accrh_halo_zmass + srvir_zmass - srvir_remz)/currenthalo_stars_zmass
frac_ism = ism_zmass/currenthalo_stars_zmass
frac_star = (srvir_zmass - srvir_remz)/currenthalo_stars_zmass
FOR i = 0, n_elements(xmass) - 1 DO BEGIN &  $
polyfill,xcoords[i,*],[0,frac_star[i],frac_star[i],0],color = colors[2],/line_fill,orientation = 45,spacing = 0.1 &  $
oplot,xcoords[i,*],[0,frac_star[i],frac_star[i],0],color = colors[2],thick = 5 &  $
polyfill,xcoords[i,*],[frac_star[i],frac_star[i]+frac_ism[i],frac_star[i]+frac_ism[i],frac_star[i]],color = colors[1],/line_fill,orientation = 60,spacing = 0.1 &  $
oplot,xcoords[i,*],[frac_star[i],frac_star[i]+frac_ism[i],frac_star[i]+frac_ism[i],frac_star[i]],color = colors[1],thick = 5 &  $
polyfill,xcoords[i,*],[frac_star[i]+frac_ism[i],frac_inhalo[i],frac_inhalo[i],frac_star[i]+frac_ism[i]],color = colors[0],/line_fill,orientation = 30,spacing = 0.1 &  $
oplot,xcoords[i,*],[frac_star[i]+frac_ism[i],frac_inhalo[i],frac_inhalo[i],frac_star[i]+frac_ism[i]],color = colors[0],thick = 5 &  $
ENDFOR
legend,["Within Rvir","Stars","Disk Gas"],color = [colors[0],colors[2],colors[1]],linestyl = [0,0,0],box = 0
IF keyword_set(redshift) THEN legend,['z = ' + strtrim(redshift,2)],box = 0,/right,/top
IF keyword_set(outplot) THEN device, /close ELSE stop
print,'Quantitative comparison to Dwarfs: '
ind = where((srvir_mass - srvir_remmass) LE 1.9e7)
print,(srvir_mass - srvir_remmass)[ind]
print,'CGM'
print,frac_inhalo[ind]
print,'Disk'
print,frac_star[ind]+frac_ism[ind]
print,'HI + H2'
print,HIH2_zmass[ind]/currenthalo_stars_zmass[ind]
print,'Star'
print,frac_star[ind]

;Fraction of metals produced by stars in current main halo, as
;calculated by peeples in that are
;still within rvir                                                                                                                           
IF keyword_set(outplot) THEN  device,filename = outplot + '_zcurrfrac_peeples2_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,(srvir_mass - srvir_remmass),(accr_zmass + s_zmass)/(zmetals_predict),/xlog,xrange = xrange,yrange = [0,1],/nodata,xtitle = xtitle,ytitle = textoidl('M_z/M_z avaliable')
oplot,(srvir_mass - srvir_remmass),(grvir_zmass + srvir_zmass - srvir_remz)/(zmetals_predict),psym = symcat(symbols[0]),color = colors[0],symsize = symsize
oplot,(srvir_mass - srvir_remmass),(grvir_zmass + srvir_zmass - srvir_remz)/(zmetals_predict),psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,(srvir_mass - srvir_remmass),ism_zmass/(zmetals_predict),psym = symcat(symbols[1]),color = colors[1],symsize = symsize
oplot,(srvir_mass - srvir_remmass),ism_zmass/(zmetals_predict),psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,(srvir_mass - srvir_remmass),(srvir_zmass - srvir_remz)/(zmetals_predict),psym = symcat(symbols[2]),color = colors[2],symsize = symsize
oplot,(srvir_mass - srvir_remmass),(srvir_zmass - srvir_remz)/(zmetals_predict),psym = symcat(sym_outline(symbols[2])),symsize = symsize
oplot,(srvir_mass - srvir_remmass),(srvir_zmass - srvir_remz + ism_zmass)/(zmetals_predict),psym = symcat(symbols[3]),color = colors[3],symsize = symsize
oplot,(srvir_mass - srvir_remmass),(srvir_zmass - srvir_remz + ism_zmass)/(zmetals_predict),psym = symcat(sym_outline(symbols[3])),symsize = symsize
legend,["Within Rvir","Stars + ISM","Stars","ISM"],color = [colors[0],colors[3],colors[2],colors[1]],psym = [symbols[0],symbols[3],symbols[2],symbols[1]],box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop
ind = where((srvir_mass - srvir_remmass) LE 1.9e7)
print,(srvir_mass - srvir_remmass)[ind]
print,'CGM'
print,frac_inhalo[ind]
print,'Disk'
print,frac_star[ind]+frac_ism[ind]
print,'HI + H2'
print,HIH2_zmass[ind]/currenthalo_stars_zmass[ind]
print,'Star'
print,frac_star[ind]

IF keyword_set(outplot) THEN  device,filename = outplot + '_zcurrfrac_peeples2_hist.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = xsize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = xsize
plot,alog10(xmass),(accr_zmass + s_zmass)/(currenthalo_stars_zmass),xrange = [9.2,11.4],yrange = [0,1],/nodata,xtitle = textoidl('Log(M_*/M')+sunsymbol()+')',ytitle = textoidl('M_z/M_z avaliable')
binsize = 0.1
xcoords = [[alog10(srvir_mass - srvir_remmass) - binsize/2],[alog10(srvir_mass - srvir_remmass) - binsize/2],[alog10(srvir_mass - srvir_remmass) + binsize/2],[alog10(srvir_mass - srvir_remmass) + binsize/2]]
frac_inhalo = (gr150_zmass + srvir_zmass - srvir_remz)/zmetals_predict ;(grvir_zmass + srvir_zmass - srvir_remz)
frac_ism = HIH2_zmass/zmetals_predict ;ism_zmass/zmetals_predict
frac_star = (srvir_zmass - srvir_remz)/zmetals_predict
FOR i = 0, n_elements(xmass) - 1 DO BEGIN &  $
IF alog10(xmass[i]) GT 9.5 THEN BEGIN & $
polyfill,xcoords[i,*],[0,frac_star[i],frac_star[i],0],color = colors[2],/line_fill,orientation = 45,spacing = 0.1 &  $
oplot,xcoords[i,*],[0,frac_star[i],frac_star[i],0],color = colors[2],thick = 5 &  $
polyfill,xcoords[i,*],[frac_star[i],frac_star[i]+frac_ism[i],frac_star[i]+frac_ism[i],frac_star[i]],color = colors[1],/line_fill,orientation = 60,spacing = 0.1 &  $
oplot,xcoords[i,*],[frac_star[i],frac_star[i]+frac_ism[i],frac_star[i]+frac_ism[i],frac_star[i]],color = colors[1],thick = 5 &  $
polyfill,xcoords[i,*],[frac_star[i]+frac_ism[i],frac_inhalo[i],frac_inhalo[i],frac_star[i]+frac_ism[i]],color = colors[0],/line_fill,orientation = 30,spacing = 0.1 &  $
oplot,xcoords[i,*],[frac_star[i]+frac_ism[i],frac_inhalo[i],frac_inhalo[i],frac_star[i]+frac_ism[i]],color = colors[0],thick = 5 &  $
  ENDIF & $
ENDFOR
legend,["Within 150 kpc","Stars",textoidl('HI + H_2')],color = [colors[0],colors[2],colors[1]],linestyl = [0,0,0],box = 0
IF keyword_set(redshift) THEN legend,['z = ' + strtrim(redshift,2)],box = 0,/right,/top
IF keyword_set(outplot) THEN device, /close ELSE stop
print,'Quantitative comparison to Peeples: '
ind = where((srvir_mass - srvir_remmass) GT 10^9.5)
print,(srvir_mass - srvir_remmass)[ind]
print,'CGM'
print,frac_inhalo[ind]
print,'Disk'
print,frac_star[ind]+frac_ism[ind]

IF redshift eq 0 THEN BEGIN
;Fraction of metals produced by stars in main progenitor in that are still within rvir
IF keyword_set(outplot) THEN  device,filename = outplot + '_zcurrfrac_mass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,xmass,(accr_zmass + s_zmass)/(primehalo_stars_zmass),/xlog,xrange = xrange,yrange = [0,1],/nodata,xtitle = xtitle,ytitle = textoidl('M_z/M_z avaliable')
oplot,xmass,(accr_halo_zmass + s_zmass)/(primehalo_stars_zmass),psym = symcat(symbols[0]),color = colors[0],symsize = symsize
oplot,xmass,(accr_halo_zmass + s_zmass)/(primehalo_stars_zmass),psym = symcat(sym_outline(symbols[0])),symsize = symsize
oplot,xmass,ism_zmass/(primehalo_stars_zmass),psym = symcat(symbols[1]),color = colors[1],symsize = symsize
oplot,xmass,ism_zmass/(primehalo_stars_zmass),psym = symcat(sym_outline(symbols[1])),symsize = symsize
oplot,xmass,s_zmass/(primehalo_stars_zmass),psym = symcat(symbols[2]),color = colors[2],symsize = symsize
oplot,xmass,s_zmass/(primehalo_stars_zmass),psym = symcat(sym_outline(symbols[2])),symsize = symsize
oplot,xmass,(s_zmass + ism_zmass)/(primehalo_stars_zmass),psym = symcat(symbols[3]),color = colors[3],symsize = symsize
oplot,xmass,(s_zmass + ism_zmass)/(primehalo_stars_zmass),psym = symcat(sym_outline(symbols[3])),symsize = symsize
legend,["Within Rvir","Stars + ISM","Stars","ISM"],color = [colors[0],colors[3],colors[2],colors[1]],psym = [symbols[0],symbols[3],symbols[2],symbols[1]],box = 0
IF keyword_set(outplot) THEN device, /close ELSE stop
ENDIF


;Fraction of metals produced in accreted material that are still within rvir
IF keyword_set(outplot) THEN  device,filename = outplot + '_zfrac_smass.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
plot,shalo_mass,(accr_zmass + s_zmass)/(accr_zmass + s_zmass),/xlog,xrange = xrange,yrange = [0,1],/nodata,xtitle = 'Stellar Mass [Msun]',ytitle = textoidl('M_z/M_z avaliable')
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
