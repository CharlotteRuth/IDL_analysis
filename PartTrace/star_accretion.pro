
; star_accretion,'h799.cosmo25cmb.3072g14HBWK','1',/molecularH
; star_accretion,'h986.cosmo50cmb.3072g14HBWK','1',/molecularH
PRO star_accretion, filename, haloid, molecularH = molecularH
units = tipsyunits(filename + '.param')
stars_struct = {haloid:0L,iord:0L,mass:0.0}
IF file_test(filename + '_2merge.starlog.fits') THEN BEGIN
   print,filename + '_2merge.starlog.fits'
   sl = mrdfits(filename + '_2merge.starlog.fits',1) 
ENDIF ELSE BEGIN
   IF file_test(filename + '.starlog.fits') THEN BEGIN
      print,filename + '.starlog.fits'
      sl = mrdfits(filename + '.starlog.fits',1)        
   ENDIF ELSE BEGIN
      print,'Starlog: ',filename + '.starlog'
      sl = rstarlog(filename + '.starlog',molecularH = molecularH)
   ENDELSE
ENDELSE
sl = sl[uniq(sl.iorderstar,sort(sl.iorderstar))]

readcol,filename + '.grp' + haloid + '.haloid.dat',files,halo,format='a,l'
spawn,'ls ' + filename + '*512/*cosmo**512',file
rtipsy,file,h,g,d,s,/justhead
spawn,'ls ' + filename + '*512/*512.iord',file_iord
read_tipsy_arr,file_iord,h,iord,type = 'long'
iord_star = iord[h.ngas + h.ndark:h.n - 1]
spawn,'ls ' + filename + '*512/*512.amiga.grp',file_grp
read_tipsy_arr,file_grp,h,grp,type = 'long'
grp_star = grp[h.ngas + h.ndark:h.n - 1]
spawn,'ls ' + filename + '*512/*512.amiga.stat',file_stat
stat = read_stat_struc_amiga(file_stat)
main = where(stat.group EQ haloid)
satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir AND stat.ngas gt 0)
inhalo = fltarr(n_elements(grp_star))
FOR j = 0, n_elements(satellites) - 1 DO BEGIN
   test = where(grp_star EQ stat[satellites[j]].group, ntest)
   IF ntest NE 0 THEN inhalo[test] = 1
ENDFOR
iord_star = iord_star[where(inhalo EQ 1)]

;halostars = mrdfits(files[0] + '.grp' + haloid + '.star.history.fits',1)
;match,sl.iorderstar,halostars.iord,keep,temp
match,sl.iorderstar,iord_star,keep,temp
sl = sl[keep] ;Only keep the stars that are ever in the halo

files = reverse(files)
halo = reverse(halo)
nsteps = n_elements(files)
halodat = mrdfits('grp' + haloid + '.alignment.fits',1)

accrz = fltarr(n_elements(sl)) + 99
accrm = fltarr(n_elements(sl))

;timestart = halodat[0].time
;timestart = 0
iordstart = 0
FOR i = 0, nsteps - 1 DO BEGIN
   rtipsy,files[i],h,g,d,s
   g = 0
   d = 0
   read_tipsy_arr,files[i]+'.amiga.grp',h,grp_star,part='star',type='long'
   read_tipsy_arr,files[i]+'.iord',h,iord_star,part='star',type='long'
   halostars = replicate(stars_struct,n_elements(grp_star))
   halostars.iord = iord_star
   halostars.haloid = grp_star
   halostars.mass = s.mass
   match,sl.iorderstar,halostars.iord,ind_sl,ind_arr ;sl[ind_sl] = halostars[formed[ind_all]]
   halostars = halostars[ind_arr] ;Select only halostars that are in the final halo at z = 0
   print,'Step: ',i,' (out of ',nsteps - 1,'); Time: ',halodat[i].time,' [Gyr]'
;   IF (where(iord_star EQ 13597670))[0] NE -1 THEN stop

   stat = read_stat_struc_amiga(files[i] + '.amiga.stat')
   main = where(halodat[i].haloid EQ stat.group)
   satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
   inhalotemp = intarr(n_elements(halostars))
   FOR j = 0, n_elements(satellites) - 1 DO BEGIN
      test = where(halostars.haloid EQ stat[satellites[j]].group, ntest)
      IF ntest NE 0 THEN inhalotemp[test] = 1
   ENDFOR
   inhalo = where(inhalotemp EQ 1, complement = outhalo)
  
   formed_outhalo = where(halostars.iord GT iordstart AND inhalotemp NE 1)
   print,'Formed: ',strtrim(n_elements(where(halostars.iord GT iordstart)),2),'; OutHalo: ',strtrim(n_elements(formed_outhalo),2),'; InHalo: ',strtrim(n_elements(where(halostars.iord GT iordstart AND inhalotemp EQ 1)),2)
;   stop
   IF formed_outhalo[0] NE -1 THEN $
      IF keyword_set(accr_iords) THEN accr_iords = [accr_iords, halostars[formed_outhalo].iord] ELSE accr_iords = [halostars[formed_outhalo].iord] 
   print,'Accr_Iords: ',strtrim(n_elements(accr_iords),2)
 
;-------------------------------------------------------------
   match,accr_iords,halostars.iord,accr_ind,halostars_ind ;Match the iords of the stars that formed outside the halo with those of the particle tracing iords
;   match,accr_iords,halostars.iord,accr_ind,halostars_ind 
   accr_inhalo = where(inhalotemp[halostars_ind] EQ 1, complement = notaccr_inhalo)
   IF accr_inhalo[0] NE -1 THEN BEGIN
      match,sl.iorderstar,halostars[halostars_ind[accr_inhalo]].iord,ind0,ind1
      accrz[ind0] = halodat[i].z
      accrm[ind0] = halostars[halostars_ind[accr_inhalo[ind1]]].mass*units.massunit
      IF (where(halostars[accr_inhalo[ind1]].mass*units.massunit EQ 0))[0] NE -1 THEN stop
   ENDIF
   print,'Accr_InHalo: ',strtrim(n_elements(accr_inhalo),2),'; NotAccr_InHalo: ',strtrim(n_elements(notaccr_inhalo),2)
;   stop

   accr_iords = accr_iords[accr_ind[notaccr_inhalo]]
   iordstart = max(iord_star)
ENDFOR
histogramp,accrz,weight = accrm,max = 5,nbins = 20
mwrfits,accrz,'grp' + haloid + '.accrstars_z.fits',/create ;Time when accreted onto the disk
mwrfits,accrm,'grp' + haloid + '.accrstars_m.fits',/create ;Time when accreted onto the disk
mwrfits,sl.iorderstar,'grp' + haloid + '.accrstars_i.fits',/create ;Time when accreted onto the disk
stop
END
