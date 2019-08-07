;find_gas_history_master,'h986.cosmo50cmb.3072g14HBWK.grp1.haloid.dat','1'

;run find_all_gas, haloidoutput ,finalid
PRO find_gas_history_master, filebase, finalid = finalid, a_dir_struc = a_dir_struc,rewrite = rewrite
IF NOT keyword_set(finalid) THEN finalid = ['1']

;Note that haloid files can vary in length based on how long the halo
;was detected for. You want the longest.
readcol,  filebase + '.grp' + finalid[0] + '.haloid.dat', stepfile, halo, format='a,l'
print,filebase + '.grp' + finalid[0] + '.haloid.dat'
endpos = strpos(stepfile[0],'/')
FOR i = 0, n_elements(stepfile) - 1 DO stepfile[i] = strmid(stepfile[i],0,endpos)

orig_iord_g = [-1]
grpgas = [-1]
FOR i_finalid = 0, n_elements(finalid) - 1 DO BEGIN
    orig_iord_temp = mrdfits('grp' + finalid[i_finalid] + '.allgas.iord.fits',0)
    orig_iord_g = [orig_iord_g,orig_iord_temp]
    grpgas = [grpgas,fltarr(n_elements(orig_iord_temp)) + fix(finalid[i_finalid])]
ENDFOR
orig_iord_g = orig_iord_g[1:n_elements(orig_iord_g) - 1]
grpgas = grpgas[1:n_elements(grpgas) - 1]

FOR i=0,n_elements(stepfile)-1 DO BEGIN
   cd,stepfile[i]
   IF NOT exist(stepfile[i]+ '.grp' + strtrim(finalid[0],2) + '.allgas.history.fits') OR keyword_set(rewrite) THEN BEGIN
      print,stepfile[i]
      IF keyword_set(a_dir_struc) THEN BEGIN
         find_history, stepfile[i], 'grp' + finalid + '.allgas', '../tipsy.units.idl', ORIG_IORD=orig_iord_g, GRP_ARR = grpgas  ;, minpart=64, /del
      ENDIF ELSE BEGIN
         find_history, stepfile[i], 'grp' + finalid + '.allgas', '../../tipsy.units.idl', ORIG_IORD=orig_iord_g, GRP_Arr = grpgas ;, minpart=64, /del
      ENDELSE
   ENDIF
   cd,'../..'
;   cd,'..'
ENDFOR

END

