;find_gas_history_master,'h986.cosmo50cmb.3072g14HBWK.grp1.haloid.dat','1'

;run find_all_gas, haloidoutput ,finalid
PRO find_gas_history_master, filebase, finalid = finalid, a_dir_struc = a_dir_struc
IF NOT keyword_set(finalid) THEN finalid = '1'

readcol,  filebase + '.grp' + finalid + '.haloid.dat', stepfile, halo, format='a,l'
endpos = strpos(stepfile[0],'/')
FOR i = 0, n_elements(stepfile) - 1 DO stepfile[i] = strmid(stepfile[i],0,endpos)

orig_iord = mrdfits('grp' + finalid + '.allgas.iord.fits',0)
FOR i=0,n_elements(stepfile)-1 DO BEGIN
   cd,stepfile[i]
;   IF NOT exist(stepfile[i]+ '.grp' + strtrim(finalid,2) + '.allgas.history.fits') THEN BEGIN
      print,stepfile[i]
      IF keyword_set(a_dir_struc) THEN BEGIN
         find_history, stepfile[i], 'grp' + finalid + '.allgas', '../tipsy.units.idl', ORIG_IORD=orig_iord ;, minpart=64, /del
      ENDIF ELSE BEGIN
         find_history, stepfile[i], 'grp' + finalid + '.allgas', '../../tipsy.units.idl', ORIG_IORD=orig_iord ;, minpart=64, /del
      ENDELSE
;   ENDIF
   cd,'../..'
;   cd,'..'
ENDFOR

END

