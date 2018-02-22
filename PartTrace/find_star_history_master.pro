;find_star_history_master,'/home/christensen/Storage2/UW/MolecH/Cosmo/h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK','h277.cosmo50cmb.3072g14HMbwK',1

;run find_all_gas, haloidoutput ,finalid
PRO find_star_history_master, dir, filebase, finalid = finalid, a_dir_struc = a_dir_struc,rewrite = rewrite
IF NOT keyword_set(finalid) THEN finalid = ['1']

;command = 'ls ' + dir + "/steps/" + filebase + "*.dir/" + filebase + ".00???.amiga.grp | grep amiga.grp | sed 's/.amiga.grp//g'"
command = 'ls ' + dir + "/" + filebase + "*/" + filebase + ".00*.amiga.grp | grep amiga.grp | sed 's/.amiga.grp//g'"
spawn,command,files
file_temp = files[0]
cut_pos = strsplit(file_temp,'.')
files_cut = strmid(file_temp,0)
steps_st = strmid(files,cut_pos[N_ELEMENTS(cut_pos) - 1])
basefile = filebase + '.' + steps_st

last_step = steps_st[n_elements(steps_st) - 1]
rtipsy,filebase + '.' + last_step + '/' + filebase + '.' + last_step, h,g,d,s
d = 0
g = 0
;sind = where(s.tform gt 0) ;Keep black holes in
sind = where(s.tform)
grp = read_lon_array(filebase + '.' + last_step + '/' + filebase + '.' + last_step + '.amiga.grp')
grp = grp[(h.ngas+h.ndark):(h.n-1)]
grp = grp[sind]
iord    = read_lon_array(filebase + '.' + last_step + '/' + filebase + '.' + last_step + '.iord')
igasord = read_lon_array(filebase + '.' + last_step + '/' + filebase + '.' + last_step + '.igasorder')

dwfstar = [-1]
grpstar = [-1]
FOR i_finalid = 0, n_elements(finalid) - 1 DO BEGIN
    dwfstar_temp = where(grp EQ  fix(finalid[i_finalid]  ))
    dwfstar = [dwfstar,dwfstar_temp]
    grpstar = [grpstar,fltarr(n_elements(dwfstar_temp)) + fix(finalid[i_finalid])]
ENDFOR
dwfstar = dwfstar[1:n_elements(dwfstar) - 1]
grpstar = grpstar[1:n_elements(grpstar) - 1]

dwf = sind[dwfstar]+(h.ngas+h.ndark)
orig_iord = iord(dwf)
orig_igasord = igasord(dwf)
FOR i = 0, n_elements(basefile) - 1 DO BEGIN
    cd,basefile[i]
    print,basefile[i]
    IF NOT exist(basefile[i]+ '.grp' + finalid[0] + '.star.history.fits') OR keyword_set(rewrite) THEN BEGIN
        IF keyword_set(a_dir_struc) THEN BEGIN
            IF NOT file_test('../tipsy.units.idl') THEN BEGIN
                units = tipsyunits('../'+filebase+'.param')
                writecol, '../tipsy.units.idl',units.lengthunit, units.massunit, units.vunit, format='(E15.7,E15.7,f)'
            ENDIF
            find_history, basefile[i], 'grp' + finalid + '.star',    '../tipsy.units.idl', IND_STARS=dwfstar, IGASORD_STARS=orig_igasord, ORIG_IORD=orig_iord, GRP_ARR = grpstar;, minpart=64 
            cd,'../'
        ENDIF ELSE BEGIN
            IF NOT file_test('../../tipsy.units.idl') THEN BEGIN
                units = tipsyunits('../../'+filebase+'.param')
                writecol, '../../tipsy.units.idl', units.lengthunit, units.massunit, units.vunit, format='(E15.7,E15.7,f)'
            ENDIF
            find_history, basefile[i], 'grp' + finalid + '.star', '../../tipsy.units.idl', IND_STARS=dwfstar, IGASORD_STARS=orig_igasord, ORIG_IORD=orig_iord, GRP_ARR = grpstar;, minpart=64 
        ENDELSE
    ENDIF
    cd,dir
ENDFOR

END
