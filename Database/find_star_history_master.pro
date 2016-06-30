;find_star_history_master,'/home/christensen/Storage2/UW/MolecH/Cosmo/h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK','h277.cosmo50cmb.3072g14HMbwK',1

;run find_all_gas, haloidoutput ,haloid
PRO find_star_history_master,dir,file,haloid,a_dir_struc = a_dir_struc
cd,dir
;command = 'ls ' + dir + "/steps/" + file + "*.dir/" + file + ".00???.amiga.grp | grep amiga.grp | sed 's/.amiga.grp//g'"
command = 'ls ' + dir + "/" + file + "*/" + file + ".00???.amiga.grp | grep amiga.grp | sed 's/.amiga.grp//g'"
spawn,command,files
file_temp = files[0]
cut_pos = strsplit(file_temp,'.')
files_cut = strmid(file_temp,0)
steps_st = strmid(files,cut_pos[N_ELEMENTS(cut_pos) - 1])
basefile = file + '.' + steps_st

last_step = steps_st[n_elements(steps_st) - 1]
rtipsy,file + '.' + last_step + '/' + file + '.' + last_step, h,g,d,s
d = 0
g = 0
sind = where(s.tform gt 0)
grp = read_lon_array(file + '.' + last_step + '/' + file + '.' + last_step + '.amiga.grp')
grp = grp[(h.ngas+h.ndark):(h.n-1)]
grp = grp[sind]
dwfstar = where(grp EQ  haloid  )
dwf = sind[dwfstar]+(h.ngas+h.ndark)
iord    = read_lon_array(file + '.' + last_step + '/' + file + '.' + last_step + '.iord')
igasord = read_lon_array(file + '.' + last_step + '/' + file + '.' + last_step + '.igasorder')
orig_iord = iord(dwf)
orig_igasord = igasord(dwf)

FOR i = 0, n_elements(basefile) - 1 DO BEGIN
   cd,basefile[i]
   IF keyword_set(a_dir_struc) THEN BEGIN
      find_history, basefile[i], 'grp' + strtrim(haloid,2) + '.star',    '../tipsy.units.idl', IND_STARS=dwfstar, IGASORD_STARS=orig_igasord, ORIG_IORD=orig_iord;, minpart=64 
      cd,'../'
   ENDIF ELSE BEGIN
      find_history, basefile[i], 'grp' + strtrim(haloid,2) + '.star', '../../tipsy.units.idl', IND_STARS=dwfstar, IGASORD_STARS=orig_igasord, ORIG_IORD=orig_iord;, minpart=64 
      cd,'../..'
   ENDELSE
ENDFOR

END
