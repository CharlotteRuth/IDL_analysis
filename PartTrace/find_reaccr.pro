PRO find_reaccr,dirs,files,halo = halo,finalstep = finalstep
n = n_elements(dirs)
IF NOT keyword_set(finalstep) THEN finalstep = '00512'

CASE 1 OF
   keyword_set(heat): ejectext = 'reheat_all'
   keyword_set(rvircut_half): ejectext = 'reeject_rvir0.5'
   keyword_set(rvircut_fifth): ejectext = 'reeject_rvir0.2'
   ELSE: ejectext = 'reeject'
ENDCASE

FOR i = 0, n_elements(dirs) - 1 DO BEGIN
   diskaccri = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_iord.fits',0) 
   diskaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdisk.fits',0)
   diskaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0)
   diskaccrt = z_to_t(diskaccrz)

   timesaccr = histogram(diskaccri, reverse_indices = r)
   accrone = where(timesaccr EQ 1)
   FOR ih = 0L, n_elements(accrone) - 1 DO IF R[accrone[ih]] NE R[accrone[ih] + 1] THEN diskaccri[R[R[accrone[ih]] : R[accrone[ih] + 1]-1]] = -1 ;Set to -1 all particles that are accreted only once
   reaccr = where(diskaccri NE -1)
   diskaccri = diskaccri[reaccr]
   diskaccrm = diskaccrm[reaccr]
   diskaccrz = diskaccrz[reaccr]
   diskaccrt = diskaccrt[reaccr]
;Sort arrays by iord values
   sorted = sort(diskaccri)
   diskaccri = diskaccri[sorted]
   diskaccrm = diskaccrm[sorted]
   diskaccrz = diskaccrz[sorted]
   diskaccrt = diskaccrt[sorted]
   uniqi = diskaccri[uniq(diskaccri)]

   reejecti  = mrdfits(dir + '/grp' + finalid + '.' + ejectext + '_iord.fits',0)
   reejectm  = mrdfits(dir + '/grp' + finalid + '.mass_at_' + ejectext + '.fits',0)
   reejectz = mrdfits(dir + '/grp' + finalid + '.' + ejectext + '_z.fits',0)
   reejectt = z_to_t(reejectz)

   ind_history_eject_r = reverse(indgen(n_elements(history_eject)))
   history_eject_r = history_eject[ind_history_eject_r]
   ind_history_accr_r = reverse(indgen(n_elements(history_accr)))
   history_accr_r = history_accr(ind_history_accr_r)

   history_eject2 = history_accr
   accriord = history_accr[uniq(history_accr.iord,sort(history_accr.iord))].iord
   FOR j = 0, n_elements(accriord) - 1 DO BEGIN
      indaccr = where(history_accr.iord EQ accriord[j])
      indeject = where(history_eject.iord EQ accriord[j])
      history_eject2[indaccr] = history_eject[indeject[0:n_elements(indaccr) - 1]]
   ENDFOR
ENDFOR

END
