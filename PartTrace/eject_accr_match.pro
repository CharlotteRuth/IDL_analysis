FUNCTION eject_accr_match,dir,halo

   history_eject =  mrdfits(dir + 'grp' + halo + '.reeject_disk.fits',1)
   reejectz = mrdfits(dir + '/grp' + halo + '.reeject_z.fits',0)
   reejectt = z_to_t(reejectz)

   sort_ejectaccr = sort(history_eject.iord)
   history_eject = history_eject[sort_ejectaccr]
   reejectt = reejectt[sort_ejectaccr]
   reejectz = reejectz[sort_ejectaccr]

   diskaccri = mrdfits(dir + '/grp' + halo + '.reaccrdisk_iord.fits',0)
   diskaccrz = mrdfits(dir + '/grp' + halo + '.reaccrdisk_z.fits',0)
   sort_diskaccr = sort(diskaccri)
   diskaccri = diskaccri[sort_diskaccr]
   diskaccrz = diskaccrz[sort_diskaccr]

   ind_uniq = uniq(diskaccri,/first) ;Delete first instance of every particle
   diskaccrz[ind_uniq] = -1
   diskaccri = diskaccri[where(diskaccrz GT -1)]
   diskaccrz = diskaccrz[where(diskaccrz GT -1)]
   diskaccrt = z_to_t(diskaccrz)

   history_accr = mrdfits(dir + 'grp' + halo + '.reaccrdisk_history.fits',1)
   sort_diskaccr = sort(history_accr.iord)
   history_accr = history_accr[sort_diskaccr]

 ;  uniq_history_eject = uniq(history_eject.iord,/first)

   match2,history_eject.iord,history_accr.iord,ind0,ind1
   history_eject = history_eject[where(ind0 NE -1)]
   reejectz = reejectz[where(ind0 NE -1)]
   reejectt = reejectt[where(ind0 NE -1)]

   histeject = histogram(history_eject.iord,locations = histiord)
   histaccr = histogram(history_accr.iord)

   iord_endedejected = histiord[where(histeject - histaccr GT 0)]

   uniq_history_eject = uniq(history_eject.iord)
   uniq_iord_eject = history_eject[uniq_history_eject].iord
   match,uniq_iord_eject,iord_endedejected,ind0,ind1

   reejectt[uniq_history_eject[ind0]] = -1
   history_eject = history_eject[where(reejectt GT -1)]


   RETURN,history_eject
END
