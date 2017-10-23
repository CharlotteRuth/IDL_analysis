PRO z_at_ejection,dir,file,haloid
   spawn,'ls ' + dir + '*512/*cosmo**512',filename
   rtipsy,filename,h,g,d,s
   spawn,'ls ' + dir + '*512/*512.iord',file_iord
   readarr,file_iord,  h,iord,/ascii,type = 'long'
   iord_gas = iord[0:h.ngas - 1]
   reejecti = mrdfits(dir + '/grp' + haloid + '.reeject_iord.fits',0,/silent)
;   reejectm = mrdfits(dir + '/grp' + haloid + '.mass_at_reeject.fits',0,/silent)
;   reejectz = mrdfits(dir + '/grp' + haloid + '.reeject_z.fits',0,/silent)
   reeject_data = mrdfits(dir + '/grp' + haloid + '.reeject_halo.fits',1,/silent)
   uniq_eject = uniq(reejecti,sort(reejecti))
   reejecti = reejecti[uniq_eject]
   reeject_data = reeject_data[uniq_eject]
   match2,iord_gas,reejecti,ind1,ind2
   metals = g.zmetal
   plot,metals[ind2(where(ind2 ne -1))],reeject_data[where(ind2 ne -1)].metallicity,psym = 3,ytitle = 'Ejected Metallicity',xtitle = 'z = 0 Metallicity'
   metals[ind2] = reeject_data.metallicity
   metals_arr = fltarr(n_elements(iord_gas))
   metals_arr[0:h.ngas - 1] = metals
   split = strsplit(file_iord,'.')
   writecol,strmid(file_iord,0,split[n_elements(split) - 1])+haloid+'.zmetals',metals_arr
END
