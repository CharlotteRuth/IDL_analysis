;Finds the amount of time between cycling

PRO time_cycling,dirs,files,halo = halo
n = n_elements(dirs)
IF NOT keyword_set(halo) THEN halo = strarr(n) + '1'

formatplot,outplot = outplot,thick = formatthick

timecycled = [0]
FOR i = 0, n - 1 DO BEGIN
;   IF file_test(dirs[i] + '/outflowtime.grp' + halo[i] + '.txt') THEN BEGIN
      diskaccri = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_iord.fits',0) 
      diskaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdisk.fits',0)
      diskaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0)
      diskaccrt = z_to_t(diskaccrz)
      
      timesaccr = histogram(diskaccri, reverse_indices = r)
;    IF R[0] NE R[1] THEN diskaccri[R[R[0] : R[1]-1]] = -1 ;Set to zero all 
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
      
      reejecti = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_iord.fits',0)
      reejectm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reeject.fits',0)
      reejectz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_z.fits',0)
      reejectt = z_to_t(reejectz)
      
      uniqi = diskaccri[uniq(diskaccri)]
;    stop
      
      outflowtime = [0]
      disktime = [0]
      outflowmass = [0]
      diskmass = [0]
      reaccrtime = [0]
      ejecttime = [0]
      reaccriord = [0]
      FOR j = 0, n_elements(uniqi) - 1 DO BEGIN
         diskaccrind = where(diskaccri EQ uniqi[j])
         reejectind  = where(reejecti EQ uniqi[j])
         diskaccrtimes = diskaccrt[diskaccrind[sort(diskaccrt[diskaccrind])]]
         reejecttimes = reejectt[reejectind[sort(reejectt[reejectind])]]
         disktime = [disktime,reejecttimes - diskaccrtimes[0:n_elements(reejectind) - 1]]
         outflowtime = [outflowtime,diskaccrtimes[1:n_elements(diskaccrind) - 1 ] - reejecttimes[0:n_elements(diskaccrind) - 2]]
         IF (where(diskaccrtimes[1:n_elements(diskaccrind) - 1 ] - reejecttimes[0:n_elements(diskaccrind) - 2] LT 0.001))[0] NE -1 THEN stop
         reaccrtime = [reaccrtime,diskaccrtimes[1:n_elements(diskaccrind) - 1 ]]
         reaccriord = [reaccriord,lonarr(n_elements(diskaccrind) - 1) + uniqi[j]]
         ejecttime = [ejecttime,reejecttimes[0:n_elements(diskaccrind) - 2]]
         diskmass = [diskmass,reejectm[reejectind[sort(reejectt[reejectind])]]]
         outflowmass = [outflowmass,(diskaccrm[diskaccrind[sort(diskaccrt[diskaccrind])]])[1:n_elements(diskaccrind) - 1 ]]
      ENDFOR
      outflowtime = outflowtime[1:n_elements(outflowtime) - 1]
      disktime = disktime[1:n_elements(disktime) - 1]
      reaccrtime = reaccrtime[1:n_elements(reaccrtime) - 1]
      ejecttime = ejecttime[1:n_elements(ejecttime) - 1]
      outflowmass = outflowmass[1:n_elements(outflowmass) - 1]
      diskmass = diskmass[1:n_elements(diskmass) - 1]
       reaccriord = reaccriord[1:n_elements(reaccriord) - 1]
      writecol,dirs[i] + '/reaccriord.grp' + halo[i] + '.txt',long(reaccriord)
      openw,1,dirs[i] + '/outflowtime.grp' + halo[i] + '.txt'
      printf,1,transpose([[outflowtime],[outflowmass]])
      close,1
      openw,1,dirs[i] + '/outflowtime.grp' + halo[i] + '.txt'
      printf,1,transpose([[outflowtime],[outflowmass]])
      close,1
      openw,1,dirs[i] + '/disktime.grp' + halo[i] + '.txt'
      printf,1,transpose([[disktime],[diskmass]])
      close,1
      openw,1,dirs[i] + '/inouttime.grp' + halo[i] + '.txt'
      printf,1,transpose([[ejecttime],[reaccrtime]])
      close,1
;    stop
;   ENDIF
ENDFOR
   
END
