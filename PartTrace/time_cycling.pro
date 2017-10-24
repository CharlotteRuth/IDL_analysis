;Finds the amount of time between cycling

PRO time_cycling,dirs,files,halo = halo,nowrite = nowrite,debug = debug
n = n_elements(dirs)
IF NOT keyword_set(halo) THEN halo = strarr(n) + '1'

formatplot,outplot = outplot,thick = formatthick
zsolar = 0.0130215
timecycled = [0]
print,'Enter 78 for the MPL purple-yellow colortable'
loadct, file='/home6/crchrist/IDL/IDL-Colorbars/mycolorbars.tbl'
openw,2,'~/metallicity_diff_inflow_outflow.txt'
FOR i = 0, n - 1 DO BEGIN
;   IF file_test(dirs[i] + '/outflowtime.grp' + halo[i] + '.txt') THEN BEGIN
      diskaccri = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_iord.fits',0) 
      diskaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdisk.fits',0)
      diskaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0)
      diskaccrt = z_to_t(diskaccrz)
      diskaccr_data = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_history.fits',1,/silent)

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
      reeject_data = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_halo.fits',1,/silent)

      uniqi = diskaccri[uniq(diskaccri)]
;    stop
      
      outflowtime = [0]
      disktime = [0]
      outflowmass = [0]
      diskmass = [0]
      reaccrtime = [0]
      ejecttime = [0]
      reaccriord = [0]
      metdiff = [0]
      reaccrmet = [0]
      ejectmet = [0]
      FOR j = 0, n_elements(uniqi) - 1 DO BEGIN
         diskaccrind = where(diskaccri EQ uniqi[j])
         reejectind  = where(reejecti EQ uniqi[j])
         IF (n_elements(diskaccrind) - 1 NE n_elements(reejectind)) AND (n_elements(diskaccrind) NE n_elements(reejectind)) THEN stop
         IF n_elements(diskaccrind) - 1 NE n_elements([where(diskaccr_data.iord EQ uniqi[j])]) THEN stop
         diskaccrtimes = diskaccrt[diskaccrind[sort(diskaccrt[diskaccrind])]]
         reejecttimes = reejectt[reejectind[sort(reejectt[reejectind])]]
         diskaccrmets = diskaccr_data[diskaccrind[sort(diskaccrt[diskaccrind])]].metallicity
         reejectmets = reeject_data[reejectind[sort(reejectt[reejectind])]].metallicity
         disktime = [disktime,reejecttimes - diskaccrtimes[0:n_elements(reejectind) - 1]]
         outflowtime = [outflowtime,diskaccrtimes[1:n_elements(diskaccrind) - 1 ] - reejecttimes[0:n_elements(diskaccrind) - 2]]
;         stop
;         print,diskaccr_data[where(diskaccr_data.iord EQ uniqi[j])]
         metdiff     = [metdiff, diskaccr_data[where(diskaccr_data.iord EQ uniqi[j])].metallicity - reejectmets[0:n_elements(diskaccrind) - 2]]
         reaccrmet = [reaccrmet,diskaccr_data[where(diskaccr_data.iord EQ uniqi[j])].metallicity]
         ejectmet = [ejectmet,reejectmets[0:n_elements(diskaccrind) - 2]]
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
      metdiff = metdiff[1:n_elements(metdiff) - 1]
      reaccrmet = reaccrmet[1:n_elements(reaccrmet) - 1]
      ejectmet = ejectmet[1:n_elements(ejectmet) - 1]
      reaccriord = reaccriord[1:n_elements(reaccriord) - 1]
      IF keyword_set(debug) THEN BEGIN
          contour_plus,outflowtime,alog10(abs(metdiff/ejectmet)),nlevels = 254,ymin = -6.0,ymax = 2.0,xmin = 0,xmax = 14.0,threshold = 1e7,xbinsize=0.2,ybinsize = 0.2,/fill,xtitle = 'Time in outflow',ytitle = textoidl('Log(|(Z_{reaccr} - Z_{eject})/Z_{eject}|)')
          stop
          contour_plus,outflowtime,metdiff/ejectmet,nlevels = 254,ymin = -2,ymax = 10,threshold = 1e7,xbinsize=0.2,ybinsize = 0.2,/fill,xtitle = 'Time in outflow',ytitle = textoidl('(Z_{reaccr} - Z_{eject})/Z_{eject}')
          stop
;      plot,outflowtime,metdiff/ejectmet,psym = 3,/ylog,yrange = [1e-5,100],ytitle = 'Fractional change in metallicity',xtitle = 'Time in outflow'
;      plot,ejectmet,reaccrmet,psym = 3,/ylog,/xlog,xrange = [1e-6,1e1],yrange = [1e-6,1e1],xtitle = 'Ejected Metallicity',ytitle = 'Reaccreted Metallicity'
          histogramp,metdiff,xtitle = textoidl('Z_{reaccr} - Z_{eject}'),weight = outflowmass
          stop
          histogramp,alog10(abs(metdiff/ejectmet)),xtitle = textoidl('Log(|(Z_{reaccr} - Z_{eject})/Z_{eject}|)'),weight = outflowmass
          stop
          hist_eject = histogram(alog10(reeject_data.metallicity/zsolar),min = -5,max = 1,locations = locations,nbins = 200)
          hist_eject_reaccr = histogram(alog10(ejectmet/zsolar),min = -5,max = 1,locations = locations,nbins = 200)
          plot,locations,hist_eject_reaccr/float(hist_eject),xtitle = 'Metallicity [Z'+sunsymbol()+']'
          oplot,locations,hist_eject/float(n_elements(ejectmet)),color = 100
          oplot,locations,hist_eject_reaccr/float(n_elements(ejectmet)),color = 200
          print,i,': ',files[i] + '.' + halo[i]
          print,'Ave(metdiff): ',strtrim(total(outflowmass*metdiff)/total(outflowmass),2),', STDEV(metdiff): ',strtrim(stdev(metdiff),2)
          print,'Ave(metdiff/metallicity): ',strtrim(total(outflowmass*metdiff/ejectmet)/total(outflowmass),2),', STDEV(metdiff/metallicity): ',strtrim(stdev(metdiff/ejectmet),2)
          print,'Total metals in ejecta [Msol]: ',strtrim(total(reejectm*reeject_data.metallicity),2)
          print,'Total metals in reaccreted ejecta [Msol]: ',strtrim(total(ejectmet*outflowmass),2)
          print,'Total metals reaccreted [Msol]: ',strtrim(total(reaccrmet*outflowmass),2)
          print,'Total metal mass loss through diffusion [Msol]: ',strtrim(total(outflowmass*metdiff),2),', fraction of reaccreted metals: ',strtrim(total(outflowmass*metdiff)/total(reaccrmet*outflowmass),2)
          print,'Average metallicity of ejected gas: [Zsol]: ',strtrim(mean(reeject_data.metallicity)/zsolar,2)
          print,'Average metallicity of reaccreted ejected gas: [Zsol]: ',strtrim(mean(ejectmet)/zsolar,2)
          print,'Relative difference: ',mean(ejectmet)/mean(reeject_data.metallicity)
          stop
      ENDIF

      printf,2,i,': ',files[i] + '.' + halo[i]
      printf,2,'Ave(metdiff): ',strtrim(total(outflowmass*metdiff)/total(outflowmass),2),', STDEV(metdiff): ',strtrim(stdev(metdiff),2)
      printf,2,'Ave(metdiff/metallicity): ',strtrim(total(outflowmass*metdiff/ejectmet)/total(outflowmass),2),', STDEV(metdiff/metallicity): ',strtrim(stdev(metdiff/ejectmet),2)
      printf,2,'Total metals in ejecta [Msol]: ',strtrim(total(reejectm*reeject_data.metallicity),2)
      printf,2,'Total metals in reaccreted ejecta [Msol]: ',strtrim(total(ejectmet*outflowmass),2)
      printf,2,'Total metals reaccreted [Msol]: ',strtrim(total(reaccrmet*outflowmass),2)
      printf,2,'Total metal mass loss through diffusion [Msol]: ',strtrim(total(outflowmass*metdiff),2),', fraction of reaccreted metals: ',strtrim(total(outflowmass*metdiff)/total(reaccrmet*outflowmass),2)
      printf,2,'Average metallicity of ejected gas: [Zsol]: ',strtrim(mean(reeject_data.metallicity)/zsolar,2)
      printf,2,'Average metallicity of reaccreted ejected gas: [Zsol]: ',strtrim(mean(ejectmet)/zsolar,2)
      printf,2,'Relative difference: ',mean(ejectmet)/mean(reeject_data.metallicity)

      openw,1,dirs[i] + '/reaccZdiff.grp' + halo[i] + '.txt'
      printf,1,transpose([[metdiff],[ejectmet],[outflowmass]])
      close,1
      IF NOT keyword_set(nowrite) THEN BEGIN
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
      ENDIF
;   ENDIF
ENDFOR
close,2   
END
