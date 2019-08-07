;Compares the metal enrichment of gas that never returns vs gas that
;cycles back. Uses output from time_cycling.pro

readcol,'unreaccZ.grp'+haloid +'.txt',ejectmet_noRet,outflowmass_noRet,ejecttime_noRet
readcol,'reaccZdiff.grp' + haloid + '.txt',metdiff,ejectmet,outflowmass,ejecttime

plot,ejecttime_noRet,ejectmet_noRet,psym = 3,/ylog
oplot,ejecttime,ejectmet,psym = 3,color = 254

tbin = 0.5

timearr = findgen(14/tbin)*tbin
ejectmet_noRet_ave = fltarr(14/tbin)
ejectmet_ave = fltarr(14/tbin)
ejectmet_ave_mass = fltarr(14/tbin)

FOR i = 0, 14/tbin - 1 DO BEGIN & $
    ind_noRet = where(ejecttime_noRet GE i*tbin AND ejecttime_noRet LT (i + 1)*tbin) & $
    IF ind_noRet[0] NE -1 THEN ejectmet_noRet_ave[i] = total(outflowmass_noRet[ind_noRet]*ejectmet_noRet[ind_noRet])/total(outflowmass_noRet[ind_noRet]) & $
    ind = where(ejecttime GE i*tbin AND ejecttime LT (i + 1)*tbin) & $
    IF ind[0] NE -1 THEN ejectmet_ave[i] = total(outflowmass[ind]*ejectmet[ind])/total(outflowmass[ind]) & $
    IF ind[0] NE -1 THEN ejectmet_ave_mass[i] = total(outflowmass[ind]) & $
  END

oplot,timearr,ejectmet_noRet_ave
oplot,timearr,ejectmet_ave,color = 254

legend,['Return','Never Return'],linestyle = [0,0],color = [254,0]
print,total(ejectmet_noRet_ave*ejectmet_ave_mass)/total(ejectmet_ave*ejectmet_ave_mass)

;799, 1:2.5; 4: 1.5; 6: 1.3

;h516, 1: 2; 2: 1.1

;h603, 1: 1.67; 2: 1.2; 3: 1.5

;h258, 1: 1.12 ; 4: 1.25

;h285, 1: 1.3; 4: 1.5, 9: 1.35

;h239, 1: 1.19
