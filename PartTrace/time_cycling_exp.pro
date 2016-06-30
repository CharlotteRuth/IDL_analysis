;Finds the amount of time between cycling

PRO time_cycling_exp,dirs,files,halo = halo
n = n_elements(dirs)
IF NOT keyword_set(halo) THEN halo = strarr(n) + '1'

formatplot,outplot = outplot,thick = formatthick

timecycled = [0]
FOR i = 0, n - 1 DO BEGIN
    reejecti = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_iord.fits',0)
    reejectm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reeject.fits',0)
    reejectz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_z.fits',0)
    reexpelli = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_iord.fits',0)
    reexpellm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reexpell.fits',0)
    reexpellz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reexpell_z.fits',0)
    reexpellm1 = fltarr(n_elements(reexpellm))
    reexpellz1 = fltarr(n_elements(reexpellz))
    FOR exi = 0, n_elements(reexpellz1) - 1 DO BEGIN
        indmatch = where(reejecti EQ reexpelli[exi])
        indmatch2 = where(reejectz[indmatch] GE reexpellz[exi])
        indmatch3 = indmatch2[n_elements(indmatch2) - 1]
        reexpellm1[exi] = reejectm[indmatch[indmatch3]]
        reexpellz1[exi] = reejectz[indmatch[indmatch3]]
    ENDFOR
    reexpellm = reexpellm1
    reexpellt = z_to_t(reexpellz1)
;    stop

    diskaccri = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_iord.fits',0) 
    diskaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdisk.fits',0)
    diskaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0)
    diskaccrt = z_to_t(diskaccrz)

    match2,reexpelli,diskaccri,rexpind,diskind
    keep = where(diskind NE -1)
    diskaccri = diskaccri[keep]
    diskaccrm = diskaccrm[keep]
    diskaccrz = diskaccrz[keep]
    diskaccrt = diskaccrt[keep]

    timesaccr = histogram(diskaccri, reverse_indices = r)
;    IF R[0] NE R[1] THEN diskaccri[R[R[0] : R[1]-1]] = -1 ;Set to zero all 
    accrone = where(timesaccr EQ 1)
    FOR ih = 0L, n_elements(accrone) - 1 DO IF R[accrone[ih]] NE R[accrone[ih] + 1] THEN diskaccri[R[R[accrone[ih]] : R[accrone[ih] + 1]-1]] = -1 ;Set to -1 all particles that are accreted only once
    reaccr = where(diskaccri NE -1)
;    stop
    diskaccri = diskaccri[reaccr]
    diskaccrm = diskaccrm[reaccr]
    diskaccrz = diskaccrz[reaccr]
    diskaccrt = diskaccrt[reaccr]

    sorted = sort(diskaccri)
    diskaccri = diskaccri[sorted]
    diskaccrm = diskaccrm[sorted]
    diskaccrz = diskaccrz[sorted]
    diskaccrt = diskaccrt[sorted]

    uniqi = diskaccri[uniq(diskaccri)]

    outflowtime = [0]
;    disktime = [0]
    outflowmass = [0]
;    diskmass = [0]
    FOR j = 0, n_elements(uniqi) - 1 DO BEGIN
        diskaccrind = where(diskaccri EQ uniqi[j])
        reexpellind  = where(reexpelli EQ uniqi[j])
        IF reexpellind[0] EQ -1 THEN stop
        diskaccrtimes = diskaccrt[diskaccrind[sort(diskaccrt[diskaccrind])]]
        diskaccrmass = diskaccrm[diskaccrind[sort(diskaccrt[diskaccrind])]]
        reexpelltimes = reexpellt[reexpellind[sort(reexpellt[reexpellind])]]
        IF max(diskaccrtimes) GT min(reexpelltimes) THEN BEGIN
            FOR k = 0, n_elements(reexpelltimes) - 1 DO BEGIN
                IF (where(diskaccrtimes GT reexpelltimes[k]))[0] NE -1 THEN BEGIN
                    diskaccrtimes = diskaccrtimes[where(diskaccrtimes GT reexpelltimes[k])]  
                    diskaccrmass = diskaccrmass[where(diskaccrtimes GT reexpelltimes[k])]  
                    outflowtime = [outflowtime, diskaccrtimes[0] - reexpelltimes[k]]
                    outflowmass = [outflowmass, diskaccrmass[0]]
                    IF diskaccrtimes[0] - reexpelltimes[k] LT 0.001 THEN stop
                ENDIF
            ENDFOR
;            stop
        ENDIF ;ELSE stop
    ENDFOR
    outflowtime = outflowtime[1:n_elements(outflowtime) - 1]
;    disktime =disktime[1:n_elements(disktime) - 1]
    outflowmass = outflowmass[1:n_elements(outflowmass) - 1]
;    diskmass =diskmass[1:n_elements(diskmass) - 1]
    openw,1,dirs[i] + '/expelltime.grp' + halo[i] + '.txt'
    printf,1,transpose([[outflowtime],[outflowmass]])
    close,1
;    openw,1,dirs[i] + '/disktime.grp' + halo[i] + '.txt'
;    printf,1,transpose([[disktime],[diskmass]])
;    close,1
;    stop
 ENDFOR

END
