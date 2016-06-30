PRO time_cycling2,dirs,files,halo = halo,finalstep = finalstep,colors = colors,outplot = outplot
n = n_elements(dirs)
maxtimes = 20
IF NOT keyword_set(halo) THEN halo = strarr(n) + '1'
IF NOT keyword_set(finalstep) THEN finalstep = '00512'

formatplot,outplot = outplot,thick = formatthick

IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 800
    ysize = 500
ENDELSE
IF keyword_set(colors) THEN BEGIN
    loadct,39
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(ctables) THEN ctables = findgen(n)
    IF NOT keyword_set(thicks) THEN $
       IF keyword_set(formatthick) THEN thicks = fltarr(n) + 4 ELSE thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2)  
ENDELSE

timecycled = [0]
smass = fltarr(n)
vmass = fltarr(n)
FOR i = 0, n - 1 DO BEGIN
    stat = read_stat_struc_amiga(dirs[i] + files[i] + '.' + finalstep + '/' + files[i] + '.' +  finalstep + '.amiga.stat')
    ind = (where(stat.group eq halo[i]))[0]
    vmass[i] = stat[ind].m_tot
    smass[i] = stat[ind].m_star
    
    diskaccri  = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_iord.fits',0) 
    diskaccrm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reaccrdisk.fits',0)
    diskaccrz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reaccrdisk_z.fits',0)
    diskaccrt = z_to_t(diskaccrz)

    timesaccr = histogram(diskaccri, reverse_indices = r)
;    IF R[0] NE R[1] THEN diskaccri[R[R[0] : R[1]-1]] = -1 ;Set to zero all 
    accrone = where(timesaccr EQ 1)
    FOR ih = 0L, n_elements(accrone) - 1 DO IF R[accrone[ih]] NE R[accrone[ih] + 1] THEN diskaccri[R[R[accrone[ih]] : R[accrone[ih] + 1]-1]] = -1 ;Set to -1 all particles that are accreted only once
    reaccr = where(diskaccri NE -1)
    stop
    diskaccri = diskaccri[reaccr]
    diskaccrm = diskaccrm[reaccm]
    diskaccrz = diskaccrz[reaccz]
    diskaccrt = diskaccrt[reacct]
    sorted = sort(diskaccri)
    diskaccri = diskaccri[sorted]
    diskaccrm = diskaccrm[sorted]
    diskaccrz = diskaccrz[sorted]
    diskaccrt = diskaccrt[sorted]

    reejecti = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_iord.fits',0)
    reejectm = mrdfits(dirs[i] + '/grp' + halo[i] + '.mass_at_reeject.fits',0)
    reejectz = mrdfits(dirs[i] + '/grp' + halo[i] + '.reeject_z.fits',0)
    reejectt = z_to_t(reejectz)

 ENDFOR

END
