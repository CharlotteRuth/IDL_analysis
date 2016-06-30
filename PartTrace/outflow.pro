;haloidoutfile = 'h516.cosmo25cmb.3072g14HBWK.grp1.haloid.dat'
;haloidoutfile = 'h603.cosmo50cmb.3072g14HBWK.haloid.dat'

;read in the amiga.stat file for each file and find the center and
;virial radius of each halo

;then, make condition sqrt((phase[j].x - cx)^2 + (phase[j].y - cy)^2 +
;                                        (phase[j].z - cz)^2) le rvir

pro outflow, haloidoutfile, haloid, phase = phase
IF NOT keyword_set(phase) THEN phase = mrdfits('grp' + haloid + '.allgas.entropy.fits',1)
readcol, haloidoutfile, files, halo, format='a,l'
files = reverse(files)
halo = reverse(halo)
nsteps = n_elements(files)
gtp_file = files + '.amiga.gtp'
stat_file = files + '.amiga.stat'
spawn,'ls h*param',pfile
units = tipsyunits(pfile)

stat = read_stat_struc_amiga(stat_file[0])
stat_array = replicate(stat[0],N_ELEMENTS(files))
ind = where(stat.group eq halo[0])
stat_array[0] = stat[ind]
FOR i=0, N_ELEMENTS(files) - 1 DO BEGIN
    stat = read_stat_struc_amiga(stat_file[i])
    stat_array[i] = stat[where(stat.group eq halo[i])]        
ENDFOR
stat_array.xc = stat_array.xc*1000.0 - units.lengthunit/2
stat_array.yc = stat_array.yc*1000.0 - units.lengthunit/2
stat_array.zc = stat_array.zc*1000.0 - units.lengthunit/2
;Need to exclude 'early' particles to compare to accrz file
early = mrdfits('early.iord.fits',0)
all = phase.iord
test = binfind(all,early)
all[test] = -1
exlude = where(all eq -1, comp=keep)
phase = phase[keep] 
ngas = n_elements(phase.iord)
z = fltarr(nsteps)
for j=0L,nsteps-1 do begin
  rheader, gtp_file[j], h
  z[j] = (1./h.time)-1.
endfor
accrz = mrdfits('grp' + haloid + '.accrz.fits',0)

print, ngas
outf = lonarr(ngas)
outz = fltarr(ngas)

FOR j=0L,ngas-1 do begin
;   ingal = where(phase[j].grp eq halo, comp=outgal) 
    ingal = where(sqrt((phase[j].x - stat_array.xc)^2 + (phase[j].y - stat_array.yc)^2 + (phase[j].z - stat_array.zc)^2) le stat_array.rvir, comp=outgal)
    IF (outgal)[0] NE -1 THEN test = where(z[outgal] lt accrz[j] and phase[j].grp[outgal] ne -1, ntest) ELSE ntest = 0 ;Find only particles that get ejected at lower z than accreted, and aren't deleted
    if ntest ne 0 then outf[j] = min(outgal[test]) else outf[j] = -1 
;    IF(phase[j].iord EQ 139718 OR phase[j].iord EQ 591774) THEN stop
ENDFOR

;To prevent a whole slew of particles appearing to have outflows at z=0, replace them with z=99
ind99 = where(outf eq -1, comp=keep)
outz[keep] = z[outf[keep]]
outz[ind99] = 99 ;If the gas is never expelled (never outside of rvir after it is accreted) then set to zero
mwrfits, outz, 'grp' + haloid + '.rvir.outflow_z.fits', /create

IF 0 THEN BEGIN
    FOR i=1,N_ELEMENTS(files) - 1 DO BEGIN
        ind = where(outz LT z[i] AND accrz GE z[i])
        indgrp = where(phase.grp[i] EQ halo[i])
        IF ind[0] NE -1 THEN plot,phase[ind].x[i] - stat_array[i].xc,phase[ind].y[i] - stat_array[i].yc,psym = 3,xrange = [-100,100],yrange = [-100,100] ELSE print,'Nothing'
        oplot,phase[indgrp].x[i] - stat_array[i].xc,phase[indgrp].y[i] - stat_array[i].yc,psym = 3,color = 240
    ENDFOR
ENDIF

massatout = fltarr(n_elements(outf))
iordout = lonarr(ngas)
FOR j=1,nsteps-1 do begin
  ind1 = where(outf eq j, nind)
  if nind ne 0 then massatout[ind1] = phase[ind1].mass[j]
  if nind ne 0 then iordout[ind1] = phase[ind1].iord
ENDFOR
;If the gas is never expelled (never outside of rvir after it is
;accreted) then iord and mass are set to sero

mwrfits, massatout, 'grp' + haloid + '.rvir.mass_at_outflow.fits', /create
mwrfits, iordout, 'grp' + haloid + '.rvir.outflow_iord.fits', /create

end

