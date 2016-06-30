pro accretion, filebase, finalid, phase, centralhalo = centralhalo

readcol,filebase + '.grp' + finalid + '.haloid.dat',files,halo,format='a,l'
files = reverse(files)
halo = reverse(halo)
nsteps = n_elements(files)
gtp_file = files+'.amiga.gtp'

;phase has been shortened to exclude particles already in the galaxy at the first step (early)
ngas = n_elements(phase.iord)

accr = lonarr(ngas)
gloc = {loc: intarr(nsteps)}
gloc = replicate(gloc, n_elements(phase)) ;set to 1 if in the halo and 0 if not

FOR i = 0, nsteps - 1 DO BEGIN
   stat = read_stat_struc_amiga(files[i] + '.amiga.stat')
   main = where(halo[i] EQ stat.group)
   satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir AND stat.ngas gt 0)   
   print,files[i],satellites
   inhalo = [0]
   IF keyword_set(centralhalo) THEN BEGIN
      test = where(gpart.grp[i] EQ stat[main].group, ntest)
      IF ntest NE 0 THEN inhalo = [inhalo,test]
   ENDIF ELSE BEGIN
 ;     print,'Step: ',i
      FOR j = 0, n_elements(satellites) - 1 DO BEGIN
         test = where(phase.grp[i] EQ stat[satellites[j]].group, ntest)
         IF ntest NE 0 THEN inhalo = [inhalo,test]
;         print,stat[satellites[j]].group,ntest
      ENDFOR
   ENDELSE
   IF n_elements(inhalo) NE 1 THEN BEGIN
      inhalo = inhalo[1:n_elements(inhalo) - 1]
      gloc[inhalo].loc[i] = 1
   ENDIF ELSE print,'PROBLEM'
;   print,n_elements(inhalo)
ENDFOR

accr     = intarr(n_elements(gloc)) - 1     ;Gas that is ever in the disk
FOR i = 0L, n_elements(accr) - 1 DO BEGIN   ;iterate through gas particles were ever in the disk
   accreted = where(gloc[i].loc EQ 1)   
   accr[i] = accreted[0]
;   IF accreted[0] EQ -1 THEN BEGIN
;        print,"PROBLEM!"
;        stop
        ;93448, 4 - 8 (halo 1), iord = 2192543, grp = 76 at 00156
;   ENDIF
ENDFOR
neverinhalo = where(accr EQ -1, complement = everinhalo)
IF neverinhalo[0] NE -1 THEN BEGIN
   print,'Unaccreted gas in "allgas" '
   stop
ENDIF

;FOR j=0L,ngas-1 do begin
;   ingal = where(phase[j].grp EQ halo)
   ;;The following ensures that the particle has been in the halo for 
   ;;two consecutive steps to count as accreted
   ;;test = lonarr(nsteps-1)
   ;;for i=0,n_elements(ingal)-2 do test[i] = ingal[i+1]-ingal[i]
   ;;ind = min(where(test eq 1))
   ;;if ind ne -1 then accr[j] = ingal[ind]
   ;;if min(ingal) eq nsteps-1 then accr[j] = nsteps-1  ;for those particles that enter at the final step 
   ;;I should test how this result differs from just 1 step: min(ingal) = accr
   ;accr[j] = min(ingal)
;ENDFOR
z = fltarr(nsteps)
for j=0L,nsteps-1 do begin
  rheader, gtp_file[j], h
  z[j] = (1./h.time)-1.
endfor
accrtime = z[accr]
mwrfits, accrtime, 'grp' + finalid + '.accrz.fits', /create

massataccr = fltarr(n_elements(accr))
iordataccr = intarr(n_elements(accr))
vxataccr   = fltarr(n_elements(accr))
vyataccr   = fltarr(n_elements(accr))
vzataccr   = fltarr(n_elements(accr))
xataccr    = fltarr(n_elements(accr))
yataccr    = fltarr(n_elements(accr))
zataccr    = fltarr(n_elements(accr))
FOR j=1,nsteps-1 do begin
  ind1 = where(accr eq j, nind)
  if nind ne 0 then begin
    massataccr[ind1] = phase[ind1].mass[j]
    iordataccr[ind1] = phase[ind1].iord
    vxataccr[ind1]   = phase[ind1].vx[j]
    vyataccr[ind1]   = phase[ind1].vy[j]
    vzataccr[ind1]   = phase[ind1].vz[j]
    xataccr[ind1]    = phase[ind1].x[j]
    yataccr[ind1]    = phase[ind1].y[j]
    zataccr[ind1]    = phase[ind1].z[j]
   endif
ENDFOR

mwrfits, massataccr, 'grp' + finalid + '.mass_at_accr.fits', /create
mwrfits, iordaraccr, 'grp' + finalid + '.iord_at_accr.fits', /create
mwrfits, vxataccr,   'grp' + finalid + '.vx_at_accr.fits',   /create
mwrfits, vyataccr,   'grp' + finalid + '.vy_at_accr.fits',   /create
mwrfits, vzataccr,   'grp' + finalid + '.vz_at_accr.fits',   /create
mwrfits, xataccr,    'grp' + finalid + '.x_at_accr.fits',    /create
mwrfits, yataccr,    'grp' + finalid + '.y_at_accr.fits',    /create
mwrfits, zataccr,    'grp' + finalid + '.z_at_accr.fits',    /create

end
