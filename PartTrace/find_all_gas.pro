pro find_all_gas, filebase, finalid = finalid, include_halo_satellites = include_halo_satellites
IF NOT keyword_set(finalid) THEN finalid = '1'

haloidoutput = filebase + '.grp' + finalid + '.haloid.dat'
;haloidoutput includes z=0 file
;order doesn't matter 
;generally will run from high z to low z

readcol, haloidoutput, basefile, halo, format='a,l'
iordfile = basefile+'.iord'
grpfile = basefile+'.amiga.grp'
statfile = basefile+'.amiga.stat'
ngas = lonarr(n_elements(basefile))
for j=0L,n_elements(basefile)-1 do begin
    print,basefile[j]
  rheader, basefile[j], h
;  rtipsy, basefile[j], h, g, d, s,/justhead
  ngas[j] = h.ngas
endfor

z0iord = read_lon_array(iordfile[0])
z0iord_gas = z0iord[0:ngas[0]-1]
z0grp = read_lon_array(grpfile[0])
z0grp_gas = z0grp[0:ngas[0]-1]
stat = read_stat_struc_amiga(statfile[0])
main = where(halo[0] EQ stat.group)
IF keyword_set(include_halo_satellites) THEN extra_halos = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir) ELSE BEGIN
    IF (where((long(stat.sat) EQ halo[0]) AND (stat.m_dark LT 2*stat.m_gas)))[0] NE -1 THEN extra_halos = [main,where((long(stat.sat) EQ halo[0]) AND (stat.m_dark LT 2*stat.m_gas))] ELSE extra_halos = main
ENDELSE
FOR j = 0, n_elements(extra_halos) - 1 DO BEGIN
   test = where(z0grp_gas EQ stat[extra_halos[j]].group, ntest)
   IF ntest NE 0 THEN BEGIN
      IF keyword_set(cum_iord) THEN cum_iord = [cum_iord,z0iord_gas[test]] ELSE cum_iord = z0iord_gas[test]
      IF j GT 0 THEN print,'Extra halo: ',strtrim(stat[extra_halos[j]].group,2),', n_gas: ',strtrim(n_elements(test),2) ELSE print,'Main halo: ',strtrim(n_elements(test),2)
  ENDIF
;   print,stat[extra_halos[j]].group,ntest
ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; The ntest = 0 was added for tracing dSphs 06/2011
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FOR i=1,n_elements(basefile)-1 DO BEGIN
    print,iordfile[i]
  hziord = read_lon_array(iordfile[i])
  hziord_gas = hziord[0:ngas[i]-1]
  grp = read_lon_array(grpfile[i])
  grpgas = grp[0:(ngas[i]-1)]
  stat = read_stat_struc_amiga(statfile[i]) 
  main = where(halo[i] EQ stat.group)
;For the g14 analysis papers, all satellites (defined as having
;centroid within the halo) were included. For future work, exculde
;satellites but include parts of the disk that amiga thinks are satellites
  IF keyword_set(include_halo_satellites) THEN extra_halos = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir) ELSE BEGIN
      IF (where((long(stat.sat) EQ halo[i]) AND (stat.m_dark LT 2*stat.m_gas)))[0] NE -1 THEN extra_halos = [main,where((long(stat.sat) EQ halo[i]) AND (stat.m_dark LT 2*stat.m_gas))] ELSE extra_halos = main
  ENDELSE
  print,'Step: ',i
  IF extra_halos[0] NE -1 THEN BEGIN
      FOR j = 0, n_elements(extra_halos) - 1 DO BEGIN
          test = where(grpgas EQ stat[extra_halos[j]].group, ntest)
          IF ntest NE 0 THEN BEGIN
            IF keyword_set(cum_iord) THEN cum_iord = [cum_iord,hziord_gas[test]] ELSE cum_iord = hziord_gas[test]
            IF j GT 0 THEN print,'Extra halo: ',stat[extra_halos[j]].group
        ENDIF
;     IF where(hziord_gas[test] eq 2192543) NE -1 THEN stop
;     print,stat[extra_halos[j]].group,ntest
      ENDFOR
  ENDIF
;  print,n_elements(cum_iord)
  ENDFOR

IF keyword_set(cum_iord) THEN BEGIN
;cum_iord = cum_iord(where(cum_iord NE -1))
    sorted = cum_iord(sort(cum_iord))
    unique = sorted(uniq(sorted))
    outfile = 'grp'+strtrim(finalid,2)+'.allgas.iord.fits'
    mwrfits, unique, outfile, /create
ENDIF ELSE print,'No gas particles ever detected in halo'

end



