pro find_all_gas, filebase, finalid = finalid
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
  rheader, basefile[j], h
  ngas[j] = h.ngas
endfor

z0iord = read_lon_array(iordfile[0])
z0iord_gas = z0iord[0:ngas[0]-1]
z0grp = read_lon_array(grpfile[0])
z0grp_gas = z0grp[0:ngas[0]-1]
stat = read_stat_struc_amiga(statfile[0])
main = where(halo[0] EQ stat.group)
satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
FOR j = 0, n_elements(satellites) - 1 DO BEGIN
   test = where(z0grp_gas EQ stat[satellites[j]].group, ntest)
   IF ntest NE 0 THEN $
      IF keyword_set(cum_iord) THEN cum_iord = [cum_iord,z0iord_gas[test]] ELSE cum_iord = z0iord_gas[test]
;   print,stat[satellites[j]].group,ntest
ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; The ntest = 0 was added for tracing dSphs 06/2011
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FOR i=1,n_elements(basefile)-1 DO BEGIN
  hziord = read_lon_array(iordfile[i])
  hziord_gas = hziord[0:ngas[i]-1]
  grp = read_lon_array(grpfile[i])
  grpgas = grp[0:(ngas[i]-1)]
  stat = read_stat_struc_amiga(statfile[i])
  main = where(halo[i] EQ stat.group)
  satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
;  print,'Step: ',i
  FOR j = 0, n_elements(satellites) - 1 DO BEGIN
     test = where(grpgas EQ stat[satellites[j]].group, ntest)
     IF ntest NE 0 THEN $
        IF keyword_set(cum_iord) THEN cum_iord = [cum_iord,hziord_gas[test]] ELSE cum_iord = hziord_gas[test]
;     IF where(hziord_gas[test] eq 2192543) NE -1 THEN stop
;     print,stat[satellites[j]].group,ntest
  ENDFOR
;  print,n_elements(cum_iord)
ENDFOR

;cum_iord = cum_iord(where(cum_iord NE -1))
sorted = cum_iord(sort(cum_iord))
unique = sorted(uniq(sorted))
outfile = 'grp'+strtrim(finalid,2)+'.allgas.iord.fits'
mwrfits, unique, outfile, /create


end



