;Tracks the mass of the halo, using information from grp?.metal.txt
;and grp?.alignment.fits

PRO track_vc,finalid = finalid
IF NOT keyword_set(finalid) THEN finalid = '1'

halodat = mrdfits('grp'+finalid+'.alignment.fits',1)
nsteps = n_elements(halodat)

vc = fltarr(nsteps)
FOR i = 0, nsteps -1 DO BEGIN
   amigastat = read_stat_struc_amiga(halodat[i].file + '.amiga.stat')
   ind = where(amigastat.group EQ halodat[i].haloid)
;   vc[i] = amigastat[ind].vc
   vc[i] = sqrt(6.67e-11*amigastat[ind].m_tot*1.9891e30/(3.08567758e19*amigastat[ind].rvir))/1000.
ENDFOR

readcol,'grp' + finalid + '.metals.txt',metal,mox,mfe,mHI,mH2,mcoldgas
mass = [[halodat.haloid],[halodat.time],[halodat.z],[vc]]
;print,transpose(mass)
format = '(I, d15.3, d15.3, d15.5)'
openw,1,'grp'+finalid+'.vc.dat'
printf,1,'Halo','Time','Redshift','vc',format='(4A20)'
printf,1,transpose(mass),format = format
close,1


END
