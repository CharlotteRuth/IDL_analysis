;Tracks the mass of the halo, using information from grp?.metal.txt
;and grp?.alignment.fits

PRO track_mass_2,finalid = finalid
IF NOT keyword_set(finalid) THEN finalid = '1'

halodat = mrdfits('grp'+finalid+'.alignment.fits',1)
nsteps = n_elements(halodat)

mtot = fltarr(nsteps)
mgas = fltarr(nsteps)
mstar = fltarr(nsteps)
mdark = fltarr(nsteps)
FOR i = 0, nsteps -1 DO BEGIN
   amigastat = read_stat_struc_amiga(halodat[i].file + '.amiga.stat')
   ind = where(amigastat.group EQ halodat[i].haloid)
   mtot[i] = amigastat[ind].m_tot
   mgas[i] = amigastat[ind].m_gas
   mstar[i] = amigastat[ind].m_star
   mdark[i] = amigastat[ind].m_dark
ENDFOR
readcol,'grp' + finalid + '.metals.txt',metal,mox,mfe,mHI,mH2,mcoldgas
mass = [[halodat.haloid],[halodat.time],[halodat.z],[mtot],[mgas],[mstar],[mdark],[metal],[mox],[mfe],[mHI],[mH2],[mcoldgas]]
;print,transpose(mass)
format = '(I,d15.3, d15.3, E15.5, E15.5, E15.5, E15.5, d15.5, d15.5, d15.5, E15.5, E15.5, E15.5)'
openw,1,'grp'+finalid+'.mass_metal_track.dat'
printf,1,'Halo','Time','Redshift','M_tot','M_gas','M_star','M_dark','Metal','Ox','Fe','HI','H2','M_coldg',format='(13A20)'
printf,1,transpose(mass),format = format
close,1


END
