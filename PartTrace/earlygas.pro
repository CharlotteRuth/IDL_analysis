;This program will output a fits file with the iords of the gas
;particles at the selected step

PRO earlygas,dir,file,step,finalid,centralhalo = centralhalo
halodat = mrdfits(dir + '/grp' + finalid + '.alignment.fits',1)
istep = where(halodat.file eq file + '.' + step + '/' + file + '.' + step)
haloid = halodat[istep].haloid

filename = dir + '/' + file + '.' + step + '/' + file + '.' + step
rtipsy,filename,h,g,d,s
readarr,filename + '.iord',h,iord,/ascii,part = 'gas'
readarr,filename + '.OxMassFrac',h,ox,/ascii,part = 'gas'
readarr,filename + '.FeMassFrac',h,fe,/ascii,part = 'gas'
stat = read_stat_struc_amiga(filename + '.amiga.stat')
readarr,filename + '.amiga.grp',h,grp,/ascii,part = 'gas'
units = tipsyunits(dir + '/' + file + '.param')
g.dens = g.dens*units.rhounit/h.time^3
g.mass = g.mass*units.massunit
g.x = (g.x*units.lengthunit/1000.0 + units.lengthunit/1000.0/2)
g.y = (g.y*units.lengthunit/1000.0 + units.lengthunit/1000.0/2)
g.z = (g.z*units.lengthunit/1000.0 + units.lengthunit/1000.0/2)

main = where(haloid EQ stat.group)
satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir AND stat.ngas gt 0)
inhalo = [0]
IF keyword_set(centralhalo) THEN BEGIN
   test = where(grp EQ stat[main].group, ntest)
   IF ntest NE 0 THEN inhalo = [inhalo,test]
ENDIF ELSE BEGIN
   FOR j = 0, n_elements(satellites) - 1 DO BEGIN
      test = where(grp EQ stat[satellites[j]].group, ntest)
      IF ntest NE 0 THEN inhalo = [inhalo,test]
   ENDFOR
ENDELSE
inhalo = inhalo[1:n_elements(inhalo)-1]
g = g[inhalo]
iord = iord[inhalo]

mwrfits,iord,dir + '/grp' + finalid + '.earlyhalo_iord.fits',/create
mwrfits,g.mass,dir + '/grp' + finalid + '.earlyhalo_mass.fits',/create
mwrfits,g.mass*(2.09*ox + 1.06*fe),dir + '/grp' + finalid + '.earlyhalo_zmetals.fits',/create

;IF keyword_set(disk) THEN BEGIN
ind = where(g.dens GE 0.1 AND g.tempg LE 1.2e4 AND abs(g.z - stat[main].zc)*h.time*1000.0 LE 3.0)
g = g[ind]
iord = iord[ind]
mwrfits,iord,dir + '/grp' + finalid + '.earlydisk_iord.fits',/create
mwrfits,g.mass,dir + '/grp' + finalid + '.earlydisk_mass.fits',/create
mwrfits,g.mass*(2.09*ox[ind] + 1.06*fe[ind]),dir + '/grp' + finalid + '.earlydisk_zmetals.fits',/create
;ENDIF

END
