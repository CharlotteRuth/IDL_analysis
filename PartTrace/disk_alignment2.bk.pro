PRO disk_alignment2,pfile,finalid
  data = mrdfits('./grp' + finalid + '.alignment.fits',1)
  units = tipsyunits(pfile) 
  halodat = {file: ' ', haloid: 0, time: 0D, z: 0D, xc: 0D, yc: 0D, zc: 0D, vx: 0D, vy: 0D, vz: 0D, xa: 0D, ya: 0D, za: 0D, rvir: 0D}
  halodat = replicate(halodat,n_elements(data))
  halodat.file = data.file
  halodat.haloid = data.haloid
  halodat.time = data.time
  halodat.z = data.z 
  halodat.xc = data.xc 
  halodat.yc = data.yc 
  halodat.zc = data.zc 
  halodat.xa = data.xa 
  halodat.ya = data.ya 
  halodat.za = data.za 
  FOR i = 0, n_elements(halodat) - 1 DO BEGIN
     amigastat = read_stat_struc_amiga(data[i].file + '.amiga.stat')
 ;    halodat.xc = halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
 ;    halodat.yc = halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
 ;    halodat.zc = halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
     halodat.vx = amigastat.xvcm ;+ units.lengthunit/2/1000.0
     halodat.vy = amigastat.yvcm ;+ units.lengthunit/2/1000.0
     halodat.vz = amigastat.zvcm ;+ units.lengthunit/2/1000.0 
  ENDFOR
END
