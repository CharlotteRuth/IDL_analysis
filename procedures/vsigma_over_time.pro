;Tracks the z-velocity dispersion of the gas in a galaxy over time using the
;aligment data from grp?.alignment.fits

PRO track_vsigma_master
  spawn,'hostname',hostname
  IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' $
  ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN prefix = '/nobackupp8/crchrist/MolecH/' $
  ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

  dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
  file799 = 'h799.cosmo25cmb.3072g14HBWK'
  key799 = 'h799'
  dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
  file516 = 'h516.cosmo25cmb.3072g14HBWK'
  key516 = 'h516'
  dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
  file986 = 'h986.cosmo50cmb.3072g14HBWK'
  key986 = 'h986'
  dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
  file603 = 'h603.cosmo50cmb.3072g14HBWK'
  key603 = 'h603'
  dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
  file258 = 'h258.cosmo50cmb.3072g14HMbwK'
  key258 = 'h258'
  dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
  file285 = 'h285.cosmo50cmb.3072g14HMbwK'
  key285 = 'h285'
  dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
  file239 = 'h239.cosmo50cmb.3072g14HMbwK'
  key239 = 'h239'
  dirs    = [ dir799, dir516, dir986, dir603, dir258, dir285, dir239]
  files   = [file799,file516,file986,file603,file258,file285,file239]
  finalstep = '00512'
  haloid =  [ '1'   , '1'   , '1'   , '1'   , '1'   , '1'   , '1'   ]
  
;  dir='/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
;  finalids = ['1','4','6']
;  units = tipsyunits('h799.cosmo25cmb.3072g14HBWK.param')
;  cd,dir
  FOR i=3, n_elements(dirs)-1 DO BEGIN
     cd,dirs[i]
     units = tipsyunits(files[i] + '.param')
     track_vsigma,files[i],units,finalid = haloid[i]
  ENDFOR
END

PRO track_vsigma,filename,units,finalid = finalid
IF NOT keyword_set(finalid) THEN finalid = '1'
print,filename

halodat = mrdfits('grp'+finalid+'.alignment.fits',1)
nsteps = n_elements(halodat)

vsigma = fltarr(nsteps)
vsigma_hi = fltarr(nsteps)
nostar = 0

FOR i = 0, nsteps -1 DO BEGIN
   amigastat = read_stat_struc_amiga(halodat[i].file + '.amiga.stat')
   rtipsy,halodat[i].file,h,g,d,s
   read_tipsy_arr,halodat[i].file + '.HI',h,hi,part = 'gas';,/ascii
   hi_frac = hi/0.758375
   read_tipsy_arr,halodat[i].file + '.amiga.grp',h,grp;,/ascii

   ind = where(grp eq halodat[i].haloid,comp=indcomp)
   IF ((where(ind ge h.ngas+h.ndark))[0] ne -1) THEN inds = ind(where(ind ge h.ngas+h.ndark)) ELSE nostar = 1
   indg = ind(where(ind lt h.ngas))
   indd = ind(where(ind ge h.ngas and ind lt h.ngas+h.ndark))
   IF NOT nostar THEN stars = s[inds-h.ngas-h.ndark]
   dark = d[indd-h.ngas] 
   gas = g[indg]
   IF NOT nostar THEN BEGIN
      foo = min(stars.phi, cmindex)
      cx = stars[cmindex].x
      cy = stars[cmindex].y
      cz = stars[cmindex].z
      
      stars.x =  stars.x - cx
      stars.y =  stars.y - cy
      stars.z =  stars.z - cz
   ENDIF ELSE BEGIN
      foo = min(dark.phi, cmindex)
      cx = dark[cmindex].x
      cy = dark[cmindex].y
      cz = dark[cmindex].z    
   ENDELSE
   dark.x = dark.x - cx
   dark.y = dark.y - cy
   dark.z = dark.z - cz
   gas.x = gas.x - cx
   gas.y = gas.y - cy
   gas.z = gas.z - cz
    IF NOT nostar THEN BEGIN
        propx=mean(stars.vx)
        propy=mean(stars.vy)
        propz=mean(stars.vz)
    ENDIF ELSE BEGIN
        propx=mean(dark.vx)
        propy=mean(dark.vy)
        propz=mean(dark.vz)
    ENDELSE

   IF NOT nostar THEN rads = sqrt(stars.x*stars.x+stars.y*stars.y+stars.z*stars.z)
   radd = sqrt(dark.x*dark.x+dark.y*dark.y+dark.z*dark.z)
   radg = sqrt(gas.x*gas.x+gas.y*gas.y+gas.z*gas.z)
   IF NOT nostar THEN maxrad = MAX([rads,radg,radd]) ELSE maxrad = MAX([radg,radd])
   g.x = g.x - cx
   g.y = g.y - cy
   g.z = g.z - cz
   g.vx = g.vx - propx
   g.vy = g.vy - propy
   g.vz = g.vz - propz
   radg = sqrt(g.x*g.x+g.y*g.y+g.z*g.z)
   indg = where(radg LT maxrad)
   g = g[indg]
   hi_frac = hi_frac[indg]

   halodat = mrdfits('grp' + finalid + '.alignment.fits',1)
   a = [[halodat.xa],[halodat.ya],[halodat.za]]
   az = reform(a[i,*])
   az0 = az[0]
   az1 = az[1]
   az2 = az[2]
   ax = [az2/sqrt(az0*az0 + az2*az2),0,-1.0*az0/sqrt(az0*az0 + az2*az2)]
   ay = crossp(az,ax)           ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
   basis = [[ax],[ay],[az]]
   gpos = [[g.x],[g.y],[g.z]] ;*scale
   gpos = transpose(transpose(basis)#transpose(gpos))
   g.x = gpos[*,0]
   g.y = gpos[*,1]
   g.z = gpos[*,2]
   gvel = [[g.vx],[g.vy],[g.vz]]
   gvel = transpose(transpose(basis)#transpose(gvel))
   g.vx = gvel[*,0]
   g.vy = gvel[*,1]
   g.vz = gvel[*,2]

   hi_frac = hi_frac[where(abs(g.z*units.lengthunit*h.time) LT 2)]
   g = g[where(abs(g.z*units.lengthunit*h.time) LT 2)] ;Select for only that gas that is within 2 kpc of the disk plane
   g.vz = g.vz*units.vunit*h.time
   gdisk = g[where(g.tempg LT 2e4 AND g.dens*units.rhounit*h.time^(-3) GT 0.1)]
   histogramp,gdisk.vz,min= -1*max(abs(gdisk.vz)),max  = max(abs(gdisk.vz))
   histogramp,g.vz,/overplot,weight = hi_frac,color = 100,min= -1*max(abs(gdisk.vz)),max  = max(abs(gdisk.vz))
   vsigma[i] = stdev(gdisk.vz)
   n = n_elements(where(hi_frac NE 0))
   vsigma_hi[i] = sqrt(total(hi_frac*(g.vz - mean(g.vz))^2)/((n - 1)*total(hi_frac)/n))
   print,halodat[i].file,vsigma[i],vsigma_hi[i]
ENDFOR

mass = [[halodat.haloid],[halodat.time],[halodat.z],[vsigma],[vsigma_hi]]
format = '(I,d15.3, d15.3,d15.3, d15.3)'
openw,1,filename + '.grp'+finalid+'.vsigma.dat'
printf,1,'Halo','Time','Redshift','sigma_coldgas','sigma_HI',format='(5A15)'
printf,1,transpose(mass),format = format
close,1


END
