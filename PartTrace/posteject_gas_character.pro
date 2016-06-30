;CC, 12/12/11

;This program determins the characteristics (density, distance from
;center etc.) of the gas at the step
;immediately after ejection.  It is based on fb_gas_character

;posteject_gas_character,'/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK'
PRO posteject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid,laststep = laststep,totalmass = totalmass

IF NOT keyword_set(finalid) THEN finalid = '1'
IF NOT keyword_set(laststep) THEN laststep = '512'

cd,dir
spawn,'ls ' + dir + '*' + laststep + '/*' + laststep + '.coolontime',file_coolon
spawn,'ls ' + dir + '*' + laststep + '/*' + laststep + '.iord',file_iord
spawn,'ls ' + dir + '*' + laststep + '/*' + laststep,file
spawn,'ls ' + dir + 'h*param',pfile

units = tipsyunits(pfile[0])
rtipsy,file,h,g,d,s,/justhead
readarr,file_coolon,h,coolon,part = 'gas',/ascii
readarr,file_iord,h,iord,part = 'gas',/ascii
timeend = 13.734377

;Reading in files
;gpart      = mrdfits(dir + 'grp' + finalid + '.allgas.entropy.fits',1)
piord      = mrdfits(dir + 'grp' + finalid + '.eject_iord.fits',0)
ejectmass  = mrdfits(dir + 'grp' + finalid + '.mass_at_eject.fits',0)
ejectz     = mrdfits(dir + 'grp' + finalid + '.eject_z.fits',0)
expellmass = mrdfits(dir + 'grp' + finalid + '.mass_at_expell.fits',0)
expellz    = mrdfits(dir + 'grp' + finalid + '.expell_z.fits',0)

;discard all gas particles that are never expelled (z = 99)
ind99 = where(ejectz  EQ 99, comp=n99)
ejectiord  = piord[n99]
ejectmass  = ejectmass[n99]
ejectz     = ejectz[n99]  ;The redshift at which they have first left the disk
ejectlb = wmap3_lookback(ejectz)

ind99 = where(expellz EQ 99, comp=n99)
expelliord = piord[n99]
expellmass = expellmass[n99]
expellz    = expellz[n99] ;The redshift at which they have first left the halo
expelllb = wmap3_lookback(expellz)

halodat = mrdfits('alignment.fits',1)
IF (max(halodat.xc) LE 0.1 AND max(halodat.xc) GT -0.1) THEN BEGIN
    halodat.xc = halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
    halodat.yc = halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
    halodat.zc = halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
ENDIF
halodat.xc = halodat.xc*1000.0 - units.lengthunit/2.0 ;"center" is in Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
halodat.yc = halodat.yc*1000.0 - units.lengthunit/2.0
halodat.zc = halodat.zc*1000.0 - units.lengthunit/2.0
a = [[halodat.xa],[halodat.ya],[halodat.za]]

match2,fix(ejectz*1000)/1000.0,fix(halodat.z*1000)/1000.0,eject_dstep,temp
;eject_dstep = eject_dstep ;The step when the ejected particles were first ejected from the dsik

match2,fix(expellz*1000)/1000.0,fix(halodat.z*1000)/1000.0,expell_hstep,temp
;expell_hstep= expell_hstep ;The step when the expelled particles were first expelled from the halo
match,ejectiord,expelliord,eject_expell_ind,temp
match2,ejectz[eject_expell_ind],halodat.z,expell_dstep,temp
;expell_dstep = expell_dstep ;The step when the expelled particles were first out of the disk
expelltime_hist = hist_2d(expell_dstep,expell_hstep)
;contour,alog10(expelltime_hist),nlevels = 254,/fill,min = 1e-5,xtitle = 'Time left disk',ytitle = 'Time left halo'
;oplot,[0,40],[0,40]

basestruct = {iord:0L, igasorder:0L, mark:0L, mass:double(0.0), x:double(0.0), y:double(0.0), z:double(0.0), vx:double(0.0), vy:double(0.0), vz:double(0.0), rho:double(0.0), temp:double(0.0), metallicity:double(0.0), haloid:0L}
ejecthistory  = replicate(basestruct,n_elements(ejectz))
expellhistory = replicate(basestruct,n_elements(expellz))

FOR i = 0, n_elements(halodat.file) - 2 DO BEGIN
   z = halodat[i].z
   scale = 1.0/(1.0 + z)
   az = [halodat[i].xa,halodat[i].ya,halodat[i].za]
   az0 = az[0]
   az1 = az[1]
   az2 = az[2]
   ax = [az2/sqrt(az0*az0 + az2*az2),0,-1.0*az0/sqrt(az0*az0 + az2*az2)]
   ay = crossp(az,ax)    ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
   basis = [[ax],[ay],[az]]

   gashistory = mrdfits(halodat[i].file + '.allgas.history.fits',1)
   eject_step_ind  = where(eject_dstep  EQ i)
   expell_step_ind = where(expell_dstep EQ i)
   IF eject_step_ind[0] NE -1 THEN BEGIN
      match,ejectiord[eject_step_ind],gashistory.iord,eject_ind,gashistory_ind
      eject_ind = eject_step_ind[eject_ind]
      ejecthistory[eject_ind] = gashistory[gashistory_ind]
      gpos = [[ejecthistory[eject_ind].x - halodat[i].xc],[ejecthistory[eject_ind].y - halodat[i].yc],[ejecthistory[eject_ind].z - halodat[i].zc]]*scale
      gpos = transpose(transpose(basis)#transpose(gpos))
      ejecthistory[eject_ind].x = reform(gpos[*,0])
      ejecthistory[eject_ind].y = reform(gpos[*,1])
      ejecthistory[eject_ind].z = reform(gpos[*,2])
      histogramp,sqrt(ejecthistory[eject_ind].x*ejecthistory[eject_ind].x + ejecthistory[eject_ind].y*ejecthistory[eject_ind].y),nbins = 100
  ENDIF
;  stop
   IF expell_step_ind[0] NE -1 THEN BEGIN
      match,expelliord[expell_step_ind],gashistory.iord,expell_ind,gashistory_ind
      expell_ind = expell_step_ind[expell_ind]
      expellhistory[expell_ind] = gashistory[gashistory_ind]
      gpos = [[expellhistory[expell_ind].x - halodat[i].xc],[expellhistory[expell_ind].y - halodat[i].yc],[expellhistory[expell_ind].z - halodat[i].zc]]*scale
      gpos = transpose(transpose(basis)#transpose(gpos))
      expellhistory[expell_ind].x = reform(gpos[*,0])
      expellhistory[expell_ind].y = reform(gpos[*,1])
      expellhistory[expell_ind].z = reform(gpos[*,2])
      histogramp,sqrt(expellhistory[expell_ind].x*expellhistory[expell_ind].x + expellhistory[expell_ind].y*expellhistory[expell_ind].y),nbins = 100,color = 240,/overplot
;      histogramp,expellhistory[expell_ind].z,nbins= 100
  ENDIF
;  stop
ENDFOR
stop
mwrfits,expellhistory,'grp' + finalid + '.postexpell_disk.fits',/create
mwrfits,ejecthistory,'grp' + finalid + '.posteject_disk.fits',/create
END

PRO eject_gas_plot,dirs

n = n_elements(dirs)

IF KEYWORD_SET(outplot) THEN BEGIN
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
IF KEYWORD_SET(colors) THEN BEGIN
    loadct,39
    IF NOT keyword_set(ctables) THEN ctables = [39,39,39]
    IF colors[0] eq 1 THEN  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN ctables = [0,0,0]
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT KEYWORD_SET(thicks) THEN thicks = (findgen(n) + 1)*6/n - 1
ENDELSE

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
maxtime = wmap3_lookback(1000)

IF KEYWORD_SET(outplot) THEN  device,filename = outplot + '_ejectchar.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize

END


PRO eject_gas_character_master,outplot = outplot
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

dir = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK','/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK']
colors = [50,240]

dir = ['/astro/store/nbody3/fabio/h986/3072g1bwK',$
       '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK',$
       '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK']
colors = [50,120,240]
thicks = [6,6,6]

dir = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.3072g1bwK',$
       '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK',$
       '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK']
colors = [50,120,240]
thicks = [6,6,6]

dir = [prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/']

formatplot,outplot=outplot,/thick
IF keyword_set(outplot) THEN device,filename = outplot,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch,/color

FOR i = 0, n_elements(dir)-1 DO BEGIN
    cd,dir[i]
    eject_gas = eject_gas_character(dir[i])
    IF i EQ 0 THEN histogramp,sqrt(eject_gas.x*eject_gas.x + eject_gas.y*eject_gas.y),xrange = [0,30],/normalize,yrange = [0,0.1],xtitle = 'Radius [kpc]',ytitle = '1/N dN/dr'
    histogramp,sqrt(eject_gas.x*eject_gas.x + eject_gas.y*eject_gas.y),/overplot,/normalize,color = colors[i],thick = thicks[i]
ENDFOR

IF keyword_set(outplot) THEN device,/close

END
