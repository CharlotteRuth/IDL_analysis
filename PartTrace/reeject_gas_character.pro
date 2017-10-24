;CC, 12/12/11

;This program determines the characteristics (density, distance from
;center etc.) of the gas at the step
;prior to ejection, each time it is ejected.  It is based on fb_gas_character and eject_gas_character

;reeject_gas_character,'/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK',ejecthistory,expellhistory,
PRO reeject_gas_character,ejecthistory,dir = dir,finalid = finalid,laststep = laststep,postdisk = postdisk,plot = plot,expell = expell,heat = heat,rvircut_half = rvircut_half,rvircut_fifth = rvircut_fifth

IF NOT keyword_set(finalid) THEN finalid = '1'
IF NOT keyword_set(laststep) THEN laststep = '512'
IF NOT keyword_set(dir) THEN dir='./'

cd,dir
spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.coolontime',file_coolon
spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord
spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file
spawn,'ls ' + dir + '/h*param',pfile

units = tipsyunits(pfile[0])
rtipsy,file,h,g,d,s,/justhead
readarr,file_coolon,h,coolon,part = 'gas',/ascii
readarr,file_iord,h,iord,part = 'gas',/ascii
timeend = 13.734377

CASE 1 OF
   keyword_set(heat): ejectext = 'reheat_all'
   keyword_set(rvircut_half): ejectext = 'reeject_rvir0.5'
   keyword_set(rvircut_fifth): ejectext = 'reeject_rvir0.2'
   ELSE: ejectext = 'reeject'
ENDCASE

;Reading in files
;gpart      = mrdfits(dir + 'grp' + finalid + '.allgas.entropy.fits',1)
ejectiord  = mrdfits(dir + '/grp' + finalid + '.' + ejectext + '_iord.fits',0)
ejectmass  = mrdfits(dir + '/grp' + finalid + '.mass_at_' + ejectext + '.fits',0)
ejectz     = mrdfits(dir + '/grp' + finalid + '.' + ejectext + '_z.fits',0)
IF keyword_set(expell) THEN BEGIN
   expelliord = mrdfits(dir + '/grp' + finalid + '.reexpell_iord.fits',0)
   expellmass = mrdfits(dir + '/grp' + finalid + '.mass_at_reexpell.fits',0)
   expellz    = mrdfits(dir + '/grp' + finalid + '.reexpell_z.fits',0)
ENDIF

halodat = mrdfits(dir + '/grp' + finalid + '.alignment.fits',1)
IF (max(halodat.xc) LE 0.1 AND max(halodat.xc) GT -0.1) THEN BEGIN
    halodat.xc = halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
    halodat.yc = halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
    halodat.zc = halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
ENDIF
halodat.xc = (halodat.xc*1000.0 - units.lengthunit/2.0) ;"center" is in Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
halodat.yc = (halodat.yc*1000.0 - units.lengthunit/2.0)
halodat.zc = (halodat.zc*1000.0 - units.lengthunit/2.0)
a = [[halodat.xa],[halodat.ya],[halodat.za]]

match2,fix(ejectz*1000)/1000.0,fix(halodat.z*1000)/1000.0,eject_dstep,temp
IF NOT keyword_set(postdisk) THEN eject_dstep = eject_dstep - 1 ;The step when the ejected particles were last in the disk

IF keyword_set(expell) THEN BEGIN
   match2,fix(expellz*1000)/1000.0,fix(halodat.z*1000)/1000.0,expell_hstep,temp
;expell_hstep: The step when the expelled particles were last in the halo
   expell_dstep = fltarr(n_elements(expell_hstep)) - 1
   eject_expell_i = lonarr(n_elements(expell_hstep)) - 1
   hist = histogram(ejectiord,min = min(ejectiord),max = max(ejectiord),nbins = max(ejectiord)-min(ejectiord)+1,locations = histiord)
   singleiord = histiord[where(hist EQ 1)]
   match,ejectiord,singleiord,eject_single_ind,temp
   match,expelliord,ejectiord[eject_single_ind],expell_single_ind,temp
   expell_dstep[expell_single_ind] = eject_dstep[eject_single_ind[temp]]
   eject_expell_i[expell_single_ind] = eject_single_ind[temp];The indice of the ejected data

   multiind = where(expell_dstep EQ -1) ;n_elements(multiind)
   FOR i = 0, n_elements(multiind) -1  DO BEGIN & $
      matchind = where(ejectiord EQ expelliord[multiind[i]] AND eject_dstep LE expell_hstep[multiind[i]]) & $
;   print,ejectiord[matchind],expelliord[multiind[i]],eject_dstep[matchind],expell_hstep[multiind[i]] & $
      expell_dstep[multiind[i]] = eject_dstep[matchind[n_elements(matchind) -1]] & $
      eject_expell_i[multiind[i]] = matchind[n_elements(matchind) -1]
   ENDFOR
   print,minmax(expelliord - ejectiord[eject_expell_i])
   mwrfits,eject_expell_i,dir+'/grp'+finalid+'.eject_expell_iord.fits',/create

   expelltime_hist = hist_2d(expell_dstep,expell_hstep)
   contour,alog10(expelltime_hist),nlevels = 254,/fill,min = 0,max = 4,xtitle = 'Time left disk',ytitle = 'Time left halo'
   oplot,expell_dstep,expell_hstep,psym = 3
   oplot,[0,140],[0,140]
ENDIF

basestruct = {iord:0L, igasorder:0L, mark:0L, mass:double(0.0), x:double(0.0), y:double(0.0), z:double(0.0), vx:double(0.0), vy:double(0.0), vz:double(0.0), rho:double(0.0), temp:double(0.0), metallicity:double(0.0), haloid:0L}
ejecthistory  = replicate(basestruct,n_elements(ejectz))
ejecthistory.iord = ejectiord

FOR i = 0, n_elements(halodat.file) - 1 DO BEGIN
   satsdata = read_stat_struc_AMIGA(halodat[i].file + '.amiga.stat')
   isats = where(satsdata.group EQ halodat[i].haloid)
   z = halodat[i].z
   scale = 1.0/(1.0 + z)
   az = [halodat[i].xa,halodat[i].ya,halodat[i].za]
   az0 = az[0]
   az1 = az[1]
   az2 = az[2]
   ax = [az2/sqrt(az0*az0 + az2*az2),0,-1.0*az0/sqrt(az0*az0 + az2*az2)]
   ay = crossp(az,ax)    ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
   basis = [[ax],[ay],[az]]

   gashistory = mrdfits(halodat[i].file + '.grp' + finalid + '.allgas.history.fits',1)
   eject_step_ind  = where(eject_dstep  EQ i) ;The ejcted particles that are in the disk at step i
;   expell_step_ind = where(expell_dstep EQ i)
   IF eject_step_ind[0] NE -1 THEN BEGIN
      match,ejectiord[eject_step_ind],gashistory.iord,eject_ind,gashistory_ind
 ;     eject_ind = eject_step_ind[eject_ind]
      gpos = [[gashistory[gashistory_ind].x - halodat[i].xc*scale],[gashistory[gashistory_ind].y - halodat[i].yc*scale],[gashistory[gashistory_ind].z - halodat[i].zc*scale]]
      gvel = [[gashistory[gashistory_ind].vx - satsdata[isats].xvcm],[gashistory[gashistory_ind].vy - satsdata[isats].yvcm],[gashistory[gashistory_ind].vz - satsdata[isats].zvcm]]
      gpos = transpose(transpose(basis)#transpose(gpos))
      gvel = transpose(transpose(basis)#transpose(gvel))

      ejecthistory[eject_step_ind[eject_ind]].x = reform(gpos[*,0])
      ejecthistory[eject_step_ind[eject_ind]].y = reform(gpos[*,1])
      ejecthistory[eject_step_ind[eject_ind]].z = reform(gpos[*,2])
      ejecthistory[eject_step_ind[eject_ind]].vx = reform(gvel[*,0])
      ejecthistory[eject_step_ind[eject_ind]].vy = reform(gvel[*,1])
      ejecthistory[eject_step_ind[eject_ind]].vz = reform(gvel[*,2])

      ejecthistory[eject_step_ind[eject_ind]].mass = gashistory[gashistory_ind].mass
      ejecthistory[eject_step_ind[eject_ind]].rho  = gashistory[gashistory_ind].rho
      ejecthistory[eject_step_ind[eject_ind]].temp = gashistory[gashistory_ind].temp
      ejecthistory[eject_step_ind[eject_ind]].metallicity = gashistory[gashistory_ind].metallicity
      ejecthistory[eject_step_ind[eject_ind]].haloid = gashistory[gashistory_ind].haloid
 ;     histogramp,sqrt(ejecthistory[eject_step_ind[eject_ind]].x*ejecthistory[eject_step_ind[eject_ind]].x + ejecthistory[eject_step_ind[eject_ind]].y*ejecthistory[eject_step_ind[eject_ind]].y),nbins = 100,xtitle = 'R'

      IF keyword_set(plot) THEN BEGIN
         v = sqrt(ejecthistory[eject_step_ind[eject_ind]].vz^2 + ejecthistory[eject_step_ind[eject_ind]].vx^2 +ejecthistory[eject_step_ind[eject_ind]].vy^2) 
         histogramp,v,nbins = 100,min = -100,max = 100,xrange = [-100,100]
         histogramp,ejecthistory[eject_step_ind[eject_ind]].vz,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 30,/overplot
         histogramp,ejecthistory[eject_step_ind[eject_ind]].vx,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 100,/overplot
         histogramp,ejecthistory[eject_step_ind[eject_ind]].vy,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 80,/overplot
      ENDIF
 ENDIF
;   IF expell_step_ind[0] NE -1 THEN BEGIN
;      match,expelliord[expell_step_ind],gashistory.iord,expell_ind,gashistory_ind
;      expell_ind = expell_step_ind[expell_ind]
;      expellhistory[expell_ind] = gashistory[gashistory_ind]
;      gpos = [[expellhistory[expell_ind].x - halodat[i].xc],[expellhistory[expell_ind].y - halodat[i].yc],[expellhistory[expell_ind].z - halodat[i].zc]]*scale
;      gvel = [[expellhistory[expell_ind].vx],[expellhistory[expell_ind].vy],[expellhistory[expell_ind].vz]]*units.vunit*scale
;      gpos = transpose(transpose(basis)#transpose(gpos))
;      gvel = transpose(transpose(basis)#transpose(gvel))
;      expellhistory[expell_ind].x = reform(gpos[*,0])
;      expellhistory[expell_ind].y = reform(gpos[*,1])
;      expellhistory[expell_ind].z = reform(gpos[*,2])
;      expellhistory[expell_ind].vx = reform(gvel[*,0])
;      expellhistory[expell_ind].vy = reform(gvel[*,1])
;      expellhistory[expell_ind].vz = reform(gvel[*,2])
;      histogramp,sqrt(expellhistory[expell_ind].x*expellhistory[expell_ind].x + expellhistory[expell_ind].y*expellhistory[expell_ind].y),nbins = 100,color = 240,/overplot
;      histogramp,expellhistory[expell_ind].z,nbins= 100
;      histogramp,expellhistory[expell_ind].vx,nbins = 100
;      stop
;  ENDIF
;  stop
ENDFOR
r = sqrt(ejecthistory.x^2 +ejecthistory.y^2+ejecthistory.z^2) 
histogramp,r,nbins = 100,min = -100,max = 100,xrange = [-100,100],xtitle = 'R'
;histogramp,r[eject_expell_i],nbins = 100,min = -100,max = 100,xrange = [-100,100],/overplot,linestyle = 2
histogramp,ejecthistory.z,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 30,/overplot
;histogramp,ejecthistory[eject_expell_i].z,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 30,/overplot,linestyle = 2
histogramp,ejecthistory.x,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 100,/overplot
;histogramp,ejecthistory[eject_expell_i].x,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 100,/overplot,linestyle = 2
histogramp,ejecthistory.y,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 80,/overplot
;histogramp,ejecthistory[eject_expell_i].y,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 80,/overplot,linestyle = 2
histogramp,r,nbins = 100,min = 0,max = 50,xrange = [0,50],xtitle = 'R'
;histogramp,r[eject_expell_i],nbins = 100,min = 0,max = 50,/overplot,linestyle = 2
v = sqrt(ejecthistory.vz^2 + ejecthistory.vx^2 +ejecthistory.vy^2) 
histogramp,v,nbins = 100,min = -100,max = 100,xrange = [-100,100],xtitle = 'V'
;histogramp,v[eject_expell_i],nbins = 100,min = -100,max = 100,xrange = [-100,100],/overplot,linestyle = 2
histogramp,ejecthistory.vz,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 30,/overplot
;histogramp,ejecthistory[eject_expell_i].vz,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 30,/overplot,linestyle = 2
histogramp,ejecthistory.vx,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 100,/overplot
;histogramp,ejecthistory[eject_expell_i].vx,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 100,/overplot,linestyle = 2
histogramp,ejecthistory.vy,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 80,/overplot
;histogramp,ejecthistory[eject_expell_i].vy,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 80,/overplot,linestyle = 2
;histogramp,v,nbins = 100,min = 0,max = 200,xrange = [0,200],xtitle = 'V'
;histogramp,v[eject_expell_i],nbins = 100,min = 0,max = 200,xrange = [0,200],/overplot,linestyle = 2

;mwrfits,expellhistory,'grp' + finalid + '.expell_disk.fits',/create
IF keyword_set(postdisk) THEN mwrfits,ejecthistory,'grp' + finalid + '.' + ejectext + '_halo.fits',/create ELSE mwrfits,ejecthistory,'grp' + finalid + '.' + ejectext + '_disk.fits',/create
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
    IF keyword_set(plot) THEN BEGIN
       IF i EQ 0 THEN histogramp,sqrt(eject_gas.x*eject_gas.x + eject_gas.y*eject_gas.y),xrange = [0,30],/normalize,yrange = [0,0.1],xtitle = 'Radius [kpc]',ytitle = '1/N dN/dr'
       histogramp,sqrt(eject_gas.x*eject_gas.x + eject_gas.y*eject_gas.y),/overplot,/normalize,color = colors[i],thick = thicks[i]
    ENDIF
ENDFOR

IF keyword_set(outplot) THEN device,/close

END


PRO reeject_gas_master

;dir = '/nobackupp2/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
;finalid = '1'
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid,/postdisk
;finalid = '2'
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid,/postdisk

;dir = '/nobackupp2/crchrist/MolecH/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
;finalid = '1'
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid,/postdisk
;finalid = '2'
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid,/postdisk
;finalid = '3'
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid,/postdisk

;dir = '/nobackupp2/crchrist/MolecH/h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
;finalid = '1'
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid,/postdisk
;finalid = '2'
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid
;reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid,/postdisk

dir = '/nobackupp2/crchrist/MolecH/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
finalid = '1'
reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid
reeject_gas_character,dir,ejecthistory,expellhistory,finalid = finalid,/postdisk

END
