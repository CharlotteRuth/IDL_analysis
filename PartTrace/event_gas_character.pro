;CC, 08/11/16

;This program determines the characteristics (density, distance from
;center etc.) of the gas at the step
;it is reaccreted, each time it is accreted.  It is based on
;fb_gas_character, eject_gas_character, reeject_character, and
;reaccr_gas_character.
;Unlike reaccr_gas_charcter, the gas does not have to have been
;ejected previously; rather, it tallys the properties of the gas
;each time it is accreted/ejected/expelled etc

;event_ext may be:
;reaccrdiskall
;relost
;reaccr

PRO event_gas_character,history,event_ext,dir = dir,finalid = finalid,laststep = laststep,plot = plot

IF NOT keyword_set(finalid) THEN finalid = '1'
IF NOT keyword_set(laststep) THEN laststep = '512'
IF NOT keyword_set(dir) THEN dir='./'

cd,dir
spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord
spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file
spawn,'ls ' + dir + '/h*param',pfile

units = tipsyunits(pfile[0])
rtipsy,file,h,g,d,s,/justhead
readarr,file_iord,h,iord,part = 'gas',/ascii
timeend = 13.734377

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

;Reading in files
accriord = mrdfits(dir + '/grp' + finalid + '.' + event_ext + '_iord.fits',0)
accrmass = mrdfits(dir + '/grp' + finalid + '.mass_at_' + event_ext + '.fits',0)
accrz    = mrdfits(dir + '/grp' + finalid + '.' + event_ext + '_z.fits',0)

match2,fix(accrz*1000)/1000.0,fix(halodat.z*1000)/1000.0,accr_step,temp
accr_step = accr_step

basestruct = {iord:0L, igasorder:0L, mark:0L, mass:double(0.0), red:double(0.0), x:double(0.0), y:double(0.0), z:double(0.0), vx:double(0.0), vy:double(0.0), vz:double(0.0), rho:double(0.0), temp:double(0.0), metallicity:double(0.0), haloid:0L}
accrhistory  = replicate(basestruct,n_elements(accrz))
accrhistory.iord = accriord

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
   accr_step_ind  = where(accr_step  EQ i) ;The ejcted particles that are in the disk at step i
   IF accr_step_ind[0] NE -1 THEN BEGIN
      match,accriord[accr_step_ind],gashistory.iord,accr_ind,gashistory_ind
 ;     accr_ind = accr_step_ind[accr_ind]
      gpos = [[gashistory[gashistory_ind].x - halodat[i].xc*scale],[gashistory[gashistory_ind].y - halodat[i].yc*scale],[gashistory[gashistory_ind].z - halodat[i].zc*scale]]
      gvel = [[gashistory[gashistory_ind].vx - satsdata[isats].xvcm],[gashistory[gashistory_ind].vy - satsdata[isats].yvcm],[gashistory[gashistory_ind].vz - satsdata[isats].zvcm]]
      gpos = transpose(transpose(basis)#transpose(gpos))
      gvel = transpose(transpose(basis)#transpose(gvel))
      accrhistory[accr_step_ind[accr_ind]].x = reform(gpos[*,0])
      accrhistory[accr_step_ind[accr_ind]].y = reform(gpos[*,1])
      accrhistory[accr_step_ind[accr_ind]].z = reform(gpos[*,2])
      accrhistory[accr_step_ind[accr_ind]].vx = reform(gvel[*,0])
      accrhistory[accr_step_ind[accr_ind]].vy = reform(gvel[*,1])
      accrhistory[accr_step_ind[accr_ind]].vz = reform(gvel[*,2])
      accrhistory[accr_step_ind[accr_ind]].red = fltarr(n_elements(accr_ind)) + z
      
      accrhistory[accr_step_ind[accr_ind]].mass = gashistory[gashistory_ind].mass
      accrhistory[accr_step_ind[accr_ind]].rho  = gashistory[gashistory_ind].rho
      accrhistory[accr_step_ind[accr_ind]].temp = gashistory[gashistory_ind].temp
      accrhistory[accr_step_ind[accr_ind]].metallicity = gashistory[gashistory_ind].metallicity
      accrhistory[accr_step_ind[accr_ind]].haloid = gashistory[gashistory_ind].haloid
;      histogramp,sqrt(accrhistory[accr_step_ind[accr_ind]].x*accrhistory[accr_step_ind[accr_ind]].x + accrhistory[accr_step_ind[accr_ind]].y*accrhistory[accr_step_ind[accr_ind]].y),nbins = 100

      IF keyword_set(plot) THEN BEGIN
;         histogramp,accrhistory[accr_step_ind[accr_ind]].vz,nbins = 100,min = -100,max = 100,xrange = [-100,100]
;         histogramp,accrhistory[accr_step_ind[accr_ind]].vx,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 100,/overplot
;         histogramp,accrhistory[accr_step_ind[accr_ind]].vy,nbins = 100,min = -100,max = 100,xrange = [-100,100],color = 80,/overplot
         histogramp,accrhistory[accr_step_ind[accr_ind]].z,nbins = 50,min = -25,max = 25,xrange = [-25,25]
         histogramp,accrhistory[accr_step_ind[accr_ind]].x,nbins = 50,min = -25,max = 25,xrange = [-25,25],color = 100,/overplot
         histogramp,accrhistory[accr_step_ind[accr_ind]].y,nbins = 50,min = -25,max = 25,xrange = [-25,25],color = 80,/overplot
      ENDIF
;      if i EQ n_elements(halodat.file) - 1 THEN stop
   ENDIF
ENDFOR

mwrfits,accrhistory,'grp' + finalid + '.' + event_ext + '_history.fits',/create
END
