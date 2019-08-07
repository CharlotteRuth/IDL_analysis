;# just input your favourite galaxy: everything else should be aok
;# lunit is in Mpc
;#

pro accrmode, filebase, finalid = finalid, skip=skip
; lunit, wmap1=wmap1,
;lunit in Mpc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This code creates a structure that finds the maximum temperature of a 
; gas particle in its history.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;if keyword_set(wmap1) then densunit = 135.98 else densunit = 147.8344
;lunit = lunit*1000.
;units = tipsyunits(pfile)
IF NOT keyword_set(finalid) THEN finalid = '1'

readcol, filebase + '.grp' + finalid + '.haloid.dat', files, halo, format='a,l'
;Output of haloid runs from low z to high z, but find_shocked requires high to low.  Fix here.
files = reverse(files)
halo = reverse(halo)
list = files+'.grp' + finalid + '.allgas.history.fits'
gtp_file = files+'.amiga.gtp'
statfile = files+'.amiga.stat'
iord_files = files+'.iord'

if keyword_set(skip) THEN GOTO, jump1
; ****************************************************
; define some stuff
;*******************************************************

alliords = mrdfits('grp' + finalid + '.allgas.iord.fits',0)
nall = n_elements(alliords)
nsteps = n_elements(halo)
readcol, statfile[0], grp, mvir, rvir0, format='l,x,x,x,x,f,f', /silent
grpind = where(grp EQ halo[0])
allstruct = replicate({iord:0L, mass:LONARR(nsteps), grp:LONARR(nsteps), temp:FLTARR(nsteps), rho:FLTARR(nsteps), entropy:FLTARR(nsteps), radius:FLTARR(nsteps), vx:FLTARR(nsteps), vy:FLTARR(nsteps), vz:FLTARR(nsteps), x:FLTARR(nsteps), y:FLTARR(nsteps), z:FLTARR(nsteps)}, nall)
allstruct.iord = alliords
alliords = 0
history = mrdfits(list[0],1)
rtipsy, gtp_file[0], h,g,d,s
allstruct.temp[0] = history.temp
allstruct.mass[0] = history.mass
;allstruct.rho[0] = history.rho*densunit/(h.time^3.)
allstruct.rho[0] = history.rho
;allstruct.entropy[0] = alog10(history.temp^1.5/(history.rho*unit.rhounit/h.time^3.))
allstruct.entropy[0] = alog10(history.temp^1.5/(history.rho))
allstruct.grp[0] = history.haloid
;allstruct.radius[0] = ((h.time*(history.x-units.lengthunit*s[0].x))^2.+(h.time*(history.y-units.lengthunit*s[0].y))^2.+(h.time*(history.z-units.lengthunit*s[0].z))^2.)^0.5 
;allstruct.radius[0] = ((h.time*(history.x-units.lengthunit*s[grpind].x))^2.+(h.time*(history.y-units.lengthunit*s[grpind].y))^2.+(h.time*(history.z-units.lengthunit*s[grpind].z))^2.)^0.5 
allstruct.radius[0] = ((history.x-s[grpind].x)^2.+(history.y-s[grpind].y)^2.+(history.z-s[grpind].z)^2.)^0.5 
allstruct.vx[0] = history.vx
allstruct.vy[0] = history.vy
allstruct.vz[0] = history.vz
allstruct.x[0] = history.x
allstruct.y[0] = history.y
allstruct.z[0] = history.z

FOR j=1L,nsteps-1 DO BEGIN
   history = mrdfits(list[j],1)
   rtipsy, gtp_file[j], h,g,d,s
   readcol, statfile[j], grp, mvir, rvir0, format='l,x,x,x,x,f,f', /silent
   grpind = where(grp EQ halo[j])
   ;iord = read_lon_array(iord_files[j])
   ;inds = binfind(iord, allstruct.iord)
   ;exist = where(inds NE -1, comp=del)
   exist = where(history.mass NE 0, comp=del)
   ninds = n_elements(exist)
   allstruct[exist].temp[j] = history[exist].temp
   allstruct[exist].mass[j] = history[exist].mass
   allstruct[exist].vx[j] = history[exist].vx
   allstruct[exist].vy[j] = history[exist].vy
   allstruct[exist].vz[j] = history[exist].vz
   allstruct[exist].x[j] = history[exist].x
   allstruct[exist].y[j] = history[exist].y
   allstruct[exist].z[j] = history[exist].z
   ;allstruct[exist].rho[j] = history[exist].rho*unit.rhounit/(h.time^3.)
   allstruct[exist].rho[j] = history[exist].rho
   ;allstruct[exist].entropy[j] = alog10(history[exist].temp^1.5/(history[exist].rho/h.time^3.))
   allstruct[exist].entropy[j] = alog10(history[exist].temp^1.5/(history[exist].rho))
   allstruct[exist].grp[j] = history[exist].haloid
   ;allstruct[exist].radius[j] = ((h.time*(history[exist].x-s[grpind].x))^2.+(h.time*(history[exist].y-s[grpind].y))^2.+(h.time*(history[exist].z-s[grpind].z))^2.)^0.5 
   allstruct[exist].radius[j] = ((history[exist].x-s[grpind].x)^2.+ (history[exist].y-s[grpind].y)^2.+(history[exist].z-s[grpind].z)^2.)^0.5 

   IF del[0] NE -1 THEN BEGIN
      allstruct[del].temp[j] = -1
      allstruct[del].mass[j] = -1
      allstruct[del].rho[j] = -1
      allstruct[del].entropy[j] = -1
      allstruct[del].grp[j] = -1
      allstruct[del].radius[j] = -1
      allstruct[del].vx[j] = -1
      allstruct[del].vy[j] = -1
      allstruct[del].vz[j] = -1
      allstruct[del].x[j] = -1
      allstruct[del].y[j] = -1
      allstruct[del].z[j] = -1
   ENDIF
ENDFOR

;Write for testing purposes.  But this file is huge.  Delete if possible.
mwrfits, allstruct, 'grp' + finalid + '.allgas.entropy.fits', /create

;Pass structure to find_smooth, to find eary and smoothly accreted gas particles
jump1: allstruct = mrdfits('grp' + finalid + '.allgas.entropy.fits',1)

;early = where(allstruct.grp[0] eq halo[0], comp=free)
;allstruct = allstruct[free]
;Get info on accretion time, mass at accretion 
;;accretion, filebase, finalid, allstruct


IF 0 THEN BEGIN
;I'm waiting on this until I've figured out the density units.  CC 12/10/12
;Now sort the results to get clumpy and unshocked
   haloidoutfile = filebase + '.grp' + finalid + '.haloid.dat'
   all = mrdfits('grp' + finalid + '.allgas.iord.fits',0)
   early = mrdfits('grp' + finalid + '.early.iord.fits',0)

   test = binfind(early,all)
   late = all(where(test eq -1))
   find_smooth, allstruct, halo
   sm = mrdfits('grp' + finalid + '.smooth.accr.iord.fits',0)
   test = binfind(sm,late)
   clumpy = late(where(test eq -1))
   mwrfits, clumpy, 'grp' + finalid + '.clumpy.accr.iord.fits', /create
;Now pass the structure to find_shocked, to find the shocked smooth particles 
   if keyword_set(wmap1) then find_shocked, allstruct, haloidoutfile, /wmap1 else $
      find_shocked, allstruct, haloidoutfile, finalid
   shock = mrdfits('grp' + finalid + '.shocked.iord.fits',0)
   test = binfind(shock,sm)
   unshock = sm(where(test eq -1))
   mwrfits, unshock, 'grp' + finalid + '.unshock.iord.fits', /create
ENDIF

end

