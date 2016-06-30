pro find_sfz,filename,firststep,WMAP3=wmap3,h=h,g=g,d=d,s=s
;# just input your favourite galaxy: everything else shoule be aok!
;# "MW1.1024g1bwK"," ", " "
;#

;This file find the height, z, above the galaxy plane when a gas particle 
;first becomes dense enough to form a star.
;phase = mrdfits('/astro/net/scratch1/abrooks/GasAccr/MW1/grp1.allgas.entropy.fits',1)
iords = mrdfits('/astro/net/scratch1/abrooks/FABIO/MW1.1024g1bwK/alyson/iord.adi.halo.gasnoprog.fits',0)

; Get a at each timestep to calculate min SF density
prefix = '/astro/net/scratch1/abrooks/FABIO/'+filename+'/'+filename
alltracefile = file_search(prefix+'.0*/*.allgas.history.fits')
liststeps = file_basename(alltracefile,'.allgas.history.fits')
allsteps = strmid(liststeps,4,5,/reverse_offset)
strsteps = allsteps(where(fix(allsteps) ge firststep))
nsteps = n_elements(strsteps)
list = prefix+'.'+strsteps+'/'+liststeps(where(fix(allsteps) ge firststep))
gtp_file = list+'.amiga.gtp'
iord_file = list+'.iord'
rot_file = list+'.999.00.01.rot'
;mu = 0.6  ; molecular weight of gas, ~0.6 for primordial gas
mh = 1.673534d-24 ;hydrogen mass, in grams

; SF min density must be both > 2.63*rhoc, and > PhysDenMin = 0.1 m_h/cc
; At redshifts back to at least z=4 (and probably much higher), PhysDenMin
; is always greater than 2.63*rhoc.

a = fltarr(nsteps)
z = fltarr(nsteps)
zh = fltarr(nsteps)
xh = fltarr(nsteps)
yh = fltarr(nsteps)
vzh = fltarr(nsteps)
vxh = fltarr(nsteps)
vyh = fltarr(nsteps)
lb = fltarr(nsteps)
FOR j=0,nsteps-1 do begin
   rtipsy, gtp_file[j], h,g,d,s
   a[j] = h.time
   z[j] = (1./a[j])-1.
   lb[j] = lookback(z[j])
   zh[j] = s[0].z
   xh[j] = s[0].x
   yh[j] = s[0].y
   vzh[j] = s[0].vz
   vxh[j] = s[0].vx
   vyh[j] = s[0].vy
ENDFOR
afterbb = (13.7e9-lb)/13.7e9
simtime = afterbb*0.3331

rtipsy, list[nsteps-1], h,g,d,s
iord = read_ascii_array(iord_file[nsteps-1])
iord = iord[h.ngas+h.ndark:h.n-1]
ind = findex(iord,iords)
tform = s[ind].tform
;az6 = where(tform gt min(simtime))
;iords = iords[az6]
;tform = tform[az6]
;ind = ind[az6]
;closest = findex(simtime,tform)
;ind2 = ceil(closest)

;rhoc = rho_crit(z)
;overden = 2.63*rhoc
;physmin = 0.1*mh*(3.08568d21)^3./2.d33

;mu = 0.6  ; molecular weight of gas, ~0.6 for primordial gas
;mh = 1.673534d-24 ;hydrogen mass, in grams
;physmin = 0.1*mh*(3.08568d21)^3./2.d33
;stop

; ****************************************************
; Find the timesteps this where particles first become dense enough for SF
;*******************************************************
;
;timestep=-1
;smooth_accr = findex(phase.iord, iords)
;smaccr = phase[smooth_accr]
;ngas = n_elements(iords)
;FOR j=0L,ngas-1 do begin
; presf = min(where(smaccr[j].rho gt physmin)) 
; timestep = [timestep, presf]
;ENDFOR
;timestep = timestep[1:n_elements(timestep)-1]

rtipsy, list[0], h,g,d,s
first = max(s.tform)
;Want the z height above the disk plane when this happens
zpart=fltarr(n_elements(ind))
radf=fltarr(n_elements(ind))
vradf=fltarr(n_elements(ind))
FOR j=1L,nsteps-1 do begin
  rtipsy, list[j], h,g,d,s
  ;rtipsy, rot_file[j], h,g,d,s
  new = max(s.tform)
  exist = where(tform gt first and tform le new)
  if exist[0] eq -1 then continue
  xvals = s[ind(exist)].x-xh[j]
  yvals = s[ind(exist)].y-yh[j]
  zvals = s[ind(exist)].z-zh[j]
  vxvals = s[ind(exist)].vx-vxh[j]
  vyvals = s[ind(exist)].vy-vyh[j]
  vzvals = s[ind(exist)].vz-vzh[j]
  rad = (xvals^2.+yvals^2.+zvals^2.)^0.5
  vrad = (vxvals^2.+vyvals^2.+vzvals^2.)^0.5
  radf[exist] = rad*a[j]
  vradf[exist] = vrad*a[j]
  zpart[exist] = abs(zvals)*a[j]
  first = new
ENDFOR

mwrfits, zpart, 'adi.halo.sf.zheight.fits', /create
mwrfits, radf, 'adi.halo.sf.rad.fits', /create
mwrfits, vradf, 'adi.halo.sf.vrad.fits', /create
mwrfits, tform, 'adi.halo.sf.tform.fits', /create

stop
end

