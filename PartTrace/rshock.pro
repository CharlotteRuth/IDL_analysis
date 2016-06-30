pro rshock,filename,firststep,dir,WMAP3=wmap3,h=h,g=g,d=d,s=s
;# just input your favourite galaxy: everything else shoule be aok!
;# "MW1.1024g1bwK"," ", " "
;#

;lengthunit = 20000./0.7
lengthunit = 50000.
shocked = mrdfits('/astro/net/scratch1/abrooks/GasAccr/'+dir+'/smooth.shocked.iord.fits',0)
accr = mrdfits('/astro/net/scratch1/abrooks/GasAccr/'+dir+'/smooth.shocked.accrtime.fits',0)
step = mrdfits('/astro/net/scratch1/abrooks/GasAccr/'+dir+'/smooth.shocked.timestep.fits',0)
test = accr-step
shocked = shocked(where(test gt 2))
;lengthunit = 1.e5;/0.7
;shocked = mrdfits('/astro/net/scratch1/abrooks/GasAccr/gal1/smooth.shocked.iord.fits',0)
nshock = n_elements(shocked)

; Get a at each timestep to calculate min SF density
;prefix = '/astro/net/scratch1/abrooks/FABIO/'+filename+'/'+filename
prefix = '/astro/net/scratch2/fabio/REPOSITORY/e12Gals/'+filename+'/'+filename
alltracefile = file_search(prefix+'.0*/*.allgas.history.fits')
liststeps = file_basename(alltracefile,'.allgas.history.fits')
allsteps = strmid(liststeps,4,5,/reverse_offset)
strsteps = allsteps(where(fix(allsteps) ge firststep))
nsteps = n_elements(strsteps)
list = prefix+'.'+strsteps+'/'+liststeps(where(fix(allsteps) ge firststep))
historyfiles = list+'.allgas.history.fits'
gtpfiles = list+'.amiga.gtp'

history = mrdfits(historyfiles[0],1)
rtipsy, gtpfiles[0], h,g,d,s
x = history.x-s[0].x*lengthunit
y = history.y-s[0].y*lengthunit
z = history.z-s[0].z*lengthunit
ind = binfind(history.iord, shocked)
keep = history[ind].iord
;radius = (history[ind].x^2.+history[ind].y^2.+history[ind].z^2.)^0.5
radius = (x[ind]^2.+y[ind]^2.+z[ind]^2.)^0.5
struct = replicate({iord:0L, shell:0L, radius:fltarr(nsteps)}, nshock)
shellbounds = indgen(101)*20.+200.
struct.radius[0] = radius
;del = where(history[ind].x eq 0.0)
;if del[0] ne -1 then struct[ind(del)].radius[0] = 0.0
struct.iord = keep
for i=0,99 do begin
  ind2 = where(radius ge shellbounds[i] and radius lt shellbounds[i+1])
  struct[ind2].shell = i+1
endfor

rvir=fltarr(nsteps)
redshift=fltarr(nsteps)
lb=fltarr(nsteps)
rvir[0] = s[0].eps*lengthunit;/h.time
redshift[0] = (1./h.time)-1.
lb[0] = lookback(redshift[0])

for i=1,nsteps-1 do begin
  history = mrdfits(historyfiles[i],1)
  rtipsy, gtpfiles[i], h,g,d,s
  x = history.x-(s[0].x*lengthunit)
  y = history.y-(s[0].y*lengthunit)
  z = history.z-(s[0].z*lengthunit)
  radius = (x[ind]^2.+y[ind]^2.+z[ind]^2.)^0.5
  struct.radius[i] = radius
  del = where(history[ind].x eq 0.0)
  if del[0] ne -1 then struct[del].radius[i] = 0.0
  rvir[i] = s[0].eps*lengthunit;/h.time
  redshift[i] = (1./h.time)-1.
  lb[i] = lookback(redshift[i]) 
endfor 

;lb = lookback(z)
struct2 = replicate({shell:0L, meanr:fltarr(nsteps)}, 101)
struct2.shell = indgen(101)
for i=0,99 do begin
  ind = where(struct.shell eq i+1)
  for j=0, nsteps-1 do begin
    good = where(struct[ind].radius[j] ne 0.0)
    struct2[i].meanr[j] = mean(struct[ind(good)].radius[j])
  endfor
endfor

lb = lb/1.e9
;x = indgen(nsteps)
plot, 13.7-lb, struct2[45].meanr, thick=2, xthick=2, ythick=2, title=filename, ytitle='Comoving kpc', charthick=2, xtitle = 'Lookback Time (Gyr)', yrange=[0,400]
oplot, 13.7-lb, struct2[35].meanr, thick=2
oplot, 13.7-lb, struct2[25].meanr, thick=2
oplot, 13.7-lb, struct2[15].meanr, thick=2
oplot, 13.7-lb, struct2[5].meanr, thick=2
oplot, 13.7-lb, struct2[0].meanr, thick=2
oplot, 13.7-lb, rvir, linestyle=2, thick=2

;plot, lb, struct2[45].meanr, thick=2, xthick=2, ythick=2, title=filename, ytitle='Comoving kpc', charthick=2, xtitle = 'Lookback Time (Gyr)', yrange=[0,400], xrange=[
;oplot, lb, struct2[35].meanr, thick=2
;oplot, lb, struct2[25].meanr, thick=2
;oplot, lb, struct2[15].meanr, thick=2
;oplot, lb, struct2[5].meanr, thick=2
;oplot, lb, struct2[0].meanr, thick=2
;oplot, lb, rvir, linestyle=2, thick=2
stop 

end

