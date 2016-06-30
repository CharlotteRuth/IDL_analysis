pro find_tvir,filename,firststep,dir,WMAP3=wmap3
;# just input your favourite galaxy: everything else shoule be aok!
;# "MW1.1024g1bwK"," ", " "
;#

;phase = mrdfits('grp1.allgas.entropy.fits',1)
phase = mrdfits('grp1.allgas.cut.entropy.fits',1)
smaccr = mrdfits('smooth.accr.iord.fits',0)
;phase = mrdfits('/astro/net/scratch1/abrooks/GasAccr/h579/grp1.allgas.entropy.fits',1)

;Concentrate only on smoothly accreted gas, to avoid complications with SF in 
;accreted halos
ind = findex(phase.iord, smaccr)
smphase = phase[ind]

; Get a at each timestep to calculate min SF density
;prefix = '/astro/net/scratch2/fabio/'+filename+'/'+filename
prefix = '/astro/net/scratch1/abrooks/FABIO/'+filename+'/'+filename
;prefix = '/astro/net/scratch2/fabio/cosmo25runs/'+filename+'/'+filename
alltracefile = file_search(prefix+'.0*/*.allgas.cut.history.fits')
liststeps = file_basename(alltracefile,'.allgas.cut.history.fits')
allsteps = strmid(liststeps,4,5,/reverse_offset)
strsteps = allsteps(where(fix(allsteps) ge firststep))
nsteps = n_elements(strsteps)
list = prefix+'.'+strsteps+'/'+liststeps(where(fix(allsteps) ge firststep))
gtp_file = list+'.amiga.gtp'

test = n_elements(smphase[0].grp)
if test ne nsteps then stop

; Assume there's an entropy floor set by the UV background temp and the 
; omega_baryon.  When gas shocks, it jumps in density by at least a factor 
; of 4, and to T_vir.  
;omega_m = rho_m(t)/rho_c(t)

IF keyword_set(wmap3) then begin
 om_m = 0.24  ;for this cosmology
 ;readcol, 'z.tgas2.dat', z1, t_uv
 z1 = [5.89,4,3,2,1,0.5,0]
 t_uv = [18781.3,24207.5,25697.6,24083.4,19720.4,16283.6,11554.3]
 readcol, '/astro/net/scratch1/abrooks/GasAccr/'+dir+'/mvir.dat', z2,mvir,r_vir,t_vir, /silent
;readcol, '/astro/net/scratch2/fabio/hz3.cosmo50cmb.2048g2bwK/alyson/mvir.dat', z2,mvir,rvir,t_vir, /silent
 a = fltarr(nsteps)
 FOR j=0,nsteps-1 do begin
   rtipsy, gtp_file[j], h,g,d,s
   a[j] = h.time
 ENDFOR
 z = (1./a)-1.
 linterp,z2,t_vir,z,tvir ;This interpolates tvir
 linterp,z2,r_vir,z,rvir ;This interpolates tvir
 linterp,z1,t_uv,z,tuv ;This interpolates mean entropy
 deltat = tvir-tuv
ENDIF ELSE BEGIN
 om_m = 0.3  ;for this cosmology
 z1 = [5.89,4,3,2,1,0.5,0]
 readcol, '/astro/net/scratch1/abrooks/GasAccr/'+dir+'/mvir.dat', z2,mvir,r_vir,t_vir, /silent
 t_uv = [18781.3,24207.5,25697.6,24083.4,19720.4,16283.6,11554.3]
 a = fltarr(nsteps)
 FOR j=0,nsteps-1 do begin
   rtipsy, gtp_file[j], h,g,d,s
   a[j] = h.time
 ENDFOR
 z = (1./a)-1.
 linterp,z2,t_vir,z,tvir ;This interpolates tvir
 linterp,z2,r_vir,z,rvir ;This interpolates tvir
 linterp,z1,t_uv,z,tuv  ;This interpolates mean entropy
 deltat = tvir-tuv
ENDELSE
plot, z, tvir, /ylog
oplot, z, tuv, linestyle=2;, /ylog
legend, linestyle=[0,2], ['tvir','tuv'], /top, /right

testmin = min(where(deltat gt 0))
if testmin lt 0 then begin
  print, 'Shocked entropy is always less than the ambient entropy.  Stopping'
  stop
endif
print, testmin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
shocked=-1
timestep=-1
radius=-1
allshock=-1
allradii=-1
alltime=-1
ngas = n_elements(smphase.iord)
accr=lonarr(ngas)+999
lowr=fltarr(ngas)
old=lonarr(ngas)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Have decided to follow particles until they reach 30 kpc 
; from the galaxy center, or rvir (whichever is smaller).
; First, find the step where the particles reach 30 kpc. If rvir 
; at that step is smaller than 30 kpc, then for those particles 
; find where particle enters rvir instead.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FOR j=0L,ngas-1 do lowr[j] = min(smphase[j].radius(where(smphase[j].radius gt 0.0)))
ind = where(lowr le 30.0, nind, comp=ind2, ncomp=nind2)
for k=0L,nind-1 do accr[ind(k)] = min(where(smphase[ind(k)].radius lt 30. and smphase[ind(k)].radius ge 0.))
FOR j=0L,nind2-1 do accr[ind2(j)] = where(smphase[ind2(j)].radius eq lowr[ind2(j)])   
change = where(rvir lt 30.0)
first = max(change)+1
ind = where(accr le first,nind)
for k=0L,nind-1 do accr[ind(k)] = min(where(smphase[ind(k)].grp eq 1))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Do all z at once
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
late = where(accr gt testmin, nlate)
IF nlate gt 0 then begin
FOR j=0L,n_elements(late)-1 do begin
  test = where(smphase[late(j)].temp[testmin+1:accr(late[j])] gt 3.*tvir[testmin+1:accr(late[j])]/8.,ntest)
  if ntest gt 0 then begin
    allshock = [allshock, late(j)]
    rtest = smphase[late(j)].radius[testmin+1:accr(late[j])]
    in = where(rtest-rvir[testmin+1:accr(late[j])] lt 0, nin)
    keep = lonarr(n_elements(test))
    for k=0L,n_elements(test)-1 do keep[k] = where(test[k] eq in)
    n = where(keep ne -1, nn)
    if nn ne 0 then begin 
;   if nin gt 0 then begin
;     keep = intersect(in,test)
;     prior = min(keep)-1
;     sudden = where(test eq prior)
;     if sudden eq -1 and prior ne -1 then begin
        keep = keep[n]
        shocked = [shocked, late(j)]
        timestep = [timestep, min(keep)+testmin+1]
        radius = [radius,smphase[late(j)].radius[min(keep)+testmin+1]]
;     endif
    endif
  endif
ENDFOR
ENDIF
    
;Eliminate any double entries that may occur in the high z and low z searches
shocked = shocked[1:n_elements(shocked)-1]
timestep = timestep[1:n_elements(timestep)-1]
radius = radius[1:n_elements(radius)-1]
allshock = allshock[1:n_elements(allshock)-1]
accrvir = accr[shocked]
accrall = accr[allshock]
print, 'Number of gas particles traced:  ', ngas
print, 'Number of particles shocked within Rvir:  ', n_elements(shocked)
print, 'Number of particles shocked anywhere:  ', n_elements(allshock)

;Sort them (just in case - they need to be sorted for other programs)
sorted = sort(shocked)
siord = shocked[sorted]
saccr = accrvir[sorted]
sstep = timestep[sorted]
srad = radius[sorted]
iords = smphase[siord].iord
mwrfits, iords, 'smooth.38tvir.iord.rvir.fits', /create
mwrfits, z[sstep], 'smooth.38tvir.zshock.rvir.fits', /create
mwrfits, sstep, 'smooth.38tvir.shockstep.rvir.fits', /create
mwrfits, saccr, 'smooth.38tvir.tracestep.rvir.fits', /create
mwrfits, srad, 'smooth.38tvir.radius.rvir.fits', /create

sorted = sort(allshock)
siord = allshock[sorted]
saccr = accrall[sorted]
;sstep = alltime[sorted]
;srad = allradii[sorted]
iords = smphase[siord].iord
;mwrfits, iords, 'smooth.tvir.iord.all.fits', /create
;mwrfits, saccr, 'smooth.tvir.tracestep.all.fits', /create
;mwrfits, z[sstep], 'smooth.shocked.zshock.all.fits', /create
;mwrfits, srad, 'smooth.shocked.radius.all.fits', /create

stop
end

