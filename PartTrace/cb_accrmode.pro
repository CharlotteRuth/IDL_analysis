;# just input your favourite galaxy: everything else shoule be aok!
;# "MW1.1024g1bwK"," ", " "
;#
pro accrmode,filename,firststep,h=h,g=g,d=d,s=s

prefix = '/astro/net/mega-1/abrooks/FABIO/'+filename+'/'+filename

alltracefile = file_search(prefix+'.0*/*.star.history.fits')
liststeps = file_basename(alltracefile,'.star.history.fits')
allsteps = strmid(liststeps,4,5,/reverse_offset)
strsteps = allsteps(where(fix(allsteps) ge firststep))
nsteps = n_elements(strsteps)
list = prefix+'.'+strsteps+'/'+liststeps(where(fix(allsteps) ge firststep))
tracefile = alltracefile(where(fix(allsteps) ge firststep))
ending= strmid(liststeps[0],5,6,/reverse_offset)
base = file_basename(liststeps[0],ending)
finalfile = prefix+'.00512/'+base+'.00512'
cmp_file = finalfile+'.cmp'

if(keyword_set(h) EQ 0 AND keyword_set(s) EQ 0) then rtipsy, finalfile, h,g,d,s

;**************************
; get steps info
;*************************

steptime=fltarr(nsteps+1)


; ****************************************************
; define some stuff
;*******************************************************

tempdata = read_ascii_array(cmp_file)
components = tempdata[h.ngas+h.ndark:h.ndark+h.nstar+h.ngas-1]

nstar = n_elements(s)
star_haloinf = replicate({haloid:intarr(nstar),Mvir:fltarr(nstar),maxhalo:fltarr(nstar),GMass:fltarr(nstar),SMass:fltarr(nstar),accrtime1:fltarr(nstar),accrtime2:fltarr(nstar),accrtime3:fltarr(nstar),delay:fltarr(nstar),accrmode:intarr(nstar),comp:intarr(nstar)},1)
star_halo = replicate({Mvir:fltarr(nstar)},nsteps)
comp = ['unbound','disk','halo','bulge','rest']
delay = fltarr(nstar)
loadct,39


;****************************************************************
;now loop over the components and find the accreted stars vs gas 
;****************************************************************

galaxy = where(components ge 1 and components le 4)


ncomp=n_elements(galaxy)
st = replicate({id:Lonarr(ncomp)},nsteps)

;********************************************************************
; do the 0th timestep first, to work out what is in sats at that time 
;********************************************************************

history = mrdfits(tracefile[0],1)
st[0].id = history.haloid
steptime[0]= fix(strsteps[0])/512.*0.3331
;free_lun,tfile

gal=1             ; the central galaxy is defined by having amiga halo 1
nonhalo0 = where(st[0].id le gal,comp =sats1)
accreted0 = where(st[0].id eq gal,comp=sats2)


;***********************************
; loop over rest of timesteps
; **********************************

for i =1,nsteps do begin

if i lt nsteps then begin
history = mrdfits(tracefile[i],1)
st[i].id = history.haloid
steptime[i]= fix(strsteps[i])/512.*0.3331

accr1 = where(st[i].id le gal,comp=tmp)

endif

if i eq nsteps then begin
steptime[nsteps] = 0.3331
accr1 = where(galaxy ge 0)
endif


accreted1 = intersect(sats1,accr1) 
sats1 = tmp

if max(accreted1) ge 0 then begin
star_haloinf.accrtime1[galaxy[accreted1]]=steptime[i-1]
star_haloinf.accrtime3[galaxy[accreted1]]=steptime[i]
endif

if i lt nsteps then begin
accr2 = where(st[i].id eq gal,comp=tmp)
endif
if i eq nsteps then begin
accr2 = where(galaxy ge 0)
endif

accreted2 = intersect(sats2,accr2) 
sats2 = tmp

if max(accreted2) ge 0  then begin
star_haloinf.accrtime2[galaxy[accreted2]]=steptime[i-1]
star_haloinf.accrtime3[galaxy[accreted2]]=steptime[i]
endif


statfile=list[i-1]+'.amiga.stat' 
readcol,statfile,haloid,Nt,Ng,Ns,Nd,Mvir,Rvir,Gmass,Smass,Dmass,Vmax,RVm,Vdisp,Xc,Yc,Zc,VXc,VYc,VZc,ID_A,FORMAT='I,I,I,I,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,I' 

if (max(accreted1 ge 0)) then begin
stat_row_number = st[i-1].id[accreted1]
loop = n_elements(accreted1)-1
for j = 0L,loop do begin
star_haloinf.haloid[galaxy[accreted1[j]]] = haloid[stat_row_number[j]-1]
star_haloinf.Mvir[galaxy[accreted1[j]]] = Mvir[stat_row_number[j]-1]
star_haloinf.Gmass[galaxy[accreted1[j]]] = Gmass[stat_row_number[j]-1]
star_haloinf.Smass[galaxy[accreted1[j]]] = Smass[stat_row_number[j]-1]
endfor
endif


;****************************************************
; find most massive halo that each star was ever in 
;*****************************************************

insats = where(st[i-1].id gt gal)
if (max(insats) ge 0) then begin
stat_row_number = st[i-1].id[insats]

for j = 0L,n_elements(insats)-1 do begin
star_halo[i-1].Mvir[galaxy[insats[j]]] = Mvir[stat_row_number[j]-1]
endfor
endif

endfor

for j = 0L,n_elements(galaxy)-1 do begin
star_haloinf.maxhalo[galaxy[j]] = max(star_halo.Mvir[galaxy[j]])
endfor

;***********************************************************************
;deteremine whether a star was gas/star at time of accretion
;******************************************************************

astars = where(s[galaxy].tform lt star_haloinf.accrtime2[galaxy] and st[0].id ne gal,comp=merge1)
agas = where(s[galaxy].tform gt star_haloinf.accrtime3[galaxy] and st[0].id ne gal,comp=merge2)
early = where(st[0].id eq gal,comp=late)
mergetmp = intersect(merge1,merge2)
merge = where(st[0].id[mergetmp] ne gal)
star_haloinf.delay[galaxy] = 13.7/0.3331*(s[galaxy].tform - star_haloinf.accrtime3[galaxy])

star_haloinf.accrmode[galaxy[early]] = 1
star_haloinf.accrmode[galaxy[astars]] = 2
star_haloinf.accrmode[galaxy[mergetmp[merge]]] = 3
star_haloinf.accrmode[galaxy[agas]] = 4

star_haloinf.comp = components

;******************************************************************
; write a binary file with all this info!
;*******************************************************************
stop
outfile=filename+'.accretion.FITS'
mwrfits,star_haloinf,outfile,/create

stop
end
