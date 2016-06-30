;# just input your favourite galaxy: everything else shoule be aok!
;# "MW1.1024g1bwK"," ", " "
;#
pro gasaccr,filename,firststep,h=h,g=g,d=d,s=s

prefix = '/astro/net/scratch1/abrooks/FABIO/'+filename+'/'+filename

alltracefile = file_search(prefix+'.0*/*.allgas.history.fits')
liststeps = file_basename(alltracefile,'.allgas.history.fits')
allsteps = strmid(liststeps,4,5,/reverse_offset)
strsteps = allsteps(where(fix(allsteps) ge firststep))
nsteps = n_elements(strsteps)
list = prefix+'.'+strsteps+'/'+liststeps(where(fix(allsteps) ge firststep))
tracefile = alltracefile(where(fix(allsteps) ge firststep))
ending= strmid(liststeps[0],5,6,/reverse_offset)
base = file_basename(liststeps[0],ending)
finalfile = prefix+'.00512/'+base+'.00512'
cmp_file = finalfile+'.cmp'
grp_file = finalfile+'.amiga.grp'
gtp_file = list+'.amiga.gtp'
iord_files = list+'.iord'

if(keyword_set(h) EQ 0 AND keyword_set(s) EQ 0) then rtipsy, finalfile, h,g,d,s

;**************************
; get steps info
;*************************

steptime=fltarr(nsteps+1)


; ****************************************************
; define some stuff
;*******************************************************

;phase contains information for each particle, and their place in the list stays constant 
;even after deletion  (but deleted particles should have grp 0)
phase = mrdfits('/astro/net/scratch1/abrooks/GasAccr/dwfh/grp1.allgas.entropy.fits',1)
nall = n_elements(phase.iord)

grp = read_ascii_array(grp_file)
grp = grp[0:h.ngas-1]

;ngas = h.ngas
star_haloinf = replicate({iord:lonarr(nall),haloid:intarr(nall),Mvir:fltarr(nall),maxhalo:fltarr(nall),GMass:fltarr(nall),SMass:fltarr(nall),accrtime1:fltarr(nall),accrtime2:fltarr(nall),accrtime3:fltarr(nall),comp:intarr(nall)},1)
star_halo = replicate({Mvir:fltarr(nall)},nsteps)
comp = ['unbound','disk','halo','bulge','rest']
loadct,39
star_haloinf.iord = phase.iord

;For all gas (including some that gets deleted by z=0), the order of the 
;particles will change due to deletion.  So let's create a structure first 
;that puts all the particles in the same place at each timestep

;allstruct = mrdfits('/astro/net/scratch1/abrooks/GasAccr/grp1.allgas.haloid.fits',1)
;allstruct = mrdfits('/astro/net/scratch1/abrooks/FABIO/MW1.1024g1bwK/grp1.allgas.haloid.fits',1)
;allstruct = replicate({iord:0L, haloid:LONARR(nsteps)}, nall)
;allstruct.iord = alliords
;history = mrdfits(tracefile[0],1)
;allstruct.haloid[0] = history.haloid
;FOR j=1L,nsteps-1 do begin
;  history = mrdfits(tracefile[j],1)
;  iord = read_ascii_array(iord_files[j])
;  inds = binfind(iord, alliords)
;  exist = where(inds NE -1, comp=del)
;  ninds = n_elements(exist)
;  print, n_elements(history.haloid)
;  print, n_elements(exist)
;  allstruct[exist].haloid[j] = history[0:ninds-1].haloid
;  allstruct[del].haloid[j] = -1
;ENDFOR

;tempdata = read_ascii_array(cmp_file)
;components = tempdata[where(grp eq 1)]
z0 = where(phase.grp[nsteps-1] EQ 1, comp=nongal)

;****************************************************************
;now loop over the components and find the accreted stars vs gas 
;****************************************************************

;********************************************************************
; do the 0th timestep first, to work out what is in sats at that time 
;********************************************************************
; These next steps depend on running from high z to low z
; For MW, assume this is timestep 35, and haloid is 1

rtipsy, gtp_file[0], h,g,d,s
z = (1./h.time)-1.
steptime[0]= lookback(z) ;[0])/512.*0.3331
gal=1             ; the central galaxy is defined by having amiga halo 1
nonhalo0 = where(phase.grp[0] le gal and phase.grp[0] ge 0,comp =sats1)
accreted0 = where(phase.grp[0] eq gal,comp=sats2)
;sats 1 is everyone grp 2 and up
;sats 2 is grp 0 and grp 2 up


;***********************************
; loop over rest of timesteps
; **********************************

for i=1,nsteps-1 do begin
if i eq 1 then gal=2 else gal=1
if i eq 2 then gal=3 else gal=1

rtipsy, gtp_file[i], h,g,d,s
z = (1./h.time)-1.
steptime[i]= lookback(z) ;fix(strsteps[i])/512.*0.3331
accr1 = where(phase.grp[i] eq gal and phase.grp[i] eq 0,comp=tmp) 
;accr1 is then all things in the main gal and in grp 0
;tmp is all other halos

accreted1 = intersect(sats1,accr1) ;find things that were in sats, now in 0 or 1
sats1 = tmp	;redefine the sats for the next timestep

;if max(accreted1) ge 0 then begin
if n_elements(accreted1) gt 1 then begin
;It was in a unique halo the timestep before this.  Make that accrtime1
star_haloinf.accrtime1[accreted1]=steptime[i-1]
endif

accr2 = where(phase.grp[i] eq gal,comp=tmp)
accreted2 = intersect(sats2,accr2) ;everything that was not in grp 1, but now is
sats2 = tmp

;if max(accreted2) ge 0  then begin
if n_elements(accreted2) gt 1  then begin
star_haloinf.accrtime2[accreted2]=steptime[i-1] 
;Could have been 0 or a halo, but is now part of grp 1
star_haloinf.accrtime3[accreted2]=steptime[i]
endif


statfile=list[i-1]+'.amiga.stat' 
readcol,statfile,haloid,Nt,Ng,Ns,Nd,Mvir,Rvir,Gmass,Smass,Dmass,Vmax,RVm,Vdisp,Xc,Yc,Zc,VXc,VYc,VZc,ID_A,FORMAT='I,I,I,I,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,I' 

;if max(accreted1) ge 0 then begin
if n_elements(accreted1) gt 1 then begin
stat_row_number = phase[accreted1].grp[i-1];.id[accreted1]
loop = n_elements(accreted1)-1
for j = 0L,loop do begin
star_haloinf.haloid[accreted1[j]] = haloid[stat_row_number[j]-1]
star_haloinf.Mvir[accreted1[j]] = Mvir[stat_row_number[j]-1]
star_haloinf.Gmass[accreted1[j]] = Gmass[stat_row_number[j]-1]
star_haloinf.Smass[accreted1[j]] = Smass[stat_row_number[j]-1]
endfor
endif

;****************************************************
; find most massive halo that each star was ever in 
;*****************************************************

insats = where(phase.grp[i-1] gt gal,nsats)
;if max(satsid) ge 0 then begin
if nsats ne 0 then begin
stat_row_number = phase[insats].grp[i-1]
for j = 0L,n_elements(insats)-1 do star_halo[i-1].Mvir[insats[j]] = Mvir[stat_row_number[j]-1]
endif

endfor

for j = 0L,nall-1 do star_haloinf.maxhalo[j] = max(star_halo.Mvir[j])

;star_haloinf.comp[z0] = components 
star_haloinf.comp[nongal] = -1

;******************************************************************
; write a binary file with all this info!
;*******************************************************************
outfile=filename+'.allgas.accretion.fits'
mwrfits,star_haloinf,outfile,/create

stop
end
