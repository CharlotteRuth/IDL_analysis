;# just input your favourite galaxy: everything else shoule be aok!
;# "MW1.1024g1bwK"," ", " "
;#
pro mergertree,gal,fraction=fraction,threshold=threshold,_extra=_extra

if n_params() eq 0 then begin
   print,"TRACE_TREE, gal"
   return
endif
if keyword_set(fraction) and keyword_set(threshold) then begin
   print,"TRACE_TREE: set only one of (threshold,fraction)"
   print,"            (default: fraction = 0.01)"
   return
endif


  mdark =  759941.              ; solar masses
  dist_units = 1000.
;  dir='/astro/net/scratch1/abrooks/FABIO/gal1c6.2304g1bwK/'
;  dir='/astro/net/scratch1/abrooks/FABIO/MW1.1024g1bwK/'
;   dir='/astro/net/scratch1/abrooks/FABIO/dwf1a.1024g1bwK/'
;dir='/net/scratch2/fabio/h986.cosmo50cmb.2304g2bwK/'
; dir='/astro/net/scratch2/fabio/cosmo25runs/h579.cosmo25.1536g1bwK/'
;  dir='/astro/net/scratch2/fabio/cosmo25runs/h590.cosmo25.1536g1bwK/'
dir='/net/scratch1/cbrook/sims/h516.2304g4/'
  command = "ls "+dir+"*/*amiga.grp | grep amiga.grp | sed 's/.amiga.grp//g'"
  spawn,command,file ; return the files without the .amiga.grp
  grpfile = file+'.amiga.grp'
  statfile = file+'.amiga.stat'
  step = fix(strmid(file,2,3,/reverse_offset))
  time = step/512.*13.66
  nsteps = n_elements(step)

  get_lun,unit
  openw,unit,'halo_step.out' ; store the information for Evan in his format
   final_halo  = read_ascii_array(grpfile[nsteps-1]) 
  rheader,file[nsteps-1],h
  final_halo = final_halo[h.ngas:h.ngas+h.ndark -1]; restrict to DM only
   readcol,statfile[nsteps-1],halo_id,ntot,ngas,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,tmp,format='l,l,l,l,l,f,f,f,f,f,f,f,f,f,f,f,f,f,f,I',numline = 1,skipline=gal ; read information about each halo

 finalgal = where(final_halo eq gal) ; restrict to virialized galaxy
  halo_distance = 0.

  if not keyword_set(fraction) then fraction = 0.0001 else fraction=fraction*1.
  if not keyword_set(threshold) then threshold = n_elements(finalgal)*fraction ; minimum halo particle number

  stepsize = 1
  old_linked = lonarr(1)+gal

  for i =nsteps-stepsize-1,0,-stepsize do begin
    if i lt 0 then break
    halo_all  = read_ascii_array(grpfile[i+stepsize]) ; get halo ids 
    rheader,file[i+stepsize],h_halo
    halo_all = halo_all[h_halo.ngas:h_halo.ngas+h_halo.ndark -1]

if  (i+stepsize eq nsteps-1 ) then begin 

readcol,statfile[i+stepsize],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,ID_A,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15,I6',/preserve_null

endif

if (i+stepsize ne nsteps-1 ) then begin 

readcol,statfile[i+stepsize],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,ID_A,clean,reality,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15,I6,A9,A9',/preserve_null

falsehalos = where(reality eq 'false',nf)
if (nf ne 0) then begin
false_id=halo_id[falsehalos[uniq(halo_id[falsehalos])]]
for j=0,n_elements(false_id)-1 do begin
false = where(halo_all eq false_id[j],nf)
if (nf ne 0) then begin
halo_all[false]= replicate(1,n_elements(false))
endif
endfor
endif

endif
    halo = halo_all[finalgal] ; only trace those in final galaxy at z=0
    progenitor_all = read_ascii_array(grpfile[i]) 
    rheader,file[i],h_prog
    progenitor_all = progenitor_all[h_prog.ngas:h_prog.ngas+h_prog.ndark -1]

readcol,statfile[i],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,ID_A,clean,reality,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15,I6,A9,A9',/silent,/preserve_null

falsehalos = where(reality eq 'false',nf)
if (nf ne 0) then begin
false_id=halo_id[falsehalos[uniq(halo_id[falsehalos])]]
for j=0,n_elements(false_id)-1 do begin
false = where(progenitor_all eq false_id[j],nf)
if (nf ne 0) then begin
progenitor_all[false]= replicate(1,n_elements(false))
endif
endfor
endif

    progenitor = progenitor_all[finalgal] ; only trace those in final galaxy at z=0
    
    max_old_halos = max(progenitor)
    num_halo = ulonarr(max_old_halos+1)
    transfer =  ulonarr(max_old_halos+1)
    
    num_prog_all1 = histogram(progenitor_all,reverse_indices=r,min=0)
    num_prog1 = histogram(progenitor,min=0)
    
    ;; find # particles from prog that transfer to halo
    for j = 1L,max_old_halos-1 do begin
       if num_prog1[j] eq 0 then continue
       if num_prog_all1[j] ge threshold then begin
          member = r[r[j]:r[j+1]-1] 
          halo_list = halo_all[member]
          h=histogram(halo_list,min=0)
          ntrans = max(h,maxbin) ; # of particles that go to new halo
          if maxbin eq 0 and n_elements(h) ne 1 then begin ; for the case of more than half the particles being stripped, but still a large number of particles remain in the halo
             if max(h[1:*]) ge threshold then begin
               ntrans =  max(h[1:*],maxbin)
                maxbin = maxbin+1
             endif
          endif          
          transfer[j] = ntrans  ; number going to new halo
          num_halo[j] = maxbin  ; number of the new halo
       endif

    endfor


    ;; find progenitor halos which have greater than "threshold" particles
    nz_all = where(progenitor_all ne 0)
    num_prog_all = histogram(progenitor_all[nz_all],min=0)
    prog_linked_all = where(num_prog_all ge threshold,nprog_all) ; 
    nz = where(progenitor ne 0)
    num_prog = histogram(progenitor[nz],min=0)
    prog_linked = where(num_prog ge threshold,nprog) ; 

    ;; get only halos who end up merged in virialized halo
    match_multi,num_halo[prog_linked],old_linked,dup=keep
    prog_linked = prog_linked[keep]

;    ;; add to file showing merger links
    nhalo = n_elements(rem_dup(num_halo[prog_linked]))
    printf,unit,'step',step[i],step[i+stepsize],n_elements(prog_linked),nhalo,format='(A4,I8,I8,I8,I8)'
     for j=0L,n_elements(prog_linked)-1 do begin
      printf,unit,prog_linked[j],num_halo[prog_linked[j]],num_prog_all[prog_linked[j]],transfer[prog_linked[j]],format='(I12,I8,I8,I8)'
    endfor


    old_linked = prog_linked

  endfor                        ;end loop over steps
  close,unit
  free_lun,unit
  stop


end
