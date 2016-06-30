;# just input your favourite galaxy: everything else shoule be aok!
;# "MW1.1024g1bwK"," ", " "
;#
pro trace_tree,gal,timer=timer,distr=distr,fraction=fraction,threshold=threshold,_extra=_extra

  if n_params() eq 0 then begin
    print,"TRACE_TREE, gal [, timer=timer, distr=distr]"
    return
  endif
  if keyword_set(fraction) and keyword_set(threshold) then begin
    print,"TRACE_TREE: set only one of (threshold,fraction)"
    print,"            (default: fraction = 0.01)"
    return
  endif
  if not keyword_set(distr) then distr = [0,799.] 
  if not keyword_set(timer) then timer = [0,13.8] 


  dist_units = 1000.
;  dist_units = 25000.0
;  dir = '/net/scratch1/cbrook/sims/h516.2304g4bwdK/'
;  prefix='h516.cosmo25cmb.2304g4bwK.'
  dir= '/net/scratch1/abrooks/FABIO/h516.cosmo25cmb.2304g2bwK/'
  prefix='h516.cosmo25cmb.2304g2bwK.'
  steps=['00048','00084','00144','00216','00288','00336','00408','00456','00512']

  dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/'
  prefix = 'h516.cosmo25cmb.1536g14HBWK'
  steps = ['00024','00048','00072','00092','00096','00100','00103','00104','00108','00112','00116','00120','00128','00132','00136','00140','00144','00148','00152','00156','00160','00164','00167','00168','00172','00176','00180','00184','00188','00192','00196','00216','00240','00264','00288','00312','00336','00360','00384','00408','00432','00456','00480','00504','00512']

  dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/steps/'
  prefix = 'h516.cosmo25cmb.1536g1MBWK'
  steps = ['00024','00037','00046','00048','00061','00072','00084','00096','00120','00128','00144','00168','00192','00216','00227','00240','00264','00271','00288','00312','00328','00336','00360','00384','00406','00408','00432','00455','00456','00480','00504','00512']


  dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/'
  prefix = 'h603.cosmo50cmb.3072g14HBWK'
  steps = ['00072','00084','00096','00108','00120','00132','00144','00168','00180','00192','00204','00216','00312','00324']
  dist_units = 50000.0

  file = dir+prefix+'.'+steps+'.dir/'+prefix+'.'+steps
  
;  dir = '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/'
;  prefix = 'h986.cosmo50cmb.3072gs1MbwK'
;  steps = ['00014','00036','00048','00060','00504']
;  dist_units = 50000.0

;  file = dir+prefix+'.'+steps+'/'+prefix+'.'+steps

  nsteps = n_elements(file)
  step = fix(strmid(steps,2,3,/reverse_offset))
  grpfile = file+'.amiga.grp'
  statfile = file+'.amiga.stat'
  time = step/512.*13.66

  
  get_lun,unit
  openw,unit,'halo_step.out'    ; store the information for Evan in his format
;  begplot,'mergertree_old_005.ps',xsize=5.5,ysize=5.5,/color
  begplot,prefix+'_mergertree_old_005.ps',xsize=5.5,ysize=5.5,/color
  final_halo  = read_ascii_array(grpfile[nsteps-1]) 
;  rheader,file[nsteps-1],h
  rtipsy,file[nsteps-1],h,g,d,s,/justhead
  final_halo = final_halo[h.ngas:h.ngas+h.ndark -1] ; restrict to DM only
;  readcol,statfile[nsteps-1],halo_id,ntot,ngas,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,tmp,format='l,l,l,l,l,f,f,f,f,f,f,f,f,f,f,f,f,f,f,I',numline = 1,skipline=gal ; read information about each halo
  readcol,statfile[nsteps-1],halo_id,ntot,ngas,nstar,ndark,mvir,rvir, $
    gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc, $
    vxc,vyc,vzc,clean, /silent, format ='(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,f,F,F,F,A)',numline = 1,skipline=gal  

  galmass=string(sigfig(mvir[0],2,/scientific))
  virial =fltarr(nsteps)
  vtime=fltarr(nsteps)
  virial[0]=rvir[0]
  vtime[0]=time[nsteps-1]
  finalgal = where(final_halo eq gal) ; restrict to virialized galaxy
  halo_distance = 0.
  
  if not keyword_set(fraction) then fraction = 0.005 else fraction=fraction*1.
  if not keyword_set(threshold) then threshold = n_elements(finalgal)*fraction ; minimum halo particle number
;  distance = sqrt((xc-xc[0])^2+(yc-yc[0])^2+(zc-zc[0])^2)
  plot,[0,0,0],[0,0,0],xrange=distr,yrange=timer,/xstyle,/ystyle,xtit='distance (kpc)',charsize=0.9,pos=[0.15,0.1,0.99,0.95],yticks=7,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ']
  
  xyouts,-distr[1]/8.,6.,"time (Gyrs)",orientation=90.
  xyouts,0.95*distr[1],12.5,"Mvir ="+galmass+'M'+sunsymbol(),alignment=1.
  xyouts,-distr[1]/14.,1.8,'2',charsize=1.4
  xyouts,-distr[1]/14.,3.8,'4',charsize=1.4
  xyouts,-distr[1]/14.,5.8,'6',charsize=1.4
  xyouts,-distr[1]/14.,7.8,'8',charsize=1.4
  xyouts,-distr[1]/12.,9.8,'10',charsize=1.4
  xyouts,-distr[1]/12.,11.8,'12',charsize=1.4
  symsize = 0.45*(alog10(mvir[0]))-2.5
  symstar = 0.45*(alog10(smass[0]))-2.5
  loadct,39,/silent
  plots,halo_distance[0],time[nsteps-1],psym=sym(1),symsize=symsize
  plots,halo_distance[0],time[nsteps-1],psym=sym(1),symsize=symstar,color=250
  
  stepsize = 1
  old_linked = lonarr(1)+gal
; halo = later timestep: progenitor = earlier timestep
  l=1
  for i =nsteps-2,1,-1 do begin
    if i lt 0 then break
    halo_all  = read_ascii_array(grpfile[i+1]) ; get halo ids 
;    rheader,file[i+1],h_halo
    rtipsy,file[i+1],h_halo,temp,temp,temp,/justhead
    halo_all = halo_all[h_halo.ngas:h_halo.ngas+h_halo.ndark -1]
    
    if  (i+stepsize eq nsteps-1 ) then begin 
      
;      readcol,statfile[i+stepsize],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,ID_A,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15,I6',/preserve_null
        readcol,statfile[i+stepsize],halo_id,ntot,ngas,nstar,ndark,mvir,rvir, $
          gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc, $
          vxc,vyc,vzc,clean, /silent, format ='(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,f,F,F,F,A)',/preserve_null
    endif
    
    if (i+stepsize ne nsteps-1 ) then begin 
      
;      readcol,statfile[i+stepsize],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,ID_A,clean,reality,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15,I6,A9,A9',/preserve_null
        readcol,statfile[i+stepsize],halo_id,ntot,ngas,nstar,ndark,mvir,rvir, $
          gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc, $
          vxc,vyc,vzc,clean,reality, /silent, format ='(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,f,F,F,F,A,A)',/preserve_null      
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
    halo = halo_all[finalgal]   ; only trace those in final galaxy at z=0
    progenitor_all = read_ascii_array(grpfile[i]) 
;    rheader,file[i],h_prog
    rtipsy,file[i],h_prog,temp,temp,temp,/justhead
    progenitor_all = progenitor_all[h_prog.ngas:h_prog.ngas+h_prog.ndark -1]
    
;    readcol,statfile[i],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,ID_A,clean,reality,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15,I6,A9,A9',/silent,/preserve_null
    readcol,statfile[i+stepsize],halo_id,ntot,ngas,nstar,ndark,mvir,rvir, $
      gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc, $
      vxc,vyc,vzc,clean,reality, /silent, format ='(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,f,F,F,F,A,A)',/preserve_null 
    
    virial[l]=rvir[0]
    vtime[l]=time[i]
;oplot,[virial[l],virial[l-1]],[vtime[l],vtime[l-1]],color=50,line=1
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
        ntrans = max(h,maxbin)  ; # of particles that go to new halo
        if maxbin eq 0 and n_elements(h) ne 1 then begin ; for the case of more than half the particles being stripped, but still a large number of particles remain in the halo
          if max(h[1:*]) ge threshold then begin
            ntrans =  max(h[1:*],maxbin)
            maxbin = maxbin+1
          endif
        endif          
        transfer[j] = ntrans    ; number going to new halo
        num_halo[j] = maxbin    ; number of the new halo
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
;    nz = where(halo ne 0)
;    g = histogram(halo[nz],min=0)
;    dummy = where(g ge threshold,nhalo)
    nhalo = n_elements(rem_dup(num_halo[prog_linked]))
    printf,unit,'step',step[i],step[i+stepsize],n_elements(prog_linked),nhalo,format='(A4,I8,I8,I8,I8)'
    for j=0L,n_elements(prog_linked)-1 do begin
      printf,unit,prog_linked[j],num_halo[prog_linked[j]],num_prog_all[prog_linked[j]],transfer[prog_linked[j]],format='(I12,I8,I8,I8)'
    endfor
    
;    readcol,statfile[i],halo_id,ntot,ngas,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,tmp,format='l,l,l,l,l,f,f,f,f,f,f,f,f,f,f,f,f,f,f,I'
    prog_distance = sqrt((xc[prog_linked-1]-xc[0])^2+(yc[prog_linked-1]-yc[0])^2+(zc[prog_linked-1]-zc[0])^2)*dist_units
    
    match_multi, old_linked, num_halo[prog_linked], dup=dist_link
    
    
    for k =0L,n_elements(prog_linked) -1 do begin
      if k gt 0 then oplot,[prog_distance[k],halo_distance[dist_link[k]]],[time[i],time[i+stepsize]],line=2
      symsize = .45*(alog10(Mvir[prog_linked[k]-1]))-2.5
      plots,prog_distance[k],time[i],psym=sym(1),symsize=symsize
      if nstar[prog_linked[k]-1] gt 1 then begin
        symstar = 0.45*(alog10(smass[prog_linked[k]-1]))-2.5
        plots,prog_distance[k],time[i],psym=sym(1),symsize=symstar,color=250
      endif
      
    endfor
    old_linked = prog_linked
    halo_distance = prog_distance 
    l=l+1
  endfor                        ;end loop over steps
  endplot
  close,unit
  free_lun,unit
  stop
  
  
end

