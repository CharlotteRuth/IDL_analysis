;# just input your favourite galaxy: everything else shoule be aok!
;# "MW1.1024g1bwK"," ", " "
;#

;dir = '/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK/'
;dist_units = 50000.
;filebase = 'h603.cosmo50cmb.2304g2bwK'
;gal = 1
;distr=[0,30000]

;dir = '/astro/net/scratch1/abrooks/FABIO/h516.cosmo25cmb.2304g2bwK/'
;dist_units = 25000.
;filebase = 'h516.cosmo25cmb.2304g2bwK'
;gal = 1
;distr=[0,15000]
;step = ['00072','00120','00168','00216','00264','00312','00360','00408','00456','00512']
;step = ['00024','00036','00048','00060','00072','00084','00096','00108','00512']

;dir = '/astro/net/nbody1/abrooks/h799.cosmo25cmb.3072g1MBWK/'
;filebase = 'h799.cosmo25cmb.3072g1MBWK'

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/'
;filebase = 'h516.cosmo25cmb.3072g14HBWK'

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/'
;filebase = 'h516.cosmo25cmb.3072g1MBWK'
;step =['00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00492']

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/'
;filebase = 'h516.cosmo25cmb.3072g14HBWK'
;step = [00024,00032,00080,00088,00092,00104,00128,00132,00144,00148,00156,00168,00180,00192,00204,00216,00228,00240,00252,00260,00272,00284,00308,00320,00332,00344,00356,00368,00380,00392,00404,00420,00432,00436,00444,00456,00468,00480,00492,00504,00512]

;dir ='/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/'
;filebase = 'h603.cosmo50cmb.3072g14HBWK'
;step = [00072,00084,00096,00108,00120,00132,00144,00156,00168,00180,00192,00204,00216,00228,00240,00252,00264,00276,00288,00312,00324]
;distr = [0,999.]
;threshold = 25000.0
;mergertree,dir,filebase,threshold = threshold,step = step,/plot

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/'
;filebase = 'h986.cosmo50cmb.3072g14HBWK'
;step = ['00012','00024','00036','00048','00060','00072','00084','00156','00168','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00396','00408','00432','00444','00456','00468','00492','00504','00512']
;step = ['00300','00324','00348','00372','00408','00444','00468','00504','00512']
;step = ['00300','00372','00444','00512']
;step = FIX(step)
;distr = [0,999.]
;threshold = 25000.0
;mergertree,dir,filebase,threshold = threshold,step = step,/plot


;dir ='/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/'
;filebase = 'h603.cosmo50cmb.3072g14HBWK'
;distr = [999.]
;mergertree,dir,filebase,/plot

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/'
;filebase = 'h603.cosmo50cmb.3072gs1MbwK'
;step = [24,36,512]

PRO mergertree,filebase,dir = dir,gal = gal,timer=timer,distr=distr,fraction=fraction,threshold=threshold,_extra=_extra,plot = plot,step = step

IF NOT keyword_set(dir) THEN dir = './'

IF n_params() EQ 0 THEN BEGIN
    print,"TRACE_TREE, gal [, timer=timer, distr=distr]"
    RETURN
ENDIF
IF keyword_set(fraction) AND keyword_set(threshold) THEN BEGIN
    print,"TRACE_TREE: set only one of (threshold,fraction)"
    print,"            (default: fraction = 0.005)"
    RETURN
ENDIF

IF NOT keyword_set(distr) THEN distr = [0,799.] 
IF NOT keyword_set(timer) THEN timer = [0,13.8] 
IF NOT keyword_set(gal) THEN gal = 1

;loadct,39
!P.thick = 1.5
!X.Charsize = 1.25
!Y.Charsize = 1.25
!P.CHARSIZE = 1.25

IF NOT KEYWORD_SET(step) THEN BEGIN
    command = "ls "+dir+filebase + "*/" + filebase + ".00???.amiga.grp | grep amiga.grp | sed 's/.amiga.grp//g'"
    spawn,command,file       ; return the files without the .amiga.grp
    grpfile = file+'.amiga.grp'
    statfile = file+'.amiga.stat'
    step_st = strmid(file,4,5,/reverse_offset)
    step = fix(step_st) ;strmid(file,2,3,/reverse_offset))
ENDIF ELSE BEGIN
    step_st = strtrim(STRING(step),2)
    IF ((where(step lt 100))[0] ne -1) THEN step_st[where(step lt 100)] = '000'+step_st[where(step lt 100)]
    IF ((where(step ge 100))[0] ne -1) THEN step_st[where(step ge 100)] = '00'+step_st[where(step ge 100)]
    file = dir+filebase+'.'+step_st+'/'+filebase+'.'+step_st
    grpfile = dir+filebase+'.'+step_st+'/'+filebase+'.'+step_st+'.amiga.grp'
    statfile = dir+filebase+'.'+step_st+'/'+filebase+'.'+step_st+'.amiga.stat'
ENDELSE
time = step/float(max(step))*13.66
nsteps = n_elements(step)
;name_steps = dir+filebase+'.'+step_st+'.dir/'+filebase+'.'+step_st
IF max(step) LT 1000 THEN name_steps = dir+filebase+'.'+step_st+'/'+filebase+'.'+step_st $ ;up to 512
ELSE name_steps = dir+filebase+'.0'+step_st+'/'+filebase+'.0'+step_st ;up to 4096
ind_steps = intarr(nsteps)
FOR i = 0, nsteps - 1 DO BEGIN
    ind_steps[i] = WHERE(strcmp(file, name_steps[i]) eq 1)
ENDFOR

final_halo  = read_ascii_array(grpfile[ind_steps[nsteps - 1]]) 
rtipsy,file[ind_steps[nsteps - 1]],h,g,d,s,/justhead
final_halo = final_halo[h.ngas:h.ngas+h.ndark -1]; restrict to DM only
halo_all = final_halo

readcol,statfile[ind_steps[nsteps - 1]],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,ID_A,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15,I6',/preserve_null,/silent
if (N_ELEMENTS(halo_id eq 0) gt 1 )then readcol,statfile[ind_steps[nsteps - 1]],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15',/preserve_null,/silent
haloind = where(halo_id eq gal)
galmass = strtrim(mvir[haloind]);*8.0) not sure why the 8.0 was there
virial =fltarr(nsteps)
vtime=fltarr(nsteps)
virial[0]=rvir[haloind]
vtime[0]=time[nsteps-1]
finalgal = where(final_halo EQ gal) ; restrict to virialized galaxy
halo = halo_all[finalgal]    ; only trace those in final galaxy at z=0
halo_distance = 0.

IF NOT keyword_set(fraction) THEN fraction = 0.005 ELSE fraction=fraction*1.
IF NOT keyword_set(threshold) THEN threshold = n_elements(finalgal)*fraction ; minimum halo particle number
;  distance = sqrt((xc-xc[0])^2+(yc-yc[0])^2+(zc-zc[0])^2)
get_lun,unit
openw,unit,dir + filebase+'.grp' + strtrim(gal,2)+'.halo_step.out' ; store the information for Evan in his format

IF keyword_set(plot) THEN BEGIN
    set_plot,'ps'
    device,filename = dir + '/mergertree_'+filebase+'.grp'+ strtrim(gal,2) + '.eps',/color,bits_per_pixel=8
ENDIF ELSE set_plot,'x'
plot,[0,0,0],[0,0,0],xrange=distr,yrange=timer,/xstyle,/ystyle,xtit='distance (kpc)',charsize=0.9,pos=[0.15,0.1,0.99,0.95],yticks=7,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ']
xyouts,-distr[1]/8.,6.,"time (Gyrs)",orientation=90.
xyouts,0.95*distr[1],12.5,textoidl('M_{vir}')+galmass+'M'+sunsymbol(),alignment=1.
xyouts,-distr[1]/14.,1.8,'2',charsize=1.4
xyouts,-distr[1]/14.,3.8,'4',charsize=1.4
xyouts,-distr[1]/14.,5.8,'6',charsize=1.4
xyouts,-distr[1]/14.,7.8,'8',charsize=1.4
xyouts,-distr[1]/12.,9.8,'10',charsize=1.4
xyouts,-distr[1]/12.,11.8,'12',charsize=1.4
symsize = 0.45*(alog10(mvir[haloind]))-2.5
symstar = 0.45*(alog10(smass[haloind]))-2.5
loadct,39,/silent
plots,halo_distance[0],time[nsteps-1],psym=sym(1),symsize=symsize
plots,halo_distance[0],time[nsteps-1],psym=sym(1),symsize=symstar,color=250

old_linked = lonarr(1)+gal
; halo = later timestep: progenitor = earlier timestep
l=1

;print,'I:  ',strtrim(i,2),'; step: ',step[nsteps-1],'; ',grpfile[ind_steps[nsteps-1]]

for i =nsteps-1,1,-1 do begin
    print,''
    print,'I:  ',strtrim(i,2),'; step: ',step[i],'; ',grpfile[ind_steps[i]]
    print,'I2: ',strtrim(i-1,2),'; step: ',step[i-1],', ',grpfile[ind_steps[i-1]]
    progenitor_all = read_ascii_array(grpfile[ind_steps[i -1]]) ;get halo ids for previous timestep
    rtipsy,file[ind_steps[i - 1]],h_prog,g_prog,d_prog,s_prog,/justhead
    progenitor_all = progenitor_all[h_prog.ngas:h_prog.ngas+h_prog.ndark -1]

    readcol,statfile[ind_steps[i - 1]],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,ID_A,clean,reality,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15,I6,A9,A9',/silent,/preserve_null
    if (N_ELEMENTS(halo_id eq 0) gt 1 )then readcol,statfile[ind_steps[i-1]],halo_id,ntot,ng,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,clean,reality,FORMAT='I5,L9,L9,L9,L9,F18,F12,F18,F18,F18,F11,F11,F11,F11,F11,F11,F15,F15,F15,A9,A9',/preserve_null,/silent

    falsehalos = where(reality eq 'false',nf)
    if (nf ne 0) then begin
        false_id=halo_id[falsehalos[uniq(halo_id[falsehalos])]]
        for j=0,n_elements(false_id)-1 do begin
            false = where(progenitor_all eq false_id[j],nf)
            if (nf ne 0) then begin
                progenitor_all[false]= replicate(1,n_elements(false))
            endif
        endfor
        print,'False Halos: ',halo_id[falsehalos]
    endif
;    progenitor = progenitor_all[finalgal] ; only trace those in final galaxy at z=0

 ;   virial[l]=rvir[0]
 ;   vtime[l]=time[i]
;oplot,[virial[l],virial[l-1]],[vtime[l],vtime[l-1]],color=50,line=1

;    max_old_halos = max(progenitor)
    max_old_halos = max(progenitor_all)
    num_halo = ulonarr(max_old_halos+1)
    transfer =  ulonarr(max_old_halos+1)
    
    num_prog_all1 = histogram(progenitor_all,reverse_indices=r,min=0) ; Halo ids for all particle
;    num_prog1 = histogram(progenitor,min=0) ;Halo ids of particles that end up in final galaxy at z=0

    ;; find # particles from prog that transfer to halo
    for j = 1L,max_old_halos-1 do begin
 ;       if num_prog1[j] eq 0 then continue
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
            transfer[j] = ntrans ; number going to new halo
            num_halo[j] = maxbin ; number of the new halo
        endif        
    endfor

    ;; find progenitor halos which have greater than "threshold" particles
    nz_all = where(progenitor_all ne 0)
    num_prog_all = histogram(progenitor_all[nz_all],min=0)
    prog_linked_all = where(num_prog_all ge threshold,nprog_all) ; 
;    nz = where(progenitor ne 0)
;    num_prog = histogram(progenitor[nz],min=0)
;    prog_linked = where(num_prog ge threshold,nprog) ;The ids of the halos with more than threshold particles that will end up in main halo
    IF (prog_linked_all[0] NE  -1 ) THEN BEGIN
        ;; get only halos who end up merged in virialized halo
        match_multi,num_halo[prog_linked_all],old_linked,dup=keep
        IF keep[0] EQ -1 THEN BREAK          
        prog_linked = prog_linked_all[keep]
        print,'Prog_linked: ',prog_linked
 ;       IF ( i eq 2 OR i eq 1) then stop
;    ;; add to file showing merger links
;    nz = where(halo ne 0)
;    g = histogram(halo[nz],min=0)
;    dummy = where(g ge threshold,nhalo)
        nhalo = n_elements(rem_dup(num_halo[prog_linked]))
        printf,unit,'step',step[i-1],step[i],n_elements(prog_linked),nhalo,format='(A4,I10,I10,I10,I10)'
        FOR j=0L,n_elements(prog_linked)-1 DO BEGIN
            printf,unit,prog_linked[j],num_halo[prog_linked[j]],num_prog_all[prog_linked[j]],transfer[prog_linked[j]],format='(I14,I10,I10,I10)'
        ENDFOR

        print,' '
        print,'step',step[i-1],step[i],n_elements(prog_linked),nhalo,format='(A4,I10,I10,I10,I10)'
        FOR j=0L,n_elements(prog_linked)-1 DO BEGIN
            print,prog_linked[j],num_halo[prog_linked[j]],num_prog_all[prog_linked[j]],transfer[prog_linked[j]],format='(I14,I10,I10,I10)'
        ENDFOR

;    readcol,statfile[i],halo_id,ntot,ngas,nstar,ndark,mvir,rvir,gmass,smass,dmass,Vmax,Rvmax,Vdisp,xc,yc,zc,vxc,vyc,vzc,tmp,format='l,l,l,l,l,f,f,f,f,f,f,f,f,f,f,f,f,f,f,I'

;        center = MIN(prog_linked) - 1
;        prog_distance = sqrt((xc[prog_linked-1]-xc[center])^2+(yc[prog_linked-1]-yc[center])^2+(zc[prog_linked-1]-zc[center])^2)*1000;dist_units
        match_multi, old_linked, num_halo[prog_linked], dup=dist_link

        centerind = where(halo_id EQ min(prog_linked))
        prog_distance = fltarr(n_elements(prog_linked))

        FOR k =0L,n_elements(prog_linked) -1 DO BEGIN
            statind  = where(halo_id EQ prog_linked[k])
            prog_distance[k] = sqrt((xc[statind]-xc[centerind])^2+(yc[statind]-yc[centerind])^2+(zc[statind]-zc[centerind])^2)*1000;dist_units
            IF dist_link[0] NE -1 THEN $
              IF k GT 0 THEN oplot,[prog_distance[k],halo_distance[dist_link[k]]],[time[i-1],time[i]],line=2
            symsize = .45*(alog10(Mvir[statind]))-2.5
            plots,prog_distance[k],time[i-1],psym=sym(1),symsize=symsize
            IF nstar[statind] GT 1 THEN BEGIN
                symstar = 0.45*(alog10(smass[statind]))-2.5
                plots,prog_distance[k],time[i-1],psym=sym(1),symsize=symstar,color=250
            ENDIF
        ENDFOR
        old_linked = prog_linked
        halo_distance = prog_distance 
    ENDIF
    IF keep[0] EQ -1 THEN BREAK
    halo_all = progenitor_all
ENDFOR                          ;end loop over steps
;  endplot

IF keyword_set(plot) THEN device,/close

close,unit
free_lun,unit

end
;file[ind_steps[nsteps - 1]]
