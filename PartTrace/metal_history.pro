;Charlotte Christensen
;5/25/12
;This program uses grp1.allgas.entropy.fits and alignment.txt to
;establish when the gas is in the main halo and when it is in the disk
;of the main halo

;run from the main simulation directory

;  dir = '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g1bwK'
;  dir = '/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK'
;  dir = '/astro/store/nbody3/christensen/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK'
  
;  dir = '/astro/store/nbody3/christensen/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK'
;  dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK'

;keyword parameter "inhalo" outputs the metals in the halo in addition
;to those in the disk

PRO metal_history,dir = dir,finalid = finalid,laststep = laststep,centralhalo = centralhalo,debug = debug,nowrite = nowrite,inhalo = inhalo,plots= plots

;-------------------- Read in simulation data for final step
  IF NOT keyword_set(finalid) THEN finalid = '1' 
  IF NOT keyword_set(laststep) THEN laststep = '512'
  IF NOT keyword_set(dir) THEN dir = '.'
  cd,dir
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.coolontime',file_coolon
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file
  spawn,'ls ' + dir + '/*.param',pfile
  units = tipsyunits(pfile[0])

;-------------------- Read in the information on the halos
  halodat = mrdfits('grp' + finalid + '.alignment.fits',1)
  nsteps = n_elements(halodat.file)
  halodat.xc = halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.yc = halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.zc = halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  center = [[halodat.time],[halodat.z],[halodat.xc],[halodat.yc],[halodat.zc]]
  center[*,2] = center[*,2]*1000.0 - units.lengthunit/2.0 ;"center" is in Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
  center[*,3] = center[*,3]*1000.0 - units.lengthunit/2.0
  center[*,4] = center[*,4]*1000.0 - units.lengthunit/2.0
  a = [[halodat.xa],[halodat.ya],[halodat.za]]
  zmetals = fltarr(nsteps)
  Fes = fltarr(nsteps)
  Oxs = fltarr(nsteps)
  coldgas=fltarr(nsteps)
  zmetals_ha = fltarr(nsteps)
  Fes_ha = fltarr(nsteps)
  Oxs_ha = fltarr(nsteps)
  zmetals_H = fltarr(nsteps)
  Fes_H = fltarr(nsteps)
  Oxs_H = fltarr(nsteps)
  HIs = fltarr(nsteps)
  H2s = fltarr(nsteps)
  Hgas=fltarr(nsteps)  

;-------------------- Read in the information about the gas particles at each timestep
  gpart_save = mrdfits('grp' + finalid + '.allgas.entropy.fits',1)

;-------------------- Need to exclude 'early' particles to compare to the outflow files
;  early = mrdfits('early.iord.fits',0)
;  all = gpart.iord
;  test = binfind(all,early)
;  early = 0           ;delete early array
;  all[test] = -1      ;mark all particles in the early array
;  exclude = where(all eq -1, comp=keep) ;find all particles in the early array
;  exclude = 0         ;delete exclude array
;  gpart = gpart[keep] ;gpart now only contains the particles that were not accreted early

;-------------------- Make sure the halo information and gpart are for equivalent outputs
  IF nsteps NE n_elements(gpart_save[0].mass) THEN BEGIN 
     print,'Files are for different time steps'
     
     stop
     stop

     indsave  = [indgen(71),indgen(51)+72] ;Skips step 308 from h239
     indsave  = [indgen(100),indgen(24)+101] ;Skips step 420 from h285
     gpart_clip = replicate({iord:0L, mass:LONARR(nsteps), grp:LONARR(nsteps), temp:FLTARR(nsteps), rho:FLTARR(nsteps), entropy:FLTARR(nsteps), radius:FLTARR(nsteps), vx:FLTARR(nsteps), vy:FLTARR(nsteps), vz:FLTARR(nsteps), x:FLTARR(nsteps), y:FLTARR(nsteps), z:FLTARR(nsteps)}, n_elements(gpart_save.iord))
     gpart_clip.iord = gpart_save.iord
     gpart_clip.mass = gpart_save.mass[indsave]
     gpart_clip.grp = gpart_save.grp[indsave]
     gpart_clip.temp = gpart_save.temp[indsave]
     gpart_clip.rho = gpart_save.rho[indsave]
     gpart_clip.entropy = gpart_save.entropy[indsave]
     gpart_clip.radius = gpart_save.radius[indsave]
     gpart_clip.x = gpart_save.x[indsave]
     gpart_clip.y = gpart_save.y[indsave]
     gpart_clip.z = gpart_save.z[indsave]
     gpart_clip.vx = gpart_save.vx[indsave]
     gpart_clip.vy = gpart_save.vy[indsave]
     gpart_clip.vz = gpart_save.vz[indsave]
     print,'Are you sure you want to rewrite grp?.allgas.entropy.fits? If so, continue.'
     stop
     mwrfits, gpart_clip, 'grp' + finalid + '.allgas.entropy.fits', /create
     gpart_save = gpart_clip
  ENDIF
  
;-------------------- Step through outputs to find the location of the
;                     particles
  FOR i = 0, nsteps - 1 DO BEGIN
     stat = read_stat_struc_amiga(halodat[i].file + '.amiga.stat')
     main = where(halodat[i].haloid EQ stat.group)
     satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir AND stat.ngas gt 0)
;     satellites2 = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
;     IF keyword_set(debug) THEN print,halodat[i].file,satellites,satellites2
     z = center[i,1] 
     scale = 1.0/(1.0 + z)

;-------------------- Read Simulation Data ------------
     IF keyword_set(debug) THEN rtipsy,halodat[i].file,h,g,d,s ELSE rtipsy,halodat[i].file,h,g,d,s,/justhead
     read_tipsy_arr,halodat[i].file + '.iord',h,iord,part = 'gas',type = 'long'
     read_tipsy_arr,halodat[i].file + '.HI',h,HI,part = 'gas',type = 'float'
     read_tipsy_arr,halodat[i].file + '.H2',h,H2,part = 'gas',type = 'float'
     read_tipsy_arr,halodat[i].file + '.FeMassFrac',h,Fe,part = 'gas',type = 'float'
     read_tipsy_arr,halodat[i].file + '.OxMassFrac',h,Ox,part = 'gas',type ='float'  
     match,gpart_save.iord,iord,indg,inda ;Only look at the gas particles that will ever be in this halo
     ;Make sure that there any gas particle that is not at this time step has a mass of zero at that timestep in gpart_save 
     IF n_elements(indg) NE n_elements(gpart_save.iord) THEN BEGIN 
        match2,gpart_save.iord,iord,indg2,inda2 ;Select for only those gas particles that are not still around at this timestep
        IF max(gpart_save[where(indg2 EQ -1)].mass[i]) GT 0 THEN stop ;Check for h603, 3
     ENDIF
     HI = HI[inda]
     H2 = 2.0*H2[inda]
     Fe = Fe[inda]
     Ox = Ox[inda]
     iord_arr = iord
     iord = iord[inda]
     gpart = gpart_save[indg]
     gloc = intarr(n_elements(gpart)) 

     HIm = HI*gpart.mass[i]
     H2m = H2*gpart.mass[i]
     Fem = Fe*gpart.mass[i]
     Oxm = Ox*gpart.mass[i]

;-------------------- Convert the position information gpart to the frame of the halo
     az = reform(a[i,*])
     az0 = az[0]
     az1 = az[1]
     az2 = az[2]
     ax = [az2/sqrt(az0*az0 + az2*az2),0,-1.0*az0/sqrt(az0*az0 + az2*az2)]
     ay = crossp(az,ax)         ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
     basis = [[ax],[ay],[az]]
     gpos = [[gpart.x[i] - center[i,2]*scale],[gpart.y[i] - center[i,3]*scale],[gpart.z[i] - center[i,4]*scale]];*scale
     gpos = transpose(transpose(basis)#transpose(gpos))
     gpart.x[i] = gpos[*,0]
     gpart.y[i] = gpos[*,1]
     gpart.z[i] = gpos[*,2]

;-------------------- Gas particles that have been deleted are set to -1
;     gone = where(gpart.grp[i] EQ -1)
;     IF gone[0] NE -1 THEN gloc[gone] = -1

;-------------------- Gas particles that are in the halo/satellites of the halo are set to 2 
;     distance = sqrt(gpart.x[i]*gpart.x[i] + gpart.y[i]*gpart.y[i] + gpart.z[i]*gpart.z[i]) ;should be the same thing to gpart.radius[i]
;     inhalo = where(distance LE halodat[i].rvir)
;     inhalo = where(gpart.grp[i] EQ halodat[i].haloid)
     inhalo = [0]
     test = where(gpart.grp[i] EQ stat[main].group, ntest)
     IF ntest NE 0 THEN inhalo = [inhalo,test]
     IF NOT keyword_set(centralhalo) THEN BEGIN
         IF satellites[0] NE -1 THEN BEGIN
             FOR j = 0, n_elements(satellites) - 1 DO BEGIN               ;& $
                 test = where(gpart.grp[i] EQ stat[satellites[j]].group, ntest) ;& $
                 IF ntest NE 0 THEN inhalo = [inhalo,test]                      ;& $
             ENDFOR                                                             ;& $
        ENDIF
     ENDIF
     IF n_elements(inhalo) NE 1 THEN BEGIN
        inhalo = inhalo[1:n_elements(inhalo) - 1]
        gloc[inhalo] = 1
        HI = HI[inhalo]
        H2 = H2[inhalo]
        Oxm = Oxm[inhalo]
        Fem = Fem[inhalo]
        iord = iord[inhalo]
        gpart = gpart[inhalo]
     ENDIF ELSE BEGIN
        print,'No Gas Particles'
        zmetals[i] = 0
        Oxs[i] = 0
        Fes[i] = 0
        HIs[i] = 0
        H2s[i] = 0
        coldgas[i] = 0
        print,i,zmetals[i],Oxs[i],Fes[i],coldgas[i],zmetals_H[i],Oxs_H[i],Fes_H[i],HIs[i],H2s[i],Hgas[i]
        CONTINUE
     ENDELSE
 ;-------------------- Gas particles that are part of the disk (in the galaxy, with density GE 0.1, temperature LE 1.2e4, and abs(z) LE 3kpc) are set to 3
     disk = where(gpart.rho[i] GE 0.1 AND gpart.temp[i] LE 1.2e4 AND abs(gpart.z[i]) LE 3.0)
;     disk = where(gpart.grp[i] EQ halodat[i].haloid AND gpart.rho[i] GE 0.1 AND gpart.temp[i] LE 2e4 AND abs(gpart.z[i]) LE 10.0)
     zmetals[i] = total(2.09*Oxm[disk] + 1.06*Fem[disk])
     Oxs[i] = total(Oxm[disk])
     Fes[i] = total(Fem[disk])
     zmetals_ha[i] = total(2.09*Oxm + 1.06*Fem)
     Oxs_ha[i] = total(Oxm)
     Fes_ha[i] = total(Fem)
     coldgas[i]=total(gpart[disk].mass[i])
     zmetals_H[i] = total((HI + H2)*(2.09*Oxm + 1.06*Fem))/max(HI+H2)
     Oxs_H[i] = total((HI + H2)*Oxm)/max(HI+H2)
     Fes_H[i] = total((HI + H2)*Fem)/max(HI+H2)
     HIs[i] = total(HIm)
     H2s[i] = total(H2m)
     Hgas[i]=total(HIm + H2m)
     print,i,zmetals[i],Oxs[i],Fes[i],coldgas[i],zmetals_H[i],Oxs_H[i],Fes_H[i],HIs[i],H2s[i],Hgas[i]

     IF keyword_set(debug) AND keyword_set(plots) THEN BEGIN
;     IF 0 THEN BEGIN
        print,halodat[i].file,minmax(gpart[inhalo].grp[i])
        window,2,xsize = 400,ysize = 400
        plot,gpart[inhalo].x[i],gpart[inhalo].z[i],psym = 3,title = halodat[i].file,xtitle = 'X',ytitle = 'Y';,yrange = [-1*stat[main].rvir,stat[main].rvir],xrange = [-1*stat[main].rvir,stat[main].rvir]
        IF disk[0] NE -1 THEN oplot,gpart[disk].x[i],gpart[disk].z[i],psym = 3,color = 245
        window,3,xsize = 600,ysize = 400
        plot,gpart[inhalo].rho[i],gpart[inhalo].temp[i],psym = 3,/xlog,/ylog,xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]',xrange = [1e-7,1e4],yrange = [10,1e8],title = halodat[i].file,xstyle = 1
        IF disk[0] NE -1 THEN oplot,gpart[disk].rho[i],gpart[disk].temp[i],psym = 3,color = 245
         stop
     ENDIF
     gpart = 0
  ENDFOR
  IF NOT keyword_set(nowrite) THEN writecol,'grp' + finalid + '.metals.txt',zmetals,Oxs,Fes,coldgas,zmetals_H,Oxs_H,Fes_H,HIs,H2s,Hgas,format='(10E)'
  IF keyword_set(inhalo) AND NOT keyword_set(nowrite) THEN writecol,'grp' + finalid + '.halometals.txt',zmetals_ha,Oxs_ha,Fes_ha,format='(3E)'
END


PRO master_metal_history

  dir = '/nobackupp2/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
  finalid = '2'
  metal_history,dir = dir,finalid = finalid
  finalid = '1'
  metal_history,dir = dir,finalid = finalid  

  dir = '/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
  finalid = '1'
  metal_history,dir = dir,finalid = finalid,/inhalo

  dir = '/nobackupp2/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
  finalid = '2'
  metal_history,dir = dir,finalid = finalid
  finalid = '1'
  metal_history,dir = dir,finalid = finalid

  dir = '/nobackupp2/crchrist/MolecH/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
  finalid = '3'
  metal_history,dir = dir,finalid = finalid
  finalid = '2'
  metal_history,dir = dir,finalid = finalid
  finalid = '1'
  metal_history,dir = dir,finalid = finalid

  dir = '/nobackupp2/crchrist/MolecH/h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
  finalid = '2'
  metal_history,dir = dir,finalid = finalid
  finalid = '1'
  metal_history,dir = dir,finalid = finalid
END
