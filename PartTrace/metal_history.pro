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

PRO metal_history,dir = dir,finalid = finalid,laststep = laststep,centralhalo = centralhalo,debug = debug

;-------------------- Read in simulation data for final step
  IF NOT keyword_set(finalid) THEN finalid = '1' 
  IF NOT keyword_set(laststep) THEN laststep = '512'
  IF NOT keyword_set(dir) THEN dir = '.'
  cd,dir
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.coolontime',file_coolon
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file
  spawn,'ls ' + dir + '/h*.param',pfile
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
  HIs = fltarr(nsteps)
  H2s = fltarr(nsteps)

;-------------------- Read in the information about the gas particles at each timestep
  gpart = mrdfits('grp' + finalid + '.allgas.entropy.fits',1)
  npart = n_elements(gpart)

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
  IF NOT nsteps EQ n_elements(gpart[0].mass) THEN print,'Files are for different time steps'
  indiskarr = intarr(npart) ; List of gas particles that are ever in the disk
  
;-------------------- Step through outputs to find the location of the particles
  FOR i = 0, nsteps - 1 DO BEGIN
     gloc = intarr(npart)  
     stat = read_stat_struc_amiga(halodat[i].file + '.amiga.stat')
     main = where(halodat[i].haloid EQ stat.group)
     satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir AND stat.ngas gt 0)
;     satellites2 = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
;     IF keyword_set(debug) THEN print,halodat[i].file,satellites,satellites2
     z = center[i,1] 
     scale = 1.0/(1.0 + z)

;-------------------- Read Simulation Data ------------
     rtipsy,halodat[i].file,h,g,d,s,/justhead
     read_tipsy_arr,halodat[i].file + '.iord',h,iord,part = 'gas',type = 'long'
     read_tipsy_arr,halodat[i].file + '.HI',h,HI,part = 'gas',type = 'float'
     read_tipsy_arr,halodat[i].file + '.H2',h,H2,part = 'gas',type = 'float'
     read_tipsy_arr,halodat[i].file + '.FeMassFrac',h,Fe,part = 'gas',type = 'float'
     read_tipsy_arr,halodat[i].file + '.OxMassFrac',h,Ox,part = 'gas',type ='float'  
     match,gpart.iord,iord,indg,inda
     HI = HI[inda];*gpart[indg].mass[i]
     H2 = 2.0*H2[inda];*gpart[indg].mass[i]
     Fe = Fe[inda];*gpart[indg].mass[i]
     Ox = Ox[inda];*gpart[indg].mass[i]

     HIm = HI*gpart[indg].mass[i]
     H2m = H2*gpart[indg].mass[i]
     Fem = Fe*gpart[indg].mass[i]
     Oxm = Ox*gpart[indg].mass[i]

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
     IF keyword_set(centralhalo) THEN BEGIN
        test = where(gpart.grp[i] EQ stat[main].group, ntest)
         IF ntest NE 0 THEN inhalo = [inhalo,test]
     ENDIF ELSE BEGIN
        FOR j = 0, n_elements(satellites) - 1 DO BEGIN
           test = where(gpart.grp[i] EQ stat[satellites[j]].group, ntest)
           IF ntest NE 0 THEN inhalo = [inhalo,test]
        ENDFOR
     ENDELSE
     IF n_elements(inhalo) NE 1 THEN BEGIN
        inhalo = inhalo[1:n_elements(inhalo) - 1]
        gloc[inhalo] = 1
     ENDIF ELSE BEGIN
        print,'No Gas Particles'
        zmetals[i] = 0
        Oxs[i] = 0
        Fes[i] = 0
        HIs[i] = 0
        H2s[i] = 0
        coldgas[i] = 0
        print,i,zmetals[i],Oxs[i],Fes[i],HIs[i],H2s[i],coldgas[i]  
        CONTINUE
     ENDELSE

;-------------------- Gas particles that are part of the disk (in the galaxy, with density GE 0.1, temperature LE 2e4, and abs(z) LE 10kpc) are set to 3
     disk = where(gloc EQ 1 AND gpart.rho[i] GE 0.1 AND gpart.temp[i] LE 1.2e4 AND abs(gpart.z[i]) LE 10.0)
;     disk = where(gpart.grp[i] EQ halodat[i].haloid AND gpart.rho[i] GE 0.1 AND gpart.temp[i] LE 2e4 AND abs(gpart.z[i]) LE 10.0)
     zmetals[i] = total((HI + H2)*(2.09*Oxm + 1.06*Fem))/max(HI+H2)
     Oxs[i] = total((HI + H2)*Oxm)/max(HI+H2)
     Fes[i] = total((HI + H2)*Fem)/max(HI+H2)
     HIs[i] = total(HIm)
     H2s[i] = total(H2m)
     coldgas[i]=total(HIm + H2m)
     print,i,zmetals[i],Oxs[i],Fes[i],HIs[i],H2s[i],coldgas[i]
;     stop

     IF keyword_set(debug) AND keyword_set(plots) THEN BEGIN
;     IF 0 THEN BEGIN
        print,halodat[i].file,minmax(gpart[inhalo].grp[i])
        window,2,xsize = 400,ysize = 400
        plot,gpart[inhalo].x[i],gpart[inhalo].z[i],psym = 3,title = halodat[i].file,xtitle = 'X',ytitle = 'Y';,yrange = [-1*stat[main].rvir,stat[main].rvir],xrange = [-1*stat[main].rvir,stat[main].rvir]
        IF disk[0] NE -1 THEN oplot,gpart[disk].x[i],gpart[disk].z[i],psym = 3,color = 245
        IF (inhalo[ind_unbound])[0] NE -1 THEN oplot,gpart[inhalo[ind_unbound]].x[i],gpart[inhalo[ind_unbound]].z[i],psym = 3,color = 100
        window,3,xsize = 600,ysize = 400
        plot,gpart[inhalo].rho[i],gpart[inhalo].temp[i],psym = 3,/xlog,/ylog,xtitle = 'Density [amu/cc]',ytitle = 'Temperature [K]',xrange = [1e-7,1e4],yrange = [10,1e8],title = halodat[i].file,xstyle = 1
        IF disk[0] NE -1 THEN oplot,gpart[disk].rho[i],gpart[disk].temp[i],psym = 3,color = 245
        IF (inhalo[ind_unbound])[0] NE -1 THEN oplot,gpart[inhalo[ind_unbound]].rho[i],gpart[inhalo[ind_unbound]].temp[i],psym = 3,color = 100
        window,0,xsize = 600,ysize = 400
        plot,magpos[ind_unbound],tanvel[ind_unbound],psym = 3

        window,1,xsize = 600,ysize = 400
        histogramp,abs(gpart[inhalo].z[i]),nbins = 100
        histogramp,abs(gpart[inhalo[ind_unbound]].z[i]),/overplot,color = 245,nbins = 100        
     ENDIF
  ENDFOR
  writecol,'grp' + finalid + '.metals.txt',zmetals,Oxs,Fes,HIs,H2s,coldgas
END


PRO master_metal_history

  dir = '/nobackupp2/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
  finalid = '2'
  metal_history,dir = dir,finalid = finalid
  finalid = '1'
  metal_history,dir = dir,finalid = finalid  

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
