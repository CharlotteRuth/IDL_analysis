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

PRO stellarmetal_history,dir = dir,finalid = finalid,laststep = laststep,centralhalo = centralhalo,debug = debug,nowrite = nowrite
  dDelta = 0.000691208 ;516, 799 0.000691088, 603, 986 (0.000691208), 239 (0.000691208)

;-------------------- Read in simulation data for final step
  IF NOT keyword_set(finalid) THEN finalid = '1' 
  IF NOT keyword_set(laststep) THEN laststep = '512'
  IF NOT keyword_set(dir) THEN dir = '.'
  cd,dir

  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.massform',file_massform
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file
  spawn,'ls ' + dir + '/h*.param',pfile
;  spawn,'ls ' + dir + '/*starlog',starlogfile
  rtipsy,file,h,g,d,s,/justhead
  read_tipsy_arr,file + '.iord',h,iordlast,part = 'star',type = 'long'
  read_tipsy_arr,file + '.massform',h,massform,part = 'star',type ='float'
  units = tipsyunits(pfile[0])
  massform = massform*units.massunit
 ; dDelta = 7.812500e-04*units.timeunit

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

;-------------------- Step through outputs to find the location of the
;                     particles
  FOR i = nsteps - 1, nsteps - 1 DO BEGIN ;For debugging purposes
     stat = read_stat_struc_amiga(halodat[i].file + '.amiga.stat')
     main = where(halodat[i].haloid EQ stat.group)
     satellites = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir AND stat.ngas gt 0)
;     satellites2 = where(sqrt((stat.xc - stat[main].xc)^2 + (stat.yc - stat[main].yc)^2 + (stat.zc - stat[main].zc)^2)*1000 LE stat[main].rvir)
;     IF keyword_set(debug) THEN print,halodat[i].file,satellites,satellites2
     z = center[i,1] 
     scale = 1.0/(1.0 + z)
     str = strsplit(halodat[i].file,'.',/extract)
     step = str[n_elements(str) - 1]
     tmax = (float(step)*dDelta + dDelta)*units.timeunit
     IF i EQ 0 THEN tmin = 0 ELSE BEGIN
        str = strsplit(halodat[i - 1].file,'.',/extract)
        step_prior = str[n_elements(str) - 1]  
        tmin = (float(step_prior)*dDelta + dDelta)*units.timeunit
     ENDELSE

     tmin = 0 ; For debugging purposes
;-------------------- Read Simulation Data ------------
     rtipsy,halodat[i].file,h,g,d,s
     read_tipsy_arr,halodat[i].file + '.iord',h,iord,part = 'star',type = 'long'
     read_tipsy_arr,halodat[i].file + '.FeMassFrac',h,Fe,part = 'star',type = 'float'
     read_tipsy_arr,halodat[i].file + '.OxMassFrac',h,Ox,part = 'star',type ='float'
;     read_tipsy_arr,halodat[i].file + '.massform',h,massform,part = 'star',type ='float'  
     read_tipsy_arr,halodat[i].file + '.amiga.grp',h,grp,part = 'star',type ='long' 
     s.x = s.x*units.lengthunit*scale
     s.y = s.y*units.lengthunit*scale
     s.z = s.z*units.lengthunit*scale
     s.tform = s.tform*units.timeunit

;-------------------- Convert the position information gpart to the frame of the halo
     az = reform(a[i,*])
     az0 = az[0]
     az1 = az[1]
     az2 = az[2]
     ax = [az2/sqrt(az0*az0 + az2*az2),0,-1.0*az0/sqrt(az0*az0 + az2*az2)]
     ay = crossp(az,ax)         ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
     basis = [[ax],[ay],[az]]
     spos = [[s.x - center[i,2]*scale],[s.y - center[i,3]*scale],[s.z - center[i,4]*scale]]
     spos = transpose(transpose(basis)#transpose(spos))
     s.x = spos[*,0]
     s.y = spos[*,1]
     s.z = spos[*,2]

;---Find all stars that are members of the halo (and, if desired, its satellites)
     inhalo = [0]
     IF keyword_set(centralhalo) THEN BEGIN
        test = where(grp[i] EQ stat[main].group, ntest)
         IF ntest NE 0 THEN inhalo = [inhalo,test]
     ENDIF ELSE BEGIN
        FOR j = 0, n_elements(satellites) - 1 DO BEGIN
           test = where(grp EQ stat[satellites[j]].group, ntest)
           IF ntest NE 0 THEN inhalo = [inhalo,test]
        ENDFOR
     ENDELSE
     IF n_elements(inhalo) NE 1 THEN BEGIN
        inhalo = inhalo[1:n_elements(inhalo) - 1]
        match,iord,iordlast,ind_step,ind_last
        IF n_elements(iord) EQ n_elements(ind_last) THEN massform_step = massform[ind_last] ELSE stop

        s = s[inhalo]
        massform_step = massform_step[inhalo]
        stop
        which_imf = 0
        z_snii = zsnovaii(tmin, tmax + dDelta, s, massform_step, which_imf)
        z_snia = zsnovaia(tmin, tmax + dDelta, s, massform_step)
        zmetals[i] = z_snii.zmassloss + z_snia.zmassloss
        Oxs[i] = z_snii.oxmassloss + z_snia.oxmassloss
        Fes[i] = z_snii.femassloss + z_snia.femassloss
        print,i,tmin,tmax,zmetals[i],Oxs[i],Fes[i]
        IF Oxs[i] LT 0 OR  Fes[i] LT 0 THEN stop
        IF NOT finite(Oxs[i]) OR NOT finite(Fes[i]) THEN STOP
;        stop
     ENDIF ELSE BEGIN
        stars = s[inhalo]
        print,'No Star Particles'
        zmetals[i] = 0
        Oxs[i] = 0
        Fes[i] = 0
        print,i,zmetals[i],Oxs[i],Fes[i]
        CONTINUE
     ENDELSE

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
  IF NOT keyword_set(nowrite) THEN writecol,'grp' + finalid + '.sfmetals.txt',zmetals,Oxs,Fes,format='(3E)'
END


PRO master_stellarmetal_history

  dir = '/nobackupp2/crchrist/MolecH/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
  finalid = '2'
  stellarmetal_history,dir = dir,finalid = finalid
  finalid = '1'
  stellarmetal_history,dir = dir,finalid = finalid  

  dir = '/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
  finalid = '1'
  stellarmetal_history,dir = dir,finalid = finalid,/nowrite

  dir = '/nobackupp2/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
  finalid = '2'
  stellarmetal_history,dir = dir,finalid = finalid
  finalid = '1'
  stellarmetal_history,dir = dir,finalid = finalid

  dir = '/nobackupp2/crchrist/MolecH/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
  finalid = '3'
  stellarmetal_history,dir = dir,finalid = finalid
  finalid = '2'
  stellarmetal_history,dir = dir,finalid = finalid
  finalid = '1'
  stellarmetal_history,dir = dir,finalid = finalid

  dir = '/nobackupp2/crchrist/MolecH/h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
  finalid = '2'
  stellarmetal_history,dir = dir,finalid = finalid
  finalid = '1'
  stellarmetal_history,dir = dir,finalid = finalid
END
