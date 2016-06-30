;Charlotte Christensen

PRO central_mass_history,dir = dir,finalid = finalid,laststep = laststep,centralhalo = centralhalo,debug = debug,outplot = outplot
  cd,dir
  formatplot,outplot = outplot
;-------------------- Read in simulation data for final step
  IF NOT keyword_set(finalid) THEN finalid = '1' 
  IF NOT keyword_set(laststep) THEN laststep = '512'
  IF NOT keyword_set(dir) THEN dir = '.'
  cd,dir
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file
;  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + 'amiga.stat',amiga_file
  amiga = read_stat_struc_amiga(file + '.amiga.stat')
  rvir = amiga[where(amiga.group EQ finalid)].rvir
  frac_rvir = 0.006;0.004,0.0058,0.0067
  spawn,'ls ' + dir + '/h*.param',pfile
  units = tipsyunits(pfile[0])
  
;-------------------- Read in the information on the halos
  halodat = mrdfits('grp' + finalid + '.alignment.fits',1)
  nsteps = n_elements(halodat.file)
  halodat.xc = halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.yc = halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.zc = halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  center = [[halodat.time],[halodat.z],[halodat.xc],[halodat.yc],[halodat.zc]]
  center[*,2] = center[*,2]*1000.0 - units.lengthunit/2.0 ;"center" is in comoving Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
  center[*,3] = center[*,3]*1000.0 - units.lengthunit/2.0
  center[*,4] = center[*,4]*1000.0 - units.lengthunit/2.0
  a = [[halodat.xa],[halodat.ya],[halodat.za]]
  gmass = fltarr(nsteps)
  smass = fltarr(nsteps)
  dmass = fltarr(nsteps)
  gmass1 = fltarr(nsteps)
  smass1 = fltarr(nsteps)
  dmass1 = fltarr(nsteps)
  slope = fltarr(nsteps)
  slope1 = fltarr(nsteps)
;-------------------- Step through outputs to find the location of the particles
  FOR i = 0, nsteps - 1 DO BEGIN
     print,halodat[i].file
     z = center[i,1] 
     scale = 1.0/(1.0 + z)

;-------------------- Read Simulation Data ------------
     rtipsy,halodat[i].file,h,g,d,s
     read_tipsy_arr,halodat[i].file + '.amiga.grp',h,grp
     
     grpg = grp[0:h.ngas - 1]
     grpd = grp[h.ngas:h.ngas + h.ndark - 1]
     grps = grp[h.ngas + h.ndark:h.n - 1]
     g.x = g.x*units.lengthunit - center[i,2]
     g.y = g.y*units.lengthunit - center[i,3]
     g.z = g.z*units.lengthunit - center[i,4]
     s.x = s.x*units.lengthunit - center[i,2]
     s.y = s.y*units.lengthunit - center[i,3]
     s.z = s.z*units.lengthunit - center[i,4]
     d.x = d.x*units.lengthunit - center[i,2]
     d.y = d.y*units.lengthunit - center[i,3]
     d.z = d.z*units.lengthunit - center[i,4]
     gr = sqrt(g.x*g.x + g.y*g.y + g.z*g.z)
     sr = sqrt(s.x*s.x + s.y*s.y + s.z*s.z)
     dr = sqrt(d.x*d.x + d.y*d.y + d.z*d.z)
     gind = where(scale*gr LE 0.5)
     sind = where(scale*sr LE 0.5)
     dind = where(scale*dr LE 0.5)
     gmass[i] = total(g[gind].mass)*units.massunit
     smass[i] = total(s[sind].mass)*units.massunit
     dmass[i] = total(d[dind].mass)*units.massunit

     IF 1 THEN BEGIN
        amiga = read_stat_struc_amiga(file + '.amiga.stat')
     IF keyword_set(outplot) THEN device,filename = halodat[i].file + '.' + finalid + '_pic.eps',/color,bits_per_pixel= 8,xsize = 12,ysize = 12,xoffset =  2,yoffset =  2 ELSE window,2,xsize = 400,ysize = 400
     plot,g.x,g.y,psym = 3,xrange = [-1*amiga[where(amiga.group EQ halodat[i].haloid)].rvir/2,amiga[where(amiga.group EQ halodat[i].haloid)].rvir/2],yrange = [-1*amiga[where(amiga.group EQ halodat[i].haloid)/2].rvir,amiga[where(amiga.group EQ halodat[i].haloid)].rvir/2],/nodata
;     ind = where(grpg EQ halodat[i].haloid)
;     plot,g[ind].x*scale,g[ind].y*scale,psym = 3
     oplot,g.x,g.y,color = 254,psym = 3
     oplot,s.x,s.y,color = 100,psym = 3
     oplot,[0,0],[-1e6,1e6]
     oplot,[-1e6,1e6],[0,0]
;     oplot,g[gind].x*scale,g[gind].y*scale,psym = 3,color = 100
;     ind = where(grps EQ halodat[i].haloid)
;     oplot,s[ind].x*scale,s[ind].y*scale,psym = 3,color = 200
;     oplot,s[sind].x*scale,s[sind].y*scale,psym = 3,color = 254
     IF keyword_set(outplot) THEN device,/close
     ENDIF

     gind = where(scale*gr LE rvir*frac_rvir)
     sind = where(scale*sr LE rvir*frac_rvir)
     dind = where(scale*dr LE rvir*frac_rvir)
     gmass1[i] = total(g[gind].mass)*units.massunit
     smass1[i] = total(s[sind].mass)*units.massunit
     dmass1[i] = total(d[dind].mass)*units.massunit

     y = weighted_histogram(scale*dr,weight = g.mass,locations = x,min = 300,max = 700)
     y = y/x^2
     fits = robust_linefit( alog10(x), alog10(y), den_fit, sigma )
;     plot,x,y,/xlog,/ylog
;     oplot,x,10^fits[0]*x^fits[1],linestyle = 2
     slope[i] = fits[1]

     y = weighted_histogram(scale*dr,weight = g.mass,locations = x,min = 800,max = 1200)
     y = y/x^2
     fits = robust_linefit( alog10(x), alog10(y), den_fit, sigma )
     slope1[i] = fits[1]
     
  ENDFOR
  writecol,dir + '/grp' + finalid + '.mass500.txt',halodat.haloid,halodat.time,halodat.z,gmass,smass,dmass,slope,gmass1,smass1,dmass1,slope1,format='(I,F,F,E,E,E,F,E,E,E,F)'
stop
END
