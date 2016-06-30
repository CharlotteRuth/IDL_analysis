PRO comp_stellar_dist,dir,file,finalid,firststep,haloid_fst,laststep,haloid_lst,outplot = outplot,multiplot = multiplot
  cd,dir
  formatplot,outplot = outplot
  spawn,'ls ' + dir + '/h*.param',pfile
  units = tipsyunits(pfile[0])
  halodat = mrdfits('grp' + finalid + '.alignment.fits',1)
  halodat.xc = halodat.xc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.yc = halodat.yc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  halodat.zc = halodat.zc*units.lengthunit/1000.0 + units.lengthunit/2/1000.0
  center = [[halodat.time],[halodat.z],[halodat.xc],[halodat.yc],[halodat.zc]]
  center[*,2] = center[*,2]*1000.0 - units.lengthunit/2.0 ;"center" is in comoving Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
  center[*,3] = center[*,3]*1000.0 - units.lengthunit/2.0
  center[*,4] = center[*,4]*1000.0 - units.lengthunit/2.0
  a = [[halodat.xa],[halodat.ya],[halodat.za]]
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord_lst
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file_lst
  spawn,'ls ' + dir + '/*' + firststep + '/*' + firststep + '.iord',file_iord_fst
  spawn,'ls ' + dir + '/*' + firststep + '/*' + firststep,file_fst

  rtipsy,file_fst,h0,g0,d0,s0
  read_tipsy_arr,file_fst + '.amiga.grp',h0,grp0
  read_tipsy_arr,file_fst + '.iord',h0,iord0
  ind_fst = -1
  FOR i = 0, n_elements(halodat) - 1 DO $
     IF dir + '/' + halodat[i].file EQ file_fst THEN ind_fst = i
  scale0 = h0.time
  grp0s = grp0[h0.ngas + h0.ndark:h0.n - 1]
  s0.x = (s0.x*units.lengthunit - center[ind_fst,2])*scale0
  s0.y = (s0.y*units.lengthunit - center[ind_fst,3])*scale0
  s0.z = (s0.z*units.lengthunit - center[ind_fst,4])*scale0
;  s0r = sqrt(s0.x*s0.x + s0.y*s0.y + s0.z*s0.z)
;  ind_s0 = where(s0r LE halodat[ind_fst].rvir)
  ind_s0 = where(grp0s EQ haloid_fst)
  s0r = sqrt(s0[ind_s0].x*s0[ind_s0].x + s0[ind_s0].y*s0[ind_s0].y + s0[ind_s0].z*s0[ind_s0].z)
  iord_s0 = (iord0[h0.ngas + h0.ndark:h0.n - 1])[ind_s0]

  rtipsy,file_lst,h1,g1,d1,s1
  read_tipsy_arr,file_lst + '.amiga.grp',h1,grp1
  read_tipsy_arr,file_lst + '.iord',h1,iord1
  ind_lst = -1
  FOR i = 0, n_elements(halodat) - 1 DO $
     IF dir + '/' + halodat[i].file EQ file_lst THEN ind_lst = i
  scale1 = h1.time
  grp1s = grp1[h1.ngas + h1.ndark:h1.n - 1]
  s1.x = (s1.x*units.lengthunit - center[ind_lst,2])*scale1
  s1.y = (s1.y*units.lengthunit - center[ind_lst,3])*scale1
  s1.z = (s1.z*units.lengthunit - center[ind_lst,4])*scale1
;  s1r = sqrt(s1.x*s1.x + s1.y*s1.y + s1.z*s1.z)
;  ind_s1 = where(s1r LE halodat[ind_lst].rvir)
  ind_s1 = where(grp1s EQ haloid_lst)
  s1r = sqrt(s1[ind_s1].x*s1[ind_s1].x + s1[ind_s1].y*s1[ind_s1].y + s1[ind_s1].z*s1[ind_s1].z)
  iord_s1 = (iord1[h1.ngas + h1.ndark:h1.n - 1])[ind_s1]
  match,iord_s0,iord_s1,inda,indb

IF keyword_set(multiplot) THEN BEGIN
   IF NOT keyword_set(outplot) THEN BEGIN
      xsize = 580
      ysize = 800
   ENDIF ELSE BEGIN
      xsize = 18
      ysize = 40
   ENDELSE
ENDIF ELSE BEGIN
   IF NOT keyword_set(outplot) THEN BEGIN
      xsize = 400
      ysize = 300
   ENDIF ELSE BEGIN
      xsize = 18
      ysize = 12
   ENDELSE
ENDELSE
IF (keyword_set(outplot)) THEN device,filename = file + finalid + '_stellar_redist0.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,0,xsize = xsize,ysize = ysize
;  multiplot,[1,3],mTitle = file + '.' + finalid
  plot,s0[ind_s0].x,s0[ind_s0].y,psym = 3
  oplot,s1[ind_s1[indb]].x,s1[ind_s1[indb]].y,psym = 3,color = 150
;  multiplot
IF (keyword_set(outplot)) THEN device,/close ELSE stop

IF (keyword_set(outplot)) THEN device,filename = file + finalid + '_stellar_redist1.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,0,xsize = xsize,ysize = ysize
  plot,s0r[inda],s1r[indb],psym = 3,/ylog,/xlog,yrange = [0.01,100],xrange = [0.01,100]
  oplot,[0.01,100],[0.01,100]
;  multiplot
IF (keyword_set(outplot)) THEN device,/close ELSE stop

IF (keyword_set(outplot)) THEN device,filename = file + finalid + '_stellar_redist2.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,0,xsize = xsize,ysize = ysize
  histogramp,s1r,nbins = 200,min = 0.00001,max = 20,thick = 4,xtitle = 'Radius'
  histogramp,s1r[indb],nbins = 200,color = 254,/overplot,min = 0.00001,max = 20,thick = 4
  histogramp,s0r[inda],nbins = 200,color = 60,/overplot,min = 0.00001,max = 20,thick = 4
  legend,['Pre-merger','Post-merger, redist','Post-merger, all'],color = [60,254,0],linestyle = [0,0,0],/right,thick = [4,4,4]
;  multiplot,/reset
  IF (keyword_set(outplot)) THEN device,/close ELSE stop

END
