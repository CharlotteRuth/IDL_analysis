;Select the stars that form in the center of the main halo and then
;determine their eventual location

PRO follow_central_sf,filename,step1,step2,dir = dir,finalid = finalid,laststep = laststep,debug = debug,outplot = outplot

  IF KEYWORD_SET(outplot) THEN BEGIN
     fgcolor = 0 
     bgcolor = 255
     xsize = 18
     ysize = 12
  ENDIF ELSE BEGIN
     fgcolor = 255
     bgcolor = 0
     xsize = 800
     ysize = 500
  ENDELSE
  loadct,39
  formatplot,outplot = outplot

 ;-------------------- Read in simulation data for final step
  IF NOT keyword_set(finalid) THEN finalid = '1' 
  IF NOT keyword_set(laststep) THEN laststep = '512'
  IF NOT keyword_set(dir) THEN dir = '.'
  cd,dir 
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + '.iord',file_iord
  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep,file
;  spawn,'ls ' + dir + '/*' + laststep + '/*' + laststep + 'amiga.stat',amiga_file
  rtipsy,file,h,g,d,s
  read_tipsy_arr,file_iord,h,iord,part = 'star',type = 'long'
  amiga = read_stat_struc_amiga(file + '.amiga.stat')
  rvir = amiga[where(amiga.group EQ finalid)].rvir
  frac_rvir = 0.006;0.004,0.0058,0.0067
  spawn,'ls ' + dir + '/h*.param',pfile
  units = tipsyunits(pfile[0])

  halodat = mrdfits('grp'+finalid+'.alignment.fits',1)
  readcol,'/nobackupp8/crchrist/MolecH/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/h239.cosmo50cmb.3072g14HMbwK.log',dTime,z,E,T,U,Eth,Lx,Ly,Lz,WallTime,dWMax,dImax,dEMax,dMultiEff,/silent
  zfix = fltarr(n_elements(halodat))
  tfix = fltarr(n_elements(halodat))
  FOR i = 0,n_elements(halodat) -1 DO BEGIN
     temp = min(abs(z - halodat[i].z),indz)
     zfix[i] = z[indz]
     tfix[i] = dtime[indz]
  END
  halodat.time = tfix           ;halodat.time*1e9/units.timeunit
  halodat.z = zfix

  IF max(halodat.time) GT 13 THEN halodat.time = halodat.time*1e9/units.timeunit
  a = 1/(1 + halodat.z)
;  halodat.xc = halodat.xc;*a ;scale center position
;  halodat.yc = halodat.yc;*a
;  halodat.zc = halodat.zc;*a
  istep1 = where(halodat.file EQ filename + '.' + step1 + '/' + filename + '.' + step1)
  time1 = halodat[istep1].time
  istep2 = where(halodat.file EQ filename + '.' + step2 + '/' + filename + '.' + step2)
  time2 = halodat[istep2].time

  s.x = s.x - halodat[n_elements(halodat) - 1].xc
  s.y = s.y - halodat[n_elements(halodat) - 1].yc
  s.z = s.z - halodat[n_elements(halodat) - 1].zc
  sr = sqrt(s.x*s.x + s.y*s.y + s.z*s.z)
  indgal = where(sr*units.lengthunit LT halodat[n_elements(halodat) - 1].rvir)
  s = s[indgal]
  sr = sr[indgal]
  iord = iord[indgal]

  IF file_test(filename + '_2merge.starlog.fits') THEN sl = mrdfits(filename + '_2merge.starlog.fits',1) ELSE $
     IF file_test(filename + '.starlog.fits') THEN sl = mrdfits(filename + '.starlog.fits',1) ELSE $
        sl = rstarlog(filename + '.starlog',/molecularH)

  match,iord,sl.iorderstar,ind1,ind2
  sl = sl[ind2]
  timeform_uniq = sl[uniq(sl.timeform,sort(sl.timeform))].timeform
  xcfit = spline(halodat.time,halodat.xc,timeform_uniq)
  ycfit = spline(halodat.time,halodat.yc,timeform_uniq)
  zcfit = spline(halodat.time,halodat.zc,timeform_uniq)
  axfit = spline(halodat.time,halodat.xa,timeform_uniq)
  ayfit = spline(halodat.time,halodat.ya,timeform_uniq)
  azfit = spline(halodat.time,halodat.za,timeform_uniq)
  afit = spline(halodat.time,a,timeform_uniq)
  rvirfit = spline(halodat.time,halodat.rvir,timeform_uniq)

  match2,sl.timeform,timeform_uniq,ind1,ind2
  slx = (sl.x - xcfit[ind1]);*afit[ind1])
  sly = (sl.y - ycfit[ind1]);*afit[ind1])
  slz = (sl.z - zcfit[ind1]);*afit[ind1])
  slx_cm = slx
  sly_cm = sly
  slz_cm = slz
  slx = slx*afit[ind1]
  sly = sly*afit[ind1]
  slz = slz*afit[ind1]

  s.x = s.x*h.time
  s.y = s.y*h.time
  s.z = s.z*h.time
  
  IF keyword_set(realign) THEN BEGIN
     pos = [[slx],[sly],[slz]]
     pos_prime = pos
     FOR i = 0, n_elements(timeform_uniq) - 1 DO BEGIN
        az_prime = [axfit[i],ayfit[i],azfit[i]]
        mag = sqrt(total(az_prime*az_prime))
        az_prime = az_prime/mag
        ax_prime = [az_prime[2]/sqrt(az_prime[0]*az_prime[0] + az_prime[2]*az_prime[2]),0,-1.0*az_prime[0]/sqrt(az_prime[0]*az_prime[0] + az_prime[2]*az_prime[2])]
        ay_prime = crossp(az_prime,ax_prime) 
        basis = [[ax_prime],[ay_prime],[az_prime]]
        ind_time = where(timeform_uniq[i] EQ sl.timeform)
        pos_prime[ind_time,*] = pos[ind_time,*]#basis
     ENDFOR
     slx = pos_prime[*,0]
     sly = pos_prime[*,1]
     slz = pos_prime[*,2]
  ENDIF

  rvirfit = rvirfit[ind1]
  rvirfit[where(rvirfit LT 0)] = 0
  slr = sqrt(slx*slx + sly*sly + slz*slz)
  indcenter = where(slr*units.lengthunit LT rvir*frac_rvir)
  slr_cm = sqrt(slx_cm*slx_cm + sly_cm*sly_cm + slz_cm*slz_cm)
  
;indicies of stars that formed within the 0.2*virial radius of the main halo
  insitu_scale = 0.2
  IF keyword_set(realign) THEN $
     slis_ind = where(sqrt(slx*slx + sly*sly)*units.lengthunit - 0.1*rvirfit LT 0 AND abs(slz)*units.lengthunit LT 0.75) ELSE $
        slis_ind = where(slr*units.lengthunit LT insitu_scale*rvirfit)
  sleps_ind = where(slr_cm LT min(s.eps))
  loadct,39
  yr_scale = 0.03

  IF keyword_set(debug) THEN BEGIN
;  IF keyword_set(outplot) THEN  device,filename = outplot + '_sfrad.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize 
     plot,sl.timeform*units.timeunit/1e9,slr*units.lengthunit,psym = 3
     oplot,sl.timeform*units.timeunit/1e9,slr*units.lengthunit,psym = 3,color = 100
     oplot,sl[slis_ind].timeform*units.timeunit/1e9,slr[slis_ind]*units.lengthunit,psym = 3,color = 60
     oplot,sl.timeform*units.timeunit/1e9,rvirfit*insitu_scale,thick = 1,color = 254
     oplot,sl.timeform*units.timeunit/1e9,rvirfit,thick = 1,color = 254,linestyle = 2
;  IF keyword_set(outplot) THEN device,/close
     stop
  ENDIF

  IF keyword_set(debug) THEN BEGIN
     window,0
     plot,sl.timeform*units.timeunit/1e9,sl.x*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_x',title = filename + ' ' + finalid,/nodata
     oplot,halodat.time/1e9*units.timeunit,halodat.xc*units.lengthunit,psym= 4,color = 254
     oplot,timeform_uniq*units.timeunit/1e9,xcfit*units.lengthunit,color = 254
     oplot,sl.timeform*units.timeunit/1e9,sl.x*units.lengthunit,psym = 3
     oplot,sl[slis_ind].timeform*units.timeunit/1e9,sl[slis_ind].x*units.lengthunit,psym = 3,color = 60
     oplot,sl.timeform*units.timeunit/1e9,xcfit[ind1]*units.lengthunit + rvirfit,thick = 1
     oplot,sl.timeform*units.timeunit/1e9,xcfit[ind1]*units.lengthunit - rvirfit,thick = 1
  
     window,2
     plot,sl.timeform*units.timeunit/1e9,slx*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_x',title = filename + ' ' + finalid,yrange = yr_scale*[-1*rvir,rvir]
     oplot,sl[slis_ind].timeform*units.timeunit/1e9,slx[slis_ind]*units.lengthunit,psym = 3,color = 60
     oplot,halodat.time/1e9*units.timeunit,0*halodat.time,psym = 4,color = 254
     oplot,[0,14],[0,0],color = 254
     oplot,[0,14],[rvir*frac_rvir,rvir*frac_rvir],color = 254
     oplot,[0,14],-1*[rvir*frac_rvir,rvir*frac_rvir],color = 254
;     oplot,sl.timeform*units.timeunit/1e9,slx*units.lengthunit,psym = 3
     stop

     window,0
     plot,sl.timeform*units.timeunit/1e9,sl.y*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_y',title = filename + ' ' + finalid,/nodata
     oplot,halodat.time/1e9*units.timeunit,halodat.yc*units.lengthunit,psym= 4,color = 254
     oplot,timeform_uniq*units.timeunit/1e9,ycfit*units.lengthunit,color = 254
     oplot,sl.timeform*units.timeunit/1e9,sl.y*units.lengthunit,psym = 3
  
     window,2
     plot,sl.timeform*units.timeunit/1e9,sly*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_y',title = filename + ' ' + finalid,yrange = yr_scale*[-1*rvir,rvir]
     oplot,sl[slis_ind].timeform*units.timeunit/1e9,sly[slis_ind]*units.lengthunit,psym = 3,color = 60
     oplot,halodat.time/1e9*units.timeunit,0*halodat.time,psym = 4,color = 254
     oplot,[0,14],[0,0],color = 254
     oplot,[0,14],[rvir*frac_rvir,rvir*frac_rvir],color = 254
     oplot,[0,14],-1*[rvir*frac_rvir,rvir*frac_rvir],color = 254
;     oplot,sl.timeform*units.timeunit/1e9,sly*units.lengthunit,psym = 3
     stop

     window,0
     plot,sl.timeform*units.timeunit/1e9,sl.z*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_z',title = filename + ' ' + finalid,/nodata
     oplot,halodat.time/1e9*units.timeunit,halodat.zc*units.lengthunit,psym= 4,color = 254
     oplot,timeform_uniq*units.timeunit/1e9,zcfit*units.lengthunit,color = 254
     oplot,sl.timeform*units.timeunit/1e9,sl.z*units.lengthunit,psym = 3
  
     window,2
     plot,sl.timeform*units.timeunit/1e9,slz*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r_z',title = filename + ' ' + finalid,yrange = yr_scale*[-1*rvir,rvir]
     oplot,sl[slis_ind].timeform*units.timeunit/1e9,slz[slis_ind]*units.lengthunit,psym = 3,color = 60
     oplot,halodat.time/1e9*units.timeunit,0*halodat.time,psym = 4,color = 254
     oplot,[0,14],[0,0],color = 254
     oplot,[0,14],[rvir*frac_rvir,rvir*frac_rvir],color = 254
     oplot,[0,14],-1*[rvir*frac_rvir,rvir*frac_rvir],color = 254
;     oplot,sl.timeform*units.timeunit/1e9,slz*units.lengthunit,psym = 3
     stop
  ENDIF
  slcenter = sl[indcenter]
  match,slcenter.iorderstar,iord,ind1,ind2
  indcenter_z0 = where(sr*units.lengthunit LT rvir*frac_rvir)
  iordcenter_z0 = iord[indcenter_z0]
  match,sl.iorderstar,iordcenter_z0,ind3,ind4
  match2,sl[sleps_ind].iorderstar,iord,ind5,sleps_ind_z0

  slis = sl[slis_ind]
  slris = slr[slis_ind]
  dt = 0.5 ;in Gyr
  nbins = 14/dt
  tbins = findgen(nbins)*dt + dt
  medr_bins = tbins*0
  FOR i = 0, n_elements(tbins) - 1 DO BEGIN
     indbin = where(slis.timeform GT (tbins[i] - dt/2)*1e9/units.timeunit AND slis.timeform LT (tbins[i] + dt/2)*1e9/units.timeunit)
     IF indbin[0] NE -1 THEN medr_bins[i] = total(slris[indbin]*slis[indbin].massform)/total(slis[indbin].massform) $
     ELSE medr_bins[i] = 0
;     indbin = where(slris*units.lengthunit LT rvir*frac_rvir AND slis.timeform GT tbins - dt/2 AND slis.timeform LT tbins + dt/2)
  ENDFOR

  IF keyword_set(outplot) THEN  device,filename = outplot + '_sfrad.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
  plot,slis.timeform*units.timeunit/1e9,slris*units.lengthunit,psym = 3,xtitle = 'Time [Gyr]',ytitle = 'r',title = filename + ' ' + finalid,yrange = yr_scale*[0,rvir]
  oplot,sl[sleps_ind].timeform*units.timeunit/1e9,slr[sleps_ind]*units.lengthunit,psym = 3,color = 130
  oplot,[0,14],[rvir*frac_rvir,rvir*frac_rvir],color = 254
  oplot,[time1,time1]*units.timeunit/1e9,[0,yr_scale*rvir],color = 254,linestyle = 2
  oplot,[time2,time2]*units.timeunit/1e9,[0,yr_scale*rvir],color = 254,linestyle = 2
  oplot,tbins,medr_bins*units.lengthunit,psym = -4,color = 60
  IF keyword_set(outplot) THEN device,/close

  IF keyword_set(outplot) THEN  device,filename = outplot + '_prof_sf_z0.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 0, xsize = xsize, ysize = ysize
  histogramp,slr*units.lengthunit,nbins = 100,max = yr_scale*rvir,xtitle = 'Radius [kpc]',title = filename + ' ' + finalid,yrange = [0,max([histogram(slr*units.lengthunit,nbins = 100,max = yr_scale*rvir),histogram(sr*units.lengthunit,nbins = 100,max = yr_scale*rvir)])]
  histogramp,slr[indcenter]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 254
;  histogramp,slr[sleps_ind]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 130
  oplot,[rvir*frac_rvir,rvir*frac_rvir],[0,1e9]
  histogramp,sr*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,linestyle = 2
  histogramp,sr[ind2]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 254,linestyle = 2
  histogramp,sr[indcenter_z0]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,linestyle = 2,color = 60
  histogramp,slr[ind3]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 60
;  histogramp,sr[sleps_ind_z0]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 130,linestyle = 2
  legend,['SF','Central SF','Stars, z = 0','Central SF, z = 0','Central stars, z = 0','Central stars, formation'],color = [0,254,0,254,60,60],linestyle = [0,0,2,2,2,0],/right,box = 0
  IF keyword_set(outplot) THEN device,/close 

  IF keyword_set(outplot) THEN  device,filename = outplot + '_prof_sf_z0_eps.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 1, xsize = xsize, ysize = ysize
  histogramp,slr*units.lengthunit,nbins = 100,max = yr_scale*rvir,xtitle = 'Radius [kpc]',title = filename + ' ' + finalid,yrange = [0,max([histogram(slr*units.lengthunit,nbins = 100,max = yr_scale*rvir),histogram(sr*units.lengthunit,nbins = 100,max = yr_scale*rvir)])],linestyle = 2
  histogramp,slr[sleps_ind]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 130,linestyle = 2
  histogramp,sr*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot
  histogramp,sr[where(sleps_ind_z0 NE -1)]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 130
  histogramp,sr[where(sleps_ind_z0 EQ -1)]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 60
  legend,['SF','SF within softening','Stars, z = 0','Stars within softening, z = 0','Stars outside softening, z = 0'],color = [0,130,0,130,60],linestyle = [0,0,2,2,2],/right,box = 0
  IF keyword_set(outplot) THEN device,/close ELSE stop
 
  IF 0 THEN  BEGIN
     window,3
     plot,s.x*units.lengthunit,s.z*units.lengthunit,psym = 3,xrange = yr_scale*[-1*rvir,rvir],yrange = yr_scale*[-1*rvir,rvir]
     oplot,s[ind2].x*units.lengthunit,s[ind2].z*units.lengthunit,psym = 3,color = 60
     oplot,slx[indcenter]*units.lengthunit,slz[indcenter]*units.lengthunit,psym = 3,color = 254

     window,1
     plot,slcenter.timeform*units.timeunit/1e9,(sr[ind2] - slr[indcenter])*units.lengthunit,psym = 3,yrange = [0,10]
     stop
  ENDIF

;-------------------------------- Comparing before and merger--------
;istep1
  spawn,'ls ' + dir + '/*' + step1 + '/*' + step1 + '.iord',file_iord1
  spawn,'ls ' + dir + '/*' + step1 + '/*' + step1,file1
  rtipsy,file1,h1,g1,d1,s1
  read_tipsy_arr,file_iord1,h1,iord1,part = 'star',type = 'long'

;istep2
  spawn,'ls ' + dir + '/*' + step2 + '/*' + step2 + '.iord',file_iord2
  spawn,'ls ' + dir + '/*' + step2 + '/*' + step2,file2
  rtipsy,file2,h2,g2,d2,s2
  read_tipsy_arr,file_iord2,h2,iord2,part = 'star',type = 'long'

  s1.x = (s1.x - halodat[istep1].xc)*h1.time
  s1.y = (s1.y - halodat[istep1].yc)*h1.time
  s1.z = (s1.z - halodat[istep1].zc)*h1.time
  sr1 = sqrt(s1.x*s1.x + s1.y*s1.y + s1.z*s1.z)*h1.time
  match2,iord1,slis.iorderstar,iord1_is,ind1
  s1is = s1[where(iord1_is NE -1)];stars formed in situ
  sr1is = sqrt(s1is.x^2 + s1is.y^2 + s1is.z^2)*h1.time
  s1es = s1[where(iord1_is EQ -1)];stars formed ex situ
  sr1es = sqrt(s1es.x^2 + s1es.y^2 + s1es.z^2)*h1.time

  s2.x = (s2.x - halodat[istep2].xc)*h2.time
  s2.y = (s2.y - halodat[istep2].yc)*h2.time
  s2.z = (s2.z - halodat[istep2].zc)*h2.time
  sr2 = sqrt(s2.x*s2.x + s2.y*s2.y + s2.z*s2.z)*h2.time
  match2,iord2,slis.iorderstar,iord2_is,ind1
  s2is = s2[where(iord2_is NE -1)];stars formed in situ
  sr2is = sqrt(s2is.x^2 + s2is.y^2 + s2is.z^2)*h2.time
  s2es = s2[where(iord2_is EQ -1)];stars formed ex situ
  sr2es = sqrt(s2es.x^2 + s2es.y^2 + s2es.z^2)*h2.time

  indcenter1 = where(sr1*units.lengthunit LT rvir*frac_rvir)
  indcenter2 = where(sr2*units.lengthunit LT rvir*frac_rvir)
  match,iord1[indcenter1],iord2,ind1,ind2
  indnew = where(s2.tform GT max(s1.tform),complement = indold)
  indisnew = where(s2is.tform GT max(s1.tform),complement = indisold)
  indesnew = where(s2es.tform GT max(s1.tform),complement = indesold)  

  IF keyword_set(outplot) THEN  device,filename = outplot + '_prof_overmerger.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 2, xsize = xsize, ysize = ysize
  histogramp,sr2*units.lengthunit,nbins = 100,max = yr_scale*rvir,xtitle = 'Radius [kpc]',tfile = filename + ' ' + finalid,yrange = [0,max([histogram(sr2*units.lengthunit,nbins = 100,max = yr_scale*rvir),histogram(sr1*units.lengthunit,nbins = 100,max = yr_scale*rvir)])],title = filename + ' ' + finalid
  histogramp,sr1*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 60;Profile before merger
  histogramp,sr1[indcenter1]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 60,linestyle = 1;Profile at center before merger
  histogramp,sr2[ind2]*units.lengthunit,nbins = 100,max = yr_scale*rvir,color = 100,/overplot;Central stars post merger
  histogramp,sr2[indnew]*units.lengthunit,nbins = 100,max = yr_scale*rvir,color = 254,/overplot
  histogramp,sr2is[indisnew]*units.lengthunit,nbins = 100,max = yr_scale*rvir,color = 254,/overplot,linestyle = 2
  IF min(sr2es[indesnew]*units.lengthunit) LT yr_scale*rvir THEN $
     histogramp,sr2es[indesnew]*units.lengthunit,nbins = 100,max = yr_scale*rvir,color = 254,/overplot,linestyle = 1
  histogramp,sr2[indold]*units.lengthunit,nbins = 100,max = yr_scale*rvir,color = 130,/overplot
  histogramp,sr2is[indisold]*units.lengthunit,nbins = 100,max = yr_scale*rvir,color = 130,/overplot,linestyle = 2
  histogramp,sr2es[indesold]*units.lengthunit,nbins = 100,max = yr_scale*rvir,color = 130,/overplot,linestyle = 1
  legend,['Profile before merger','Profile at center before merger','Profile post merger','Central stars post merger','Stars formed during merger','Stars formed before merger'],linestyle = [0,1,0,0,0,0],color = [60,60,0,100,254,130],/right,box = 0
  IF keyword_set(outplot) THEN device,/close

  match2,iord,slis.iorderstar,iord_is,ind1
  sis = s[where(iord_is NE -1)];stars formed in situ
  sris = sqrt(sis.x^2 + sis.y^2 + sis.z^2)
  ses = s[where(iord_is EQ -1)];stars formed ex situ
  sres = sqrt(ses.x^2 + ses.y^2 + ses.z^2)
  indpremerger = where(sis.tform LT max(s1.tform))
  indmerger = where(sis.tform GT max(s1.tform) AND sis.tform LT max(s2.tform))
  indpostmerger = where(sis.tform GT max(s2.tform))

  IF keyword_set(outplot) THEN  device,filename = outplot + '_prof_z0.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window, 0, xsize = xsize, ysize = ysize
  histogramp,sr*units.lengthunit,nbins = 100,max = yr_scale*rvir,xtitle = 'Radius [kpc]',tfile = filename + ' ' + finalid,yrange = [0,max(histogram(sr*units.lengthunit,nbins = 100,max = yr_scale*rvir))],title = filename + ' ' + finalid
  histogramp,sris*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot;in situ stars
  histogramp,sres*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 254 ;ex situ stars
  histogramp,sris[indpremerger]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 100
  histogramp,sris[indmerger]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 130
  histogramp,sris[indpostmerger]*units.lengthunit,nbins = 100,max = yr_scale*rvir,/overplot,color = 190
  legend,['All stars','All in situ stars','In situ stars pre merger','In situ stars during merger','In situ stars post merger','Ex situ stars'],linestyle = [0,0,0,0,0,0],color = [0,0,100,130,190,254],/right,box = 0
  IF keyword_set(outplot) THEN device,/close ELSE stop

END
