pro checkSurfaceDen_obs_master 
;filename = 'h603.cosmo50cmb.2304g5bwK.00512'
;cd,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG'

cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir'
filename = 'h516.cosmo25cmb.3072g1MBWK.00492'
pfile = 'h516.cosmo25cmb.3072g1MBWK.param'
;outplot = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1_' + strtrim(res,2)
;schmidtlaw_res_obs,filename,pfile,res = res,outplot = outplot

cd,'/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse_1e6/Disk_Collapse_1e6_sol_H2'
filename = 'Disk_Collapse.H2.solmet.ch0.00100.00020'
pfile = 'Disk_Collapse.H2.solmet.ch0.param'
;outplot = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1_' + strtrim(res,2)
checkSurfaceDen_obs,filename,pfile,res = 0.10;,outplot = outplot

end

pro checkSurfaceDen_obs,filename,pfile,res = res,outplot = outplot,overplot = overplot,useH2 = useH2
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
speed_o_light = 299792458 ;m per sec
molec_weight = (0.76*1 + 0.24*4.0)
units = tipsyunits(pfile)
massunit = units.massunit
lengthunit = units.lengthunit
convert2msol = 1/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18

!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps' 
    nbins=100.0
    linestyles = [0,2]
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25 
    !X.MARGIN = [12,5]
    !Y.MARGIN = [6,4]
    cb_charsize = 0.75
    black = 0
    white = 255
    fgcolor = 0
    l_charsize = 0.75
ENDIF ELSE BEGIN
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    set_plot,'x'
    nbins=100.0
    linestyles = [0,2]
    !X.MARGIN = [12,5]
    !Y.MARGIN = [6,4]
    cb_charsize = 0.75
    white = 255
    black = 0
    l_charsize = 1.0
    fgcolor = 255
ENDELSE
!p.multi = 0

;************** Obs Data**********************
readcol,'/astro/users/christensen/code/MolecH/Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
readcol,'/astro/users/christensen/code/MolecH/Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'

;*************** Tipsy Data ***********************
rtipsy,filename,h,g,d,s
g.x = g.x*units.lengthunit
g.y = g.y*units.lengthunit
g.z = g.z*units.lengthunit
g.mass = g.mass*massunit
g.dens = g.dens*units.rhounit/h.time/h.time/h.time
readarr,filename + '.HI',h,HI,/ascii,part = 'gas'
readarr,filename + '.H2',h,H2,/ascii,part = 'gas'
readarr,filename+'.lw',h,lw,/ascii,part = 'gas'
readarr,filename+'.smoothlength',h,hl,/ascii,part = 'gas'
hl = hl*units.lengthunit*cm_per_kpc
IF (FILE_TEST(filename+".correL")) THEN BEGIN
    readarr,filename+'.correL',h,correL,/ascii,part = 'gas'
    correL = correL*units.lengthunit*cm_per_kpc
ENDIF
totalH = HI + 2.0*H2
length = hl
H2 = H2*2


;*************** Rotation ****************
angle = 1.0*!PI/4.0
angle0 = ATAN(g.z/g.y)
angle0[where(g.y lt 0)] = angle0[where(g.y lt 0)]+!PI
g1 = g
Ryz = SQRT(g.y*g.y + g.z*g.z)
g1.y = Ryz*COS(angle0 + angle)
g1.z = Ryz*SIN(angle0 + angle)
g1.x = g1.x
radiusg1 = sqrt(g1.x*g1.x + g1.y*g1.y) 

;***************** Gas Cubes **********************
;cubeHI_surface_den = read_dencube_fits(filename + '.halo.1.cube.smoothed.mom0.fits',headerHI,/noscale);/amu_per_gm/gm_per_msol*cm_per_kpc*cm_per_kpc/1d6
;cubeHI_surface_den = read_dencube_fits(filename + '.HI.fo.fits',headerHI,/noscale);/amu_per_gm/gm_per_msol*cm_per_kpc*cm_per_kpc/1d6
;cubeHI_surface_den_prof = cubeHI_surface_den
;IF NOT (where(~finite(cubeHI_surface_den)))[0] EQ -1 THEN cubeHI_surface_den[where(~finite(cubeHI_surface_den))] = 0

cubeHI = read_dencube_fits(filename+'.HI.fo.fits',headerHI)
pix_areaHI = headerHI.CDELT1*headerHI.CDELT2*3.08568021d21*3.08568021d21
cubeHI_surface_den = TRANSPOSE(cubeHI)/pix_areaHI*gm_per_msol*amu_per_gm
vcubeHI = read_dencube_fits(filename+'.cubeHI.fits')/pix_areaHI*gm_per_msol*amu_per_gm

cubeH2 = read_dencube_fits(filename+'.H2.fo.fits',headerH2)
pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*3.08568021d21*3.08568021d21
cubeH2_surface_den = TRANSPOSE(cubeH2)/pix_areaH2*gm_per_msol*amu_per_gm*2.0
vcubeH2 = read_dencube_fits(filename+'.cubeH2.fits')/pix_areaHI*gm_per_msol*amu_per_gm

;Axes for the HI/H2 surface density
xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1 + headerHI.CDELT1/2
yaxes = (findgen(headerHI.NAXIS2) - headerHI.NAXIS1/2.0)*headerHI.CDELT2 + headerHI.CDELT1/2
dxcell = headerHI.cdelt1
areapc = 1e6*dxcell*dxcell
areakpc = dxcell*dxcell
minrange = headerHI.CRVAL1

;*************** Plotting **************
loadct,0
nlevels = 60
range = 8.0
;!p.multi = 0
;position = [0.42, 0.6, 0.45, 0.9]
position =  [0.85, 0.15, 0.9, 0.9]
IF (KEYWORD_SET(outplot)) THEN device,filename=outplot+'_H2sd.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 800,ysize = 800
loadct,0
mingas = 1e-2
maxgas = 1e2
userlevels = indgen(nlevels)*(alog10(maxgas) - alog10(mingas))/nlevels + alog10(mingas)
cubeH2_surface_den_plot = cubeH2_surface_den*convert2msol
ind = where(cubeH2_surface_den_plot lt mingas)
IF ind[0] NE -1 THEN cubeH2_surface_den_plot[ind] = mingas*1.01
;;cubeH2_surface_den_plot[where(NOT FINITE(cubeH2_surface_den_plot))] = 1e-2
contour,alog10(cubeH2_surface_den_plot),xaxes,yaxes,$
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),levels = userlevels ;,title = textoidl('H_2 Surface Density [M')+sunsymbol()+textoidl(' pc^{-2}]')
colorbar,range = [alog10(mingas),alog10(maxgas)],/vertical,divisions = 4,position = position,color = white
IF KEYWORD_SET(outplot) THEN device,/close

IF (KEYWORD_SET(outplot)) THEN device,filename=outplot+'_HIsd.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,1,xsize = 800,ysize = 800
mingas = 1e-2
maxgas = 1e2
userlevels = indgen(nlevels)*(alog10(maxgas) - alog10(mingas))/nlevels + alog10(mingas)
cubeHI_surface_den_plot = cubeHI_surface_den*convert2msol
ind = where(cubeHI_surface_den_plot lt mingas)
IF ind[0] NE -1 THEN cubeHI_surface_den_plot[ind] = mingas*1.01
contour,alog10(cubeHI_surface_den_plot),xaxes,yaxes,$
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),levels = userlevels ;,title = textoidl('HI Surface Density [M')+sunsymbol()+textoidl(' pc^{-2}]')
colorbar,range = [alog10(mingas),alog10(maxgas)],/vertical,divisions = 4,position = position,color = white
IF KEYWORD_SET(outplot) THEN device,/close


;************************** Profile **********************
minr = 0
maxr = range
range = maxr - minr
nbins = 40.0
r_bin = findgen(nbins + 1)*range/(nbins)
dr = range/nbins
area = ((r_bin + dr)*(r_bin + dr) - r_bin*r_bin)*!PI*1e6
prof_H2 = fltarr(nbins)
prof_HI = fltarr(nbins)
prof_H2_part = fltarr(nbins)
prof_HI_part = fltarr(nbins)
xaxes_grid = REBIN(xaxes + 0.0234375,N_ELEMENTS(xaxes),N_ELEMENTS(xaxes))
yaxes_grid = TRANSPOSE(REBIN(xaxes + 0.0234375,N_ELEMENTS(xaxes),N_ELEMENTS(xaxes)))
distance = SQRT((xaxes_grid)*(xaxes_grid) + (yaxes_grid)* (yaxes_grid))
FOR ir = 0, nbins - 1 DO BEGIN
    ind = where(distance le r_bin(ir + 1) AND distance gt r_bin(ir))
    prof_H2[ir] = MEAN(cubeH2_surface_den[ind]*convert2msol)*cos(angle) 
    prof_HI[ir] = MEAN(cubeHI_surface_den[ind]*convert2msol)*cos(angle) 
    indg = where(radiusg1 le r_bin(ir + 1) AND radiusg1 gt r_bin(ir))
    prof_H2_part[ir] = TOTAL(g1[indg].mass*H2[indg])/area[ir]*cos(angle) 
    prof_HI_part[ir] = TOTAL(g1[indg].mass*HI[indg])/area[ir]*cos(angle)
ENDFOR

formatplot,outplot = outplot
!p.multi = 0
IF KEYWORD_SET(outplot) THEN  device,filename = outplot + 'gasProf.eps' ELSE BEGIN
    window,2,xsize = 400,ysize = 400
    loadct,39
ENDELSE
plot,r_bin[0:nbins - 1] + dr/2.0,alog10(prof_HI),yrange = [-1,2],xrange = [0,range],xtitle = 'Radius [kpc]',ytitle = textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",psym = 10
oplot,r_bin[0:nbins - 1] + dr/2.0,alog10(prof_H2),linestyle = 1,psym = 10
oplot,r_bin[0:nbins - 1] + dr/2.0,alog10(prof_HI_part),color = 150,psym = 10
oplot,r_bin[0:nbins - 1] + dr/2.0,alog10(prof_H2_part),linestyle = 1,color = 150,psym = 10
legend,['HI',textoidl('H_2')],linestyle = [0,1],/right
IF KEYWORD_SET(outplot) THEN device,/close
stop

formatplot,outplot = outplot
!p.multi = 0
IF KEYWORD_SET(outplot) THEN  device,filename = outplot + 'gasProf.eps' ELSE BEGIN
    window,3,xsize = 400,ysize = 400
    loadct,39
ENDELSE
;indnozero = where(vcubeHI ne 0)
;vcubeHI = vcubeHI[indnozero]
;vcubeH2 = vcubeH2[indnozero]
indplot = where(vcubeH2 gt 1e19)
;indplot = where(vcubeH gt 1e19)
vcubeH = vcubeHI + 2.0*vcubeH2
plot,vcubeH[indplot],2.0*vcubeH2[indplot]/vcubeH[indplot],psym = 3,/xlog,/ylog,yrange=[1e-6,1],xtitle = textoidl('\Sigma_{HI + H_{2}} [amu/cm^2]'),ytitle = textoidl('H_2/H'),xrange = [1e19,1e23]
stop

min1 = 19
max1 = 23
min2 = -6
max2 = 0.4
bin1 = 0.1 ;(max1 - min1)/4.0/10.0
bin2 = 0.1 ;(max2 - min2)/6.0/10.0
nxbins=FLOOR((max1-min1) / bin1)
nybins=FLOOR((max2-min2) / bin2)
xbins=(FINDGEN(nxbins)*bin1)+min1+bin1/2.0
ybins=(FINDGEN(nybins)*bin2)+min2+bin2/2.0
histv=HIST_2D(alog10(vcubeH),alog10(2.0*vcubeH2/(vcubeH)),min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)
histd=HIST_2D(alog10((cubeH2_surface_den + cubeHI_surface_den)*cos(angle)),alog10(cubeH2_surface_den/(cubeH2_surface_den + cubeHI_surface_den)),min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)
indcool = where(g.tempg lt 1e4)
histp=HIST_2D(alog10(g[indcool].dens*length[indcool]),alog10(H2[indcool]/totalH[indcool]*2),min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)

loadct,0
device,filename = 'H2abund.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  
contour,histv,xbins,ybins,nlevels = 4,xtitle = textoidl('LOG(\Sigma_{HI + H_{2}}) [amu/cm^2]'),ytitle = textoidl('LOG(H_2/H)'),xrange = [19,23],yrange = [-6,0]
contour,histp,xbins,ybins,nlevels = 20,xtitle = textoidl('LOG(\Sigma_{HI + H_{2}}) [amu/cm^2]'),ytitle = textoidl('LOG(H_2/H)'),xrange = [19,23],yrange = [-6,0],color = 100,/overplot
contour,histv,xbins,ybins,nlevels = 4,xtitle = textoidl('LOG(\Sigma_{HI + H_{2}}) [amu/cm^2]'),ytitle = textoidl('LOG(H_2/H)'),xrange = [19,23],yrange = [-6,0],/overplot
;oplot,alog10(g.dens*length),alog10(H2/totalH*2),psym = 3
oplot,nhw,logfH2w,psym = 7,color = 50
oplot,nhg,logfH2g,psym = 1,color = 50

;plot,(cubeH2_surface_den + cubeHI_surface_den)*cos(angle),cubeH2_surface_den/(cubeH2_surface_den + cubeHI_surface_den),xtitle = textoidl('\Sigma_{HI + H_{2}} [amu/cm^2]'),ytitle = textoidl('H_2/H'),psym = 3,/ylog,/xlog,xrange = [1e19,1e23],yrange = [1e-6,1]
;oplot,10^nhw,10^logfH2w,psym = 4,color = 150
;oplot,10^nhg,10^logfH2g,psym = 6,color = 150
legend,['Simulated Observation','Particle Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [-3,-3,7,1],/right,/bottom,charsize = l_charsize,color = [fgcolor,150,50,50]
device,/close

stop

;**************** Rebinning and Spatial Axes
;Axes rebinned
IF NOT KEYWORD_SET(res) THEN res = 0.75 ;kpc
sr_range = abs(headerH2.CRVAL1*2.0)
ncell = ROUND(sr_range/res)
res = 1.0*sr_range/ncell
xaxes_lr = (findgen(ncell)*res - sr_range/2.0 + res/2.0)
yaxes_lr = (findgen(ncell)*res - sr_range/2.0 + res/2.0)
HIrebin = fltarr(ncell,ncell)
H2rebin = fltarr(ncell,ncell)
gas_sd = fltarr(ncell,ncell)

FOR ix = 0, ncell - 1 DO BEGIN
    FOR iy = 0,ncell - 1 DO BEGIN 
        indx = where((xaxes gt (ix*res + minrange)) AND (xaxes le ((ix+1)*res + minrange)))
        indy = where((yaxes gt (iy*res + minrange)) AND (yaxes le ((iy+1)*res + minrange)))
        IF ((indx[0] eq -1) OR (indy[0] eq -1)) THEN HIrebin[ix,iy] = 0 ELSE HIrebin[ix,iy] = MEAN((cubeHI_surface_den[indx,*])[*,indy])
        IF ((indx[0] eq -1) OR (indy[0] eq -1)) THEN H2rebin[ix,iy] = 0 ELSE H2rebin[ix,iy] = MEAN((cubeH2_surface_den[indx,*])[*,indy])
;        HIrebin[ix,iy] = MEAN(cubeHI_surface_den[ix*dcell : ix*dcell + dcell - 1, iy*dcell :iy*dcell + dcell - 1])
;        IF keyword_set( useH2) THEN H2rebin[ix,iy] = MEAN(cubeH2_surface_den[ix*dcell : ix*dcell + dcell - 1, iy*dcell :iy*dcell + dcell - 1])

        indg = where((g1.x gt (res*ix + minrange)) AND (g1.x le (res*(ix+1) + minrange)) AND (g1.y gt (res*iy + minrange)) AND (g1.y lt (res*(iy+1) + minrange)))
        IF indg[0] NE -1 THEN gas_sd[ix,iy] = TOTAL(g1[indg].mass*(HI[indg] + H2[indg]))/res/res/1000.0/1000.0 ELSE gas_sd[ix,iy] = 0
     ENDFOR
 ENDFOR


;********************* Plotting *****************************
mingas = 1e18
maxgas = 1e23
cubeHI_surface_den_cut = cubeHI_surface_den
ind = where(cubeHI_surface_den_cut le mingas)
IF ind[0] NE -1 THEN cubeHI_surface_den_cut[ind] = 1.01*mingas
HIrebin_cut = HIrebin
ind = where(HIrebin_cut le mingas)
IF ind[0] NE -1 THEN HIrebin_cut[ind] = 1.01*mingas
userlevels = indgen(nlevels)*(alog10(maxgas) - alog10(mingas))/nlevels + alog10(mingas)
cubeH2_surface_den_cut = cubeH2_surface_den
ind = where(cubeH2_surface_den_cut le mingas) 
IF ind[0] NE -1 THEN cubeH2_surface_den_cut[ind] = 1.01*mingas
H2rebin_cut = H2rebin
ind = where(H2rebin_cut le mingas)
IF ind[0] NE -1 THEN H2rebin_cut[ind] = 1.01*mingas
position = [0.415, 0.59, 0.44, 0.90]
!X.MARGIN = [4,1]
!Y.MARGIN = [2,1]
IF (KEYWORD_SET(outplot)) THEN device,filename=outplot+'_gas.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,2,xsize = 800,ysize = 800
multiplot,[2,2],gap=0.06,/square,/doxaxis,/doyaxis

contour,alog10(cubeHI_surface_den_cut + cubeH2_surface_den_cut),xaxes,yaxes,title = textoidl('Log \Sigma_{HI + H_2} [amu cm^{-2}]'), $
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels 
colorbar,range = [alog10(mingas),alog10(maxgas)],/vertical,divisions = 4,position = position,charsize=cb_charsize,color = white
multiplot,/doxaxis,/doyaxis

contour,alog10(HIrebin_cut + H2rebin_cut),xaxes_lr,yaxes_lr,title = textoidl('Log \Sigma_{HI + H_2} [amu cm^{-2}]'),$                
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels 
;colorbar,range = [19,23],/vertical,divisions = 3
multiplot,/doxaxis,/doyaxis
contour,alog10(cubeH2_surface_den_cut),xaxes,yaxes,title = textoidl('Log \Sigma_{H_2} [amu cm^{-2}]'), $                        
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels 
;colorbar,range = [12,23],/vertical,divisions = 3
multiplot,/doxaxis,/doyaxis
contour,alog10(cubeHI_surface_den_cut),xaxes,yaxes,title = textoidl('Log \Sigma_{HI} [amu cm^{-2}]'),$                         
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels 
;colorbar,range = [19,23],/vertical,divisions = 3
IF KEYWORD_SET(outplot) THEN device,/close
multiplot,/reset

finite = where(FINITE(H2rebin),complement = nfinite)
if ((nfinite[0]) ne -1) then H2rebin[nfinite] = 0
finite = where(FINITE(HIrebin),complement = nfinite)
if ((nfinite[0]) ne -1) then HIrebin[nfinite] = 0
finite = where(FINITE(gas_sd),complement = nfinite)
if ((nfinite[0]) ne -1) then gas_sd[nfinite] = 0

;************************** Profile **********************
minr = 0
maxr = range
range = maxr - minr
nbins = 40.0
r_bin = findgen(nbins + 1)*range/(nbins)
dr = range/nbins
area = ((r_bin + dr)*(r_bin + dr) - r_bin*r_bin)*1e6
prof_H2 = fltarr(nbins)
prof_HI = fltarr(nbins)
xaxes_grid = REBIN(xaxes + 0.0234375,N_ELEMENTS(xaxes),N_ELEMENTS(xaxes))
yaxes_grid = TRANSPOSE(REBIN(xaxes + 0.0234375,N_ELEMENTS(xaxes),N_ELEMENTS(xaxes)))
distance = SQRT((xaxes_grid)*(xaxes_grid) + (yaxes_grid)* (yaxes_grid))
FOR ir = 0, nbins - 1 DO BEGIN
    ind = where(distance le r_bin(ir + 1) AND distance gt r_bin(ir))
    prof_H2[ir] = MEAN(cubeH2_surface_den[ind]*convert2msol) 
    prof_HI[ir] = MEAN(cubeHI_surface_den[ind]*convert2msol) 
ENDFOR

!X.MARGIN = [12,5]
!Y.MARGIN = [6,4]
!p.multi = 0
IF KEYWORD_SET(outplot) THEN  device,filename = outplot + 'gasProf.eps' ELSE BEGIN
    window,3,xsize = 400,ysize = 400
    loadct,39
ENDELSE
plot,r_bin[0:nbins - 1] + dr/2.0,alog10(prof_HI),yrange = [-1,2],xrange = [0,5],xtitle = 'Radius [kpc]',ytitle = textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"
oplot,r_bin[0:nbins - 1] + dr/2.0,alog10(prof_H2),linestyle = 1
legend,['HI',textoidl('H_2')],linestyle = [0,1],/right
IF KEYWORD_SET(outplot) THEN device,/close

IF KEYWORD_SET(outplot) THEN  device,filename = outplot + 'gasProf.eps' ELSE BEGIN
    window,4,xsize = 400,ysize = 400
    loadct,39
ENDELSE
plot,(H2rebin + HIrebin)*cos(angle),2.0*H2rebin/(H2rebin + HIrebin),xtitle = textoidl('\Sigma_{HI + H_{2}} [amu/cm^2]'),ytitle = textoidl('H_2/H'),psym = 1,/ylog,/xlog,xrange = [1e19,1e23],yrange = [1e-6,1]
readcol,'/astro/users/christensen/code/MolecH/Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
readcol,'/astro/users/christensen/code/MolecH/Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'
oplot,10^nhw,10^logfH2w,psym = 4,color = 150
oplot,10^nhg,10^logfH2g,psym = 6,color = 150
legend,['Simulated Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [1,4,6],/right,/bottom,charsize = l_charsize,color = [fgcolor,150,150]
stop 

END

