PRO z_rho,outplot = outplot

formatplot,outplot = outplot
loadct,39

rtipsy,'h986.cosmo50cmb.3072g14HBWK.00512.halo.1',h,g,d,s
readarr,'h986.cosmo50cmb.3072g14HBWK.00512.halo.1.decomp',h,deco,part = 'star',/ascii
readarr,'h986.cosmo50cmb.3072g14HBWK.00512.halo.1.iord',h,iord,part = 'star',/ascii
readarr,'h986.cosmo50cmb.3072g14HBWK.00512.halo.1.OxMassFrac',h,ox,part = 'star',/ascii
readarr,'h986.cosmo50cmb.3072g14HBWK.00512.halo.1.FeMassFrac',h,fe,part = 'star',/ascii
;rtipsy,'h986.cosmo50cmb.3072g14HBWK.00512',h,g,d,s
;readarr,'h986.cosmo50cmb.3072g14HBWK.00512.decomp',h,deco,part = 'star',/ascii
;readarr,'h986.cosmo50cmb.3072g14HBWK.00512.iord',h,iord,part = 'star',/ascii
;readarr,'h986.cosmo50cmb.3072g14HBWK.00512.OxMassFrac',h,ox,part = 'star',/ascii
;readarr,'h986.cosmo50cmb.3072g14HBWK.00512.FeMassFrac',h,fe,part = 'star',/ascii

units = tipsyunits('../../h986.cosmo50cmb.3072g14HBWK.param')


sl = rstarlog('../../h986.cosmo50cmb.3072g14HBWK.starlog',/molecularH)
match,iord,sl.iorderstar,indt,indsl
sl = sl[indsl]
s = s[indt]
deco = deco[indt]
iord = iord[indt]
ox = ox[indt]
fe = fe[indt]
z =  1.06*fe + 2.09*ox
zsolar = 0.0130215
z = z/zsolar
sl.rhoform = sl.rhoform*units.rhounit
;plot,sl.rhoform,z,/xlog,/ylog,psym = 3,xtitle = 'Density [amu/cc]',ytitle = 'Metallicity [Z'+sunsymbol()+']'
;stop
;plot,sl.rhoform,z,/xlog,/ylog,psym = 3,xtitle = 'Density [amu/cc]',ytitle = 'Metallicity [Z'+sunsymbol()+']',xrange = [1,1e4],yrange = [1e-2,1]

binsize1 = 0.1
binsize2 = 0.1
min1 = 1
min2 = -2
max1 = 4
max2 = 0

rho_log = alog10(sl.rhoform)
z_log = alog10(z)
;plot,rho_log,z_log,psym = 3,xtitle = 'Log Density [amu/cc]',ytitle = 'Log Metallicity [Z'+sunsymbol()+']',xrange = [0,4],yrange = [-2,0]

IF keyword_set(outplot) THEN device,filename='~/h986_z_rho_rho.eps', /color,xsize = 14, ysize = 21,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 600,ysize = 900
multiplot,[1,2]

hist = hist2d(rho_log, z_log, min = [min1,min2], max = [max1,max2], bin = [binsize1,binsize2],out = outarray)
;stop
;hist = better_hist2d(rho_log, z_log, sl.massform, min1 = min1, max1 = max1, min2 = min2, max2 = max2, binsize1 = binsize1,binsize2 = binsize2)
;xarray = (findgen((size(hist))[1]) + 0.5)*binsize1 + min1
;yarray = (findgen((size(hist))[2]) + 0.5)*binsize2 + min2
contour,hist,outarray.xax,outarray.yax,ytitle = 'Log Metallicity [Z'+sunsymbol()+']',nlevels = 20,xrange = [min1,max1],yrange = [min2,max2]

inddisk = where(deco EQ 1)
indbulge = where(deco EQ 3)

rho_log_disk = alog10(sl[inddisk].rhoform)
z_log_disk = alog10(z[inddisk])
rho_log_bulge = alog10(sl[indbulge].rhoform)
z_log_bulge = alog10(z[indbulge])

;hist_disk = better_hist2d(rho_log_disk, z_log_disk, sl[inddisk].massform, min1 = min1, max1 = 4, min2 = min2, max2 = 0, binsize1 = binsize1,binsize2 = binsize2)
;xarray = (findgen((size(hist_disk))[1]) + 0.5)*binsize1 + min1
;yarray = (findgen((size(hist_disk))[2]) + 0.5)*binsize2 + min2
hist_disk = hist2d(rho_log_disk, z_log_disk, min = [min1,min2], max = [max1,max2], bin = [binsize1,binsize2],out = outarray)
contour,hist_disk,outarray.xax,outarray.yax,nlevels = 20,color= 60,/overplot
;hist_bulge = better_hist2d(rho_log_bulge, z_log_bulge, sl[indbulge].massform, min1 = min1, max1 = 4, min2 = min2, max2 = 0, binsize1 = binsize1,binsize2 = binsize2)
;xarray = (findgen((size(hist_bulge))[1]) + 0.5)*binsize1 + min1
;yarray = (findgen((size(hist_bulge))[2]) + 0.5)*binsize2 + min2
hist_bulge = hist2d(rho_log_bulge, z_log_bulge, min = [min1,min2], max = [max1,max2], bin = [binsize1,binsize2],out = outarray)
contour,hist_bulge,outarray.xax,outarray.yax,nlevels = 20,color= 254,/overplot
legend,['All Stars','Disk','Bulge'],color = [0,60,254],linestyle = [0,0,0],/right
multiplot

thick = 4
histogramp,rho_log,min = min1,max = max1,nbins = 100,xtitle = 'Log Density [amu/cc]',xrange = [min1,max1],thick = thick
histogramp,rho_log_disk,min = min1,max = max1,nbins = 100,/overplot,color = 60,thick = thick
histogramp,rho_log_bulge,min = min1,max = max1,nbins = 100,/overplot,color = 254,thick = thick
multiplot,/reset
IF keyword_set(outplot) THEN device,/close ELSE stop

IF keyword_set(outplot) THEN device,filename='~/h986_z_rho.eps', /color,xsize = 18, ysize = 15,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 600,ysize = 400
multiplot,[1,1]
hist = hist2d(rho_log, z_log, min = [min1,min2], max = [max1,max2], bin = [binsize1,binsize2],out = outarray)
;stop
;hist = better_hist2d(rho_log, z_log, sl.massform, min1 = min1, max1 = max1, min2 = min2, max2 = max2, binsize1 = binsize1,binsize2 = binsize2)
;xarray = (findgen((size(hist))[1]) + 0.5)*binsize1 + min1
;yarray = (findgen((size(hist))[2]) + 0.5)*binsize2 + min2
contour,hist,outarray.xax,outarray.yax,ytitle = 'Log Metallicity [Z'+sunsymbol()+']',nlevels = 20,xrange = [min1,max1],yrange = [min2,max2],xtitle = 'Log Density [amu/cc]'

hist_disk = hist2d(rho_log_disk, z_log_disk, min = [min1,min2], max = [max1,max2], bin = [binsize1,binsize2],out = outarray)
contour,hist_disk,outarray.xax,outarray.yax,nlevels = 20,color= 60,/overplot

hist_bulge = hist2d(rho_log_bulge, z_log_bulge, min = [min1,min2], max = [max1,max2], bin = [binsize1,binsize2],out = outarray)
contour,hist_bulge,outarray.xax,outarray.yax,nlevels = 20,color= 254,/overplot
legend,['All Stars','Disk','Bulge'],color = [0,60,254],linestyle = [0,0,0],/right
multiplot,/reset
IF keyword_set(outplot) THEN device,/close ELSE stop

IF keyword_set(outplot) THEN device,filename='~/h986_z_rho2.eps', /color,xsize = 18, ysize = 15,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 600,ysize = 400
multiplot,[1,1]
;hist = hist2d(rho_log, z_log, min = [min1,min2], max = [max1,max2], bin = [binsize1,binsize2],out = outarray)
;stop
;hist = better_hist2d(rho_log, z_log, sl.massform, min1 = min1, max1 = max1, min2 = min2, max2 = max2, binsize1 = binsize1,binsize2 = binsize2)
;xarray = (findgen((size(hist))[1]) + 0.5)*binsize1 + min1
;yarray = (findgen((size(hist))[2]) + 0.5)*binsize2 + min2
;contour,hist,outarray.xax,outarray.yax,ytitle = 'Log Metallicity [Z'+sunsymbol()+']',nlevels = 20,xrange = [min1,max1],yrange = [min2,max2],xtitle = 'Log Density [amu/cc]'

max2 = -0.1
min2 = -1.5
max1 = 3.7

hist_disk = hist2d(rho_log_disk, z_log_disk, min = [min1,min2], max = [max1,max2], bin = [binsize1,binsize2],out = outarray)
contour,hist_disk,outarray.xax,outarray.yax,ytitle = 'Log Metallicity [Z'+sunsymbol()+']',nlevels = 20,xrange = [min1,max1],yrange = [min2,max2],xtitle = 'Log Density [amu/cc]',/nodata
contour,hist_disk,outarray.xax,outarray.yax,nlevels = 15,color= 60,/overplot,thick = thick

hist_bulge = hist2d(rho_log_bulge, z_log_bulge, min = [min1,min2], max = [max1,max2], bin = [binsize1,binsize2],out = outarray)
contour,hist_bulge,outarray.xax,outarray.yax,nlevels = 15,color= 254,/overplot,thick = thick
legend,['Disk','Bulge'],color = [60,254],linestyle = [0,0],/right
multiplot,/reset
IF keyword_set(outplot) THEN device,/close ELSE stop

END
