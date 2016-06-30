pro surfaceden_contour,outplot = outplot
loadct,39
rtipsy,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00512.dir/h516.cosmo25cmb.3072g1MBWK.00512.halo.1.std',h,g,d,s
rtipsy,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1',h1,g1,d1,s1
hist0 = HIST_2D( alog10(g.dens*units.rhounit), alog10(g.tempg),min1 = -2, max1 = 3, min2 = 2, max2 = 5,bin1 = 5.0/100,bin2 = 3.0/100)
hist1 = HIST_2D( alog10(g1.dens*units.rhounit), alog10(g1.tempg),min1 = -2, max1 = 3, min2 = 2, max2 = 5,bin1 = 5.0/100,bin2 = 3.0/100)
xarray = findgen(101)/101*5 - 2
yarray = findgen(101)/101*3 + 2

if keyword_set(outplot) then device,filename = outplot,/color,bits_per_pixel= 8,/times,ysize=18,xsize=18,xoffset =  2,yoffset =  2
contour,hist0,xarray,yarray,levels = [5,50,500,5000,50000],xtitle = 'log Density [amu/cc]',ytitle = 'log Temp [K]'
contour,hist1,xarray,yarray,/overplot,color = 240,levels = [5,50,500,5000,50000]
if keyword_set(outplot) then device,/close
end
