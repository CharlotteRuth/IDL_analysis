PRO HIsurfaceden
rtipsy,'h277.cosmo50cmb.3072g14HMbwK.00512.1.std',h,g,d,s
readarr,'h277.cosmo50cmb.3072g14HMbwK.00512.1.HI',h,HI,part = 'gas',/ascii

units = tipsyunits('../../h277.cosmo50cmb.3072g14HMbwK.param')
HI = HI*g.mass*units.massunit
g.x = g.x*units.lengthunit
g.y = g.y*units.lengthunit
g.z = g.z*units.lengthunit

range = 450.0
nbins = 450.0
xaxis = findgen(nbins)*range/nbins - range/2.0
yaxis = findgen(nbins)*range/nbins - range/2.0

result = better_hist2d(g.x,g.y,HI,min1 = -1*range/2,max1 = range/2,binsize1 = range/nbins,min2 = -1*range/2,max2 = range/2,binsize2 = range/nbins)
window,0,xsize = 800,ysize = 800
contour,alog10(result),xaxis,yaxis,/fill,nlevels = 254,min = -4,max = 7

result = better_hist2d(g.x,g.z,HI,min1 = -1*range/2,max1 = range/2,binsize1 = range/nbins,min2 = -1*range/2,max2 = range/2,binsize2 = range/nbins)
window,0,xsize = 800,ysize = 800
contour,alog10(result),xaxis,yaxis,/fill,nlevels = 254,min = -4,max = 7

END
