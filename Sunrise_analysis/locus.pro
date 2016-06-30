FUNCTION locus, xcolor, ycolor, minx, maxx, miny, maxy, binsize, xtitle1, ytitle1,scale

nxbin=FIX((maxx-minx)/binsize)+1
nybin=FIX((maxy-miny)/binsize)+1
hess=hist_2d(xcolor,ycolor,bin1=binsize,bin2=binsize,min1=minx,min2=miny,max1=maxx,max2=maxy)
peak = MAX(hess)
theselevels = scale*[peak/500.,peak/200.,peak/100.,peak/50.,peak/20.,peak/10.,peak/5.,peak/2.,peak]
hessx=FINDGEN(nxbin)*binsize+minx
hessy=FINDGEN(nybin)*binsize+miny
contour,hess,hessx,hessy,levels=theselevels,/fill,xtitle=xtitle1,ytitle=ytitle1,xstyle=1,ystyle=1,xthick=2,ythick=2,charsize=2.5,xmargin=[8,3],ymargin=[3,2],charthick=3

RETURN, 1

END

