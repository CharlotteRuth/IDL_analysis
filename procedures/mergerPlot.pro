;2/28/14
;This plot reads in some files listing the alignment and position of
;the two halos at each timestep and creates a plot showing the orbit

PRO mergerPlot,alignFile1,alignFile2,pfile,nstep = nstep

units = tipsyunits(pfile)
align1 = mrdfits(alignFile1,1)
align2 = mrdfits(alignFile2,1)

IF NOT keyword_set(nstep) THEN nstep = n_elements(align2)
align1 = align1[0:nstep - 1]
align2 = align2[0:nstep - 1]
align2.xc = align2.xc - align1.xc
align2.yc = align2.yc - align1.yc
align2.zc = align2.zc - align1.zc
spinvectorz = [align1[nstep - 1].xa,align1[nstep - 1].za,align1[nstep - 1].ya]
x = spinvectorz[0]
y = spinvectorz[1]
z = spinvectorz[2]
spinvectorx = [z/sqrt(x*x + z*z),0,-1.0*x/sqrt(x*x + z*z)]
spinvectory = crossp(spinvectorz,spinvectorx);[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
basis = [[spinvectorx],[spinvectory],[spinvectorz]]
position = basis#transpose([[align2.xc],[align2.yc],[align2.zc]])

;device,filename = 'mergerplot.eps',/color,bits_per_pixel= 8,xsize = 8,ysize = 4,/inch
window,0,xsize = 800,ysize = 400
multiplot,/reset
multiplot,[2,1],/square
plot,position[0,*]*units.lengthunit,position[1,*]*units.lengthunit,xtitle = 'X',ytitle = 'Y'
oplot,[0,0],[-1e4,1e4],linestyle = 1
oplot,[-1e4,1e4],[0,0],linestyle = 1
multiplot
plot,position[2,*]*units.lengthunit,position[1,*]*units.lengthunit,xtitle = 'Z'
oplot,[0,0],[-1e4,1e4],linestyle = 1
oplot,[-1e4,1e4],[0,0],linestyle = 1
multiplot,/reset
;device,/close
END
