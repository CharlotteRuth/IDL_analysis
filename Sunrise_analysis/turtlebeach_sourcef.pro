PRO turtlebeach_sourcef,dir,file,unitfile,outplot = outplot
formatplot,outplot = outplot
restore,dir+'/sunrise_grid.1.save'

;for CO 1-0
freq = 115.d9
AUL = 7.203d-8
a_h = 6.6260755d-27
a_c = 3d10
a_k = 1.380658d-16 	

;popfile1=turtledir+'/COdumps/co'+strcompress(SNAP,/remove_all)+'_1pop1dump.dat'
;popfile0=turtledir+'/COdumps/co'+strcompress(SNAP,/remove_all)+'_1pop0dump.dat'
popfile1 = dir+'/co1_1pop1dump.dat'
popfile0 = dir+'/co1_1pop0dump.dat'

print,'reading in '+popfile0
readcol,popfile0,pop0
print,'reading in '+popfile1
readcol,popfile1,pop1


;code saves in log
pop0=10D^pop0
pop1=10D^pop1


g1 = (2.*1+1)
g0 = (2.*0+1)


sfunc = 2.*a_h*(freq)^3D/a_c^2D
sfunc /= ((g1*pop0)/(g0*pop1)-1)
tb_ico = sfunc * (a_c^2.)/(2.*freq^2.*a_k)

IF keyword_set(outplot) THEN thick = 4 ELSE thick = 1
IF keyword_set(outplot) THEN device,filename=outplot + '_co10line.eps', /color,xsize = 15, ysize =  12
histogramp,vx*1e-5,weight = tb_ico,min = -400,max = 400,xtitle = 'Velocity [km/s]',ytitle = 'Source Function',title = 'CO 1-0 Emission, z = 1.9',thick = thick,yrange = [0,3e6]
IF keyword_set(outplot) THEN device,/close ELSE stop

;units = tipsyunits(unitfile)
;rtipsy,dir+'/'+file,h,g,d,s
;readarr,dir+'/'+file+'.HI',h,hi,part = 'gas',/ascii
;histogramp,g.vx*units.vunit*h.time,weight = HI*g.mass*units.massunit,min = -400,max = 400
;stop
END
