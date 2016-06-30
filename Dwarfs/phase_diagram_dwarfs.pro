; .r /astro/users/christensen/code/MolecH/hist_2d_weighted.pro
; .r /astro/users/christensen/code/MolecH/contour_plus.pro

pro phase_diagram_dwarfs,infile,outfile
;phase_diagram,'1E5R/10M/o10M_1.00300','10M_phaseD30.ps
;phase_diagram,'1E5R/10M_original/o10M_1.00300','10M_phaseD30_original.ps
;This will plot the phase diagram for the gas in a given simulation
lengthunit = 3.0857d21 ;system length unit in cm (=1kpc)
massunit = 2.362d5 ;system mass unit in solar masses
solarmass = 1.989d33
k = 1.38d-16
mp = 1.67d-24
loadct,39

set_plot,'x'
rtipsy,'../'+infile,h,g,d,s
density = g.dens*massunit*solarmass/lengthunit^3
;plot,g.dens,g.tempg,psym = 1,xtitle = 'Density/Critical Density',ytitle = 'Temperature (K)'
;stop
;plot,g.tempg,g.tempg*density*k/mp,psym = 1,xtitle = 'Temperature',ytitle = 'Pressure',title = 'Phase Diagram'
;stop
x = g.dens
y = g.tempg
;xrange = MAX(x) - MIN(X)
;yrange = MAX(y) - MIN(y)
;contour_plus,x,y,xbinsize = xrange/30.,ybinsize = yrange/30.,threshold = 5,nlevels = 5,xtitle = 'Density/Critical Density',ytitle = 'Temperature (K)',title = 'Phase Diagram',psym = 1
;stop
x = ALOG10(x)
y = ALOG10(y)
xrange = MAX(x) - MIN(X)
yrange = MAX(y) - MIN(y)
contour_plus,x,y,xbinsize = xrange/30.,ybinsize = yrange/30.,threshold = 15,nlevels = 5,xtitle = 'LOG(Density/Critical Density)',ytitle = 'LOG(Temperature K)',title = 'Phase Diagram',psym = 3,yrange = [0,7],xrange = [-8,8]
set_plot,'ps'
device,filename=outfile,/color,bits_per_pixel=8
contour_plus,x,y,xbinsize = xrange/30.,ybinsize = yrange/30.,threshold = 15,nlevels = 5,xtitle = 'LOG(Density/Critical Density)',ytitle = 'LOG(Temperature K)',title = 'Phase Diagram'+infile,psym = 3,yrange = [0,7],xrange = [-8,8]
device,/close
set_plot,'x'
;contour_plus,x,y,xbinsize = xrange/30.,ybinsize = yrange/30.,threshold = 15,nlevels = 5,xtitle = 'LOG(Density/Critical Density)',ytitle = 'LOG(Temperature K)',title = 'Phase Diagram',psym = 3,/loglevel
;stop
END
