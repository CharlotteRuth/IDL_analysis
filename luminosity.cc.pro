pro luminosities

;laroy chase 2/24/10
;read in broadband.fits and determine colors and luminosities at
;different locations in the galaxy


;Make a histogram of the surface brightnesses for each cell for the different filters
loadct,39

dir ="/astro/users/lacowboy/sunrise/h277.cosmo50cmb.1536g2bwdK.00504.1/"
dir ="/astro/users/lacowboy/sunrise/h277.cosmo50cmb.1536g2bwdK.00512.1/"

file = "broadband.fits"



filename = dir + file
data=mrdfits(filename,13,header)
help,data,/struct
print,data.filter
data1=mrdfits(filename,14,header)
help,data1

;*********************************************************

u = data1[*,*,2]  ;Get the layer of the cube that corresponds to the u_SDSS filter and put it in the variable 'u'
;data in u is in units of Watts/m/m^2/sterradians

num_cells_SB = histogram(u, locations = SB,binsize=4) ;returns the number of cells of a given surface brightness, SB is now a variable containing the surface brightness
 
;********************************************************

device,decomposed = 0
;set_plot,'ps' ;These two lines save a copy of the graph.  To plot to the screen comment them out
;device,filename = 'surfacebrightness.ps',/color
window,0

plot,SB,num_cells_SB,psym = 10,yrange = [1,1e6],xrange = [0,400],ytitle = 'num_cells_SB',xtitle = 'surface brightness[SB]',/ylog
oplot,SB,num_cells_SB,psym = 10,color = 30
;*********************************************************
g = data1[*,*,3]  ;Get the layer of the cube that corresponds to the u_SDSS filter and put it in the variable 'u'
;data in u is in units of Watts/m/m^2/sterradians
num_cells_SB = histogram(g, locations = SB,binsize=4) ;returns the number of cells of a given surface brightness, SB is now a variable containing the surface brightness
;********************************************************
oplot,SB,num_cells_SB,psym = 10,color = 140
;*********************************************************
r = data1[*,*,4]  ;Get the layer of the cube that corresponds to the u_SDSS filter and put it in the variable 'u'
;data in u is in units of Watts/m/m^2/sterradians
num_cells_SB = histogram(r, locations = SB,binsize=4) ;returns the number of cells of a given surface brightness, SB is now a variable containing the surface brightness
;********************************************************
oplot,SB,num_cells_SB,psym = 10,color = 220
;*********************************************************
i = data1[*,*,5]  ;Get the layer of the cube that corresponds to the u_SDSS filter and put it in the variable 'u'
;data in u is in units of Watts/m/m^2/sterradians
num_cells_SB = histogram(i, locations = SB,binsize=4) ;returns the number of cells of a given surface brightness, SB is now a variable containing the surface brightness
;********************************************************
oplot,SB,num_cells_SB,psym = 10,color = 240
legend,['u','g','r','i'],linestyle=[0,0,0,0],color=[30,140,220,240],/right
;device,/close;This line saves a copy of the graph.  To plot to the screen comment it out
set_plot,'x'

stop
window,2
filters=mrdfits(filename,13,headerfilter) ;Read in the filter extension where the luminosities in each filter are found
u_lum_0 = filters[2].L_LAMBDA_EFF_NONSCATTER0 ;Face-on non-scattered luminosity in the SDSS g band
sfr_lum = u_lum_0*3.04e-43 ;Equation from Gilbank et al 2010
print,'SFR_u = ',STRTRIM(sfr_lum,2),' M_sol yr^-1'
tipsyfile = dir + 'galaxy.std'
massunit = 1.84793e16
lengthunit = 50000
massform = 3.4645e-12
timeunit = SQRT((lengthunit*3.086d21)^3/(6.67d-8*massunit*1.99d33))/(3600.*24.*365.24)
rtipsy,tipsyfile,h,g,d,s
sfr,s,massunit = massunit,timeunit =timeunit,massform = massform
oplot,[0,15],[sfr_lum,sfr_lum],color = 240
stop
nlevels = 60
colors = findgen(nlevels)/nlevels*240
window,1,xsize = 600,ysize = 600
contour,alog10(u),nlevels = nlevels,/fill

end
