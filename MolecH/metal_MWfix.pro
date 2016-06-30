;file = 'Disk_Collapse.met.ch0.00100'
;outfile = 'Disk_Collapse.scalez0.3.ch0.00100'
pro metal_MWfix,file,outfile,z = z
;This written to scale Disk_Collapse.met.ch0.00100 up to a more
;reasonable metallicity 
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
!P.CHARTHICK=4
!X.THICK=4
!Y.THICK=4
!p.charsize=1.0
!x.charsize=1.5                 ;2.25
!y.charsize=1.5                 ;2.25
!X.MARGIN = [12,3]
!Y.MARGIN = [6,2]
!p.multi = 0
formatplot
set_plot,'x'
;set_plot,'ps'
window,0
zsolar = 0.0130215
rtipsy,file,h,g,d,s
gnew = g
cutz = 0.0005 ;This is the minimum metallicty for the disk gas
IF NOT KEYWORD_SET(z) THEN z = 1
indhigh = where(g.zmetal gt cutz)
disk_ind = where(g.zmetal gt cutz)
indlow = where(g.zmetal le cutz)
meanz = mean(g[indhigh].zmetal)
offset = z*zsolar - meanz
gnew.zmetal = g.zmetal + offset
IF (where(gnew.zmetal lt 0))[0] ne -1 THEN gnew[where(gnew.zmetal lt 0)].zmetal = 0
;gnew[indhigh].zmetal = g[indhigh].zmetal + offset
;gnew[indlow].zmetal = g[indlow].zmetal*(cutz + offset)/cutz 
y = histogram(g.zmetal,locations = x,min = 0, max = 0.028,nbins = 500)
ynew = histogram(gnew.zmetal,locations = x,min = 0, max = 0.028,nbins = 500)
;device,filename = '~/diskmetallicity.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize= 12,xoffset =  2,yoffset =  2
plot,x,y,/ylog,yrange = [1,1e4],xtitle = 'Metallicity',ytitle = 'Distribution'
oplot,x,y,thick = 2
oplot,x,ynew,thick = 1
oplot,[0.00471488,0.00471488],[1,1e4]
oplot,[0.0130215,0.0130215],[1,1e4],linestyle = 2
;device,/close
wtipsy,outfile,h,gnew,d,s,/standard

spawn,'ls *param',result
units = tipsyunits(result[0])
r = sqrt(g.x*g.x + g.y*g.y)
;disk_ind = where(abs(g.z) lt 0.5)
window,1
zmetalhist = weighted_histogram(r[disk_ind],input = g[disk_ind].zmetal*g.zmetal,max = 10,locations = locations)
gmasshist = weighted_histogram(r[disk_ind],input = g.zmetal,max = 10,locations = locations)
zmetalnewhist = weighted_histogram(r[disk_ind],input = gnew[disk_ind].zmetal*gnew.zmetal,max = 10,locations = locations)
gmassnewhist = weighted_histogram(r[disk_ind],input = gnew.zmetal,max = 10,locations = locations)
plot,locations,zmetalhist/gmasshist,yrange = [0,0.03]
oplot,locations,zmetalnewhist/gmassnewhist,linestyle = 2
oplot,[0,10],[zsolar,zsolar],linestyle = 1
;histogramp,r[disk_ind],weight = g[disk_ind].zmetal,max = 10

end
