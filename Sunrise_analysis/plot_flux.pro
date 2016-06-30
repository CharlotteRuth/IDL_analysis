pro plot_flux
loadct,39
set_plot,'x'
set_plot,'ps'

file1 = '~/Scratch2/Sunrise/MW1lr_v2/broadband.fits'
grid1 = mrdfits(file1,13,h) ;For filter i, that is layer 5
print,file1,h
;plot,grid1.mass_gas,grid1.mass_metals,psym = 3
;stop
file2 = '~/Scratch2/Sunrise/MW1lr_v3/broadband.fits'
grid2 = mrdfits(file2,12,h)
;plot,grid2.mass_gas,grid2.mass_metals,psym = 3
;stop
print,file2,h
iflux1 = REFORM(grid1[*,*,5])
iflux2 = REFORM(grid2[*,*,5])
itotal1 = total(iflux1)
itotal2 = total(iflux2)
iflux1 = iflux1/itotal1
iflux2 = iflux2/itotal2
db = MAX([iflux1,iflux2])/50.
y1 = histogram(iflux1,binsize = db)
x1 = findgen(N_elements(y1))*db
y2 = histogram(iflux2,binsize = db)
x2 = findgen(N_elements(y2))*db
device,filename='images/MW1lr_surfacebrightness_smallstar.eps',/color,bits_per_pixel=8,/times
!X.style = 1
plot,x1,y1,xtitle = '% Total Surface Brightness of Grid Cells',ytitle = 'Number of Cells',title = 'MW1lr -- Surface Brightness Distribution',yrange = [1,MAX([y1,y2])],/ylog,xrange = [0,MAX([x1,x2])]
oplot,x1,y1,color = 240,psym = -2
oplot,x2,y2,color = 50,psym = -4
legend,['MW1lr v2','MW1lr v3'],psym=[-2,-4],color=[240,50],/right
device,/close
stop

file1 = '/astro/net/scratch1/fabio/h258/h258.cosmo50cmb.1536g2bwK.00264.1/broadband.fits'
grid1 = mrdfits(file1,13,h)
print,file1,h
;file2 ='/astro/net/scratch2/cbrook/analysis/sunrise/h258/h258.cosmo50cmb.1536g2bwK.00264.1/broadband.fits'
file2 = '/astro/net/scratch1/fabio/h258/h258.cosmo50cmb.1536g2bwK.00264.1.v3/broadband.fits'

grid2 = mrdfits(file2,12,h)
print,file2,h
iflux1 = REFORM(grid1[*,*,4])
;iflux2 = REFORM(grid2[*,*,5])
iflux2 = REFORM(grid2[*,*,4])
itotal1 = total(iflux1)
itotal2 = total(iflux2)
iflux1 = iflux1/itotal1
iflux2 = iflux2/itotal2
db = MAX([MAX(iflux1),MAX(iflux2)])/50.
y1 = histogram(iflux1,binsize = db)
x1 = findgen(N_elements(y1))*db
y2 = histogram(iflux2,binsize = db)
x2 = findgen(N_elements(y2))*db
device,filename='h258_surfacebrightness_smallstar.eps',/color,bits_per_pixel=8,/times
plot,x1,y1,xtitle = '% Total Surface Brightness of Grid Cells',ytitle = 'Number of Cells',title = 'h258 -- Surface Brightness Distribution',yrange = [1,MAX([y1,y2])],/ylog,xrange = [0,MAX([x1,x2])]
oplot,x1,y1,color = 240,psym = -2
oplot,x2,y2,color = 50,psym = -4
legend,['h258 v2','h258 v3'],psym=[-2,-4],color=[240,50],/right
device,/close

end
