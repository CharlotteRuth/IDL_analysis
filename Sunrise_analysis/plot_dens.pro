pro plot_dens
loadct,39
set_plot,'x'
set_plot,'ps'

file1 = '~/Scratch2/Sunrise/MW1lr_v3/sfrhist.fits'
grid1 = mrdfits(file1,3)
file2 = '~/Scratch2/Sunrise/MW1lr_v2/grid_MW1lr_001.fits'
grid2 = mrdfits(file2,8)
density1 = grid1.mass_gas/grid1.cell_volume
temp = where(Finite(density1),complement = nan)
IF (nan[0] ne -1) THEN density1[nan] = 0
density2 = grid2.mass_gas/grid2.cell_volume
temp = where(Finite(density2),complement = nan)
IF (nan[0] ne -1) THEN density2[nan] = 0

db = MAX([density1,density2])/50.
y1 = histogram(density1,binsize = db)
x1 = findgen(N_elements(y1))*db
y2 = histogram(density2,binsize = db)
x2 = findgen(N_elements(y2))*db

mass1 = x1
for i = 0, N_elements(y1)-1 DO BEGIN
    ind = where(density1 ge x1[i] AND density1 lt x1[i]+db)
    IF (ind[0] ne -1) THEN mass1[i] = TOTAL(grid1[ind].mass_gas) 
ENDFOR
mass2 = x2
for i = 0, N_elements(y2)-1 DO BEGIN
    ind = where(density2 ge x2[i] AND density2 lt x2[i]+db)
    IF (ind[0] ne -1) THEN mass2[i] = TOTAL(grid2[ind].mass_gas)     
ENDFOR

device,filename='MW1lr_density.eps',/color,bits_per_pixel=8,/times
plot,x1,mass1/Total(mass1),xtitle = 'Gas Density',ytitle = 'Mass Fraction of Total Gas',title = 'MW1lr -- Density of Gas',yrange = [1e-5,MAX([y1/Total(y1),y2/Total(y2)])],/ylog,xrange = [0,50*db]
oplot,x1,mass1/Total(mass1),color = 50
oplot,x2,mass2/Total(mass2),color = 240
legend,['MW1lr v2','MW1lr v3'],color=[240,50],linestyle=[0,0],/right
device,/close
stop

file1 = '/astro/net/scratch2/cbrook/analysis/sunrise/h258/h258.cosmo50cmb.1536g2bwK.00264.1/grid.fits'
grid1 = mrdfits(file1,3)
file2 = '/astro/net/scratch1/fabio/h258/h258.cosmo50cmb.1536g2bwK.00264.1/makegrid.fits'
grid2 = mrdfits(file2,8)
density1 = grid1.mass_gas/grid1.cell_volume
temp = where(Finite(density1),complement = nan)
IF (nan[0] ne -1) THEN density1[nan] = 0
density2 = grid2.mass_gas/grid2.cell_volume
temp = where(Finite(density2),complement = nan)
IF (nan[0] ne -1) THEN density2[nan] = 0

db = MAX([density1,density2])/50.
y1 = histogram(density1,binsize = db)
x1 = findgen(N_elements(y1))*db
y2 = histogram(density2,binsize = db)
x2 = findgen(N_elements(y2))*db

mass1 = x1
for i = 0, N_elements(y1)-1 DO BEGIN
    ind = where(density1 ge x1[i] AND density1 lt x1[i]+db)
    IF (ind[0] ne -1) THEN mass1[i] = TOTAL(grid1[ind].mass_gas) 
ENDFOR
mass2 = x2
for i = 0, N_elements(y2)-1 DO BEGIN
    ind = where(density2 ge x2[i] AND density2 lt x2[i]+db)
    IF (ind[0] ne -1) THEN mass2[i] = TOTAL(grid2[ind].mass_gas)     
ENDFOR

device,filename='h258_density.eps',/color,bits_per_pixel=8,/times
plot,x1,mass1/Total(mass1),xtitle = 'Gas Density',ytitle = 'Mass Fraction of Total Gas',title = 'h258 -- Density of Gas',yrange = [1e-10,MAX([y1/Total(y1),y2/Total(y2)])],/ylog,xrange = [0,50*db]
oplot,x1,mass1/Total(mass1),color = 50
oplot,x2,mass2/Total(mass2),color = 240
legend,['h258 v2','h258 v3'],color=[240,50],linestyle=[0,0],/right
device,/close
stop
end
