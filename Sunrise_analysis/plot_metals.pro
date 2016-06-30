pro plot_metals
loadct,39
;set_plot,'x'
set_plot,'ps'

file1 = '~/Scratch2/Sunrise/MW1lr_v3/sfrhist.fits'
grid1 = mrdfits(file1,3)
;plot,grid1.mass_gas,grid1.mass_metals,psym = 3
;stop
file2 = '~/Scratch2/Sunrise/MW1lr_v2/grid_MW1lr_001.fits'
grid2 = mrdfits(file2,8)
;plot,grid2.mass_gas,grid2.mass_metals,psym = 3
;stop
metallicity1 = grid1.mass_metals/grid1.mass_gas
temp = where(Finite(metallicity1),complement = nan)
metallicity1[nan] = 0
db = 0.0005
y1 = histogram(metallicity1,binsize = db)
x1 = findgen(N_elements(y1))*db
metallicity2 = grid2.mass_metals/grid2.mass_gas
temp = where(Finite(metallicity2),complement = nan)
metallicity2[nan] = 0
y2 = histogram(metallicity2,binsize = db)
x2 = findgen(N_elements(y2))*db

;plot,x1,y1,xtitle = 'Metallicity of Grid Cells',ytitle = 'Number of Cells',title = 'MW1lr -- Metallicity of Gas',yrange = [0,MAX([y1,y2])]
;oplot,x1,y1,color = 50
;oplot,x2,y2,color = 240
;legend,['MW1lr v2','MW1lr v3'],color=[240,50]
;stop
mass1 = x1
for i = 0, N_elements(y1)-1 DO BEGIN
    ind = where(metallicity1 ge x1[i] AND metallicity1 lt x1[i]+db)
    IF (ind[0] ne -1) THEN mass1[i] = TOTAL(grid1[ind].mass_gas) 
ENDFOR
mass2 = x2
for i = 0, N_elements(y2)-1 DO BEGIN
    ind = where(metallicity2 ge x2[i] AND metallicity2 lt x2[i]+db)
    IF (ind[0] ne -1) THEN mass2[i] = TOTAL(grid2[ind].mass_gas)     
ENDFOR
device,filename='MW1lr_metals.eps',/color,bits_per_pixel=8,/times
plot,x1,mass1,xtitle = 'Metallicity of Grid Cells',ytitle = 'Mass of Gas at that Metallicity',title = 'MW1lr -- Metallicity of Gas',yrange = [1e2,MAX([mass1,mass2])],/ylog
oplot,x1,mass1,color = 50
oplot,x2,mass2,color = 240
legend,['MW1lr v2','MW1lr v3'],color=[240,50]
device,/close
stop


file1 = '/astro/net/scratch2/cbrook/analysis/sunrise/h258/h258.cosmo50cmb.1536g2bwK.00264.1/grid.fits'
grid1 = mrdfits(file1,3)
;plot,grid1.mass_gas,grid1.mass_metals,psym = 3
;stop
file2 = '/astro/net/scratch1/fabio/h258/h258.cosmo50cmb.1536g2bwK.00264.1/makegrid.fits'
grid2 = mrdfits(file2,8)
;plot,grid2.mass_gas,grid2.mass_metals,psym = 3
;stop
metallicity1 = grid1.mass_metals/grid1.mass_gas
db = 0.0005
y1 = histogram(metallicity1,binsize = db)
x1 = findgen(N_elements(y1))*db
metallicity2 = grid2.mass_metals/grid2.mass_gas
y2 = histogram(metallicity2,binsize = db)
x2 = findgen(N_elements(y2))*db

;plot,x1,y1,xtitle = 'Metallicity of Grid Cells',ytitle = 'Number of Cells',title = 'h258.264 -- Metallicity of Gas',yrange = [0,MAX([y1,y2])]
;oplot,x1,y1,color = 50
;oplot,x2,y2,color = 240
;legend,['h258.264 v2','h258.264 v3'],color=[240,50]
;stop
mass1 = x1
for i = 0, N_elements(y1)-1 DO BEGIN
    ind = where(metallicity1 ge x1[i] AND metallicity1 lt x1[i]+db)
    IF (ind[0] ne -1) THEN mass1[i] = TOTAL(grid1[ind].mass_gas) 
ENDFOR
mass2 = x2
for i = 0, N_elements(y2)-1 DO BEGIN
    ind = where(metallicity2 ge x2[i] AND metallicity2 lt x2[i]+db)
    IF (ind[0] ne -1) THEN mass2[i] = TOTAL(grid2[ind].mass_gas)     
ENDFOR
device,filename='h258_metals.eps',/color,bits_per_pixel=8,/times
plot,x1,mass1,xtitle = 'Metallicity of Grid Cells',ytitle = 'Mass of Gas at that Metallicity',title = 'h258.264 -- Metallicity of Gas',yrange = [1e4,MAX([mass1,mass2])],/ylog
oplot,x1,mass1,color = 50
oplot,x2,mass2,color = 240
legend,['h258.264 v2','h258.264 v3'],color=[240,50]
device,/close
;stop



end
