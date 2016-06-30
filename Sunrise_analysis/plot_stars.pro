pro plot_stars
loadct,39
set_plot,'x'
;set_plot,'ps'
solarlum = 3.839e26
file1 = '/astro/net/scratch2/christensen/Sunrise/MW1lr_v2/MW1lr_001.std.stars.star99_ab.Mv.fits'
file1 = '/astro/net/scratch2/christensen/Sunrise/MW1lr_v2_oldstars/MW1lr_oldstars.std.stars.star99_ab.Mv.fits'
;file1 = '/astro/net/scratch2/christensen/Sunrise/onestar_test/galaxy.std.stars.star99_ab.Mv.fits'
grid1 = MRDFITS(file1,1)
print,MINMAX(alog10(grid1.age))
print,MINMAX(grid1.metallicity)
file2 = '~/Scratch2/Sunrise/MW1lr_v2/grid.fits'
file2 = '/astro/net/scratch2/christensen/Sunrise/MW1lr_v2_oldstars/makegrid.fits'
;file2 = '/astro/net/scratch2/christensen/Sunrise/onestar_test/testV2/grid_000.fits'
grid2 = mrdfits(file2,9)
file3 = '~/Scratch2/Sunrise/MW1lr_v3_mc0_multi0/sfrhist.fits'
;file3 = '/astro/net/scratch2/christensen/Sunrise/onestar_test/testV3/sfrhist.fits'
grid3 = mrdfits(file3,6)

age2 = alog10(grid2.age_m)
age3 = alog10(grid3.age)
age1 = alog10(grid1.age)
db = (MAX([age2,age3]) - MIN([age2,age3]))/50.
mina = MIN([age1,age2,age3])
age2_y = histogram(age2,binsize = db, min = mina)
age2_x = findgen(N_elements(age2_y))*db + mina
age3_y = histogram(age3,binsize = db, min = mina)
age3_x = findgen(N_elements(age3_y))*db + mina
age1_y = histogram(age1,binsize = db, min = mina)
age1_x = findgen(N_elements(age1_y))*db + mina

metallicity2 = grid2.mass_stellar_metals/grid2.mass_stars
metallicity3 = grid3.metallicity
metallicity1 = grid1.metallicity
db = MAX([metallicity2,metallicity3])/50.
minm = MIN([metallicity1,metallicity2,metallicity3])
metallicity2_y = histogram(metallicity2,binsize = db,min = minm)
metallicity2_x = findgen(N_elements(metallicity2_y))*db + minm
metallicity3_y = histogram(metallicity3,binsize = db,min = minm)
metallicity3_x = findgen(N_elements(metallicity3_y))*db + minm
metallicity1_y = histogram(metallicity1,binsize = db,min = minm)
metallicity1_x = findgen(N_elements(metallicity1_y))*db + minm

lum1 = alog10(10.^(-0.4*grid1.mb))
lum2 = alog10(grid2.l_bol/3.839e26)-1.9
lum3 = alog10(grid3.l_bol/3.839e26)-1.9
db = MAX([lum2,lum3])/400.
minl = MIN([lum1,lum2,lum3])
lum1_y = histogram(lum1,binsize = db,min = minl)
lum1_x = findgen(N_elements(lum1_y))*db + minl
lum2_y = histogram(lum2,binsize = db,min = minl)
lum2_x = findgen(N_elements(lum2_y))*db + minl
lum3_y = histogram(lum3,binsize = db,min = minl)
lum3_x = findgen(N_elements(lum3_y))*db + minl

window,0
;device,filename='MW1lr_stellar_age.eps',/color,bits_per_pixel=8,/times
plot,age2_x,age2_y,xtitle = 'Stellar Age',ytitle = 'Number of Stars',yrange = [1,MAX([age2_y,age3_y])],/ylog,psym = 10,charsize = 2
oplot,age3_x,age3_y,linestyle = 1,psym = 10,color = 240
oplot,age1_x,age1_y,linestyle = 2,psym = 10,color = 50
legend,['IDL','v2','v3'],linestyle=[0,1,2],/left,charsize = 2,color=[0,240,50]
;device,/close

;device,filename='MW1lr_stellar_luminosity.eps',/color,bits_per_pixel=8,/times
;plot,mass2_x,mass2_y,xtitle = 'Stellar Mass',ytitle = 'Number of Stars',xrange = [1,1e6],/ylog,psym = 10,charsize = 2,yrange = [1,MAX([mass2_y,mass3_y])]
;oplot,mass3_x,mass3_y,linestyle = 1,psym = 10
;legend,['v2','v3'],linestyle=[0,1],/right,charsize = 2
;device,/close
;stop

window,1
;device,filename='MW1lr_stellar_luminosity.eps',/color,bits_per_pixel=8,/times
plot,metallicity1_x,metallicity2_y,xtitle = 'Stellar Metallicity',ytitle = 'Number of Stars',yrange = [1,MAX([metallicity2_y,metallicity3_y])],/ylog,psym = 10,charsize = 2
oplot,metallicity2_x,metallicity2_y,linestyle = 1,psym = 10,color = 240
oplot,metallicity3_x,metallicity3_y,linestyle = 2,psym = 10,color = 50
legend,['IDL','v2','v3'],linestyle=[0,1,2],/right,charsize = 2,color=[0,240,50]
;device,/close

window,2
plot,lum1_x,lum1_y,xtitle = 'Log Stellar Luminosity',ytitle = 'Number of Stars',yrange = [1,MAX([lum2_y,lum3_y])],/ylog,psym = 10,charsize = 2
oplot,lum2_x,lum2_y,linestyle = 1,psym = 10,color = 240
oplot,lum3_x,lum3_y,linestyle = 2,psym = 10,color = 50
legend,['Star99','v2','v3'],linestyle=[0,1,2],/right,charsize = 2,color=[0,240,50]

set_plot,'ps'
device,filename='/astro/net/scratch2/christensen/Sunrise/images/MW1lr_stellar_luminosities.eps',/color,bits_per_pixel=8,/times
plot,lum1_x,lum1_y,xtitle = 'Log Stellar Luminosity',ytitle = 'Number of Stars',yrange = [1,MAX([lum2_y,lum3_y])],/ylog,psym = 10,charsize = 2
oplot,lum2_x,lum2_y,linestyle = 1,psym = 10,color = 240
oplot,lum3_x,lum3_y,linestyle = 2,psym = 10,color = 50
legend,['Star99','v2','v3'],linestyle=[0,1,2],/right,charsize = 2,color=[0,240,50]
device,/close

print,'Startburst99 Luminosity: ',Total(lum1),format='(A25,I10)'
print,'Sunrise v2   Luminosity: ',Total(lum2),format='(A25,I10)'
print,'Sunrise v3   Luminosity: ',Total(lum3),format='(A25,I10)'

print,'Startburst99 Magnitude: ',-2.5*alog10(Total(lum1))
print,'Sunrise v2   Magnitude: ',-2.5*alog10(Total(lum2))
print,'Sunrise v3   Magnitude: ',-2.5*alog10(Total(lum3))

lum1s = lum1[SORT(lum1)]
lum2s = lum2[SORT(lum2)]
lum3s = lum3[SORT(lum3)]
stop
end
