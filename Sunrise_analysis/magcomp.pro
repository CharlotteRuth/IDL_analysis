pro magcomp
loadct,39
set_plot,'x'
set_plot,'ps'

file2 = '~/Scratch2/Sunrise/MW1lr_v3/broadband.fits'
filters2 = mrdfits(file2,11)
file1 = '~/Scratch2/Sunrise/MW1lr_v2/broadband.fits'
filters1 = mrdfits(file1,12)
file3 = '~/Scratch2/Sunrise/MW1lr_v3_mcl/broadband.fits'
filters3 = mrdfits(file3,11)
;rint,filters1[2:6].filter
;rint,filters2[2:6].filter
unscatmag1 = filters1[2:6].AB_MAG_NONSCATTER0
scatmag1 =   filters1[2:6].AB_MAG0
unscatmag2 = filters2[2:6].AB_MAG_NONSCATTER0
scatmag2 =   filters2[2:6].AB_MAG0
unscatmag3 = filters3[2:6].AB_MAG_NONSCATTER0
scatmag3 =   filters3[2:6].AB_MAG0
print,'MW1lr'
print,'Version 2'
print,filters1[2].filter,filters1[2].AB_MAG_NONSCATTER0,filters1[2].AB_MAG0
print,filters1[3].filter,filters1[3].AB_MAG_NONSCATTER0,filters1[3].AB_MAG0
print,filters1[4].filter,filters1[4].AB_MAG_NONSCATTER0,filters1[4].AB_MAG0
print,filters1[5].filter,filters1[5].AB_MAG_NONSCATTER0,filters1[5].AB_MAG0
print,filters2[6].filter,filters1[6].AB_MAG_NONSCATTER0,filters1[6].AB_MAG0
print,''
print,'Version 3, mappings_mcl = 5e6'
print,filters2[2].filter,filters2[2].AB_MAG_NONSCATTER0,filters2[2].AB_MAG0
print,filters2[3].filter,filters2[3].AB_MAG_NONSCATTER0,filters2[3].AB_MAG0
print,filters2[4].filter,filters2[4].AB_MAG_NONSCATTER0,filters2[4].AB_MAG0
print,filters2[5].filter,filters2[5].AB_MAG_NONSCATTER0,filters2[5].AB_MAG0
print,filters2[6].filter,filters2[6].AB_MAG_NONSCATTER0,filters2[6].AB_MAG0
print,''
print,'Version 3, mappings_mcl = 1e6'
print,filters3[2].filter,filters3[2].AB_MAG_NONSCATTER0,filters3[2].AB_MAG0
print,filters3[3].filter,filters3[3].AB_MAG_NONSCATTER0,filters3[3].AB_MAG0
print,filters3[4].filter,filters3[4].AB_MAG_NONSCATTER0,filters3[4].AB_MAG0
print,filters3[5].filter,filters3[5].AB_MAG_NONSCATTER0,filters3[5].AB_MAG0
print,filters3[6].filter,filters3[6].AB_MAG_NONSCATTER0,filters3[6].AB_MAG0

!X.style = 1
!Y.style = 1
device,filename='MW1lr_faceon_filters_mcl.eps',/color,bits_per_pixel=8,/times
plot,[-23,-20],[-23,-20],xrange=[-23,-20],yrange = [-23,-20],xtitle = 'Magnitude without Dust',ytitle = 'Magnitude with Dust',title = 'Face on Magnitudes'
oplot,unscatmag1,scatmag1,psym = 2,color=240
oplot,unscatmag2,scatmag2,psym = 4,color=50
oplot,unscatmag3,scatmag3,psym = 5,color=50
legend,['MW1lr v2','MW1lr v3, map_mcl = 5e6', 'MW1lr v3, map_mcl = 1e6'],psym=[2,4,5],color=[240,50,50]
device,/close
;stop

unscatmag1 = filters1[2:6].AB_MAG_NONSCATTER1
scatmag1 =   filters1[2:6].AB_MAG1
unscatmag2 = filters2[2:6].AB_MAG_NONSCATTER1
scatmag2 =   filters2[2:6].AB_MAG1
unscatmag3 = filters3[2:6].AB_MAG_NONSCATTER0
scatmag3 =   filters3[2:6].AB_MAG0
print,'Version 2'
print,filters1[2].filter,filters1[2].AB_MAG_NONSCATTER1,filters1[2].AB_MAG1
print,filters1[3].filter,filters1[3].AB_MAG_NONSCATTER1,filters1[3].AB_MAG1
print,filters1[4].filter,filters1[4].AB_MAG_NONSCATTER1,filters1[4].AB_MAG1
print,filters1[5].filter,filters1[5].AB_MAG_NONSCATTER1,filters1[5].AB_MAG1
print,filters2[6].filter,filters1[6].AB_MAG_NONSCATTER1,filters1[6].AB_MAG1
print,''
print,'Version 3, mappings_mcl = 5e6'
print,filters2[2].filter,filters2[2].AB_MAG_NONSCATTER1,filters2[2].AB_MAG1
print,filters2[3].filter,filters2[3].AB_MAG_NONSCATTER1,filters2[3].AB_MAG1
print,filters2[4].filter,filters2[4].AB_MAG_NONSCATTER1,filters2[4].AB_MAG1
print,filters2[5].filter,filters2[5].AB_MAG_NONSCATTER1,filters2[5].AB_MAG1
print,filters2[6].filter,filters2[6].AB_MAG_NONSCATTER1,filters2[6].AB_MAG1
print,''
print,'Version 3, mappings_mcl = 1e6'
print,filters3[2].filter,filters3[2].AB_MAG_NONSCATTER0,filters3[2].AB_MAG0
print,filters3[3].filter,filters3[3].AB_MAG_NONSCATTER0,filters3[3].AB_MAG0
print,filters3[4].filter,filters3[4].AB_MAG_NONSCATTER0,filters3[4].AB_MAG0
print,filters3[5].filter,filters3[5].AB_MAG_NONSCATTER0,filters3[5].AB_MAG0
print,filters3[6].filter,filters3[6].AB_MAG_NONSCATTER0,filters3[6].AB_MAG0
device,filename='MW1lr_edgeon_filters_mcl.eps',/color,bits_per_pixel=8,/times
plot,[-23,-19.5],[-23,-19.5],xrange=[-23,-19.5],yrange = [-23,-19.5],xtitle = 'Magnitude without Dust',ytitle = 'Magnitude with Dust',title = 'Edge on Magnitudes'
oplot,unscatmag1,scatmag1,psym = 2,color=240
oplot,unscatmag2,scatmag2,psym = 4,color=50
oplot,unscatmag3,scatmag3,psym = 5,color=50
;legend,['MW1lr v2','MW1lr v3'],psym=[2,4],color=[240,50]
legend,['MW1lr v2','MW1lr v3, map_mcl = 5e6', 'MW1lr v3, map_mcl = 1e6'],psym=[2,4,5],color=[240,50,50]
device,/close

stop
set_plot,'ps'
file2 = '/astro/net/scratch2/cbrook/analysis/sunrise/h258/h258.cosmo50cmb.1536g2bwK.00264.1/broadband.fits'
file2 = '/astro/net/scratch1/fabio/h258/h258.cosmo50cmb.1536g2bwK.00264.1.v3/broadband.fits'
filters2 = mrdfits(file2,10)
file1 = '/astro/net/scratch1/fabio/h258/h258.cosmo50cmb.1536g2bwK.00264.1/broadband.fits'
filters1 = mrdfits(file1,12)
;rint,filters1[2:6].filter
;rint,filters2[2:6].filter
unscatmag1 = filters1[1:4].AB_MAG_NONSCATTER0
scatmag1 =   filters1[1:4].AB_MAG0
unscatmag2 = filters2[1:4].AB_MAG_NONSCATTER0
scatmag2 =   filters2[1:4].AB_MAG0
print,'h258.264'
print,'Version 2'
print,filters1[1].filter,filters1[1].AB_MAG_NONSCATTER0,filters1[1].AB_MAG0
print,filters1[2].filter,filters1[2].AB_MAG_NONSCATTER0,filters1[2].AB_MAG0
print,filters1[3].filter,filters1[3].AB_MAG_NONSCATTER0,filters1[3].AB_MAG0
print,filters1[4].filter,filters1[4].AB_MAG_NONSCATTER0,filters1[4].AB_MAG0
print,''
print,'Version 3'
print,filters2[1].filter,filters2[1].AB_MAG_NONSCATTER0,filters2[1].AB_MAG0
print,filters2[2].filter,filters2[2].AB_MAG_NONSCATTER0,filters2[2].AB_MAG0
print,filters2[3].filter,filters2[3].AB_MAG_NONSCATTER0,filters2[3].AB_MAG0
print,filters2[4].filter,filters2[4].AB_MAG_NONSCATTER0,filters2[4].AB_MAG0

device,filename='~/Scratch2/Sunrise/h258.264_filters.eps',/color,bits_per_pixel=8,/times
plot,[-23.5,-21],[-23.5,-21],xrange=[-23.5,-21],yrange = [-23.5,-21],xtitle = 'Magnitude without Dust',ytitle = 'Magnitude with Dust',title = 'Face on Magnitudes'
oplot,unscatmag1,scatmag1,psym = 2,color=240
oplot,unscatmag2,scatmag2,psym = 4,color=50
legend,['h258.264 v2','h258.264 v3'],psym=[2,4],color=[240,50]
device,/close

END
