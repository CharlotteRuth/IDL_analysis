pro two_star_grid
loadct,39
set_plot,'x'
;set_plot,'ps'
base = '/astro/net/scratch2/christensen/Sunrise/onestar_test/test1e'
f1SR = base +'7_0.0004/broadband.fits'
f1IDL = base +'7_0.0004/galaxy.star99_ab.Mv.fits'
filters1SR = (mrdfits(f1SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f1IDL,1)
filters1IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f2SR = base +'7_0.004/broadband.fits'
f2IDL = base +'7_0.004/galaxy.star99_ab.Mv.fits'
filters2SR = (mrdfits(f2SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f2IDL,1)
filters2IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f3SR = base +'7_0.019/broadband.fits'
f3IDL = base +'7_0.019/galaxy.star99_ab.Mv.fits'
filters3SR = (mrdfits(f3SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f3IDL,1)
filters3IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f4SR = base +'8_0.0004/broadband.fits'
f4IDL = base +'8_0.0004/galaxy.star99_ab.Mv.fits'
filters4SR = (mrdfits(f4SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f4IDL,1)
filters4IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f5SR = base +'8_0.004/broadband.fits'
f5IDL = base +'8_0.004/galaxy.star99_ab.Mv.fits'
filters5SR = (mrdfits(f5SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f5IDL,1)
filters5IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f6SR = base +'8_0.019/broadband.fits'
f6IDL = base +'8_0.019/galaxy.star99_ab.Mv.fits'
filters6SR = (mrdfits(f6SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f6IDL,1)
filters6IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]


f7SR = base +'9_0.0004/broadband.fits'
f7IDL = base +'9_0.0004/galaxy.star99_ab.Mv.fits'
filters7SR = (mrdfits(f7SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f7IDL,1)
filters7IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f8SR = base +'9_0.004/broadband.fits'
f8IDL = base +'9_0.004/galaxy.star99_ab.Mv.fits'
filters8SR = (mrdfits(f8SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f8IDL,1)
filters8IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f9SR = base +'9_0.019/broadband.fits'
f9IDL = base +'9_0.019/galaxy.star99_ab.Mv.fits'
filters9SR = (mrdfits(f9SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f9IDL,1)
filters9IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f10SR = base +'10_0.0004/broadband.fits'
f10IDL = base +'10_0.0004/galaxy.star99_ab.Mv.fits'
filters10SR = (mrdfits(f10SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f10IDL,1)
filters10IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f11SR = base +'10_0.004/broadband.fits'
f11IDL = base +'10_0.004/galaxy.star99_ab.Mv.fits'
filters11SR = (mrdfits(f11SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f11IDL,1)
filters11IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]

f12SR = base +'10_0.019/broadband.fits'
f12IDL = base +'10_0.019/galaxy.star99_ab.Mv.fits'
filters12SR = (mrdfits(f12SR,10))[7:13].AB_mag_nonscatter0
temp = mrdfits(f12IDL,1)
filters12IDL = [temp.u,temp.b,temp.v,temp.r,temp.i,temp.k]



filter = 2 ;v
print,'V'
print,'----------------------------------'
print,filters1IDL[filter],'-',filters1SR[filter],'     ',filters2IDL[filter],'-',filters2SR[filter],'     ',filters3IDL[filter],'-',filters3SR[filter]
print,filters4IDL[filter],'-',filters4SR[filter],'     ',filters5IDL[filter],'-',filters5SR[filter],'     ',filters6IDL[filter],'-',filters6SR[filter]
print,filters7IDL[filter],'-',filters7SR[filter],'     ',filters8IDL[filter],'-',filters8SR[filter],'     ',filters9IDL[filter],'-',filters9SR[filter]
print,filters10IDL[filter],'-',filters10SR[filter],'     ',filters11IDL[filter],'-',filters11SR[filter],'     ',filters12IDL[filter],'-',filters12SR[filter]

print,''
filter = 1 ;b
print,'B'
print,'----------------------------------'
print,filters1IDL[filter],'-',filters1SR[filter],'     ',filters2IDL[filter],'-',filters2SR[filter],'     ',filters3IDL[filter],'-',filters3SR[filter]
print,filters4IDL[filter],'-',filters4SR[filter],'     ',filters5IDL[filter],'-',filters5SR[filter],'     ',filters6IDL[filter],'-',filters6SR[filter]
print,filters7IDL[filter],'-',filters7SR[filter],'     ',filters8IDL[filter],'-',filters8SR[filter],'     ',filters9IDL[filter],'-',filters9SR[filter]
print,filters10IDL[filter],'-',filters10SR[filter],'     ',filters11IDL[filter],'-',filters11SR[filter],'     ',filters12IDL[filter],'-',filters12SR[filter]


print,''
filter = 1 ;b
print,'B, Vega'
print,'----------------------------------'
print,filters1IDL[filter]+ 0.163,'-',filters1SR[filter]+ 0.163,'     ',filters2IDL[filter]+ 0.163,'-',filters2SR[filter]+ 0.163,'     ',filters3IDL[filter]+ 0.163,'-',filters3SR[filter]+ 0.163
print,filters4IDL[filter]+ 0.163,'-',filters4SR[filter]+ 0.163,'     ',filters5IDL[filter]+ 0.163,'-',filters5SR[filter]+ 0.163,'     ',filters6IDL[filter]+ 0.163,'-',filters6SR[filter]+ 0.163
print,filters7IDL[filter]+ 0.163,'-',filters7SR[filter]+ 0.163,'     ',filters8IDL[filter]+ 0.163,'-',filters8SR[filter]+ 0.163,'     ',filters9IDL[filter]+ 0.163,'-',filters9SR[filter]+ 0.163
print,filters10IDL[filter]+ 0.163,'-',filters10SR[filter]+ 0.163,'     ',filters11IDL[filter+ 0.163],'-',filters11SR[filter]+ 0.163,'     ',filters12IDL[filter]+ 0.163,'-',filters12SR[filter]+ 0.163
stop

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
