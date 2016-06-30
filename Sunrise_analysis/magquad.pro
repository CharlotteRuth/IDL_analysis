pro magquad
loadct,39
set_plot,'x'
;set_plot,'ps'

file1 = '~/Scratch2/Sunrise/MW1lr_v3_mc1_multi1/broadband.fits'
filters1 = mrdfits(file1,11)
file2 = '~/Scratch2/Sunrise/MW1lr_v3_mc1_multi0/broadband.fits'
filters2 = mrdfits(file2,11)
file3 = '~/Scratch2/Sunrise/MW1lr_v3_mc0_multi1/broadband_uncommented.fits'
filters3 = mrdfits(file3,11)
file4 = '~/Scratch2/Sunrise/MW1lr_v3_mc0_multi0/broadband.fits'
filters4 = mrdfits(file4,11)

file5 = '~/Scratch2/Sunrise/MW1lr_v2/broadband.fits'
filters5 = mrdfits(file5,12)

unscatmag1 = filters1[2:6].AB_MAG_NONSCATTER0
scatmag1 =   filters1[2:6].AB_MAG0
unscatmag2 = filters2[2:6].AB_MAG_NONSCATTER0
scatmag2 =   filters2[2:6].AB_MAG0
unscatmag3 = filters3[2:6].AB_MAG_NONSCATTER0
scatmag3 =   filters3[2:6].AB_MAG0
unscatmag4 = filters4[2:6].AB_MAG_NONSCATTER0
;scatmag4 =   filters4[2:6].AB_MAG0
unscatmag5 = filters5[2:6].AB_MAG_NONSCATTER0
scatmag5 =   filters5[2:6].AB_MAG0

x = indgen(N_ELEMENTS(unscatmag1))
print,''
print,'MW1lr','                            NonScatttered','       Scattered'
print,'mc1_multi1'
print,filters1[2].filter,filters1[2].AB_MAG_NONSCATTER0,filters1[2].AB_MAG0
print,filters1[3].filter,filters1[3].AB_MAG_NONSCATTER0,filters1[3].AB_MAG0
print,filters1[4].filter,filters1[4].AB_MAG_NONSCATTER0,filters1[4].AB_MAG0
print,filters1[5].filter,filters1[5].AB_MAG_NONSCATTER0,filters1[5].AB_MAG0
print,filters2[6].filter,filters1[6].AB_MAG_NONSCATTER0,filters1[6].AB_MAG0
print,''
print,'mc1_multi0'
print,filters2[2].filter,filters2[2].AB_MAG_NONSCATTER0,filters2[2].AB_MAG0
print,filters2[3].filter,filters2[3].AB_MAG_NONSCATTER0,filters2[3].AB_MAG0
print,filters2[4].filter,filters2[4].AB_MAG_NONSCATTER0,filters2[4].AB_MAG0
print,filters2[5].filter,filters2[5].AB_MAG_NONSCATTER0,filters2[5].AB_MAG0
print,filters2[6].filter,filters2[6].AB_MAG_NONSCATTER0,filters2[6].AB_MAG0
print,''
print,'mc0_multi1'
print,filters3[2].filter,filters3[2].AB_MAG_NONSCATTER0,filters3[2].AB_MAG0
print,filters3[3].filter,filters3[3].AB_MAG_NONSCATTER0,filters3[3].AB_MAG0
print,filters3[4].filter,filters3[4].AB_MAG_NONSCATTER0,filters3[4].AB_MAG0
print,filters3[5].filter,filters3[5].AB_MAG_NONSCATTER0,filters3[5].AB_MAG0
print,filters3[6].filter,filters3[6].AB_MAG_NONSCATTER0,filters3[6].AB_MAG0
print,''
print,'mc0_multi0'
print,filters4[2].filter,filters4[2].AB_MAG_NONSCATTER0;,filters4[2].AB_MAG0
print,filters4[3].filter,filters4[3].AB_MAG_NONSCATTER0;,filters4[3].AB_MAG0
print,filters4[4].filter,filters4[4].AB_MAG_NONSCATTER0;,filters4[4].AB_MAG0
print,filters4[5].filter,filters4[5].AB_MAG_NONSCATTER0;,filters4[5].AB_MAG0
print,filters4[6].filter,filters4[6].AB_MAG_NONSCATTER0;,filters4[6].AB_MAG0
print,''
print,'Version 2'
print,filters5[2].filter,filters5[2].AB_MAG_NONSCATTER0,filters5[2].AB_MAG0
print,filters5[3].filter,filters5[3].AB_MAG_NONSCATTER0,filters5[3].AB_MAG0
print,filters5[4].filter,filters5[4].AB_MAG_NONSCATTER0,filters5[4].AB_MAG0
print,filters5[5].filter,filters5[5].AB_MAG_NONSCATTER0,filters5[5].AB_MAG0
print,filters5[6].filter,filters5[6].AB_MAG_NONSCATTER0,filters5[6].AB_MAG0


!X.style = 1
!Y.style = 1
brightrange = -23
faintrange = -20;-19

;window,0
set_plot,'ps'
device,filename='MW1lr_faceon_SP.eps',/color,bits_per_pixel=8,/times
plot,x,unscatmag4,xrange=[-0.5,4.5],xtitle = 'SDSS Filters',ytitle = 'Magnitudes',title = 'Single-Phase ISM',xtickname = ['u','g','r','i','z'],psym = 1,charsize = 2,symsize = 2,yrange = [faintrange,brightrange]
oplot,x,unscatmag4,psym = 1,color=20,symsize = 2
oplot,x,scatmag2,psym = 7,color=240,symsize = 2
oplot,x,unscatmag5,psym = 4,color=150,symsize = 2
oplot,x,scatmag5,psym = 4,color=200,symsize = 2

;legend,['MC1, Multi1','MC1, Multi0','MC0, Multi1','MC0, Multi0'],psym=[7,4,7,4],color=[240,240,50,50],/right,charsize = 1.5
legend,['v2 Unscat','v2 Scat','v3 no MC, unscat','v3, MC, scat'],psym=[4,4,1,7],color=[150,200,20,240],/right,charsize = 1.5,/bottom
device,/close
;stop

;window,1
device,filename='MW1lr_faceon_v2VSv3.eps',/color,bits_per_pixel=8,/times
plot,x,unscatmag4,xrange=[-0.5,4.5],xtitle = 'SDSS Filters',ytitle = 'Magnitudes',title = 'Comparing v2 to v3',xtickname = ['u','g','r','i','z'],psym = 1,charsize = 2,symsize = 2,yrange = [faintrange,brightrange]
oplot,x,unscatmag4,psym = 1,color=20,symsize = 2
;oplot,x,scatmag4,psym = 1,color=240,symsize = 2
oplot,x,unscatmag5,psym = 4,color=150,symsize = 2
oplot,x,scatmag5,psym = 4,color=200,symsize = 2
legend,['v2 Unscat','v2 Scat','v3, no MC, SP, Unscat','v3 no MC, SP, Scat'],psym=[4,4,1,1],color=[150,200,20,240],/right,charsize = 1.5,/bottom
device,/close
;stop


;window,2
device,filename='MW1lr_faceon_v3ScatUnscat.eps',/color,bits_per_pixel=8,/times
plot,x,unscatmag4,xrange=[-0.5,4.5],xtitle = 'SDSS Filters',ytitle = 'Magnitudes',title = 'v3 Unscat to Scat',xtickname = ['u','g','r','i','z'],psym = 1,charsize = 2,symsize = 2,yrange = [faintrange,brightrange]
oplot,x,unscatmag4,psym = 1,color=20,symsize = 2
oplot,x,scatmag1,psym = 7,color=50,symsize = 2
;oplot,x,unscatmag5,psym = 7,color=120,symsize = 2
;oplot,x,scatmag5,psym = 1,color=120,symsize = 2
legend,['v3 Unscat','v3 Scat (MC, MP)'],psym=[1,7],color=[20,50],/right,charsize = 1.5,/bottom
;legend,['v2 Unscat','v2 Scat','v3 Unscat','v3 Scat'],psym=[7,6,7,1],color=[120,120,240,50],/right,charsize =1.5,/bottom
device,/close
;stop

;window,3
device,filename='MW1lr_faceon_compMC.eps',/color,bits_per_pixel=8,/times
plot,x,unscatmag4,xrange=[-0.5,4.5],xtitle = 'SDSS Filters',ytitle = 'Magnitudes',title = 'Comparing MC and No MC',xtickname = ['u','g','r','i','z'],psym = 1,charsize = 2,symsize = 2,yrange = [faintrange,brightrange]
oplot,x,unscatmag3,psym = 1,color=20,symsize = 2
oplot,x,unscatmag1,psym = 7,color=20,symsize = 2
;stop
;oplot,x,unscatmag4,psym = 1,color=240,symsize = 2
;oplot,x,unscatmag2,psym = 6,color=240,symsize = 2
oplot,x,scatmag3,psym = 1,color=50,symsize = 2
oplot,x,scatmag1,psym = 7,color=50,symsize = 2
;legend,['No MC (unscat)','MC (unscat)','No MC (scat SP)','MC (scat SP)','No MC (scat MP)','MC (scat MP)'],psym=[7,4,1,6,1,6],color=[50,50,240,240,50,50],/right,charsize = 1.5,/bottom
legend,['No MC (unscat)','MC (unscat)','No MC (scat MP)','MC (scat MP)'],psym=[1,7,1,7],color=[20,20,50,50],/right,charsize = 1.5,/bottom
;legend,['No MC (unscat)','MC (unscat)'],psym=[7,4],color=[50,50],/right,charsize = 1.5,/bottom
device,/close
;stop

;window,4
device,filename='MW1lr_faceon_compMedia.eps',/color,bits_per_pixel=8,/times
plot,x,unscatmag1,xrange=[-0.5,4.5],xtitle = 'SDSS Filters',ytitle = 'Magnitudes',title = 'Comparing Single-Phase and Multi-Phase Media',xtickname = ['u','g','r','i','z'],psym = 7,charsize = 2,symsize = 2,yrange = [faintrange,brightrange]
oplot,x,unscatmag1,psym = 7,color=20,symsize = 2
oplot,x,scatmag2,psym = 7,color=240,symsize = 2
oplot,x,scatmag1,psym = 7,color=50,symsize = 2
legend,['Unscatered','SP','MP'],psym=[7,7,7],color=[20,240,50],/right,charsize = 1.5,/bottom
device,/close
;stop

END
