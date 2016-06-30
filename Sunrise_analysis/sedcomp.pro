pro sedcomp
loadct,39
set_plot,'x'
;set_plot,'ps'

file2 = '~/Scratch2/Sunrise/MW1lr_v2/mcrx.fits'
data2 = mrdfits(file2,21,h)
file3 = '~/Scratch2/Sunrise/MW1lr_v3/mcrx.fits'
data3 = mrdfits(file3,27,h)
plot,data2.lambda,data2.L_LAMBDA,/ylog,/xlog,xtitle = 'lambda',ytitle='luminosity',xrange=[1e-8,1]
oplot,data3.lambda,data3.L_LAMBDA_absorbed,color = 240,linestyle = 2
oplot,data3.lambda,data3.L_LAMBDA_nonscatter0,color = 100,linestyle = 2
oplot,data3.lambda,data3.L_LAMBDA_scatter0,color = 50,linestyle = 2

oplot,data2.lambda,data2.L_LAMBDA_absorbed,color = 240
oplot,data2.lambda,data2.L_LAMBDA_nonscatter0,color = 100
oplot,data2.lambda,data2.L_LAMBDA_scatter0,color = 50
legend,['v2, Original SED','v2, Absorbed','v2, Scatter','v3, Original SED','v3, Absorbed','v3, Scatter'],color = [100,240,50,100,240,50],linestyle = [0,0,0,2,2,2],/bottom

plot,data2.lambda,data2.L_LAMBDA,/ylog,/xlog,xtitle = 'lambda',ytitle='luminosity',xrange=[1e-7,1e-6],yrange=[1d42,3d44]
oplot,data3.lambda,data3.L_LAMBDA_absorbed,color = 240,linestyle = 2
oplot,data3.lambda,data3.L_LAMBDA_nonscatter0,color = 100,linestyle = 2
oplot,data3.lambda,data3.L_LAMBDA_scatter0,color = 50,linestyle = 2

oplot,data2.lambda,data2.L_LAMBDA_absorbed,color = 240
oplot,data2.lambda,data2.L_LAMBDA_nonscatter0,color = 100
oplot,data2.lambda,data2.L_LAMBDA_scatter0,color = 50

oplot,[3720.0d-10,3720.0d-10],[1.0d30,1.0d45],linestyle = 1,color = 140;SDSS g
oplot,[5720.0d-10,5720.0d-10],[1.0d30,1.0d45],linestyle = 1,color = 140

oplot,[5220.0d-10,5220.0d-10],[1.0d30,1.0d45],linestyle = 1,color = 220;SDSS r
oplot,[7220.0d-10,7220.0d-10],[1.0d30,1.0d45],linestyle = 1,color = 220

oplot,[3020.00d-10,3020.00d-10],[1.0d30,1.0d45],linestyle = 1,color = 20;SDSS u
oplot,[4120.0d-10,4120.0d-10],[1.0d30,1.0d45],linestyle = 1,color = 20
legend,['v2, Original SED','v2, Absorbed','v2, Scatter','v3, Original SED','v3, Absorbed','v3, Scatter'],color = [100,240,50,100,240,50],linestyle = [0,0,0,2,2,2],/bottom


set_plot,'ps'
device,filename='absvsScat.eps',/color,bits_per_pixel=8,/times
plot,[1e-7,1e-6],[1,1],/ylog,/xlog,xtitle = 'lambda [m]',ytitle='Fraction of Total Luminosity',xrange=[1e-7,1e-6],yrange=[0.001,1],charsize = 2
oplot,data3.lambda,data3.L_LAMBDA_absorbed/data3.L_LAMBDA_nonscatter0,color = 240,linestyle = 2
;oplot,data3.lambda,data3.L_LAMBDA_nonscatter0,color = 100,linestyle = 2
oplot,data3.lambda,(data3.L_LAMBDA_nonscatter0 - data3.L_LAMBDA_scatter0)/data3.L_LAMBDA_nonscatter0,color = 50,linestyle = 2
oplot,data2.lambda,data2.L_LAMBDA_absorbed/data2.L_LAMBDA_nonscatter0,color = 240
;oplot,data2.lambda,data2.L_LAMBDA_nonscatter0,color = 100
oplot,data2.lambda,(data2.L_LAMBDA_nonscatter0 - data2.L_LAMBDA_scatter0)/data2.L_LAMBDA_nonscatter0,color = 50
oplot,[3720.0d-10,3720.0d-10],[0.001,1],linestyle = 1,color = 140;SDSS g
oplot,[5720.0d-10,5720.0d-10],[0.001,1],linestyle = 1,color = 140

oplot,[5220.0d-10,5220.0d-10],[0.001,1],linestyle = 1,color = 220;SDSS r
oplot,[7220.0d-10,7220.0d-10],[0.001,1],linestyle = 1,color = 220

oplot,[3020.00d-10,3020.00d-10],[0.001,1],linestyle = 1,color = 20;SDSS u
oplot,[4120.0d-10,4120.0d-10],[0.001,1],linestyle = 1,color = 20
legend,['v2: Absorbed/Total','v2: Scatter/Total','v3: Absorbed/Total','v3: Scatter/Total'],color = [240,50,240,50],linestyle = [0,0,2,2],/bottom,charsize = 1.5
device,/close
stop
end
