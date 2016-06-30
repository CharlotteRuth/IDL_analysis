rtipsy,'h516.cosmo25cmb.1536g3HBWK.00408',h,g,d,s
units = tipsyunits('h516.cosmo25cmb.1536g7HBWK.noFB.param')
readarr,'h516.cosmo25cmb.1536g8HBWK.00408.test.eCool',h,eCool,/ascii
plot,g.tempg,(eCool*g.dens*units.rhounit/6.06912e+23),/ylog,/xlog,psym = 3,xrange = [1e1,1e7],yrange = [1e-30,1e-20]
