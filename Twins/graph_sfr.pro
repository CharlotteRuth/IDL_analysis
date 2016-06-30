pro graph_sfr

datafile = '/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512/h603.cosmo50cmb.2304g2bwK.00512.1.std'
;outdir = '/astro/net/scratch1/christensen/Twins/h603.cosmo50cmb.2304g2bwK'
outdir = '/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512/h603.cosmo50cmb.2304g2bwK.00512.1'

massunit = 1.84793e16
timeunit = 3.8785614e+10
rtipsy,datafile,h,g,d,s
set_plot,'x'
sfr,s,avesfr,massu=massunit,time=timeunit,OVERPLOT=0,xtitle = 'Time (Billion Years)', ytitle = 'SFR (Solar Masses / yr)',title = 'Star Formation History of h603'
stop
set_plot,'ps'
device,filename = outdir + '/sfh.eps',/color,bits_per_pixel=8
sfr,s,massu=massunit,time=timeunit,OVERPLOT=0,xtitle = 'Time (Billion Years)', ytitle = 'SFR (Solar Masses / yr)',title = 'Star Formation History of h603'
device,/close
print,'Average SFR: ',avesfr
stop
end
