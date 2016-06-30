pro darkprofile

files = ['../1E5R/10M/dark0.profile','../1E5R/10M/dark3.profile']

dataread0 = dblarr(12,196)
dataread1 = dblarr(12,196)
close,1
openr,1,files[0]
readf,1,dataread0
close,1
openr,1,files[1]
readf,1,dataread1
close,1

plot,dataread0[0,*],dataread0[1,*],xtitle = 'Radial Bin',ytitle='Dark Matter Density',title = 'Dark Matter Profile',linestyle=0

oplot,dataread1[0,*],dataread1[1,*],linestyle=1

END


