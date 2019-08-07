filebase1='/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g1HBWK/h516.cosmo25cmb.2304g1HBWK'
filebase2='/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g1HBWKS/h516.cosmo25cmb.2304g1HBWKS'
step = '000148'
rtipsy,filebase1+'.'+step,h1,g1,d1,s1
rtipsy,filebase2+'.'+step,h2,g2,d2,s2
sl1 = rstarlog(filebase1+'.starlog',/molecularH,/big)
sl2 = rstarlog(filebase2+'.starlog',/molecularH,/big)
filebase1a='/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1HBWK/h516.cosmo25cmb.3072g1HBWK'
filebase2a='/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1HBWKS/h516.cosmo25cmb.3072g1HBWKS'
stepa = '000084'
rtipsy,filebase1a+'.'+stepa,h1a,g1a,d1a,s1a
rtipsy,filebase1a+'.'+'000116',h2a,g2a,d2a,s2a
filebase0 = '/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK_2.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512'
rtipsy,filebase0,h0,g0,d0,s0
units = tipsyunits(filebase1+'.param')

formatplot,/outplot
device,filename = '~/h516.shield.SFH.eps',bits_per_pixel= 8,/times,/color
histogramp,s2.tform,weight = s2.mass
histogramp,s1.tform,weight = s1.mass,/overplot,linestyle = '2'
histogramp,s2a.tform,weight = s2a.mass,/overplot,color = 60
histogramp,s1a.tform,weight = s1a.mass,/overplot,color = 60,linestyle = 2
histogramp,s0.tform,weight = s0.mass,/overplot,color = 100,linestyle = 1
legend,['2306, H2','2306, Shield','3072, H2','3072, Shield','3072, g14'],linestyle = [2,0,2,0,1],color = [0,0,60,60,100]
device,/close

maxtime = max(sl1a.timeform)
window,0
plot,sl1[where(sl1.timeform LT maxtime)].rhoform*units.rhounit,sl1[where(sl1.timeform LT maxtime)].tempform,psym = 3,/ylog,/xlog
oplot,sl2[where(sl2.timeform LT maxtime)].rhoform*units.rhounit,sl2[where(sl2.timeform LT maxtime)].tempform,psym = 3,color = 100
oplot,sl1[where(sl1.timeform LT maxtime)].rhoform*units.rhounit,sl1[where(sl1.timeform LT maxtime)].tempform,psym = 3

window,1
plot,sl1a[where(sl1a.timeform LT maxtime)].rhoform*units.rhounit,sl1a[where(sl1a.timeform LT maxtime)].tempform,psym = 3,/ylog,/xlog
oplot,sl2a[where(sl2a.timeform LT maxtime)].rhoform*units.rhounit,sl2a[where(sl2a.timeform LT maxtime)].tempform,psym = 3,color = 100
oplot,sl1a[where(sl1a.timeform LT maxtime)].rhoform*units.rhounit,sl1a[where(sl1a.timeform LT maxtime)].tempform,psym = 3
