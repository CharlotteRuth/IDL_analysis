maxdistance = 4.0
linestyles = [0,0]
ctables = [39,39]
colors = [50,254]
baseplots = '/home6/crchrist/Plots/'
baseplots = '/home/christensen/Plots/'

;-----------------------------------------
dir = '/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
base = 'h799.cosmo25cmb.3072g14HBWK'

dist_units = 25000
mass_units = 2.310e15

step1 = '00152'
halo1a = '1'
halo1b = '4'
step2 = '00172'
halo2 = '1'

yrange = [0,50]

;------------------------------------------
dir = '/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
base = 'h516.cosmo25cmb.3072g14HBWK'

dist_units = 25000
mass_units = 2.310e15

step1 = '00224'
halo1a = '1'
halo1b = '3'
step2 = '00260'
halo2 = '1'

yrange = [0,55]

;----------------------------------------------------------------------
dir1 = dir + base + '.' + step1 + '/'
file1 = dir + base + '.' + step1 + '/' + base + '.' + step1
dir2 = dir + base + '.' + step2 + '/'
file2 = dir + base + '.' + step2 + '/' + base + '.' + step2

cd,dir1
tipsysatshi,base + '.' + step1,halo1a,dist_units,mass_units,/cutout_rad
cd,dir2
tipsysatshi,base + '.' + step2,halo2,dist_units,mass_units,/cutout_rad

vcirc,[file2 + '.halo.' + halo2 + '.std',file1 + '.halo.' + halo1a + '.std'],[mass_units,mass_units],[dist_units,dist_units],keys = [step2,step1],linestyle = linestyles,color = colors,maxdistance = maxdistance,type = 'all',yrange = yrange,ctables = ctables,outfile = baseplots + base + '.' + halo2 + '_vcirc.eps'

outfile = baseplots + base + '.' + halo2 + '_sprof.eps'
formatplot,outplot = outfile
loadct,39
IF keyword_set(outfile) THEN device,filename=outfile,/color,bits_per_pixel= 8 ELSE window,0
rtipsy,file2 + '.halo.' + halo2 + '.std',h,g,d,s
s_prof = prof(s, 'star', h.time, nbins = nbins,rmax = maxdistance/dist_units)
s_prof.rbins = s_prof.rbins*dist_units
s_prof.rho = s_prof.rho*mass_units/dist_units/dist_units  ; Msol/kpc^2
plot,s_prof.rbins,s_prof.rho,/ylog,xtitle = 'Radius [kpc]',ytitle = 'Sigma_star [Msol/kpc^2]',/nodata
oplot,s_prof.rbins,s_prof.rho,color = colors[0],linestyle = linestyles[0]
rtipsy,file1 + '.halo.' + halo1a + '.std',h,g,d,s
s_prof = prof(s, 'star', h.time, nbins = nbins,rmax = maxdistance/dist_units)
s_prof.rbins = s_prof.rbins*dist_units
s_prof.rho = s_prof.rho*mass_units/dist_units/dist_units  ; Msol/kpc^2
oplot,s_prof.rbins,s_prof.rho,color = colors[1],linestyle = linestyles[1]
legend,[step2,step1],color = colors,linestyle = linestyles,/right
IF keyword_set(outfile) THEN device,/close

