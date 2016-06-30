dist_units = 50000
mass_units = 1.84793e16

step1 = '00152'
halo1a = '1'
halo1b = '4'
step2 = '00172'
halo2 = '1'

dir = '/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
base = 'h799.cosmo25cmb.3072g14HBWK'
file = dir + base + '.' + step1 + '/' + base + '.' + step1

tipsysatshi,file,halo1a,dist_units,mass_units,/cutout_rad
tipsysatshi,file,halo2,dist_units,mass_units,/cutout_rad
