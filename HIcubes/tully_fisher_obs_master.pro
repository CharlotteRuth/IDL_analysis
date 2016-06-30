
PRO tully_fisher_obs_master
prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

broadband = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.1/broadband.fits',$
             'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.1/broadband.fits']
tfiles = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
         'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
msol_per_sysmass = [2.310e15,2.310e15]
key = ['no '+textoidl('H_2'), textoidl('H_2')]
key = ['DnoH2','DH2']
outfile = '~/plots/h516.cosmo25cmb.paper'
velocities_true = [53.2281,56.8234]
velocities = velocities_true ;taken from flat part of rotation curve
velocitiesHI = [53.406370,54.337259] ;from HI line widths

IF 0 THEN BEGIN
    broadband = ['h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00324.dir/h603.cosmo50cmb.3072g14HBWK.00324.1/broadband.fits']
    tfiles =    ['h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00324.dir/h603.cosmo50cmb.3072g14HBWK.00324.halo.1']
    msol_per_sysmass = [1.84793e16]
    key = ['h603 H2']
    outfile = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00324.dir/h603.cosmo50cmb.3072g14HBWK.00324'
    velocities_true = [110]
    velocities = velocities_true ;taken from flat part of rotation curve
    velocitiesHI = [121.4]      ;from HI line widths
ENDIF

gmass = tully_fisher_obs_gasmass(prefix + tfiles, msol_per_sysmass,smass = smass,gmassall = gmassall)
print,'True HI and H2 Gas Mass:     ',gmass
print,'Total Gas Mass:              ',gmassall
print,'True Stellar Mass:           ',smass
tully_fisher_obs,prefix + broadband, velocitiesHI, gmass, key = key,outfile = outfile

tully_fisher_obs_btf, prefix + broadband, velocities, gmass, key = key,smass_true = smass,velocities_true = velocities_true,outfile = outfile
END
