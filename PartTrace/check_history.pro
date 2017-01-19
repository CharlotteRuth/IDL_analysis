;Check to make sure that the allgas.history files are consistent with
;the outputs in the directories

dir = '/nobackupp8/crchrist/MolecH/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
filebase = 'h799.cosmo25cmb.3072g14HBWK.00512/h799.cosmo25cmb.3072g14HBWK.00512'
dir = '/nobackupp8/crchrist/MolecH/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
filebase = 'h516.cosmo25cmb.3072g14HBWK.00512/h516.cosmo25cmb.3072g14HBWK.00512'
dir = '/nobackupp8/crchrist/MolecH/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
filebase = 'h603.cosmo50cmb.3072g14HBWK.00512/h603.cosmo50cmb.3072g14HBWK.00512'
dir = '/nobackupp8/crchrist/MolecH/h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
filebase = 'h285.cosmo50cmb.3072g14HMbwK.00512/h285.cosmo50cmb.3072g14HMbwK.00512'
dir = '/nobackupp8/crchrist/MolecH/h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
filebase = 'h258.cosmo50cmb.3072g14HMbwK.00512/h258.cosmo50cmb.3072g14HMbwK.00512'
dir = '/nobackupp8/crchrist/MolecH/h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
filebase = 'h239.cosmo50cmb.3072g14HMbwK.00512/h239.cosmo50cmb.3072g14HMbwK.00512'
haloid = '1'
cd,dir
rtipsy,dir + filebase,h,g,d,s
readarr,dir + filebase + '.iord',h,iord,part = 'gas',/ascii
history = mrdfits(dir + filebase + '.grp' + haloid + '.allgas.history.fits',1)
match,iord,history.iord,indt,indh
IF n_elements(where(history.mass NE 0)) NE n_elements(indt) THEN print,'Missing particles'
readarr,dir + filebase + '.FeMassFrac',h,fe,part = 'gas',/ascii
readarr,dir + filebase + '.OxMassFrac',h,ox,part = 'gas',/ascii
metals = 2.09*ox + 1.06*fe
print,'Difference between history metallicity and true metallicity: ',max(metals[indt]-history[indh].metallicity)
