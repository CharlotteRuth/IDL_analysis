;dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir'
;file = 'h516.cosmo25cmb.1536g3HBWK.00512.halo.1'
;pfile = '../../h516.cosmo25cmb.1536g3HBWK.param'

;After using cube on a galaxy, this will smooth it, remove gas below
;the observable limit and make the moments

;Best run on Elektra

pro make_smooth_mom,dir,file,pfile
cd,dir
units = tipsyunits(pfile)
cube = read_cube_fits(file + '.cube.fits',header,kpcunit = 1.0,munit = 1.0)
smoothed = smooth_cube(cube,header,telres = 10.0,distance = 5.0,outfile = file + '.cube.smoothed.fits')
cube = [0]
smoothed[where(smoothed lt 135.218)] = 0
print,'Smoothed Calculated'
moments,smoothed,header,mom0,mom1,mom2,outplot=file+'.cube.smoothed',/fits
print,'Moments Calculated'
end
