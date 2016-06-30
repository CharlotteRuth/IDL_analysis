;CC 12/16/12
;Where in the disk to stars form?

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK'
;base = 'h603.cosmo50cmb.3072g14HBWK'

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK'
;base = 'h516.cosmo25cmb.3072g14HBWK'
;maxr =4.0

PRO starform_radius,dir,base,maxr,molecularH = molecularH,outplot = outplot
ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)
z = REVERSE(findgen(100)*10.0/100.0)
t = ageUniverse - wmap3_lookback(z)

loadct,39
formatplot,outplot = outplot

cd,dir
units = tipsyunits(base+'.param')

readcol,base+'.haloid.dat',files,haloid,format= '(A,F)'
files = REVERSE(files)          ;going forward in time
haloid = REVERSE(haloid)
center = read_center(files = files + '.amiga.stat',halos = haloid)
center[*,2] = center[*,2]*1000.0 - units.lengthunit/2.0 ;"center" is in Mpc for a box that goes from [0,0,0] to [kpcunit/1000,kpcunit/1000,kpcunit/1000]
center[*,3] = center[*,3]*1000.0 - units.lengthunit/2.0
center[*,4] = center[*,4]*1000.0 - units.lengthunit/2.0

rtipsy,files(N_ELEMENTS(files) - 1) + '.halo.1',h,g,d,s
s.tform = s.tform*units.timeunit
s.x = s.x*units.lengthunit
s.y = s.y*units.lengthunit
s.z = s.z*units.lengthunit

starform = rstarlog(base + '.starlog',molecularH = molecularH)
starform.timeform = starform.timeform*units.timeunit
starform.x = starform.x*units.lengthunit
starform.y = starform.y*units.lengthunit
starform.z = starform.z*units.lengthunit
starform.massform = starform.massform*units.massunit

sf_center = extrap_center(starform.timeform/1e9,center)

starform.x = starform.x - sf_center[*,0]
starform.y = starform.y - sf_center[*,1]
starform.z = starform.z - sf_center[*,2]
zstarform = spline(t,z,starform.timeform)
astarform = 1/(1 + zstarform)

radiusform = sqrt(starform.x*starform.x + starform.y*starform.y + starform.z*starform.z)*astarform
radiusnow = sqrt(s.x*s.x + s.y*s.y + s.z*s.z)

dt  = 2.0*1e9
tbins = findgen(7)*dt
colors = (indgen(N_ELEMENTS(tbins)) + 1)*254/N_ELEMENTS(tbins)
if (KEYWORD_SET(outplot)) then device,filename = outplot + '_rform.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
histogramp,radiusform,min = 0, max = maxr
FOR i = 0, N_ELEMENTS(tbins) - 1 DO histogramp,radiusform[where(starform.timeform GE tbins[i] AND starform.timeform LT tbins[i] + dt)],/overplot,color = colors[i]
if (KEYWORD_SET(outplot)) then device,/close ELSE stop

if (KEYWORD_SET(outplot)) then device,filename = outplot + '_rformnorm.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
histogramp,radiusform,min = 0, max = maxr,/normalize,yrange = [0,0.75]
FOR i = 0, N_ELEMENTS(tbins) - 1 DO histogramp,radiusform[where(starform.timeform GE tbins[i] AND starform.timeform LT tbins[i] + dt)],/overplot,color = colors[i],/normalize
if (KEYWORD_SET(outplot)) then device,/close ELSE stop

if (KEYWORD_SET(outplot)) then device,filename = outplot + '_rnow.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
histogramp,radiusnow,min = 0, max = maxr
FOR i = 0, N_ELEMENTS(tbins) - 1 DO histogramp,radiusnow[where(s.tform GE tbins[i] AND s.tform LT tbins[i] + dt)],/overplot,color = colors[i]
if (KEYWORD_SET(outplot)) then device,/close ELSE stop

if (KEYWORD_SET(outplot)) then device,filename = outplot + '_rnownorm.eps',/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2
histogramp,radiusnow,min = 0, max = maxr,/normalize,yrange = [0,0.5]
FOR i = 0, N_ELEMENTS(tbins) - 1 DO histogramp,radiusnow[where(s.tform GE tbins[i] AND s.tform LT tbins[i] + dt)],/overplot,color = colors[i],/normalize
if (KEYWORD_SET(outplot)) then device,/close ELSE stop

END
