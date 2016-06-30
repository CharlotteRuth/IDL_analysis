FUNCTION lyalphaforest,lambda,z
;does the actual absorption
nlines = 4
Ly_lambda = [1216, 1026, 973, 950]
Ly_A = [0.0036, 1.7e-3, 1.2e-3, 9.3e-4]
Ly_limit = 912
tau=fltarr(N_ELEMENTS(lambda))
l_rest=lambda/(1 + z)
xem = 1 + z

FOR i = 0, nlines-1 DO BEGIN
    ind = where(l_rest lt Ly_lambda[i])
    tau[ind] = tau[ind] + Ly_A[i]*(lambda/Ly_lambda[i])^(3.46)
ENDFOR
xc = lambda/Ly_limit
xc[where(xc lt 1)] = 1

ind = where(xem gt xc)
tau[ind] = tau[ind] + 0.25*xc[ind]^3*(xem^0.46 - xc[ind]^0.46)
tau[ind] = tau[ind] + 9.40*xc[ind]^1.5*(xem^0.18 - xc[ind]^0.18)
tau[ind] = tau[ind]  -0.70*xc[ind]^3*(xc[ind]^(-0.32) - xem^(-0.32))
tau[ind] = tau[ind]  -0.0023*(xem^1.68 - xc[ind]^1.68)

absorb = exp(-1.0*tau)
RETURN,absorb
END

PRO lyman_absorb,filename,z,outfile
;lyman_absorb,'mcrx_hz3.fits',3,'lyman_absorb.eps'
;This does post-processing to add in lyman absorption and redshift
loadct,39
angstrom = '!6!sA!r!u!9 %!6!n' 
Ly_limit = 912*(1+z)
Ly_lambda = [1216, 1026, 973, 950]*(1+z)
unit_conv = 1.0/(3.827d26*(3.24077649d-17)^2)*1d-10
;Convert from Watts/meter/m^2/sterradian to Solar Luminosities per
;per parsec^2 per angstoms
lambda = MRDFITS(filename,4,header)            ;This will change, depending on your run!
sed_scat = MRDFITS(filename,28,header)            ;This will change, depending on your run!
print,'Data File: ',filename
help,sed_scat,/struct
length = (size(sed_scat.compressed_data))[2]
lengthgrid = LONG(N_ELEMENTS(sed_scat.COMPRESSED_DATA))
flattened_sed_scat = fltarr(length)
FOR ct = 0, length-1 DO BEGIN
    flattened_sed_scat[ct] = total(sed_scat[ct].compressed_data)
ENDFOR
flattened_sed_scat = flattened_sed_scat*unit_conv
lambda.lambda = lambda.lambda*1d10 ;put in terms of angstroms
set_plot,'x'
plot,lambda.lambda,flattened_sed_scat,xtitle = 'Wavelength ('+angstrom+')',/xlog,xrange=[100,1e5],ytitle = 'Flux (arbitrary units)'
;'Flux (L'+sunsymbol()+' '+angstrom+TeXtoIDL("^{-1} parsecs^{-2}")+')'

lambda_z = lambda.lambda*(1+z)
absorb = lyalphaforest(lambda_z,z)
 
flux =  flattened_sed_scat*absorb

oplot,lambda_z,flattened_sed_scat,color = 50
oplot,lambda_z,flux,color = 240
;oplot,[Ly_limit,Ly_limit],[1,1e10],linestyle = 2
;FOR ct = 0, N_ELEMENTS(Ly_lambda) - 1 DO oplot,[Ly_lambda[ct],Ly_lambda[ct]],[1,1e10],linestyle = 1
Legend,['Rest Frame','Redshifted, z = 3','Lyman Absorption, z = 3'],color=[0,50,240],linestyle=[0,0,0],/bottom,/right

set_plot,'x'
set_plot,'ps'
device,filename=outfile,/color,bits_per_pixel=8
plot,lambda.lambda,flattened_sed_scat,xtitle = 'Wavelength ('+angstrom+')',/xlog,xrange=[100,1e5],ytitle = 'Flux (arbitrary units)'
;'Flux (L'+sunsymbol()+' '+angstrom+TeXtoIDL("^{-1} parsecs^{-2}")+')'
oplot,lambda_z,flattened_sed_scat,color = 50
oplot,lambda_z,flux  ,color = 240
;oplot,[Ly_limit,Ly_limit],[1,1e10],linestyle = 2
;FOR ct = 0, N_ELEMENTS(Ly_lambda) - 1 DO oplot,[Ly_lambda[ct],Ly_lambda[ct]],[1,1e10],linestyle = 1
Legend,['Rest Frame','Redshifted, z = 3','Lyman Absorption, z = 3'],color=[0,50,240],linestyle=[0,0,0],/bottom,/right

device,/close

END
