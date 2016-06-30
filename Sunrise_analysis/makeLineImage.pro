;Halpha
;lambda1 = 6.555e-7 6.54e-7
;lambda2 = 6.58e-7
;outfile = 'Halpha.fits'

;NII
;lambda1 = 6.70e-7
;lambda2 = 6.781e-7
;outfile = 'NII.fits'

;fitsfile = 'mcrx.fits'
;extno = 20
;extno = 21
;extno_lambda = 4
pro makeLineImage,fitsfile,extno,lambda1,lambda2,outfile,extno_lambda = extno_lambda,outplot = outplot
loadct,39
formatplot,outplot = outplot

IF NOT keyword_set(extno_lambda) THEN extno_lambda = 4
lc1 = 6.4e-7                    ;Continuum is defined at these wavelengths in Kennicutt 89, page 1097
lc2 = 6.59e-7

sedcube = mrdfits(fitsfile,extno,header)
sed = mrdfits(fitsfile,extno_lambda)
resx = double((strsplit(header[19],' ',/extract))[2])
resy = double((strsplit(header[22],' ',/extract))[2])

temp = MIN(ABS(lc1 - sed.lambda),c1)
temp = MIN(ABS(lc2 - sed.lambda),c2)
cont = (sedcube[*,*,c1] + sedcube[*,*,c2])/2.0
;cont = (sedcube[*,*,47] + sedcube[*,*,51])/2.0

indlambda = where(sed.lambda gt lambda1 AND sed.lambda lt lambda2)
indlambda = [indlambda[0] - 1,indlambda,indlambda[N_ELEMENTS(indlambda) - 1] + 1]
;47:51

;Calculate the intensity to find the sed of the entire simulation
sedintensity = FLTARR(n_elements(sed))
FOR i = 0, n_elements(sed) - 1 DO sedintensity[i] = total(sedcube[100:400,100:400,i])
sedfit = gaussfit(sed[indlambda].lambda,sedintensity[indlambda],nterms = 3)

;mid = (SIZE(sedcube))[1]/2
plot,sed.lambda,sedintensity,/ylog,xrange = [6.5e-7,6.6e-7],xtitle = 'Lambda [m]', ytitle = 'Intensity',xstyle = 1
;plot,sed.lambda,sedcube[250,250,*],/xlog,/ylog,yrange = [1e-5,1e4],xtitle = 'Lambda [m]', ytitle = 'Intensity'
;plot,sed.lambda,sedcube[mid,mid,*],/ylog,xrange = [6.5e-7,6.6e-7],xtitle = 'Lambda [m]', ytitle = 'Intensity',xstyle = 1;,yrange = [1e3,4e3]
;plot,sed.lambda,sedcube[250,250,*],/ylog,xrange = [6.0e-7,7.5e-7],xtitle = 'Lambda [m]', ytitle = 'Intensity',yrange = [1e3,4e3]
oplot,[lambda1,lambda1],[1,1e6],linestyle = 1
oplot,[lambda2,lambda2],[1,1e6],linestyle = 1
oplot,sed[indlambda].lambda,sedintensity[indlambda],color = 100
oplot,sed.lambda,gaussian(sed.lambda,sedfit),color = 100
oplot,[sed[c1].lambda,sed[c2].lambda],[sedintensity[c1],sedintensity[c2]], color = 50, linestyle = 3
;oplot,[sed[c1].lambda,sed[c2].lambda],[sedcube[mid,mid,c1],sedcube[mid,mid,c2]], color = 50, linestyle = 3
;oplot,sed.lambda,sed.lambda*0 + cont[mid,mid],color = 240, linestyle = 2
;oplot,sed[indlambda].lambda,sedcube[mid,mid,indlambda],color = 100
stop

linelambda = sed[indlambda].lambda
linecube = sedcube[*,*,indlambda]
lineimage = FLTARR((size(linecube))[1],(size(linecube))[2])
FOR ix = 0, (size(linecube))[1] - 1 DO $
   FOR iy = 0, (size(linecube))[2] - 1 DO $
      lineimage[ix,iy] = tsum(linelambda, reform(linecube[ix,iy,*] - cont[ix,iy]));, lambda1, lambda2 )

meanimage = FLTARR((size(linecube))[1],(size(linecube))[2])
FOR ix = 0, (size(linecube))[1] - 1 DO $
   FOR iy = 0, (size(linecube))[2] - 1 DO $
      IF (lineimage[ix,iy] GT 1e-10) THEN meanimage[ix,iy] = (gaussfit(sed[indlambda].lambda,reform(sedcube[ix,iy,indlambda]),coeff,nterms = 3))[1] - sedfit[1] ELSE meanimage[ix,iy] = !VALUES.F_NAN

stdevimage = FLTARR((size(linecube))[1],(size(linecube))[2])
FOR ix = 0, (size(linecube))[1] - 1 DO $
   FOR iy = 0, (size(linecube))[2] - 1 DO $
      IF (lineimage[ix,iy] GT 1e-10) THEN stdevimage[ix,iy] = (gaussfit(sed[indlambda].lambda,reform(sedcube[ix,iy,indlambda]),coeff,nterms = 3))[2] ELSE stdevimage[ix,iy] = !VALUES.F_NAN

ind0 = where(lineimage lt 0)
IF NOT ind0[0] eq -1 THEN lineimage[ind0] = 0
print,minmax(meanimage[where(meanimage NE 0)]) 
IF keyword_set(outplot) THEN device,filename = 'intenseimage.eps',/color,ysize=18,xsize=18 ELSE window,0,xsize = 800,ysize = 800
contour,alog10(lineimage),(findgen((size(linecube))[1]) - ((size(linecube))[1])/2)*resx,(findgen((size(linecube))[2]) - ((size(linecube))[1])/2)*resy,nlevels = 245,/cell_fill,min = -11,max = -6.5,/isotropic,xtitle = 'x [kpc]',ytitle = 'y [kpc]'  ;-12
IF keyword_set(outplot) THEN device,/close
mwrfits,lineimage,outfile
stop

range = min(abs([sedfit[1] - lambda1,lambda2 - sedfit[1]]))/sedfit[1]*3e8
range = 4e5
IF keyword_set(outplot) THEN device,filename = 'meanimage.eps',/color,ysize=18,xsize=18 ELSE window,0,xsize = 800,ysize = 800
contour,meanimage/sedfit[1]*3e8,(findgen((size(linecube))[1]) - ((size(linecube))[1])/2)*resx,(findgen((size(linecube))[2]) - ((size(linecube))[1])/2)*resy,nlevels = 245,/cell_fill,min_value = -range,max_value = range,/isotropic,xtitle = 'x [kpc]',ytitle = 'y [kpc]'
colorbar,range = [-range/1000,range/1000],/vertical,tickinterval=1e5,divisions = 8
IF keyword_set(outplot) THEN device,/close

;print,stdevimage[where(finite(stdevimage))]/sedfit[1]*3e8
range = (lambda2 - lambda1)/sedfit[1]*3e8/2
range = 6e5
IF keyword_set(outplot) THEN device,filename = 'stdevimage.eps',/color,ysize=18,xsize=18 ELSE window,0,xsize = 800,ysize = 800
contour,stdevimage/sedfit[1]*3e8,nlevels = 245,/cell_fill,min_value = 0,max_value = range,/isotropic,xtitle = 'x [kpc]',ytitle = 'y [kpc]'
colorbar,range = [0,range/1000],/vertical,tickinterval=1e5,divisions = 6
IF keyword_set(outplot) THEN device,/close


end
