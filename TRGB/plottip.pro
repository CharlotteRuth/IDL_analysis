;This procedure will plot all of the isochrones with their fitted
;tip.  It will also plot the difference between the observed and
;predicted tips.

;Compiled moduales needed:
;fitfunc.pro


;plottip,/fitcurve

PRO plottip,fitcurve=fitcurve
obsdata_file = '/astro/net/scratch2/christensen/TRGBDust/NGC0253-WIDE1/'
iso_file = '/astro/net/scratch2/christensen/TRGBDust/isochrones/metalgrid/'
galaxy = 'NGC0253-WIDE1'
infile=obsdata_file + galaxy + '.vi2.fits'
columns = ['mag1_acs', 'mag1_err','mag2_acs','mag2_err']
data = mrdfits(infile,1,columns=columns)
mI = data.mag2_acs
mColor = data.mag1_acs - data.mag2_acs
addy = 26.8 + 0.88

infile = obsdata_file+'iso_tips.dat'
infile = 'iso_tips.dat'
READCOL,infile,iso1,iso2,color,color_error,magnitude,magnitude_error,FORMAT='A,A,F,F,F,F',/SILENT
elements = 150
points = N_ELEMENTS(color)

IF (KEYWORD_SET(fitcurve)) THEN BEGIN
    IF points - 1 gt 5 THEN BEGIN
        a = [1.,1.,1.,1.,1.,1.] 
        tipcurve_x = findgen(20)*(MAX(color)-MIN(color))*1.1/20.0 + MIN(COLOR)+0.01
        tipfit = curvefit(color,magnitude,w,a,FUNCTION_NAME = 'fitfunc')
        fitfunc,tipcurve_x,a,tipcurve_y
    ENDIF ELSE BEGIN
        tipfit = poly_fit(color,magnitude,4)
        tipcurve_y = tipfit[0] + tipcurve_x*tipfit[1] + tipcurve_x^2*tipfit[2] + tipcurve_x^3*tipfit[3] + tipcurve_x^4*tipfit[4]
    ENDELSE
ENDIF

loadct,39
set_plot,'x'
window,5
set_plot,'ps'
device,filename='iso_plot.eps',/color,bits_per_pixel=8

plot,mColor,mI,xtitle = 'V - I', ytitle = 'I', title = 'Isochrone CMDs',yrange=[25.5,22],xrange=[-0.5,4.0],psym=3
oploterr,color,magnitude,magnitude_error
oplot,color,magnitude,psym=2,color=220

IF (KEYWORD_SET(fitcurve)) THEN oplot,tipcurve_x,tipcurve_y,color=220

isochrones = dblarr(3,elements,N_ELEMENTS(color));[quality, points, isochrones]
;stop
FOR i = 0, points-1 DO BEGIN
    file = extrap_iso(iso_file + iso1[i],elements)
    isochrones[0,*,i] = file.z                      ;metalicity
    isochrones[1,*,i] = file.f606mag - file.f814mag ;color
    isochrones[2,*,i] = file.f814mag + addy                ;imagnitude
    oplot,file.f606mag - file.f814mag,file.f814mag + addy,color =  220 - i*(220./(points+1))
ENDFOR
file = extrap_iso(iso_file + iso2[points - 1],elements)
oplot,file.f606mag - file.f814mag,file.f814mag + addy,color =  220 - points*(220./(points+1))

filelength = 300 ;Find using wc on it (300)or decide what range is desired
isos_name_init = strarr(filelength)
infile_iso = iso_file + 'index.txt'
close,1
openr,1,infile_iso
readf,1,isos_name_init
close,1
isochrones = dblarr(3,filelength);[quality, points, isochrones]

FOR i = 0, filelength - 1 DO BEGIN
    file=iso_file + isos_name_init[i]
    READCOL,file,z,f475mag,f606mag,f814mag,FORMAT='F,X,X,X,X,X,X,X,X,F,X,X,F,X,X,X,X,F',/SILENT
    isochrones[0,i] = z[0]                      ;metalicity
    isochrones[1,i] = f606mag[N_ELEMENTS(f606mag)-1] - f814mag[N_ELEMENTS(f606mag)-1] ;color
    isochrones[2,i] = f814mag[N_ELEMENTS(f814mag)-1] + addy  ;imagnitude
ENDFOR
oplot,isochrones[1,*],isochrones[2,*],color=30
device,/close

IF (KEYWORD_SET(fitcurve)) THEN BEGIN
    print,"Fitting Curve"
    set_plot,'ps'
;    stop
device,filename='iso_plot2.eps',/color,bits_per_pixel=8
    plot,color,magnitude,psym=2,xtitle = 'V - I', ytitle = 'I', title = 'Isochrone CMDs',yrange=[25.5,22],xrange=[0.5,4.0]
    oploterr,color,magnitude,magnitude_error
    oplot,color,magnitude,psym=2,color=220
    oplot,tipcurve_x,tipcurve_y,color=220
    oplot,isochrones[1,*],isochrones[2,*],color=30

    a = [1.,1.,1.,1.,1.,1.] 
    tipfit = curvefit(isochrones[1,*],isochrones[2,*],w,a,FUNCTION_NAME = 'fitfunc')
    fitfunc,color,a,isotippoint_y
    fitfunc,tipcurve_x,a,isotipcurve_y
device,/close
device,filename='iso_diff.eps',/color,bits_per_pixel=8    
;    stop
    plot,color,magnitude-isotippoint_y,psym=2,xtitle = 'Color',ytitle = 'Magnitude',title='Difference Between Predicted and Actual TRGB'
    oploterr,color,magnitude-isotippoint_y,magnitude_error
    oplot,tipcurve_x,tipcurve_y-isotipcurve_y
device,/close
;    stop
ENDIF
END
