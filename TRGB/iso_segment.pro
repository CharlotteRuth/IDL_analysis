;Needed functions:
;iso.pro
;trgb.pro
;extrap_iso.pro
;interpol_color.pro

PRO LUM_FUNC,imag,ierr,imin,imax,binsize,binval,lumfunc
;Finds the smoothed luminosity function.  This will be fit to several
;powerlaws to give a good first guess at the parameters for the broken
;powerlaw function
nbin=FIX((imax-imin)/binsize)+1 ;Number of bins
binval=FINDGEN(nbin)*binsize+imin 
lumfunc=FLTARR(nbin)
edgedet_mendez=FLTARR(nbin)
FOR i=0,nbin-1 DO BEGIN
    lumfunc[i]=TOTAL(1./(SQRT(2.*!PI)*ABS(ierr))*exp(-(imag-binval[i])^2./(2.*ierr^2.)))
ENDFOR
END

PRO prepare_tip,galaxy,vmag,verr,imag,ierr,lefttip,righttip,leftbound,rightbound, iteration, USEVI=usevi,NOPLOT=noplot
;This code is used to prepare the sections of the CMD to be sent to
;trgb.cpp so the tip can be found

;;;;;;;;;Add in Reddening and select stars;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Allows control over the number of stars
;num_stars = [30148,24297,17089,13081,8440,8262,8020,6465,3967,3299,2921]
;num_stars = [23147,14953,10963,7148,6627,11316,5424];Fake_NGC0253-WIDE1
;num_stars = [23852,18855,17652,11538,11174,17439,8537];Fake_NGC0253-WIDE2
num_stars = [7425,6738,7237,5376,6000,10695,6458];Fake_NGC0253-WIDE3
 IF(N_ELEMENTS(vmag) GT num_stars[iteration]) THEN BEGIN;
    vmag = vmag[0:num_stars[iteration]-1]
    verr = verr[0:num_stars[iteration]-1]
    imag = imag[0:num_stars[iteration]-1]
    ierr = ierr[0:num_stars[iteration]-1]
     print,'Number of Stars in ',iteration,': ',N_ELEMENTS(vmag)
ENDIF ELSE print,'Fewer stars in current isochrones: ',N_ELEMENTS(vmag),', ',num_stars[iteration]

base = '/astro/net/scratch2/christensen/TRGBDust/'
outfile_stars = base + 'Datafiles/trgbstars'+galaxy+'_' + strtrim(string(iteration),2) +'.dat' 
outfile_mags = base + 'Datafiles/magANDerr'+galaxy+'_' + strtrim(string(iteration),2) +'.dat' 
vminusi=vmag-imag
nind = (size(vmag))[1]
;Set of points within 0.8 from the theoretical tip
small_set = WHERE(SQRT(((imag-(lefttip[1]+righttip[1])/2.0))^2 + ((vmag - imag)-(lefttip[0]+righttip[0])/2.0)^2) le 0.8) 
;oplot,vmag[small_set]-imag[small_set],imag[small_set],psym = 3, color = 200
;stop
;;;;;;;;;;;;;;;;;;Rotates Reference Frame;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;I want to rotate the frame of reference such that the average slope
;of the line through the two tips of the bounding isochrones is one of the set of two orthonganal basis
;Rotates all the points and the errors
points = trans_rotate(lefttip[0],righttip[0],leftbound[0,*],leftbound[1,*],rightbound[0,*],rightbound[1,*],[[vmag-imag],[imag]])
lowerlim_b = MIN(points[small_set,1])
upperlim_b = MAX(points[small_set,1])
errors = ABS(trans_rotate(lefttip[0],righttip[0],leftbound[0,*],leftbound[1,*],rightbound[0,*],rightbound[1,*],[[(verr+ierr)/2.0],[ierr]],/error))

;Rotates the tips
tips_rot = trans_rotate(lefttip[0],righttip[0],leftbound[0,*],leftbound[1,*],rightbound[0,*],rightbound[1,*],Transpose([[lefttip],[righttip]]))
lefttip_rot = tips_rot[0,*]
righttip_rot = tips_rot[1,*]
upperlim = [MAX([lefttip_rot[0],righttip_rot[0]]),MIN([lefttip_rot[1],righttip_rot[1]])]
lowerlim = [MIN([lefttip_rot[0],righttip_rot[0]]),MAX([lefttip_rot[1],righttip_rot[1]])]

mtrgb = (upperlim[1]+lowerlim[1])/2
print,base+galaxy+'/rotate_' + strtrim(iteration,2)+'.eps'
set_plot,'ps'
device,filename = base+galaxy+'/'+galaxy+'rotate_' + strtrim(iteration,2)+'.eps';,/color,bits_per_pixel = 8
plot,points[*,0],points[*,1],xtitle = 'F606W - F814W [VEGAmag] (Rotated)', ytitle = 'F814W [VEGAmag] (Rotated)',yrange=[25.5,22],xrange=[-0.5,4.0],psym = 3,charsize = 1.5
oplot,[-1,5],[mtrgb,mtrgb]
device,/close
;stop
set_plot,'x'
;;;;;;;;;;;;;;;;;;;;Find and Initial Fit to the Luminosity Function;;;;;;;;;;;;;;;;;;;;;;
binsize=0.01                    ;Width of bin in magnitude
lum_func,points[*,1],errors[*,1],MIN(points[*,1]),MAX(points[*,1]),binsize,binval,lumfunc
window,1
log_lum = ALOG10(lumfunc)
plot,binval,log_lum,yrange = [-2,4],xrange = [26,22]
mtrgb = (upperlim[1]+lowerlim[1])/2
lowerlim[1] = mtrgb - (mtrgb - lowerlim_b)*0.50
upperlim[1] = mtrgb - (mtrgb - upperlim_b)*0.50
uppertrgb = mtrgb - (mtrgb - upperlim_b)*0.1
rgb_i  = where(binval ge uppertrgb AND binval lt upperlim_b)
agb_i  = where(binval le uppertrgb AND binval gt lowerlim_b)
rgb_fit = poly_fit(binval[rgb_i],log_lum[rgb_i],1)
agb_fit = poly_fit(binval[agb_i],log_lum[agb_i],1)
oplot,binval[where(binval lt upperlim_b and binval gt lowerlim_b)],log_lum[where(binval lt upperlim_b and binval gt lowerlim_b)],color = 240
oplot,binval[rgb_i],log_lum[rgb_i],color = 220
oplot,binval[agb_i],log_lum[agb_i],color = 200
IF(finite(agb_fit[0]) eq 0) THEN agb_fit[0] = 1.5
IF(finite(agb_fit[1]) eq 0) THEN agb_fit[1] = 0
IF(finite(rgb_fit[0]) eq 0) THEN last_fit[0] = 0.7
IF(finite(rgb_fit[1]) eq 0) THEN last_fit[1] = 0
oplot,binval,agb_fit[0] + agb_fit[1]*binval,linestyle = 2, color = 50
oplot,binval,rgb_fit[0] + rgb_fit[1]*binval,linestyle = 2, color = 50
oplot,[upperlim[1],upperlim[1]],[-2,5],linestyle = 1, color = 120
oplot,[lowerlim[1],lowerlim[1]],[-2,5],linestyle = 1, color = 120

print,"RGB Fit:    ",rgb_fit[1],rgb_fit[0] + rgb_fit[1]*mtrgb
print,"AGB Fit:    ",agb_fit[1],agb_fit[0] + agb_fit[1]*mtrgb
print,'Tip Limits: ',lowerlim[1],upperlim[1]
print,'Mag Limits: ',lowerlim_b,upperlim_b
print,'Number of Stars: ',N_ELEMENTS(where(points[*,1] gt lowerlim_b AND points[*,1] lt upperlim_b))

;*****************************Switch to C*****************************************
;My tip determination cannot be done quickly in idl so I am doing it in c
;This program will output a file of stars and stop
;When this happens, run trgb.c (make sure the output files are appropriate)
;Fitting the error to a function of magnitude
sigmafit = poly_fit(points[*,1],errors[*,1],4)
merr = sigmafit[0] + sigmafit[1]*binval + sigmafit[2]*binval^2.0+sigmafit[3]*binval^3 + sigmafit[4]*binval^4 
print,outfile_stars
OPENW,3,outfile_stars,WIDTH=1000
printf,3,TRANSPOSE([[points[*,1]],[errors[*,1]]]),FORMAT='(F,F)'
close,3

lowertemp = MIN(ABS(binval - lowerlim[1]),lowerpos)
uppertemp = MIN(ABS(binval - upperlim[1]),upperpos)
lowertemp = MIN(ABS(binval - lowerlim_b),lowerpos_b)
uppertemp = MIN(ABS(binval - upperlim_b),upperpos_b)
print,outfile_mags
OPENW,3,outfile_mags,WIDTH=1000
;width=0.05
printf,3,rgb_fit[1],rgb_fit[0] + rgb_fit[1]*mtrgb,FORMAT='(F,F,F)'
;parameters of fit: a, d
printf,3,lowerpos,upperpos,FORMAT='(I,I,I)'
;range of possible mrtgb magnitudes
printf,3,lowerpos_b,upperpos_b,N_ELEMENTS(where(points[*,1] gt lowerlim_b AND points[*,1] lt upperlim_b)),FORMAT='(I,I,I)'
;range the function should be fit over
printf,3,TRANSPOSE([[binval],[merr]]),FORMAT='(F,F)'
close,3
;stop
END


;iso_segment,'M81-DEEP',adjust = 27.79
;iso_segment,'NGC0253-WIDE1',adjust = 27.60
;iso_segment,'Fake3',/Fakefile,adjust = 27.62
pro iso_segment,galaxy,Adjust = adjust, FAKEFILE=fakefile, USEVI=usevi
;Charlotte Christensen, 7/17/07
;The program reads in isochrones of different metallicities
;It is used to divide a CMD into different sections to be used for
;calculating the trgb
close,/all
;loadct,39
set_plot,'x'
;set_plot,'ps'
base = '/astro/net/scratch2/christensen/TRGBDust/'
iso_file = base + 'Isochrones/metalgrid/'
;strmetals = ['0.00010','0.00110','0.00210','0.00310','0.00410','0.00510','0.00610','0.00710','0.00810','0.00910','0.01010','0.01110']
strmetals = ['0.00110','0.00210','0.00310','0.00410','0.00510','0.00610','0.00810','0.01010']
;strmetals = ['0.00110','0.00310','0.00510','0.00710','0.00910','0.01110']
;strmetals = ['0.00210','0.00810']
metals = FLOAT(strmetals)
IF(KEYWORD_SET(adjust)) then addy = adjust else addy = 0
ebv = 0.019;  Dunno about this.  Let's wait
ebv = 0.08
;;;;;;;;;;;;;;;;;;;;;;;;READ IN STELLAR DATA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;psfphotfile=galaxy+"/"+galaxy+".vi"
;read in data - use vega2vi to transform to Johnson if necessary
;IF (KEYWORD_SET(usevi)) THEN BEGIN
    ;The next like needs to be altered
;    RDFLOAT, psfphotfile,x,dx,y,dy,f606,verr,f814,ierr 
;    vega2vi, f606,f814,vmag,imag
;ENDIF ELSE RDFLOAT, psfphotfile,x,y,vmag,Transmag_1,verr,chi_psf_1, sharp_1, crowd_1, err_flag_1, imag,Transmag_2,ierr,chi_psf_2, sharp_2, crowd_2, err_flag_2

;now determine reddening 
;READCOL, 'procedures/reddenings.dat',rgalaxy,ebv,FORMAT='A,F',/SILENT
;match=WHERE(rgalaxy EQ galaxy)
;ebv=ebv[match[0]]
;commented out about above because galaxy names don't match and I'm
;not sure what to use, set the reddening below, from data in the file

IF (KEYWORD_SET(usevi)) THEN BEGIN
    vred=ebv*3.1
    ired=ebv*3.1*0.482 ;B&M p.136-7
ENDIF ELSE BEGIN
    vred=ebv*2.716     
    ired=ebv*1.796 ;from Roelof 16 June 2004 email, M0 star
ENDELSE

IF(KEYWORD_SET(fakefile)) THEN BEGIN  
    readcol,base + galaxy + '/' + galaxy + '.vi2',fakev,v_err,fakei,i_err
    fakev = fakev - vred
    fakei = fakei - ired
    mI = fakei
    mColor = fakev - fakei
ENDIF ELSE BEGIN
    infile=base + galaxy + '/' + galaxy + '.vi2.fits'
    columns = ['mag1_acs', 'mag1_err','mag2_acs','mag2_err']
    data = mrdfits(infile,1,columns=columns)
    data.mag1_acs=data.mag1_acs -vred
    data.mag2_acs=data.mag2_acs-ired    
    mI = data.mag2_acs
    mColor = data.mag1_acs - data.mag2_acs
ENDELSE


;;;;;;;;;;;;;;;;;;;;;;READ IN ISOCHRONES;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filelength = 120 ;Find using wc on it (300)or decide what range is desired
num_isos = N_ELEMENTS(metals)
openw,1,base+"Datafiles/metals"+galaxy+".txt"
printf,1,addy
printf,1,TRANSPOSE(strmetals)
close,1
isos_names = "isoc_"+STRING(strmetals)+".dat"
elements = 150
isochrones = dblarr(3,elements,num_isos);[quality, points, isochrones]
FOR i = 0, (num_isos - 1) DO BEGIN
    file = extrap_iso(iso_file + isos_names[i],elements)
    isochrones[0,*,i] = file.z                      ;metalicity
    isochrones[1,*,i] = file.f606mag - file.f814mag ;color
    isochrones[2,*,i] = file.f814mag + addy                ;imagnitude
ENDFOR

READCOL,iso_file+"Tips.dat",z,f475mag,f606mag,f814mag,FORMAT='F,X,X,X,X,X,X,X,X,F,X,X,F,X,X,X,X,F',/SILENT
tips = fltarr(3,filelength)
tips[0,*] = z[0:filelength-1]
tips[1,*] = f606mag[0:filelength-1] - f814mag[0:filelength-1]
tips[2,*] = f814mag[0:filelength-1] + addy

;;;;;;;;;;;;;;; FIND TIPS OF EACH SEGMENT ;;;;;;;;;;;;;;;;;;;
length = N_ELEMENTS(mI)
outfile_lumlist = base + "Datafiles/lum_func_list"+galaxy+".txt"
outfile_tipslist = base + "Datafiles/found_tips_list"+galaxy+".txt"
openw,7,outfile_lumlist
openw,8,outfile_tipslist
trgb_result = fltarr(num_isos-1,6) ; Structures to hold results

set_plot,'ps'
device,filename = base+galaxy+'/'+galaxy+'isoCMD.eps' ;,/color,bits_per_pixel = 8
plot,mColor,mI,xtitle = 'F606W - F814W [VEGAmag]', ytitle = 'F814W [VEGAmag]',yrange=[25.5,22],xrange=[-0.5,4.0],psym = 3,charsize = 1.5
FOR i = 0, (num_isos - 1) DO BEGIN
    oplot,isochrones[1,*,i],isochrones[2,*,i] ;,color = 100;220 - i*(220./num_isos)
ENDFOR
device,/close

FOR ct = 0, num_isos-2 DO BEGIN
    print,"Metallicity: ",isos_names[ct]
    length = N_ELEMENTS(mI)
    between = replicate(0,length)
    FOR i = 0L, length - 1 DO BEGIN
            IF ((below(isochrones[1,*,ct],isochrones[2,*,ct],mColor[i],mI[i]) EQ 1) AND (above(isochrones[1,*,ct+1],isochrones[2,*,ct+1],mColor[i],mI[i]) EQ 1)) THEN between[i] = 1 ;Determins whether a star is between two isochrones
    ENDFOR
    between = Where(between EQ 1)

;    window,0
    set_plot,'ps'
    device,filename = base+galaxy+'/'+galaxy+'selectiso_' + strtrim(ct,2)+'.eps';,/color,bits_per_pixel = 8
    plot,mColor,mI,xtitle = 'F606W - F814W [VEGAmag]', ytitle = 'F814W [VEGAmag]',yrange=[25.5,22],xrange=[-0.5,4.0],psym = 3,charsize = 1.5
    oplot,mColor,mI,psym = 3,color = 100
    FOR i = 0, (num_isos - 1) DO BEGIN
       oplot,isochrones[1,*,i],isochrones[2,*,i];,color = 100;220 - i*(220./num_isos)
    ENDFOR
    oplot,mColor[between], mI[between],psym = 3;,color=230
device,/close
;stop
set_plot,'x'
;    stop
   IF (KEYWORD_SET(fakefile)) THEN prepare_tip,galaxy,fakev[between],v_err[between],fakei[between],i_err[between],tips[1:2,where(tips[0,*] eq metals[ct])],tips[1:2,where(tips[0,*] eq metals[ct+1])],isochrones[1:2,*,ct],isochrones[1:2,*,ct+1],ct,/usevi ELSE prepare_tip,galaxy,(data.mag1_acs)[between],(data.mag1_err)[between],(data.mag2_acs)[between],(data.mag2_err)[between],tips[1:2,where(tips[0,*] eq metals[ct])],tips[1:2,where(tips[0,*] eq metals[ct+1])],isochrones[1:2,*,ct],isochrones[1:2,*,ct+1],ct,/usevi
   printf,7,"../Datafiles/lum_func"+galaxy+"_"+STRTRIM(STRING(ct),2)+".dat"
   printf,8,"../Datafiles/found_tips"+galaxy+"_"+STRTRIM(STRING(ct),2)+".dat"
ENDFOR
close,/all
print,"Now run trgb"
END
