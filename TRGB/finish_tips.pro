;Needed modulues
;gbasis.pro
;iso.pro
;fitfunc.pro

FUNCTION BROKEN_POWERLAW, m, a, b, c, d, mtrgb
;This function gives the functional form of the tip of the RGB
  g = m
  g[WHERE(m GT mtrgb)] = 10.0^(a*(m[WHERE(m GT mtrgb)] - mtrgb)+ d)
  g[WHERE(m LE mtrgb)] = 10.0^(b*(m[WHERE(m LE mtrgb)] - mtrgb)+ c) 
  RETURN,g
END


;This procedure will plot all of the isochrones with their fitted
;tip.  It will also plot the difference between the observed and
;predicted tips.
PRO plottip,galaxy,fitcurve=fitcurve,FAKEFILE=fakefile,USEVI = USEVI, OUTPLOT = OUTPLOT,adjust = adjust
;loadct,39
home = '/astro/net/scratch2/christensen/TRGBDust/'
iso_file = '/astro/net/scratch2/christensen/TRGBDust/Isochrones/metalgrid/'
ebv = 0.019;  Dunno about this.  Let's wait
ebv = 0.08
IF (keyword_set(adjust)) THEN adjust = adjust ELSE adjust = 0
IF (KEYWORD_SET(usevi)) THEN BEGIN 
    vred=ebv*3.1
    ired=ebv*3.1*0.482 ;B&M p.136-7
    colortitle = 'V - I [Johnson]'
    ititle = 'I [Johnson]'
    difftitle = 'I_Observed - I_Predicted [Johnson]'
ENDIF ELSE BEGIN
    vred=ebv*2.716     
    ired=ebv*1.796 ;from Roelof 16 June 2004 email, M0 star
    colortitle = 'F606W - F814W [VEGAmag]'
    ititle = 'F814W [VEGAmag]'
    difftitle = textoidl('F814W_{Obs}')+' - '+ textoidl('F814W_{Theory} [VEGAmag]')
ENDELSE

IF(KEYWORD_SET(fakefile)) THEN BEGIN  
    infile = home+galaxy+'/' + galaxy + '.vi2'
    readcol,infile,fakev,v_err,fakei,i_err,/silent
    fakev = fakev - vred
    fakei = fakei - ired
    mI = fakei
    mI_err = i_err
    mV = fakev   
    mV_err = v_err
    mColor = fakev - fakei
ENDIF ELSE BEGIN
    infile= home+galaxy+'/'+ galaxy + '.vi2.fits'
    columns = ['mag1_acs', 'mag1_err','mag2_acs','mag2_err']
    data = mrdfits(infile,1,columns=columns)
    data.mag1_acs=data.mag1_acs - vred
    data.mag2_acs=data.mag2_acs - ired
    mI = data.mag2_acs          ;i
    mI_err = data.mag2_err      ;i
    mV = data.mag1_acs          ;v
    mV_err = data.mag1_err      ;i
    mColor = data.mag1_acs - data.mag2_acs ;v - i  
ENDELSE

infile = home +galaxy+'/iso_tips'+galaxy+'.dat'
READCOL,infile,iso1,iso2,color,color_error,magnitude,magnitude_error,model_mag,addy,diff, FORMAT='A,A,F,F,F,F,F,F,F',/SILENT
addy = addy[0]
elements = 150
points = N_ELEMENTS(color)


;*********Fit a curve to tip**************************************************
IF (KEYWORD_SET(fitcurve)) THEN BEGIN
    tipcurve_x = findgen(20)*(MAX(color)-MIN(color))*1.1/20.0 + MIN(COLOR)+0.01
    IF points - 1 gt 5 THEN BEGIN
        a = [1.,1.,1.,1.,1.,1.] 
        tipfit = curvefit(color,magnitude,1/magnitude_error^2,a,sigma,FUNCTION_NAME = 'fitfunc')
        fitfunc,tipcurve_x,a,tipcurve_y
        print,'Function fitted: ',a,sigma
    ENDIF ELSE BEGIN
        IF (points - 1 gt 0) THEN BEGIN
            tipfit = poly_fit(color,magnitude,4)
            tipcurve_y = tipfit[0] + tipcurve_x*tipfit[1] + tipcurve_x^2*tipfit[2] + tipcurve_x^3*tipfit[3] + tipcurve_x^4*tipfit[4]
        ENDIF ELSE BEGIN
            tipcurve_x = [-2,tipcurve_x,5]
            tipcurve_y = tipcurve_x*0 + magnitude[0] 
        ENDELSE
    ENDELSE
ENDIF
IF (KEYWORD_SET(OUTPLOT)) THEN BEGIN 
    set_plot,'ps'
    device,filename=home + galaxy+'/iso_plot'+galaxy+'.eps',/color,bits_per_pixel=8
ENDIF ELSE BEGIN
    set_plot,'x'
    window,4
ENDELSE

;*********Plot CMD and the tip****************************************************
plot,mColor,mI,xtitle = colortitle, ytitle = ititle,yrange=[25.5,22],xrange=[-0.5,4.0],psym=3,charsize = 1.5; title = 'Isochrone CMDs, '+galaxy,
oplot,mColor,mI,psym=3,color = 100
oploterror,color,magnitude,color_error,magnitude_error,psym = 3
oplot,color,magnitude,psym=2;,color=220
;IF (KEYWORD_SET(fitcurve)) THEN oplot,tipcurve_x,tipcurve_y,color=220
;*********Read in the isochrones used to divide the data***************************
isochrones = dblarr(3,elements,N_ELEMENTS(color));[quality, points, isochrones]
FOR i = 0, points-1 DO BEGIN
    file = extrap_iso(iso_file + iso1[i],elements)
    isochrones[0,*,i] = file.z                      ;metalicity
    isochrones[1,*,i] = file.f606mag - file.f814mag ;color
    isochrones[2,*,i] = file.f814mag + addy                ;imagnitude
    oplot,file.f606mag - file.f814mag,file.f814mag + addy;,color =  220 - i*(220./(points+1))
ENDFOR
file = extrap_iso(iso_file + iso2[points - 1],elements)
oplot,file.f606mag - file.f814mag,file.f814mag + addy;,color =  220 - points*(220./(points+1))
;Read in all the isochronses to find their tip of the red giant branch
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
oplot,isochrones[1,*],isochrones[2,*];,color=30
IF (KEYWORD_SET(OUTPLOT)) THEN device,/close

;*************Plot difference between observed and predicted tip**********************************
print,"Fitting Curve"
IF (KEYWORD_SET(OUTPLOT)) THEN BEGIN 
    set_plot,'ps'
    device,filename=home + galaxy+'/iso_diff'+galaxy+'.eps',/color,bits_per_pixel=8
ENDIF ELSE BEGIN
    set_plot,'x'
    window,5
ENDELSE
a_pred = [1.,1.,1.,1.,1.,1.] ;Fit a curve to the set of modeled isochrone tips
tipfit = curvefit(isochrones[1,*],isochrones[2,*],w,a_pred,FUNCTION_NAME = 'fitfunc')
fitfunc,color,a_pred,isotippoint_y
weighted_mean = Total((magnitude-isotippoint_y)/magnitude_error^2)/Total(1/magnitude_error^2)
weighted_sigma = SQRT(1/TOTAL(1/magnitude_error^2))

plot,color,magnitude-isotippoint_y - weighted_mean,psym=2,xtitle = colortitle,ytitle = difftitle,yrange = [MIN(magnitude-isotippoint_y-weighted_mean)-MAX(magnitude_error),MAX(magnitude-isotippoint_y-weighted_mean)+MAX(magnitude_error)],charsize = 1.5;,title='Difference Between Predicted and Actual TRGB, '+galaxy,
oploterror,color,magnitude-isotippoint_y-weighted_mean,color_error,magnitude_error,psym = 3

IF (KEYWORD_SET(fitcurve)) THEN BEGIN
    fitfunc,tipcurve_x,a_pred,isotipcurve_y
;    oplot,tipcurve_x,tipcurve_y-isotipcurve_y
ENDIF

oplot,[0,4],[0,0],linestyle=2

meanstring = strmid(strtrim(round(1e4*(weighted_mean+adjust))*1e-4,2),0,6)
errorstring = strmid(strtrim(round(1e4*weighted_sigma)*1e-4,2),0,6)
xyouts,1.1,-0.1,meanstring+' '+textoidl('\pm')+' '+errorstring,charsize=2
IF (KEYWORD_SET(OUTPLOT)) THEN device,/close
END

;*********************************************************************************************************

;This program will read in a list of tip files (generated by trgb.cpp)
;and finish finding the magnitude and plotting
PRO finish_tips,galaxy,tipfiles_name=tipfiles_name,lumfiles_name = lumfiles_name,FAKEFILE=fakefile, USEVI = USEVI, PLOTTIP = PLOTTIP, OUTPLOT = OUTPLOT
close,/all
;loadct,39
IF KEYWORD_SET(OUTPLOT) THEN OUTPLOT = 1 ELSE OUTPLOT = 0
IF outplot THEN set_plot,'ps' ELSE set_plot,'x'

IF NOT KEYWORD_SET(tipfiles_name) THEN tipfiles_name = "../Datafiles/found_tips_list"+galaxy+".txt"
IF NOT KEYWORD_SET(lumfiles_name) THEN lumfiles_name = "../Datafiles/lum_func_list"+galaxy+".txt"
readcol,tipfiles_name,format = 'A60',tipfiles,/silent
readcol,lumfiles_name,format = 'A60',lumfiles,/silent
base = '/astro/net/scratch2/christensen/TRGBDust/'+galaxy+'/'
iso_file = '/astro/net/scratch2/christensen/TRGBDust/Isochrones/metalgrid/'

;;;;;;;;;;;;;;;;;;;;;;;;READ IN STELLAR DATA;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ebv = 0.019;  Dunno about this.  Let's wait
ebv = 0.08
IF (KEYWORD_SET(usevi)) THEN BEGIN 
    vred=ebv*3.1
    ired=ebv*3.1*0.482 ;B&M p.136-7
    outfile=base+galaxy+"_trgb_vi_reg.dat" 
    OPENW,2,outfile,WIDTH=1000
    colortitle = "V - I [Johnson]"
    ititle = "I [Johnson]"
ENDIF ELSE BEGIN
    vred=ebv*2.716     
    ired=ebv*1.796 ;from Roelof 16 June 2004 email, M0 star
    outfile=base+galaxy+"_trgb_reg.dat"
    OPENW,2,outfile,WIDTH=1000
    colortitle="F606W - F814W [VEGAmag]"        
    ititle="F814W [VEGAmag]"
ENDELSE

IF(KEYWORD_SET(fakefile)) THEN BEGIN  
    infile = base + galaxy + '.vi2'
    readcol,infile,fakev,v_err,fakei,i_err,/silent
    fakev = fakev - vred
    fakei = fakei - ired
    mI = fakei
    mI_err = i_err
    mV = fakev   
    mV_err = v_err
    mColor = fakev - fakei
ENDIF ELSE BEGIN
    infile=base + galaxy + '.vi2.fits'
    columns = ['mag1_acs', 'mag1_err','mag2_acs','mag2_err']
    data = mrdfits(infile,1,columns=columns)
    data.mag1_acs=data.mag1_acs - vred
    data.mag2_acs=data.mag2_acs - ired
    mI = data.mag2_acs          ;i
    mI_err = data.mag2_err      ;i
    mV = data.mag1_acs          ;v
    mV_err = data.mag1_err      ;i
    mColor = data.mag1_acs - data.mag2_acs ;v - i
ENDELSE

;;;;;;;;;;;;;;;;;;;;;;READ IN ISOCHRONES;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filelength = 120 ;Find using wc on it (300)or decide what range is desired
readcol,'../Datafiles/metals' + galaxy + '.txt',addy,FORMAT='F',/silent
readcol,'../Datafiles/metals' + galaxy + '.txt',strmetals,FORMAT='A',/silent
addy = addy[0]
strmetals = strmetals[1:N_ELEMENTS(strmetals)-1]
metals = FLOAT(strmetals)
num_isos = N_ELEMENTS(metals)
trgb_result = fltarr(num_isos-1,6) ; Structures to hold results
isos_name_init = strarr(filelength)
infile_iso = iso_file + 'index.txt'
openr,1,infile_iso
readf,1,isos_name_init
close,1
openw,1,'../'+galaxy+'/iso_tips'+galaxy+'.dat'
printf,1,'First Isochrone ','Second Isochrone ','TRGB Color','Error in Color ','Obs Magnitude ','Error in Magnitude',' Model Magnitude','Added Mag','Net Difference',FORMAT='(9A20)'

isos_names = "isoc_"+STRTRIM(STRING(strmetals),2)+".dat"
elements = 150
isochrones = dblarr(3,elements,num_isos);[quality, points, isochrones]
FOR i = 0, (num_isos - 1) DO BEGIN
    file = extrap_iso(iso_file + isos_names[i],elements)
    isochrones[0,*,i] = file.z                      ;metalicity
    isochrones[1,*,i] = file.f606mag - file.f814mag ;color
    isochrones[2,*,i] = file.f814mag + addy                ;imagnitude
ENDFOR

READCOL,iso_file+"Tips.dat",z,f475mag,f606mag,f814mag,FORMAT='F,X,X,X,X,X,X,X,X,F,X,X,F,X,X,X,X,F',/SILENT
tips = dblarr(3,filelength)
tips[0,*] = z[0:filelength-1]
tips[1,*] = f606mag[0:filelength-1] - f814mag[0:filelength-1]
tips[2,*] = f814mag[0:filelength-1]+addy

;;;;;;;;;;;;;;; Finding the tip for each of the segments;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FOR ct = 0, num_isos-2 DO BEGIN
;Looping through the isochrones
    limit_file = base + '../Datafiles/magANDerr'+galaxy+'_' + strtrim(string(ct),2) +'.dat'
    openr,3,limit_file
    readf,3, trash1,trash2,trash3,trash4,lowerlim,upperlim,trash5
    close,3
    print,"Metallicity: ",isos_names[ct]
    length = N_ELEMENTS(mI)
    between = replicate(0,length)
    FOR i = 0L, length - 1 DO BEGIN
        IF ((below(isochrones[1,*,ct],isochrones[2,*,ct],mColor[i],mI[i]) EQ 1) AND (above(isochrones[1,*,ct+1],isochrones[2,*,ct+1],mColor[i],mI[i]) EQ 1)) THEN between[i] = 1 
;Determins whether a star is between two isochrones
    ENDFOR
    between = Where(between EQ 1)
    vmag = mV[between]
    verr = mV_err[between]
    imag = mI[between]
    ierr = mI_err[between]
    lefttip = tips[1:2,where(tips[0,*] eq metals[ct])]
    righttip = tips[1:2,where(tips[0,*] eq metals[ct+1])]
    leftbound = isochrones[1:2,*,ct]
    rightbound = isochrones[1:2,*,ct+1]
    vminusi=vmag-imag
    nind = (size(vmag))[1]
    print,"__________________________________________"
    print,"Run Number: ",ct
    print, "Using ", STRTRIM(nind,2), " Stars"

;I want to rotate the frame of reference such that the average slope
;of the bounding isochrones is one of the set of two orthonganal basis
    set = Where(imag gt 23.0 AND imag lt 25.0)
    leftbound_first = leftbound
    rightbound_first = rightbound
    points = trans_rotate(lefttip[0],righttip[0],leftbound_first[0,*],leftbound_first[1,*],rightbound_first[0,*],rightbound_first[1,*],[[vmag[set]-imag[set]],[imag[set]]])
 
    errors = ABS(trans_rotate(lefttip[0],righttip[0],leftbound_first[0,*],leftbound_first[1,*],rightbound_first[0,*],rightbound_first[1,*],[[(verr[set]+ierr[set])/2.0],[ierr[set]]]))
    rightbound = trans_rotate(lefttip[0],righttip[0],leftbound_first[0,*],leftbound_first[1,*],rightbound_first[0,*],rightbound_first[1,*],TRANSPOSE(rightbound))
    leftbound = trans_rotate(lefttip[0],righttip[0],leftbound_first[0,*],leftbound_first[1,*],rightbound_first[0,*],rightbound_first[1,*],TRANSPOSE(leftbound))
;    window,5
;    plot,points[*,0],points[*,1],psym = 3,yrange=[max(points[*,1]),Min(points[*,1])]
;    plot,vmag[set]-imag[set],imag[set],psym = 3,yrange = [MAX(imag[set]),MIN(imag[set])]
;    oplot,points[*,0],points[*,1],psym = 3,color = 100
;    points2 = trans_rotate(lefttip[0],righttip[0],leftbound_first[0,*],leftbound_first[1,*],rightbound_first[0,*],rightbound_first[1,*],points,/inverse)
;    oplot,points2[*,0],points2[*,1],psym = 3,color = 50

    readcol,tipfiles[ct],mendezerr,FORMAT='F',/silent
    readcol,lumfiles[ct],binval,lumfunc,FORMAT='F,F',/silent
    a_param = binval[N_ELEMENTS(binval) - 2]
    c = binval[N_ELEMENTS(binval) - 1]  
    binval = binval[0:N_ELEMENTS(binval) - 3]
    b = lumfunc[N_ELEMENTS(lumfunc) - 2]
    d = lumfunc[N_ELEMENTS(lumfunc) - 1]    
    lumfunc = lumfunc[0:N_ELEMENTS(lumfunc) - 3]
    mendez = mendezerr[0]
    mendezerr = mendezerr[1:N_ELEMENTS(mendezerr)-1]

    IF (KEYWORD_SET(OUTPLOT)) THEN BEGIN 
        set_plot,'ps'
        device, filename = base + galaxy+'error_'+strtrim(ct,2)+'.eps'
    ENDIF ELSE BEGIN
        set_plot,'x'
        window,0
    ENDELSE
    !P.MULTI=[0,1,2]

;now fit monte carlo results to gaussian
    nruns = N_ELEMENTS(mendezerr)
    estimates=[FLOAT(nruns)/10.,mendez,0.1]
    plothist,mendezerr,xhist,yhist,bin=0.01,/xsty,ytitle="Number",xtitle=ititle,yrange = [0,N_ELEMENTS(mendezerr)+N_ELEMENTS(mendezerr)/10.],ystyle = 1,charsize=1.5,xrange = [MIN(mendezerr) - (MAX(mendezerr) - MIN(mendezerr)),MAX(mendezerr) + (MAX(mendezerr) - MIN(mendezerr))];,xrange = [binval[lowerlim],binval[upperlim]] ,charsize=1.5;,xrange=[imin,imax],title = "Error for Run "+strtrim(ct,2)+', '+galaxy,
    IF (MIN(mendezerr) ne MAX(mendezerr)) THEN BEGIN
        IF (N_Elements(xhist) gt 3) THEN BEGIN 
            yhist_cut = yhist
            result=gaussfit(xhist,yhist_cut,a,nterms=3) 
            mcwidth=a[2]
            oplot,xhist,result
        ENDIF ELSE mcwidth = stddev(mendezerr)
        plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2
        IF (mcwidth gt (MAX(xhist) - MIN(xhist))/2) THEN  mcwidth = stddev(mendezerr)
        oplot,[mendez - mcwidth, mendez + mcwidth],[MAX(yhist),MAX(yhist)],thick=2,linestyle = 2
     ENDIF ELSE BEGIN
        bin = 0.01
        xmin = MIN([mendezerr,mendez]) - 5.0*bin
        xmax = MAX([mendezerr,mendez]) + 5.0*bin
        plot,[mendezerr[0]-bin/2.,mendezerr[0]-bin/2.,mendezerr[0]+bin/2.,mendezerr[0]+bin/2.],[0,N_ELEMENTS(mendezerr),N_ELEMENTS(mendezerr),0],ytitle="Number",xtitle=ititle,yrange = [0,N_ELEMENTS(mendezerr)+N_ELEMENTS(mendezerr)/10.],ystyle = 1,charsize=1.5;,xrange=[xmin,xmax],xrange=[imin,imax],title = "Error for Run "+strtrim(ct,2)+', '+galaxy,
        plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2 
        mcwidth = bin/2.
        oplot,[mendez - mcwidth, mendez + mcwidth],[10,10],thick=2,linestyle = 2        
    ENDELSE

;see Seth et al. 2005 for derivation of M_F814W=-4.08 as TRGB magnitude
;;;;;;;;;;;;;;;;now determine colors at tip;;;;;;;;;;;;;;;;;;;;
    tipstars=WHERE((points[*,1] GE mendez) AND (points[*,1] LE (mendez+0.15)),ntipstars) ;Selects out stars right below the tip
    yhist = histogram(points[tipstars,0],binsize = 0.02)
    ymax = MAX([1.1*MAX(yhist)])
;    stop
    plothist,points[tipstars,0],xhist,yhist,bin=0.02, colortitle = 'V-I', ytitle = 'Number',charsize=1.5,yrange = [0,ymax];,title = 'Color of TRGB',
 ;   stop
    colortrgb = MEAN(points[tipstars,0])
    widthtrgb = STDDEV(points[tipstars,0])
    IF (colortrgb lt MIN(points[tipstars,0]) OR colortrgb GT MAX(points[tipstars,0])) THEN BEGIN
        a = moment(points[tipstars,0])
        colortrgb=a[0]
        widthtrgb=(MAX(points[tipstars,0]) - MIN(points[tipstars,0]))/2.0
    ENDIF
    colortrgb = MEAN(points[tipstars,0])
 ;   plots,[colortrgb,colortrgb],!Y.CRANGE,thick=2,linestyle=2 
    oplot,[colortrgb,colortrgb],[0,100],thick=2,linestyle=2 
    oplot,[colortrgb-widthtrgb, colortrgb+widthtrgb],[ymax/2,ymax/2],thick=2,linestyle=2
    IF (KEYWORD_SET(OUTPLOT)) THEN device,/close
    print, galaxy, " TRGB Magnitude (wrong basis)  ", mendez,mcwidth
    print, galaxy," TRGB Color (wrong basis)       ", colortrgb,widthtrgb
    file=galaxy+'error.eps'

;;;;;;;;;;;;;;;;;;;;;;Using the color and magnitude of the tip stars, extrapolate the color of the tip;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    IF (KEYWORD_SET(OUTPLOT)) THEN BEGIN 
        set_plot,'ps'
        file=base+galaxy+'lum_' + strtrim(ct,2)+'.eps'
        device,filename=file,/color,bits_per_pixel = 8
    ENDIF ELSE BEGIN
        set_plot,'x'
        window,1
    ENDELSE
    !P.MULTI=[0,1,1]
    plot,binval,lumfunc,ytitle="Number Density",thick=2,/ylog,yrange=[1,MAX(lumfunc)*1.1],charsize=1.5,xtitle = ititle,xrange=[26,22.5];,title = "TRGB Magnitude for Run "+strtrim(ct,2)+', '+galaxy,
    oplot,binval,broken_powerlaw(binval, a_param, b, c, d, mendez),linestyle = 3;,color = 55
    print,'Parameters -- a: ',strtrim(a_param,2),'; b: ',strtrim(b,2),'; c: ',strtrim(c,2),'; d: ',strtrim(d,2),'; mtrgb: ',strtrim(mendez,2),'; c/d: ',c/d
    !Y.CRANGE=10.^!Y.CRANGE
    plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2
    readcol,'../Datafiles/lumtest_'+galaxy+'_'+strtrim(ct,2)+'.dat',x,y,/silent
    oplot,x,y,color = 100;240
    legend,['Luminosity Function','Simulated Luminosity Function','Best Fit'],linestyle = [0,0,3],thick = [2,1,1],color = [0,100,0],/left,/bottom,charsize = 1
    IF (KEYWORD_SET(OUTPLOT)) THEN device,/close

;    window,4
;    plot,points[*,0],points[*,1],psym = 3,yrange=[max(points[*,1]),Min(points[*,1])]
;    oplot,points[tipstars,0],points[tipstars,1],psym = 3, color = 240
;    oplot,[-30,30],[mendez,mendez],thick=2
;    oplot,[colortrgb,colortrgb],[-30,30],thick=2  
;    print,colortrgb

;Now, results rotate back
    tippoint=[colortrgb,mendez]
    IF (num_isos-2 eq 0 OR widthtrgb eq 0) THEN widthtrgb = MEAN(errors[tipstars,0]) 
    tiperror=[widthtrgb,mcwidth]
    tip_array = trans_rotate(lefttip[0],righttip[0],leftbound_first[0,*],leftbound_first[1,*],rightbound_first[0,*],rightbound_first[1,*],TRANSPOSE([tippoint]),/inverse)
    error_array = trans_rotate(lefttip[0],righttip[0],leftbound_first[0,*],leftbound_first[1,*],rightbound_first[0,*],rightbound_first[1,*],TRANSPOSE([tiperror]),/inverse,/error)
    colortrgb = tip_array[0]
    mendez =    tip_array[1]
    widthtrgb = ABS(error_array[0])
    mcwidth =   ABS(error_array[1])
    print, galaxy, " TRGB Magnitude  ", mendez,mcwidth
    print, galaxy, " TRGB Color      ", colortrgb,widthtrgb

;;;;;;;;Print Results to a file;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    printf,2,galaxy,colortrgb,widthtrgb,mendez,mcwidth,nind,FORMAT='(A20,4F7.3,I8)'
    print,galaxy,colortrgb,widthtrgb,mendez,mcwidth,nind

;Plots the CMB with the tip
;    IF NOT KEYWORD_SET(noplot) THEN BEGIN
;        set_plot,'x'
;        window,2
;        titlestring = "CMD for Run "+strtrim(ct,2)+', '+galaxy
;        plot,vminusi,imag,xrange=[-1,5],yrange=[28,22],xstyle=1,ystyle=1,psym=3,title=titlestring,xtitle=colortitle,ytitle=ititle
;        oplot,[-1,5],[mendez,mendez],thick=2
;        oplot,[colortrgb,colortrgb],[28,22],thick=2
;        IF (KEYWORD_SET(usevi)) THEN savefile=galaxy+"_trgb_vi.idl" ELSE savefile=+galaxy+"_trgb.idl"
;        SAVE,file=savefile
;    ENDIF

    temp = MIN(ABS(tips[1,*] - colortrgb),tip_ind)
    printf,1,isos_names[ct],isos_names[ct+1],colortrgb,widthtrgb,mendez,ABS(mcwidth),tips[2,tip_ind]- addy,addy,mendez - (tips[2,tip_ind]-addy),FORMAT='(2A20,7F20)'
ENDFOR
close,1

IF KEYWORD_SET(plottip) THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN 
        IF KEYWORD_SET(FAKEFILE) THEN plottip,galaxy,/fitcurve,/Fakefile,/OUTPLOT,adjust = addy ELSE plottip,galaxy,/fitcurve,/OUTPLOT, adjust = addy
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(FAKEFILE) THEN plottip,galaxy,/fitcurve,/Fakefile, adjust = addy ELSE plottip,galaxy,/fitcurve, adjust = addy
    ENDELSE
ENDIF
END
