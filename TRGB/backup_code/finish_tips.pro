;Needed modulues
;gbasis.pro
;iso.pro
;fitfunc.pro

FUNCTION BROKEN_POWERLAW, m, a, b, c, d, e, mtrgb
;This function gives the functional form of the tip of the RGB
  g = m
  m1= m*(-1.0)
  mtrgb1= mtrgb*(-1.0)
  y1 = a*(e/2.0) + d
  y2 = b*(-1*e/2.0) + c
  slope = (y2 - y1)/e
  g[WHERE(m1 GT (mtrgb1 + e/2))] = 10.0^(a*(m1[WHERE(m1 GT (mtrgb1 + e/2))] - mtrgb1)+ d)
  g[WHERE(m1 LT (mtrgb1 - e/2))] = 10.0^(b*(m1[WHERE(m1 LE (mtrgb1 - e/2))] - mtrgb1 )+ c) 
  g[WHERE(m1 LE (mtrgb1 + e/2) AND m1 GE (mtrgb1 - e/2))] = 10.0^(-1*slope*(m1[WHERE(m1 LE (mtrgb1 + e/2) AND m1 GE (mtrgb1 - e/2))] - (mtrgb1 +e/2.0 )) + y1)
  RETURN,g
END


;This procedure will plot all of the isochrones with their fitted
;tip.  It will also plot the difference between the observed and
;predicted tips.
PRO plottip,galaxy,fitcurve=fitcurve,FAKEFILE=fakefile,USEVI = USEVI, OUTPLOT = OUTPLOT
loadct,39
home = '/astro/net/scratch2/christensen/TRGBDust/'
iso_file = '/astro/net/scratch2/christensen/TRGBDust/Isochrones/metalgrid/'
ebv = 0.019;  Dunno about this.  Let's wait
ebv = 0.08
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
    difftitle = 'F814W_Observed - F814W_Predicted [VEGAmag]'
ENDELSE

IF(KEYWORD_SET(fakefile)) THEN BEGIN  
    infile = home+galaxy+'/' + galaxy + '.vi2'
    readcol,infile,fakev,v_err,fakei,i_err
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
        tipfit = curvefit(color,magnitude,w,a,FUNCTION_NAME = 'fitfunc')
        fitfunc,tipcurve_x,a,tipcurve_y
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
plot,mColor,mI,xtitle = colortitle, ytitle = ititle, title = 'Isochrone CMDs, '+galaxy,yrange=[25.5,22],xrange=[-0.5,4.0],psym=3
oploterror,color,magnitude,color_error,magnitude_error,psym = 3
oplot,color,magnitude,psym=2,color=220
IF (KEYWORD_SET(fitcurve)) THEN oplot,tipcurve_x,tipcurve_y,color=220
;*********Read in the isochrones used to divide the data***************************
isochrones = dblarr(3,elements,N_ELEMENTS(color));[quality, points, isochrones]
FOR i = 0, points-1 DO BEGIN
    file = extrap_iso(iso_file + iso1[i],elements)
    isochrones[0,*,i] = file.z                      ;metalicity
    isochrones[1,*,i] = file.f606mag - file.f814mag ;color
    isochrones[2,*,i] = file.f814mag + addy                ;imagnitude
    oplot,file.f606mag - file.f814mag,file.f814mag + addy,color =  220 - i*(220./(points+1))
ENDFOR
file = extrap_iso(iso_file + iso2[points - 1],elements)
oplot,file.f606mag - file.f814mag,file.f814mag + addy,color =  220 - points*(220./(points+1))
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
oplot,isochrones[1,*],isochrones[2,*],color=30
IF (KEYWORD_SET(OUTPLOT)) THEN device,/close

;*************Plot difference between observed and predicted tip**********************************
print,"Fitting Curve"
IF (KEYWORD_SET(OUTPLOT)) THEN BEGIN 
    set_plot,'ps'
    device,filename=home + galaxy+'/iso_diff'+galaxy+'.eps',/color,bits_per_pixel=8
ENDIF ELSE BEGIN
    stop
    set_plot,'x'
    window,4
ENDELSE
a_pred = [1.,1.,1.,1.,1.,1.] ;Fit a curve to the set of modeled isochrone tips
tipfit = curvefit(isochrones[1,*],isochrones[2,*],w,a_pred,FUNCTION_NAME = 'fitfunc')
fitfunc,color,a_pred,isotippoint_y
plot,color,magnitude-isotippoint_y,psym=2,xtitle = colortitle,ytitle = difftitle,title='Difference Between Predicted and Actual TRGB, '+galaxy,yrange = [MIN(magnitude-isotippoint_y)-MAX(magnitude_error),MAX(magnitude-isotippoint_y)+MAX(magnitude_error)]
oploterror,color,magnitude-isotippoint_y,color_error,magnitude_error,psym = 3

IF (KEYWORD_SET(fitcurve)) THEN BEGIN
    fitfunc,tipcurve_x,a_pred,isotipcurve_y
    oplot,tipcurve_x,tipcurve_y-isotipcurve_y
ENDIF
oplot,[0,4],[0,0],linestyle=2
IF (KEYWORD_SET(OUTPLOT)) THEN device,/close
END

;*********************************************************************************************************

;This program will read in a list of tip files (generated by trgb.cpp)
;and finish finding the magnitude and plotting
PRO finish_tips,galaxy,tipfiles_name=tipfiles_name,lumfiles_name = lumfiles_name,FAKEFILE=fakefile, USEVI = USEVI, PLOTTIP = PLOTTIP, OUTPLOT = OUTPLOT
close,/all
loadct,39
IF KEYWORD_SET(OUTPLOT) THEN OUTPLOT = 1 ELSE OUTPLOT = 0
IF outplot THEN set_plot,'ps' ELSE set_plot,'x'

IF NOT KEYWORD_SET(tipfiles_name) THEN tipfiles_name = "../Datafiles/found_tips_list"+galaxy+".txt"
IF NOT KEYWORD_SET(lumfiles_name) THEN lumfiles_name = "../Datafiles/lum_func_list"+galaxy+".txt"
readcol,tipfiles_name,format = 'A60',tipfiles
readcol,lumfiles_name,format = 'A60',lumfiles
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
    colortitle="F814W [VEGAmag]"        
    ititle="F814W [VEGAmag]"
ENDELSE

IF(KEYWORD_SET(fakefile)) THEN BEGIN  
    infile = base + galaxy + '.vi2'
    readcol,infile,fakev,v_err,fakei,i_err
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
readcol,'../Datafiles/metals' + galaxy + '.txt',addy,FORMAT='F'
readcol,'../Datafiles/metals' + galaxy + '.txt',strmetals,FORMAT='A'
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
    basis = GBASIS_tips(lefttip[0],righttip[0],leftbound[0,*],leftbound[1,*],rightbound[0,*],rightbound[1,*])
    set = Where(imag gt 23.0 AND imag lt 25.0)
    points = [[vmag[set]-imag[set]],[imag[set]]]#basis
    errors = ABS([[(verr[set]+ierr[set])/2.0],[ierr[set]]]#basis)
    rightbound = transpose(rightbound)#basis
    leftbound = transpose(leftbound)#basis

    readcol,tipfiles[ct],mendezerr,FORMAT='F'
    readcol,lumfiles[ct],binval,lumfunc,FORMAT='F,F'
    a_param = binval[N_ELEMENTS(binval) - 3]
    c = binval[N_ELEMENTS(binval) - 2]  
    e = binval[N_ELEMENTS(binval) - 1] 
    binval = binval[0:N_ELEMENTS(binval) - 4]
    b = lumfunc[N_ELEMENTS(lumfunc) - 3]
    d = lumfunc[N_ELEMENTS(lumfunc) - 2]    
    lumfunc = lumfunc[0:N_ELEMENTS(lumfunc) - 4]
    mendez = mendezerr[0]
    mendezerr = mendezerr[1:N_ELEMENTS(mendezerr)-1]

    IF (KEYWORD_SET(OUTPLOT)) THEN BEGIN 
 ;       device,/close
        set_plot,'ps'
        device, filename = galaxy+'error_'+strtrim(ct,2)+'.eps'
    ENDIF ELSE BEGIN
        set_plot,'x'
        window,0
    ENDELSE
    !P.MULTI=[0,1,2]
    IF (MIN(mendezerr) ne MAX(mendezerr)) THEN BEGIN
        plothist,mendezerr,xhist,yhist,bin=0.01,/xsty,ytitle="Number",xtitle=ititle,title = "Error for Run "+strtrim(ct,2)+', '+galaxy;,charsize=2,xrange=[imin,imax]
        plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2
    ENDIF

;now fit monte carlo results to gaussian
    nruns = N_ELEMENTS(mendezerr)
    estimates=[FLOAT(nruns)/10.,mendez,0.1]
    IF (N_Elements(xhist) gt 3) THEN BEGIN 
        yhist_cut = yhist
        result=gaussfit(xhist,yhist_cut,a,nterms=3) 
        mcpeak=a[1]
        mcwidth=a[2]
        oplot,xhist,result
    ENDIF ELSE BEGIN
        result=MOMENT(mendezerr)
        mcpeak = result[0]
        mcwidth = SQRT(result[1])
    ENDELSE
    IF (mcwidth gt (MAX(xhist) - MIN(xhist))/2) THEN mcwidth = stddev(mendezerr)

;see Seth et al. 2005 for derivation of M_F814W=-4.08 as TRGB magnitude
;;;;;;;;;;;;;;;;now determine colors at tip;;;;;;;;;;;;;;;;;;;;
    tipstars=WHERE((points[*,0] LE mendez) AND (points[*,0] GE (mendez-0.15)),ntipstars) ;Selects out stars right below the tip
    If(N_ELEMENTS(tipstars) gt 25)  THEN BEGIN
        plothist,points[tipstars,1],xhist,yhist,bin=0.02, colortitle = 'V-I', ytitle = 'Number',title = 'Color of TRGB',charsize=1
        result=gaussfit(xhist,yhist,a,nterms=3)
        oplot,xhist,result
        colortrgb=a[1]
        plots,[colortrgb,colortrgb],!Y.CRANGE,thick=2,linestyle=2
        widthtrgb=MIN([a[2],(MAX(points[tipstars,1]) - MIN(points[tipstars,1]))/2.0])
    ENDIF ELSE BEGIN
        colortrgb = MEAN(points[tipstars,1])
        widthtrgb =(MAX(points[tipstars,1]) - MIN(points[tipstars,1]))/2.0
    ENDELSE 
    IF (colortrgb lt MIN(points[tipstars,1]) OR colortrgb GT MAX(points[tipstars,1])) THEN BEGIN
        a = moment(points[tipstars,1])
        colortrgb=a[0]
        widthtrgb=(MAX(points[tipstars,1]) - MIN(points[tipstars,1]))/2.0
    ENDIF
    colortrgb = MEAN(points[tipstars,1])
    print,'Inital Width: ', widthtrgb
  ;;  widthtrgb = MEAN([mI_err, mV_err])
    print,'Modified Width: ',MEAN([mI_err])
    print,' '
    print, galaxy, " TRGB Magnitude (wrong basis)  ", mendez,mcwidth
    print, galaxy," TRGB Color (wrong basis)       ", colortrgb,widthtrgb
    file=galaxy+'error.eps'

;;;;;;;;;;;;;;;;;;;;;;Using the color and magnitude of the tip stars, extrapolate the color of the tip;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    IF (KEYWORD_SET(OUTPLOT)) THEN BEGIN 
        device,/close
        set_plot,'ps'
        file=base+galaxy+'lum_' + strtrim(ct,2)+'.eps'
        device,filename=file,/color,bits_per_pixel = 8
    ENDIF ELSE BEGIN
        set_plot,'x'
        window,1
    ENDELSE
    !P.MULTI=[0,1,1]
    plot,binval,lumfunc,ytitle="Number Density",thick=2,/ylog,yrange=[1,MAX(lumfunc)*1.1],charsize=1,title = "TRGB Magnitude for Run "+strtrim(ct,2)+', '+galaxy,/ysty,xtitle = ititle
    oplot,binval,broken_powerlaw(binval,a_param, b, c, d,e, mendez),linestyle = 2,color = 55
    !Y.CRANGE=10.^!Y.CRANGE
    plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2
    readcol,'../Datafiles/lumtest_'+galaxy+'_'+strtrim(ct,2)+'.dat',x,y
    oplot,x,y,color = 240

;    window,4
;    plot,points[*,0],points[*,1],psym = 3
;    oplot,points[tipstars,0],points[tipstars,1],psym = 3, color = 240
;    oplot,[mendez,mendez],[-30,30],thick=2
;    oplot,[-30,30],[colortrgb,colortrgb],thick=2    
;stop

    IF (KEYWORD_SET(OUTPLOT)) THEN device,/close
    colortrgb = interpol_color(colortrgb,MEAN(points[tipstars,0]),mendez,rightbound[*,0],rightbound[*,1],leftbound[*,0],leftbound[*,1])
;Now, results rotate back
    inverse_basis = -1.0*basis
    tippoint=[mendez,colortrgb]
    IF (num_isos-2 eq 0 OR widthtrgb eq 0) THEN widthtrgb = MEAN(errors[tipstars,0]) 
    tiperror=[mcwidth,widthtrgb]
    colortrgb = -1.0*(tippoint#inverse_basis)[0]
    mendez =    -1.0*(tippoint#inverse_basis)[1]
    widthtrgb = (tiperror#inverse_basis)[0]
    mcwidth =   (tiperror#inverse_basis)[1]
    print, galaxy, " TRGB Magnitude  ", mendez,mcwidth
    print, galaxy, " TRGB Color      ", colortrgb,widthtrgb

;;;;;;;;Print Results to a file;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    printf,2,galaxy,colortrgb,widthtrgb,mendez,mcwidth,nind,FORMAT='(A20,4F7.3,I8)'
    print,galaxy,colortrgb,widthtrgb,mendez,mcwidth,nind

;Plots the CMB with the tip
    IF NOT KEYWORD_SET(noplot) THEN BEGIN
        set_plot,'x'
        window,2
        titlestring = "CMD for Run "+strtrim(ct,2)+', '+galaxy
        plot,vminusi,imag,xrange=[-1,5],yrange=[28,22],xstyle=1,ystyle=1,psym=3,title=titlestring,xtitle=colortitle,ytitle=ititle
        oplot,[-1,5],[mendez,mendez],thick=2
        oplot,[colortrgb,colortrgb],[28,22],thick=2
        IF (KEYWORD_SET(usevi)) THEN savefile=galaxy+"_trgb_vi.idl" ELSE savefile=+galaxy+"_trgb.idl"
        SAVE,file=savefile
    ENDIF
    temp = MIN(ABS(tips[1,*] - colortrgb),tip_ind)
    printf,1,isos_names[ct],isos_names[ct+1],colortrgb,widthtrgb,mendez,ABS(mcwidth),tips[2,tip_ind]- addy,addy,mendez - (tips[2,tip_ind]-addy),FORMAT='(2A20,7F20)'
ENDFOR
close,1

IF KEYWORD_SET(plottip) THEN BEGIN
    IF KEYWORD_SET(outplot) THEN BEGIN 
        IF KEYWORD_SET(FAKEFILE) THEN plottip,galaxy,/fitcurve,/Fakefile,/OUTPLOT ELSE plottip,galaxy,/fitcurve,/OUTPLOT 
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(FAKEFILE) THEN plottip,galaxy,/fitcurve,/Fakefile ELSE plottip,galaxy,/fitcurve
    ENDELSE
ENDIF
END
