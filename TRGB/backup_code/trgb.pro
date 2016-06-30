;FUNCTION interpol_color,x1,y1,y2,yline1,xline1,yline2,xline2
;segment1 = where(yline1 LE y2 AND yline1 GE y1)
;IF(segment1[0] EQ -1) THEN diff = MIN(ABS(yline1-y1),segment1)
;segment1 = [MIN(segment1)-1,segment1,MAX(segment1)+1];This brackets off the points below and above the magnitdues we are interested in
;segment2 = where(yline2 LE y2 AND yline2 GE y1)
;IF(segment2[0] EQ -1) THEN diff = MIN(ABS(yline2-y1),segment2)
;segment2 = [MIN(segment2)-1,segment2,MAX(segment2)+1]
;xs1 = SPLINE(yline1[segment1],xline1[segment1],[y1,y2])
;xs2 = SPLINE(yline2[segment2],xline2[segment2],[y1,y2])
;weight1 = (x1-xs1[0])/(xs2[0]-xs1[0]) ;weight by how close to each boundery the point is
;x2  = xs1[1] + weight1*(xs2[1]-xs1[1])
;RETURN,x2
;END

;FUNCTION GBASIS,x1,x2,xline1,yline1,xline2,yline2
;Generates an othonormal basis function for which one of the basis
;vectors is the average slope of the two bounding isochrones near the
;tip
;xline1 = reform(xline1)
;yline1 = reform(yline1)
;xline2 = reform(xline2)
;yline2 = reform(yline2)

;temp = MIN(ABS(xline1-x1),tip1)
;temp = MIN(ABS(xline2-x2),tip2) 
;x1 = ((xline1[tip1+1] - xline1[tip1-1]) + (xline2[tip2+1] - xline2[tip2-1]))/2
;y1 = ((yline1[tip1+1] - yline1[tip1-1]) + (yline2[tip2+1] - yline2[tip2-1]))/2
;mag = SQRT(x1*x1+y1*y1)
;basis = [[x1/mag,y1/mag],[y1/mag,-1.0*x1/mag]]
;RETURN,basis
;END

FUNCTION TRGB,galaxy,vmag,verr,imag,ierr,lefttip,righttip,leftbound,rightbound,iter,USEVI=usevi,NOPLOT=noplot
;This code was written by Anil Seth and can be found in /net/mega-1/seth/halo/data/programs
;I (Charlotte Christensen) have made numerous modifications

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
;commedted out about above because galaxy names don't match and I'm
;not sure what to use, set the reddening below, from data in the file

;;;;;;;;;Add in Reddening and select stars;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ebv = 0.019;  Dunno about this.  Let's wait
outfile_stars = 'trgbstars.dat'
outfile_mags = 'magANDerr.dat'
infile_tip = 'found_tips_'+ STRTRIM(iter,2)+'.dat'
infile_lumfunc = 'lum_func_'+strtrim(iter,2)+'.dat'

IF (KEYWORD_SET(usevi)) THEN BEGIN
    vred=ebv*3.1
    ired=ebv*3.1*0.482 ;B&M p.136-7
ENDIF ELSE BEGIN
    vred=ebv*2.716     
    ired=ebv*1.796 ;from Roelof 16 June 2004 email, M0 star
ENDELSE
vmag=vmag-vred
imag=imag-ired
vminusi=vmag-imag
nind = (size(vmag))[1]
print,""
print,"__________________________________________"
print, "Using ", STRTRIM(nind,2), " Stars"

;I want to rotate the frame of reference such that the average slope
;of the bounding isochrones is one of the set of two orthonganal basis
basis = GBASIS(lefttip[0],righttip[0],leftbound[0,*],leftbound[1,*],rightbound[0,*],rightbound[1,*])
lowerlim = [(lefttip[0]+righttip[0])/2.0,MIN([lefttip[1],righttip[1]])-0.5]
upperlim = [(lefttip[0]+righttip[0])/2.0,MAX([lefttip[1],righttip[1]])+0.5]
;Rotates each point in an array of points by the provided basis and
;return the values
set = Where(imag gt 23.0 AND imag lt 25.0)
points = [[vmag[set]-imag[set]],[imag[set]]]#basis
errors = ABS([[(verr[set]+ierr[set])/2.0],[ierr[set]]]#basis)
upperlim = upperlim#basis
lowerlim = lowerlim#basis
rightbound = transpose(rightbound)#basis
leftbound = transpose(leftbound)#basis

;;;;;;;;;;this is the function that detects the mag of the  TRGB;;;;;;;;;;;
;Changing this to reflect that we want stars of other metallicites to
;be sampled -- Charlotte 3/9/07
imin = MIN(points[*,0])
imax = MAX(points[*,0])
;DET_TRGB2,points[*,0],errors[*,0],imin,imax,binval,lumfunc,edgedet_mendez 
;possible = Where(binval LT lowerlim[0] AND binval GT upperlim[0]) 
;I'm tired of having unreasonable edge detections.  I'm going to hardwire this in.
;max=MAX(edgedet_mendez[possible],maxpos) ;Find the maximum of the edge detection
;mendez=binval[possible[maxpos]] ;Set this to the magnitude at the
;edge detection

;*****************************Switch to C*****************************************
;mendez = ML_DET_TRGB(imag,ierr,imin,imax,binval,lumfunc)
;My tip determination cannot be done quickly in idl so I am doing it
;in c
;This program will output a file of stars and stop
;When this happens, run trgb.c
;That will also output a file.  Continue the idl program which will
;read in the dat
dm = 0.01 ; width of magnitude bin
numm = FIX((imax - imin)/dm) + 1
m = FINDGEN(numm)*dm + imin ;array of magnitude
sigmafit = poly_fit(points[*,0],errors[*,0],4)
merr = sigmafit[0] + sigmafit[1]*m + sigmafit[2]*m^2.0+sigmafit[3]*m^3 + sigmafit[4]*m^4 
;grid_erf2 = erf_m(m,m,merr)

OPENW,3,outfile_stars,WIDTH=1000
printf,3,TRANSPOSE([[points[*,0]],[errors[*,0]]]),FORMAT='(F,F)'
close,3

lowertemp = MIN(ABS(m - lowerlim[0]),lowerpos)
uppertemp = MIN(ABS(m - upperlim[0]),upperpos)

OPENW,3,outfile_mags,WIDTH=1000
printf,3,upperpos,lowerpos,FORMAT='(F,F)'
printf,3,TRANSPOSE([[m],[merr]]),FORMAT='(F,F)'
close,3

print,'Now run the c program, trgb.c'

readcol,infile_tip,mendezerr,FORMAT='F'
readcol,infile_lumfunc,binval,lumfunc,FORMAT='F,F'

mendez = mendezerr[0]
mendezerr = mendezerr[0:N_ELEMENTS(mendezerr)-2]


;;;;;;;;;;Finding Errors in mag ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;now determine errors ala mendez02, do monte carlo detection with
;random resampling of stars.
;nruns=5
;mendezarr=FLTARR(nruns)
;smfarr=FLTARR(nruns)
;print, "Commencing Error Calc"
;FOR i=0,nruns-1 DO BEGIN
;    print,i,FORMAT='(I,$)'
;    random=LONG(RANDOMU(seed,nind)*FLOAT(nind))
;    inewerr=points[random,0] (shouldn't it be errors[random,0]?
;    inew=points[random,0]+RANDOMN(seed,nind)*inewerr
;    DET_TRGB2,inew,inewerr,imin,imax,binval,lumfunc2,edgedet_mendez2
;    possible = Where(binval LT lowerlim[0] AND binval GT upperlim[0]) 
;I'm tired of having unreasonable edge detections.  I'm going to hardwire this in.
;    max=MAX(edgedet_mendez2[possible],maxpos)  
;    mendezarr[i]=binval[possible[maxpos]]
;     mendezarr[i] = ML_DET_TRGB(inew,inewerr,imin,imax,binval,lumfunc2)
;ENDFOR
;Plot errors

window,4
IF (KEYWORD_SET(usevi)) THEN BEGIN
    xtitle="I [Johnson]"
ENDIF ELSE BEGIN
    xtitle="F814W [VEGAmag]"
ENDELSE
IF (MIN(mendezerr) ne MAX(mendezerr)) THEN BEGIN
     plothist,mendezerr,xhist,yhist,bin=0.01,xrange=[imin,imax],/xsty,ytitle="Number",xtitle=xtitle,charsize=2   
    plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2
ENDIF

;now fit monte carlo results to gaussian
nruns = N_ELEMENTS(mendezerr)
estimates=[FLOAT(nruns)/10.,mendez,0.1]
IF (N_Elements(xhist) gt 3) THEN BEGIN 
    result=gaussfit(xhist,yhist,a,nterms=3,estimates=estimates) 
    mcpeak=a[1]
    mcwidth=a[2]
    oplot,xhist,result
ENDIF ELSE BEGIN
    result=MOMENT(mendezerr)
    mcpeak = result[0]
    mcwidth = SQRT(result[1])
ENDELSE
print,' '
print,'MC Results',mcpeak,mcwidth
!P.MULTI=0



;***************
print, galaxy, " TRGB Magnitude (wrong basis)  ", mendez
;see Seth et al. 2005 for derivation of M_F814W=-4.08 as TRGB magnitude
;;;;;;;;;;;;;;;;now determine colors at tip;;;;;;;;;;;;;;;;;;;;;
window,3
!P.MULTI=[0,1,2]
tipstars=WHERE((points[*,0] LE mendez) AND (points[*,0] GE (mendez-0.5)),ntipstars) ;Selects out stars right below the tip
If(N_ELEMENTS(tipstars) gt 1)  THEN BEGIN
    plothist,points[tipstars,1],xhist,yhist,bin=0.02, xtitle = 'V-I', ytitle = 'Number',title = 'Color of TRGB',charsize=1
    result=gaussfit(xhist,yhist,a,nterms=3)
    oplot,xhist,result
    colortrgb=a[1]
    widthtrgb=a[2]
ENDIF ELSE BEGIN
    colortrgb = MEAN(points[tipstars,1])
    widthtrgb = 0
ENDELSE    
print, galaxy," TRGB Color (wrong basis)     ", colortrgb,widthtrgb
;DEVICE,/CLOSE
;file=galaxy+'error.eps'
;Device,filename=file

;;;;;;;;;;;;;;;;;;;;;;Using the color and magnitude of the tip stars,
;extrapolate the color of the tip;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

set_plot,'x'
window,2
;set_plot,'ps'
;file=galaxy+'lum.eps'
;device,filename=file
;IF NOT KEYWORD_SET(noplot) THEN BEGIN
    !P.MULTI=[0,1,2]
    plot,binval,lumfunc,ytitle="Number Density",thick=2,/xsty,ymargin=[0,2],/ylog,yrange=[1,MAX(lumfunc)*1.1],charsize=1,title = "TRGB Magnitude",xcharsize=0.0001,/ysty
    !Y.CRANGE=10.^!Y.CRANGE
    plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2
;ENDIF

colortrgb = interpol_color(colortrgb,MEAN(points[tipstars,0]),mendez,rightbound[*,0],rightbound[*,1],leftbound[*,0],leftbound[*,1])
;Now, results rotate back
inverse_basis = -1.0*basis
tippoint=[mendez,colortrgb]
tiperror=[mcwidth,widthtrgb]
colortrgb = -1.0*(tippoint#inverse_basis)[0]
mendez = -1.0*(tippoint#inverse_basis)[1]
widthtrgb = (tiperror#inverse_basis)[0]
mcwidth = (tiperror#inverse_basis)[1]
print, galaxy, " TRGB Magnitude  ", mendez,mcwidth
print, galaxy, " TRGB Color      ", colortrgb,widthtrgb

;IF NOT KEYWORD_SET(noplot) THEN BEGIN
    plot,vminusi,imag,xrange=[-1,5],yrange=[28,22],xstyle=1,ystyle=1,psym=3,title=titlestring,xtitle='V - I',ytitle='I'
    oplot,[-1,5],[mendez,mendez],thick=2
    oplot,[colortrgb,colortrgb],[28,22],thick=2
;ENDIF
;device,/close
;set_plot,'x'

;;;;;;;;Print Results to a file;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF (KEYWORD_SET(usevi)) THEN BEGIN 
    outfile=galaxy+"_trgb_vi.dat" 
    OPENW,2,outfile,WIDTH=1000
    printf,2,galaxy,colortrgb,widthtrgb,mendez,mcwidth,nind,FORMAT='(A20,4F7.3,I8)'
ENDIF ELSE BEGIN
    outfile=galaxy+"_trgb.dat"
    OPENW,2,outfile,WIDTH=1000
    printf,2,galaxy,colortrgb,widthtrgb,mendez,mcwidth,nind,FORMAT='(A20,4F7.3,I8)'
ENDELSE
CLOSE,2
;FREE_LUN,2
;print,galaxy,colortrgb,widthtrgb,mendez,mcwidth,nind

;Plots the CMB with the tip
;!P.MULTI=[0,1,1]
IF NOT KEYWORD_SET(noplot) THEN BEGIN
    ;DEVICE,/CLOSE
    file=galaxy+'CMD.eps'
    ;Device,filename=file
    ;SET_PLOT,'ps'
;    !P.MULTI=[0,1,1]
    IF (KEYWORD_SET(usevi)) THEN BEGIN
        psfile=galaxy+"_trgbcmd_vi.ps"
        xtitle="V-I [Johnson]"
        ytitle="I [Johnson]"
    ENDIF ELSE BEGIN
        psfile=galaxy+"_trgbcmd.ps"
        xtitle="F606W-F814W [VEGAmag]"
        ytitle="F814W [VEGAmag]"
    ENDELSE
;    plot,vminusi,imag,xrange=[-1,5],yrange=[28,22],xstyle=1,ystyle=1,psym=3,title=titlestring,xtitle=xtitle,ytitle=ytitle
;    oplot,[-1,5],[mendez,mendez],thick=2
;    oplot,[colortrgb,colortrgb],[28,22],thick=2
    ;DEVICE,/CLOSE
    SET_PLOT,'x'
    IF (KEYWORD_SET(usevi)) THEN savefile=galaxy+"_trgb_vi.idl" ELSE savefile=+galaxy+"_trgb.idl"
    SAVE,file=savefile
ENDIF
return_data = fltarr(4)
return_data[0] = colortrgb
return_data[1] = widthtrgb
return_data[2] = mendez
return_data[3] = mcwidth
return,return_data
END

FUNCTION erf_m,m,m2,sigma
;simple error function
;sigma is the error for m2 so they must have the same dimensions
step1 = [[m2],[sigma]]
step2 = TRANSPOSE(REBIN(step1,N_ELEMENTS(m2),2,N_ELEMENTS(m)),[2,0,1])
m_grid = REBIN(m,N_ELEMENTS(m),N_ELEMENTS(m2))
variable_grid = [[[m_grid]],[[step2[*,*,0]]],[[step2[*,*,1]]]]
;,(1./(SQRT(2+!PI)*sigma)*exp(-(m-m2)^2./(2*sigma^2.)))
RETURN,(1./(SQRT(2+!PI)*variable_grid[*,*,2])*exp(-(variable_grid[*,*,0]-variable_grid[*,*,1])^2./(2*variable_grid[*,*,2]^2.)))
;grid has dimensions (m,m2)
END

;FUNCTION broken_powerlaw,m,b,c,mtrgb
;broken power law to fit to luminosity function
;a = -0.3
;g = fltarr(N_ELEMENTS(m))
;g[where(m GT mtrgb)] = 10.0^(a*(m[where(m GT mtrgb)] - mtrgb))

;g[where(m LE mtrgb)] = 10.0^(b*(m[where(m LE mtrgb)] - mtrgb[where(m LE mtrgb)])-c)
;g has dimensions (m)
;RETURN,g
;END

FUNCTION ML_DET_TRGB,imag,ierr,imin,imax,m,lumfunc
;From Mendez 2002, maximum likelyhood method
dm = 0.01 ; width of magnitude bin
dc = 0.01 
db = 0.01 
cmin = 0
cmax = 1.2
bmin = 0.5
bmax = 1.0 
numb = (bmax - bmin)/db
numc = (cmax - cmin)/dc
imin = MIN(imag)
imax = MAX(imag)
numm = FIX((imax - imin)/dm) + 1
numi = N_ELEMENTS(imag)
b = FINDGEN(numb)*db + bmin
c = FINDGEN(numc)*dc + cmin
m = FINDGEN(numm)*dm + imin ;array of magnitudes
phi = FLTARR(numc)
phi2 = FLTARR(numb)
ln_likely = FLTARR(numm)

lumfunc=TOTAL(TOTAL(erf_m(m,imag,ierr),1),1)
grid_erf1 = erf_m(m,imag,ierr)

sigmafit = poly_fit(imag,ierr,4)
merr = sigmafit[0] + sigmafit[1]*m + sigmafit[2]*m^2.0+sigmafit[3]*m^3 + sigmafit[4]*m^4 
grid_erf2 = erf_m(m,m,merr)

FOR mct = 0, numm-1 DO BEGIN
    FOR bct = 0, numb-1 DO BEGIN
        FOR cct = 0, numc-1 DO BEGIN
            grid_g = broken_powerlaw(m,b[bct],c[cct],m[mct])
            grid_g1 = TRANSPOSE(REBIN(grid_g,numm,numi),[1,0])
            grid_g2 = TRANSPOSE(REBIN(grid_g,numm,numm),[1,0])
            phi[cct] = TOTAL(ALOG(TOTAL(grid_g1*grid_erf1*dm,1)),1) -  N_ELEMENTS(imag) * ALOG(TOTAL(TOTAL(grid_g2*grid_erf2*dm,1)*dm,1))
        ENDFOR
        phi2[bct] = TOTAL(phi*dc)
    ENDFOR
    ln_likely[mct] = TOTAL(phi2*dc)
ENDFOR   
max_l = MAX(ln_likely,maxct)
stop
RETURN,m[maxct]
END


PRO DET_TRGB2,imag,ierr,imin,imax,binval,lumfunc,edgedet_mendez
binsize=0.01 ;Width of bin in magnitude
nbin=FIX((imax-imin)/binsize)+1 ;Number of bins
binval=FINDGEN(nbin)*binsize+imin ;Creates the initial magnitude for each bin
;For some reason in the file the ierr is zero.  Not sure why and it
;will hvae to be fixed.  For now, I'll going to add .0001 to ever
;error (10% of minimum error)
;ierr=ierr+.0001
;stop
lumfunc=FLTARR(nbin)
edgedet_mendez=FLTARR(nbin)
print,'Number of bins: '
FOR i=0,nbin-1 DO BEGIN
    lumfunc[i]=TOTAL(1./(SQRT(2.*!PI)*ierr)*exp(-(imag-binval[i])^2./(2.*ierr^2.)))
    nearby=WHERE(imag GT binval[i]-0.05 AND imag LT binval[i]+0.05,nnearby)
    IF (nnearby GT 2.0) THEN BEGIN
;        shift = findgen(10)*0.1 - 0.5
;        sigma=shift*MEAN(ierr[nearby])
        sigma = MEAN(ierr[nearby])
;        print,sigma
        edgedet_mendez[i]=-(ALOG10(TOTAL(1./(SQRT(2.*!PI)*ierr)*exp(-(imag+sigma-binval[i])^2./(2.*ierr^2.))))-ALOG10(TOTAL(1./(SQRT(2.*!PI)*ierr)*exp(-(imag-sigma-binval[i])^2./(2.*ierr^2.)))))*SQRT(lumfunc[i])
;        high[i] = ALOG10(TOTAL(1./(SQRT(2.*!PI)*ierr)*exp(-(imag+sigma-binval[i])^2./(2.*ierr^2.)))
    ENDIF ELSE BEGIN
        edgedet_mendez[i]=0.0
    ENDELSE
    print,i
ENDFOR
;stop
END


