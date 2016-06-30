function TRGB2,galaxy,vmag,verr,imag,ierr,leftbound = leftbound, rightbound = rightbound, USEVI=usevi,DOERROR=doerror,NOPLOT=noplot
;This code was written by Anil Seth and can be found in /net/mega-1/seth/halo/data/programs
cp 
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
ebv = 0

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
;now set limits on which stars to use
imin=20.5
imax=MAX(imag)
;ind=WHERE(vminusi GT 0.5 AND vminusi LT 1.9,nind) ;Discarding this to reflect that we want stars of other metallicites to
;be sampled -- Charlotte 3/9/07
lowerlim = MAX([MIN(imag)-0.2,23.6])
upperlim = 24.4
nind = (size(vmag))[1]
print,""
print,"__________________________________________"
print, "Using ", STRTRIM(nind,2), " Stars"

;;;;;;;;;;this is the function that detects the mag of the  TRGB;;;;;;;;;;;
;DET_TRGB,imag[ind],ierr[ind],imin,imax,binval,lumfunc,edgedet_mendez
;Changing this to reflect that we want stars of other metallicites to
;be sampled -- Charlotte 3/9/07
IF (KEYWORD_SET(rightbound) AND (KEYWORD_SET(leftbound))) THEN DET_TRGB,imag,ierr,imin,imax,binval,lumfunc,edgedet_mendez,leftbound = leftbound, rightbound = rightbound $
ELSE DET_TRGB,imag,ierr,imin,imax,binval,lumfunc,edgedet_mendez 
;This if statement allows for normalization for unequal binsize

possible = Where(binval GT lowerlim AND binval LT upperlim) 
;I'm tired of having unreasonable edge detections.  I'm going to hardwire this in.
max=MAX(edgedet_mendez[possible],maxpos) ;Find the maximum of the edge detection
mendez=binval[possible[maxpos]] ;Set this to the magnitude at the edge detection
print, galaxy, " TRGB magnitude", mendez,mendez+4.08
;see Seth et al. 2005 for derivation of M_F814W=-4.08 as TRGB magnitude

set_plot,'x'
window,2
;file='/net/mega-4/christensen/TRGBDust/'+galaxy+'/'+galaxy+'edge.eps'
;device,filename=file
xtitle='Luminosity Function'
;IF NOT KEYWORD_SET(noplot) THEN BEGIN
    !P.MULTI=[0,1,2]
    plot,binval,lumfunc,ytitle="Number Density",thick=2,/xsty,ymargin=[0,2],/ylog,yrange=[1,MAX(lumfunc)*1.1],/ysty,charsize=1,title = "TRGB Magnitude",xrange = [28,22],xcharsize=0.0001
    !Y.CRANGE=10.^!Y.CRANGE
    plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2

    plot,binval,edgedet_mendez,thick=2,xtitle='Magnitude',ytitle="ED Response",/xsty,ymargin=[2,0],xcharsize=1,charsize=1,yrange=[1.2*MIN(edgedet_mendez),1.2*MAX(edgedet_mendez)],/ysty,xrange = [28,22]
    plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2
;ENDIF


;;;;;;;;;;;;;;;;now determine colors at tip and tip-3.5;;;;;;;;;;;;;;;;;;;;;
window,3
!P.MULTI=[0,1,2]
tipstars=WHERE((imag GE mendez) AND (imag LE (mendez+0.2)),ntipstars) ;Selects out stars inbetween tip (mendez - 0.2)
IF (KEYWORD_SET(usevi)) THEN mag35=0.56  ELSE mag35=0.58
min35=mag35-0.1 & max35=mag35+0.1
stars35=WHERE(imag GE mendez+min35 AND imag LE (mendez+max35),ntipstars)
mincolor = MIN(vmag[stars35]-imag[stars35])
maxcolor = MAX(vmag[tipstars]-imag[tipstars])

If(N_ELEMENTS(tipstars) gt 1)  THEN BEGIN
    plothist,vmag[tipstars]-imag[tipstars],xhist,yhist,bin=0.02, xtitle = 'V-I', ytitle = 'Number',title = 'Color of TRGB',charsize=1,xrange=[mincolor,maxcolor]
    result=gaussfit(xhist,yhist,a,nterms=3)
    oplot,xhist,result
    colortrgb=a[1]
    widthtrgb=a[2]
ENDIF ELSE BEGIN
    colortrgb = vmag[tipstars]-imag[tipstars]
    widthtrgb = 0
ENDELSE    
print, galaxy," TRGB Color     ", colortrgb,widthtrgb

plothist,vmag[stars35]-imag[stars35],xhist,yhist,bin=0.02,xtitle='V-I',ytitle='Number',title='Color of Stars with Mag = TRGB + 0.56',charsize=1,xrange=[mincolor,maxcolor]
result=gaussfit(xhist,yhist,a,nterms=3)
oplot,xhist,result
color35=a[1]
width35=a[2]
print, galaxy," Color at -3.5  ", color35,width35


!P.MULTI = 0
;DEVICE,/CLOSE
;file=galaxy+'error.eps'
;Device,filename=file

;I'm not sure what this bit does.  I should find outcalculate Fe/H if using VI
IF (KEYWORD_SET(usevi)) THEN BEGIN
    bolocorr=0.881-0.243*colortrgb
    feh=-12.64+(12.6*color35)-(3.3*color35^2)
    print, galaxy," [Fe/H]",feh
    mitrgb=-(0.19*feh)-3.81-bolocorr
    print, "M_I,TRGB",mitrgb
ENDIF    


;;;;;;;;;;Finding Errors in mag ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;now determine errors ala mendez02, do monte carlo detection with
;random resampling of stars.
IF KEYWORD_SET(doerror) THEN BEGIN
    nruns=30
    mendezarr=FLTARR(nruns)
    smfarr=FLTARR(nruns)
    print, "Commencing Error Calc"
    FOR i=0,nruns-1 DO BEGIN
        print,i,FORMAT='(I,$)'
        random=LONG(RANDOMU(seed,nind)*FLOAT(nind))
;        inewerr=ierr[ind[random]]
        inewerr=ierr[random]
;        inew=imag[ind[random]]+RANDOMN(seed,nind)*inewerr
        inew=imag[random]+RANDOMN(seed,nind)*inewerr
        DET_TRGB,inew,inewerr,imin,imax,binval,lumfunc2,edgedet_mendez2
        possible = Where(binval GT lowerlim AND binval LT upperlim) 
;I'm tired of having unreasonable edge detections.  I'm going to hardwire this in.
        max=MAX(edgedet_mendez2[possible],maxpos)  
        mendezarr[i]=binval[possible[maxpos]]
    ENDFOR


    IF (KEYWORD_SET(usevi)) THEN BEGIN
        xtitle="I [Johnson]"
    ENDIF ELSE BEGIN
        xtitle="F814W [VEGAmag]"
    ENDELSE
    IF (MIN(mendezarr) ne MAX(mendezarr)) THEN BEGIN
        plothist,mendezarr,xhist,yhist,bin=0.01,xrange=[imin,imax],/xsty,ytitle="Number",xtitle=xtitle,ymargin=[4,-2],charsize=2
        plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2
    ENDIF

;now fit monte carlo results to gaussian
    estimates=[FLOAT(nruns)/10.,mendez,0.1]
    IF (N_Elements(xhist) gt 3) THEN BEGIN result=gaussfit(xhist,yhist,a,nterms=3,estimates=estimates) 
        mcpeak=a[1]
        mcwidth=a[2]
        oplot,xhist,result
    ENDIF ELSE BEGIN
        result=MOMENT(mendezarr)
        mcpeak = result[0]
        mcwidth = SQRT(result[1])
    ENDELSE

    print, 'MC Results',mcpeak,mcwidth
    IF (KEYWORD_SET(usevi)) THEN BEGIN 
        outfile=galaxy+"_trgb_vi.dat" 
        OPENW,2,outfile,WIDTH=1000
        printf,2,galaxy,mendez,mcpeak,mcwidth,colortrgb,widthtrgb,color35,width35,feh,mitrgb,nind,FORMAT='(A20,9F7.3,I8)'
    ENDIF ELSE BEGIN
        outfile=galaxy+"_trgb.dat"
        OPENW,2,outfile,WIDTH=1000
        printf,2,galaxy,mendez,mcpeak,mcwidth,colortrgb,widthtrgb,color35,width35,nind,FORMAT='(A20,7F7.3,I8)'
    ENDELSE
    CLOSE,2
    FREE_LUN,2
    doerror = mcwidth
    print,2,galaxy,mendez,mcpeak,mcwidth,colortrgb,widthtrgb,color35,width35,feh,mitrgb,nind
ENDIF ELSE mcwidth = 0

!P.MULTI=[0,1,1]
IF NOT KEYWORD_SET(noplot) THEN BEGIN
;DEVICE,/CLOSE
file=galaxy+'CMD.eps'
window,4
;Device,filename=file
;SET_PLOT,'ps'
!P.MULTI=[0,1,1]
IF (KEYWORD_SET(usevi)) THEN psfile=galaxy+"_trgbcmd_vi.ps" ELSE psfile=galaxy+"_trgbcmd.ps"
;myplot,file=psfile,ysize=5,yoff=1,/inch,xsize=6

IF (KEYWORD_SET(usevi)) THEN BEGIN
    xtitle="V-I [Johnson]"
    ytitle="I [Johnson]"
ENDIF ELSE BEGIN
    xtitle="F606W-F814W [VEGAmag]"
    ytitle="F814W [VEGAmag]"
ENDELSE

plot,vminusi,imag,xrange=[-1,5],yrange=[28,22],xstyle=1,ystyle=1,psym=3,title=titlestring,xtitle=xtitle,ytitle=ytitle
plots,[-1,5],[mendez,mendez],thick=2

;DEVICE,/CLOSE
SET_PLOT,'x'

IF (KEYWORD_SET(usevi)) THEN savefile=galaxy+"_trgb_vi.idl" ELSE savefile=+galaxy+"_trgb.idl"
SAVE,file=savefile
ENDIF
return_data = fltarr(7)
return_data[0] = mendez
return_data[1] = colortrgb
return_data[2] = widthtrgb
return_data[3] = mcwidth
return_data[4] = MEAN(imag[tipstars]) ;Returns the mean magnitude of the stars used to calculate the tip
return_data[5] = color35
return_data[6] = width35
return,return_data
END


PRO DET_TRGB,imag,ierr,imin,imax,binval,lumfunc,edgedet_mendez,leftbound = leftbound, rightbound = rightbound
binsize=0.01 ;Width of bin in magnitude
nbin=FIX((imax-imin)/binsize)+1 ;Number of bins
binval=FINDGEN(nbin)*binsize+imin ;Creates the initial magnitude for each bin

IF(KEYWORD_SET(rightbound) AND (KEYWORD_SET(leftbound))) THEN BEGIN
    length = N_ELEMENTS(rightbound[1,*])
    steps = (leftbound[1,length-1] - rightbound[1,length - 1])/binsize ;Steps in magnitude needed
    IF (steps gt 0) THEN BEGIN
        add=[rightbound[1,length-1]*fltarr(steps),findgen(steps)*binsize + rightbound[2,length-1]] ;Color remains the same for each of these but magnitude increases
        rightbound = [rightbound + add]
        print,'leftbound: ',leftbound
        print,'rightbound: ',rightbound
    ENDIF
ENDIF

;For some reason in the file the ierr is zero.  Not sure why and it
;will hvae to be fixed.  For now, I'll going to add .0001 to ever
;error (10% of minimum error)
ierr=ierr+.0001
lumfunc=FLTARR(nbin)
edgedet_mendez=FLTARR(nbin)
FOR i=0,nbin-1 DO BEGIN
    IF(KEYWORD_SET(rightbound) AND (KEYWORD_SET(leftbound))) THEN BEGIN ;Finds area of each bin and normalizes by it; if keyword is not set, no normalization
        x = MIN(ABS(rightbound[1,*] - binval[i]),index)
        rightcolor = rightbound[0,index]
        x = MIN(ABS(leftbound[1,*] - binval[i]),index)        
        leftcolor = leftbound[0,index]
        binarea = (rightcolor - leftcolor)*binsize
    ENDIF ELSE binarea = 1
;    print,binarea

    lumfunc[i]=TOTAL(1./(SQRT(2.*!PI)*ierr)*exp(-(imag-binval[i])^2./(2.*ierr^2.)))/binarea
    nearby=WHERE(imag GT binval[i]-0.05 AND imag LT binval[i]+0.05,nnearby)
    IF (nnearby GT 2.0) THEN BEGIN
        sigma=MEAN(ierr[nearby])
        edgedet_mendez[i]=-(ALOG10(TOTAL(1./(SQRT(2.*!PI)*ierr)*exp(-(imag+sigma-binval[i])^2./(2.*ierr^2.))))-ALOG10(TOTAL(1./(SQRT(2.*!PI)*ierr)*exp(-(imag-sigma-binval[i])^2./(2.*ierr^2.)))))*SQRT(lumfunc[i])
    ENDIF ELSE BEGIN
        edgedet_mendez[i]=0.0
    ENDELSE
ENDFOR
END
