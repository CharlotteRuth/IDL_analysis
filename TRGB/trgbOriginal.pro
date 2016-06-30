;Anil's program to determine the TRGB

 
PRO TRGBOriginal,galaxy,USEVI=usevi,DOERROR=doerror,NOPLOT=noplot,REDDENING = reddening

;read in data - use vega2vi to transform to Johnson if necessary
;psfphotfile=galaxy+"_psfphot.dat";Comentend out by CC

;IF (KEYWORD_SET(usevi)) THEN BEGIN;Comentend out by CC
;    RDFLOAT, psfphotfile,x,dx,y,dy,f606,verr,f814,ierr;Comentend out by CC
;    vega2vi, f606,f814,vmag,imag;Comentend out by CC
;ENDIF ELSE RDFLOAT, psfphotfile,x,dx,y,dy,vmag,verr,imag,ierr;Comentend out by CC
base = '/astro/net/scratch2/christensen/TRGBDust/'

infile=base + galaxy + '/' + galaxy + '.vi2.fits'
columns = ['mag1_acs', 'mag1_err','mag2_acs','mag2_err']
data = mrdfits(infile,1,columns=columns)
imag = data.mag2_acs
ierr = data.mag2_err
vmag = data.mag1_acs
verr = data.mag1_err

;now determine reddening 
IF(KEYWORD_SET(REDDENING)) THEN READCOL, 'programs/reddenings.dat',rgalaxy,ebv,FORMAT='A,F',/SILENT ELSE ebv = 0.019

;match=WHERE(rgalaxy EQ galaxy)  ;Comentend out by CC
;ebv=ebv[match[0]];Comentend out by CC
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
;get_fitdata,galaxy,xcen,ycen,pa,sl,z0
;zhalf=0.5493*z0
;galactic_coord,x,y,xgal,ygal,xcen,ycen,pa
imin=20.5 & imax=26.0

;IF (galaxy EQ 'ngc0055-disk' OR galaxy EQ 'ngc0055') THEN BEGIN
;    imin=20.5 & imax=24.0 
;ENDIF 

;scaleheight=3.0
;IF (galaxy EQ 'ngc0055-disk') THEN scaleheight=3.0
;IF (galaxy EQ 'ngc4631' OR galaxy EQ 'ngc4631-disk') THEN scaleheight=7.0;


;Gilmore & Reid 1983 suggest thick disk dominates above ~5 scale heights.
;IF (KEYWORD_SET(usevi)) THEN $
;  ind=WHERE(ABS(ygal) GT scaleheight*zhalf AND vminusi GT 0.5 AND vminusi LT 1.9,nind) ELSE $
;  ind=WHERE(ABS(ygal) GT scaleheight*zhalf AND vminusi GT 0.3 AND vminusi LT 1.6,nind);

;IF (galaxy EQ 'ngc4631') THEN BEGIN
;    ind2=WHERE(ABS(ygal[ind]) LT 10.*zhalf,nind)
;    ind=ind[ind2]
;ENDIF

ind=WHERE(vminusi GT 0.5 AND vminusi LT 1.9,nind)
print, "Using ", STRTRIM(nind,2), " Stars"


;imfile=galaxy+'_f814_drz.fits';Comentend out by CC
;im=READFITS(imfile,ext=1,/SILENT);Comentend out by CC
;tvgal,im,0.005,10.,outim,shrink=10;Comentend out by CC
;imsize=SIZE(im);Comentend out by CC
;xsize=imsize[1];Comentend out by CC
;ysize=imsize[2];Comentend out by CC
;plot,findgen(10),/nodata,xrange=[0,xsize],yrange=[0,ysize],xstyle=1,ystyle=1;Comentend out by CC
;imgunder,outim;Comentend out by CC
;loadct,13,/SILENT;Comentend out by CC
;oplot, x[ind],y[ind],psym=3,color=120;Comentend out by CC
;loadct,0,/SILENT;Comentend out by CC

DET_TRGB,imag[ind],ierr[ind],imin,imax,binval,lumfunc,edgedet_mendez
max=MAX(edgedet_mendez,maxpos)
mendez=binval[maxpos]
print, galaxy, " TRGB magnitude", mendez,mendez+4.08

IF NOT KEYWORD_SET(noplot) THEN BEGIN
!P.MULTI=[0,1,3]
IF (KEYWORD_SET(usevi)) THEN psfile=galaxy+'/'+galaxy+"_figtrgb_vi.ps" ELSE psfile=galaxy+'/'+galaxy+"_figtrgb.ps"
;myplot,file=psfile,ysize=5,yoff=1,/inch,xsize=5
set_plot,'x'

plot,binval,lumfunc,xtitle=xtitle,ytitle="Number Density",xthick=2,ythick=2,thick=2,/xsty,ymargin=[0,2],/ylog,yrange=[0.7,MAX(lumfunc)*1.1],/ysty,xcharsize=0.0001,charsize=2
!Y.CRANGE=10.^!Y.CRANGE
plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2

plot,binval,edgedet_mendez,thick=2,xtitle=xtitle,ytitle="ED Response",/xsty,title=titlestring,ymargin=[2,0],xcharsize=0.0001,charsize=2,yrange=[1.2*MIN(edgedet_mendez),1.2*MAX(edgedet_mendez)],/ysty
plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2
ENDIF

;now determine colors at tip and -3.5
tipstars=WHERE(imag[ind] GE mendez AND (imag[ind]-0.1) LE (mendez+0.1),ntipstars)
plothist,vmag[ind[tipstars]]-imag[ind[tipstars]],xhist,yhist,bin=0.02,/noplot
result=gaussfit(xhist,yhist,a,nterms=3)
colortrgb=a[1]
widthtrgb=a[2]
print, "TRGB color", colortrgb, widthtrgb

IF (KEYWORD_SET(usevi)) THEN mag35=0.56  ELSE mag35=0.58
min35=mag35-0.1 & max35=mag35+0.1
stars35=WHERE(imag[ind] GE mendez+min35 AND imag[ind] LE (mendez+max35),ntipstars)
plothist,vmag[ind[stars35]]-imag[ind[stars35]],xhist,yhist,bin=0.02,/noplot
result=gaussfit(xhist,yhist,a,nterms=3)
;oplot,xhist,result
color35=a[1]
width35=a[2]
print, galaxy," Color at -3.5", color35,width35

;calculate Fe/H if using VI
IF (KEYWORD_SET(usevi)) THEN BEGIN
bolocorr=0.881-0.243*colortrgb
feh=-12.64+(12.6*color35)-(3.3*color35^2)
print, galaxy," [Fe/H]",feh
mitrgb=-(0.19*feh)-3.81-bolocorr
print, "M_I,TRGB",mitrgb
ENDIF    


;now determine errors ala mendez02, do monte carlo detection with
;random resampling of stars.
IF KEYWORD_SET(doerror) THEN BEGIN
nruns=500
mendezarr=FLTARR(nruns)
smfarr=FLTARR(nruns)
print, "Commencing Error Calc"
FOR i=0,nruns-1 DO BEGIN
    print,i,FORMAT='(I,$)'
    random=LONG(RANDOMU(seed,nind)*FLOAT(nind))
    inewerr=ierr[ind[random]]
    inew=imag[ind[random]]+RANDOMN(seed,nind)*inewerr
    DET_TRGB,inew,inewerr,imin,imax,binval,lumfunc2,edgedet_mendez2
    max=MAX(edgedet_mendez2,maxpos)
    mendezarr[i]=binval[maxpos]
ENDFOR


IF (KEYWORD_SET(usevi)) THEN BEGIN
    xtitle="I [Johnson]"
ENDIF ELSE BEGIN
    xtitle="F814W [VEGAmag]"
ENDELSE
plothist,mendezarr,xhist,yhist,bin=0.01,xrange=[imin,imax],/xsty,ytitle="Number",xtitle=xtitle,ymargin=[4,-2],charsize=2
plots,[mendez,mendez],!Y.CRANGE,thick=2,linestyle=2

;now fit monte carlo results to gaussian
estimates=[FLOAT(nruns)/10.,mendez,0.1]
result=gaussfit(xhist,yhist,a,nterms=3,estimates=estimates)
mcpeak=a[1]
mcwidth=a[2]
oplot,xhist,result

print, 'MC Results',mcpeak,mcwidth
IF (KEYWORD_SET(usevi)) THEN BEGIN 
    outfile=galaxy+'/'+galaxy+"_trgb_vi.dat" 
    OPENW,1,outfile,WIDTH=1000
    printf,1,galaxy,mendez,mcpeak,mcwidth,colortrgb,widthtrgb,color35,width35,feh,mitrgb,nind,FORMAT='(A20,9F7.3,I8)'
ENDIF ELSE BEGIN
    outfile=galaxy+'/'+galaxy+"_trgb.dat"
    OPENW,1,outfile,WIDTH=1000
    printf,1,galaxy,mendez,mcpeak,mcwidth,colortrgb,widthtrgb,color35,width35,nind,FORMAT='(A20,7F7.3,I8)'
ENDELSE
CLOSE,1
FREE_LUN,1
ENDIF ;doerror

IF NOT KEYWORD_SET(noplot) THEN BEGIN
;DEVICE,/CLOSE
SET_PLOT,'X'
!P.MULTI=[0,1,1]

IF (KEYWORD_SET(usevi)) THEN psfile=galaxy+'/'+galaxy+"_trgbcmd_vi.ps" ELSE psfile=galaxy+'/'+galaxy+"_trgbcmd.ps"
;myplot,file=psfile,ysize=5,yoff=1,/inch,xsize=6

IF (KEYWORD_SET(usevi)) THEN BEGIN
    xtitle="V-I [Johnson]"
    ytitle="I [Johnson]"
ENDIF ELSE BEGIN
    xtitle="F606W-F814W [VEGAmag]"
    ytitle="F814W [VEGAmag]"
ENDELSE

plot,vminusi[ind],imag[ind],xrange=[0,2],yrange=[28,18],xstyle=1,ystyle=1,psym=4,title=titlestring,xtitle=xtitle,ytitle=ytitle,symsize=0.5
plots,[0,2],[mendez,mendez],thick=2

;DEVICE,/CLOSE
SET_PLOT,'X'

IF (KEYWORD_SET(usevi)) THEN savefile=galaxy+'/'+galaxy+"_trgb_vi.idl" ELSE savefile=galaxy+'/'+galaxy+"_trgb.idl"
SAVE,file=savefile
ENDIF

;STOP
END


PRO DET_TRGB,imag,ierr,imin,imax,binval,lumfunc,edgedet_mendez
binsize=0.01
nbin=FIX((imax-imin)/binsize)+1
binval=FINDGEN(nbin)*binsize+imin
lumfunc=FLTARR(nbin)
edgedet_mendez=FLTARR(nbin)
;edgedet_smf=FLTARR(nbin)
FOR i=0,nbin-1 DO BEGIN
    lumfunc[i]=TOTAL(1./(SQRT(2*!PI)*ierr)*exp(-(imag-binval[i])^2/(2*ierr^2)))
    nearby=WHERE(imag GT binval[i]-0.05 AND imag LT binval[i]+0.05,nnearby)
    IF (nnearby GT 2) THEN BEGIN
        sigma=MEAN(ierr[nearby])
        edgedet_mendez[i]=-(ALOG10(TOTAL(1./(SQRT(2*!PI)*ierr)*exp(-(imag+sigma-binval[i])^2/(2*ierr^2))))-ALOG10(TOTAL(1./(SQRT(2*!PI)*ierr)*exp(-(imag-sigma-binval[i])^2/(2*ierr^2)))))*SQRT(lumfunc[i])
;        edgedet_smf[i]=-TOTAL(1./(SQRT(2*!PI)*ierr)*exp(-(imag+sigma-binval[i])^2/(2*ierr^2)))+TOTAL(1./(SQRT(2*!PI)*ierr)*exp(-(imag-sigma-binval[i])^2/(2*ierr^2)))
    ENDIF ELSE BEGIN
;        edgedet_smf[i]=0.0
        edgedet_mendez[i]=0.0
    ENDELSE
ENDFOR

END
