@locus.pro
@median_fit.pro
PRO FIT_LOCUS



; READCOLMORE, 'star8band_000752.dat', ra,dec,primary,a_r,run,rerun,camcol,field,id,saturflag,blendflag,childflag,crflag,interflag,edgeflag,doublemvflag,psfwidth,umag,uerr,gmag,gerr,rmag,rerr,imag,ierr,zmag,zerr,dist_2mass,jmag,jerr,hmag,herr,kmag,kerr,rdflag,blflag,ccflag,FORMAT='D,D,I,F,X,X,I,I,I,I,I,X,X,X,X,I,I,I,I,I,I,I,X,X,X,X,X,F,X,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,X,A,A,A,X,X,X'

restore, file='goodstars.idl'
;ind=WHERE(umag GT 14.0 AND umag LT 20.5 AND saturflag EQ 0 AND
;blendflag EQ 0 AND doublemvflag EQ 0 AND kmag LT 14.3 AND rdflag EQ
;'222' AND blflag EQ '111' AND ccflag EQ '000' AND rmag LT 98.0 AND
;umag LT 98.0 AND gmag LT 98.0 AND imag LT 98.0 AND zmag LT 98.0 AND
;a_r LT 0.1,nind)
Kcor=kmag - a_r*0.133   ; 0.133 is the k-band coefficient
ind=WHERE(Kcor LT 14.0 ,nind)
print, "Using ",nind," data points with A_r less than 0.05"

ur=umag[ind]-rmag[ind]
urerr=SQRT(uerr[ind]^2+rerr[ind]^2)
gr=gmag[ind]-rmag[ind]
grerr=SQRT(gerr[ind]^2+rerr[ind]^2)
ri=rmag[ind]-imag[ind]
rierr=SQRT(rerr[ind]^2+ierr[ind]^2)
rz=rmag[ind]-zmag[ind]
rzerr=SQRT(rerr[ind]^2+zerr[ind]^2)
rj=rmag[ind]-jmag[ind]
rjerr=SQRT(rerr[ind]^2+jerr[ind]^2)
rh=rmag[ind]-hmag[ind]
rherr=SQRT(rerr[ind]^2+herr[ind]^2)
rk=rmag[ind]-kmag[ind]
rkerr=SQRT(rerr[ind]^2+kerr[ind]^2)


;resultall=LINFIT(rk[use],ur[use],measure_errors=urerr[use],SIGMA=sigma,YFIT=yfit)
OPENW, 1, "pcfits.dat",WIDTH=1000
plotter=0
if plotter eq 1 then begin
 SET_PLOT, 'PS'
 DEVICE,file="pcfits_newest_1.ps",/color,bits_per_pixel=8,ysize=25,yoffset=2,encapsulated=0;
;,ysize=9,yoff=1,/inch
 ;myplot,file="pcfits2.ps",ysize=9,yoff=1,/inch
endif else set_plot,'x'
!P.MULTI=[0,2,3]
!p.charsize=2
!p.charthick=2
;these are extinction values (in terms of A_r) for each band
;OLD
;extinction=[1.87, 1.38, 1.0, 0.76, 0.54, 0.327, 0.209, 0.133]
;NEW
extinction=[1.84, 1.40, 1.0, 0.78, 0.57, 0.311, 0.19, 0.133]
urext=extinction[0]-extinction[2]
grext=extinction[1]-extinction[2]
riext=extinction[2]-extinction[3]
rzext=extinction[2]-extinction[4]
rhext=extinction[2]-extinction[5]
rjext=extinction[2]-extinction[6]
rkext=extinction[2]-extinction[7]


;fit each stellar locus between r-K of 1.2 and 2.7
;UR
median_fit,rk,ur,1.2,2.7,0.1,coeffs,errors,bincenter,yfit
a=locus(rk,ur,0.5,5,0.5,5,0.05,'r-K','u-r',1.0)
oplot, bincenter,yfit,thick=3
blah=n_elements(yfit)
arrow,bincenter[0],yfit[0],bincenter[0]+2*rkext,yfit[0]+2*urext,/data,thick=3,hsize=-.15,color=50;Ar arrow
arrow,bincenter[0],yfit[0],3.2,4.62575,thick=3,hsize=-.1,/data ;PC1 arrow
arrow,bincenter[0],yfit[0],0.7,1.840537,/data,thick=3,hsize=-.3;,color=100 ; PC2 arrow
plots,[bincenter[0],bincenter[0]+2*rkext],[yfit[0],yfit[0]+2*urext],linestyle=2,thick=3;,color=255

xyouts,3.3,4.7,'PC!i1!n',charsize=1
xyouts,3.2,2.9,'A!ir!n=1',charsize=2
xyouts,0.6,2.1,'PC!i2!n',charsize=1

;legend,["PC!i1!n","PC!i2!n","A!ir!n=1"],lines=[0,0,2],colors=[0,100,0],/right,/bottom,charsize=0.8

;stop
angle=ATAN(coeffs[1])
;print, angle*180./!PI
p1=rk*cos(angle)+ur*sin(angle)
p2=-rk*sin(angle)+ur*cos(angle)
modep1=MODE(p1)
modep2=MODE(p2)
p1=p1-modep1[0]
p2=p2-modep2[0]
printf,1, "u-r",angle,modep1[0],modep2[0]

;GR
median_fit,rk,gr,1.2,2.7,0.1,coeffs,errors,bincenter,yfit
a=locus(rk,gr,0.5,5,0,2.0,0.05,'r-K','g-r',1.0)
oplot, bincenter,yfit,thick=3
arrow,bincenter[0],yfit[0],bincenter[0]+2*rkext,yfit[0]+2*grext,/data,thick=3,hsize=-.15;,color=255
plots,[bincenter[0],bincenter[0]+2*rkext],[yfit[0],yfit[0]+2*grext],linestyle=2,thick=3,color=255
angle=ATAN(coeffs[1])
;print, angle*180./!PI
p1=rk*cos(angle)+gr*sin(angle)
p2=-rk*sin(angle)+gr*cos(angle)
modep1=MODE(p1)
modep2=MODE(p2)
p1=p1-modep1[0]
p2=p2-modep2[0]
printf,1, "g-r",angle,modep1[0],modep2[0]

;RI
median_fit,rk,ri,1.2,2.7,0.1,coeffs,errors,bincenter,yfit
a=locus(rk,ri,0.5,5,-0.5,2,0.05,'r-K','r-i',1.0)
oplot, bincenter,yfit,thick=3
arrow,bincenter[0],yfit[0],bincenter[0]+2*rkext,yfit[0]+2*riext,/data,thick=3,hsize=-.15;,color=255
plots,[bincenter[0],bincenter[0]+2*rkext],[yfit[0],yfit[0]+2*riext],linestyle=2,thick=3,color=255
angle=ATAN(coeffs[1])
;print, angle*180./!PI
p1=rk*cos(angle)+ri*sin(angle)
p2=-rk*sin(angle)+ri*cos(angle)
modep1=MODE(p1)
modep2=MODE(p2)
p1=p1-modep1[0]
p2=p2-modep2[0]
printf,1, "r-i",angle,modep1[0],modep2[0]

;RZ
median_fit,rk,rz,1.2,2.7,0.1,coeffs,errors,bincenter,yfit
a=locus(rk,rz,0.5,5,-0.5,2.5,0.05,'r-K','r-z',1.0)
oplot, bincenter,yfit,thick=3
arrow,bincenter[0],yfit[0],bincenter[0]+2*rkext,yfit[0]+2*rzext,/data,thick=3,hsize=-.15;,color=255
plots,[bincenter[0],bincenter[0]+2*rkext],[yfit[0],yfit[0]+2*rzext],linestyle=2,thick=3,color=255
angle=ATAN(coeffs[1])
;print, angle*180./!PI
p1=rk*cos(angle)+rz*sin(angle)
p2=-rk*sin(angle)+rz*cos(angle)
modep1=MODE(p1)
modep2=MODE(p2)
p1=p1-modep1[0]
p2=p2-modep2[0]
printf,1, "r-z",angle,modep1[0],modep2[0]

;RJ
median_fit,rk,rj,1.2,2.7,0.1,coeffs,errors,bincenter,yfit
a=locus(rk,rj,0.5,5,0,4.5,0.05,'r-K','r-J',1.0)
oplot, bincenter,yfit,thick=3
arrow,bincenter[0],yfit[0],bincenter[0]+2*rkext,yfit[0]+2*rjext,/data,thick=3,hsize=-.15;,color=255
plots,[bincenter[0],bincenter[0]+2*rkext],[yfit[0],yfit[0]+2*rjext],linestyle=2,thick=3,color=255
angle=ATAN(coeffs[1])
;print, angle*180./!PI
p1=rk*cos(angle)+rj*sin(angle)
p2=-rk*sin(angle)+rj*cos(angle)
modep1=MODE(p1)
modep2=MODE(p2)
p1=p1-modep1[0]
p2=p2-modep2[0]
printf,1, "r-J",angle,modep1[0],modep2[0]

;RH
median_fit,rk,rh,1.2,2.7,0.1,coeffs,errors,bincenter,yfit
a=locus(rk,rh,0.5,5,0,4.5,0.05,'r-K','r-H',1.0)
oplot, bincenter,yfit,thick=3
arrow,bincenter[0],yfit[0],bincenter[0]+2*rkext,yfit[0]+2*rhext,/data,thick=3,hsize=-.15;,color=255
plots,[bincenter[0],bincenter[0]+2*rkext],[yfit[0],yfit[0]+2*rhext],linestyle=2,thick=3,color=255
angle=ATAN(coeffs[1])
;print, angle*180./!PI
p1=rk*cos(angle)+rh*sin(angle)
p2=-rk*sin(angle)+rh*cos(angle)
modep1=MODE(p1)
modep2=MODE(p2)
p1=p1-modep1[0]
p2=p2-modep2[0]
printf,1, "r-H",angle,modep1[0],modep2[0]



if plotter eq 1 then DEVICE,/CLOSE
SET_PLOT,'x'
CLOSE,1
FREE_LUN,1




; p1red=0.32*(rkext*cos(angle)+urext*sin(angle))

; higha=WHERE(a_r[ind] GT 0.3,nhigha)
; lowa=WHERE(a_r[ind] LT 0.08,nlowa)

; SET_PLOT, 'PS'
; DEVICE, file="first_red.ps"
; plothist,p1[higha],bin=0.05,peak=1.0,yrange=[0,1.1],thick=1,xthick=2,ythick=2,title='High vs. Low Reddening data, (r-K,u-r) data',xtitle='P1 (r-K,u-r)',ytitle='Normalized Number'
; plothist,p1[lowa],bin=0.05,/over,peak=1.0,thick=2
; legend,['Ar < 0.08','Ar > 0.3'],linestyle=[0,0],thick=[2,1],/right

; plots,[-.2,(-.2+p1red)],[1.05,1.05],thick=3
; DEVICE, /CLOSE
; SET_PLOT, 'X'



END
