FUNCTION rhoNFW,r
; from pkdNFWSpheroid in pkd.c
G = 1
M200 = 5.379e-6
R200 = 0.00153
c = 20.648
rc = R200/c
rho200 = M200/(4*!PI*(R200^3)/3)
rho0 = rho200*c*(1+c)*(1+c)

return,rho0/(r/rc*(1+r/rc)^2.0)
END

PRO plottoomreq, tipsyfile, deltar=deltar,OVER=over,VERBOSE=verbose,COLOR=color,_EXTRA=_extra
; plot Toomre Q parameter vs. radius

KpcUnit = 1.e5
MsolUnit = 13.6e16
KPCCM=3.09e21
KBOLTZ=1.38e-16
GCGS=6.67e-8
MSOLG=2e33
MU=0.58
EOSGAMMA=5/3
MHYD=1.67e-24
CsFactor=EOSGAMMA*KpcUnit*KPCCM*KBOLTZ/GCGS/MsolUnit/MSOLG/MU/MHYD

rtipsy,tipsyfile,h,g,d,s
rgas=SQRT(g.x^2+g.y^2)
rstar=SQRT(s.x^2+s.y^2)

IF (KEYWORD_SET(deltar) EQ 0) THEN deltar=5e-6 ;in sysunits

rmin=0.0
rmax=0.00015 ;max radius in kpc
nbin=FIX(((rmax-rmin)/deltar)+1)

rarray=FINDGEN(nbin+1)*deltar+rmin
rmid=FINDGEN(nbin)*deltar+rmin+deltar/2.0
sigmasfr=FINDGEN(nbin) ;in Msol/yr/kpc^2
omega2=FLTARR(nbin+1)
Q=FINDGEN(nbin) ;in Msol/pc^2

FOR i=0,nbin-1 DO BEGIN
    area=!PI*rarray[i+1]^2-!PI*rarray[i]^2

    indstar=WHERE(rstar GT rarray[i] AND rstar LT rarray[i+1],nindstar)
    indgas=WHERE(rgas GT rarray[i] AND rgas LT rarray[i+1],nindgas)

    IF (nindstar GT 0) THEN BEGIN
      ;Lzstar = fltarr(nindstar) 
      Lzstar = s[indstar].x*s[indstar].vy - s[indstar].y*s[indstar].vx 
    ENDIF ELSE Lzstar = 0.0
    IF (nindgas GT 0) THEN BEGIN
      sigmagas=(TOTAL(g[indgas].mass)+total(s[indstar].mass))/area 
      soundspeed=MEAN(SQRT(CsFactor*g[indgas].tempg))
      Lz=MEAN([g[indgas].x*g[indgas].vy - g[indgas].y*g[indgas].vx,Lzstar])
      kappa=MEAN(SQRT(4*!PI*rhoNFW(rmid[i])+3*Lz*Lz/(rmid[i]^4.0)))
      omega2[i]=(Lz/rmid[i]/rmid[i])^2
      IF (i GT 0) THEN BEGIN 
        kappa=SQRT(rmid[i]*(omega2[i]-omega2[i-1])/(rmid[i]-rmid[i-1])+4*omega2[i])
        Q[i]=soundspeed*kappa/!PI/sigmagas
      ENDIF
    ENDIF ELSE BEGIN 
      Q[i]=0
    ENDELSE
    IF (keyword_set(verbose)) THEN print, rarray[i],areakpc,sigmasfr[i],sigmagas[i]

ENDFOR

IF (keyword_set(over)) THEN oplot,rmid*1e5,Q,color=color,_EXTRA=_extra ELSE $
  plot, rmid*1e5,Q,xtitle="R (kpc)",ytitle="Toomre Q",_EXTRA=_extra
;  plot, sigmagas,sigmasfr,/xlog,/ylog,psym=1,xrange=[1e-1,5e5],yrange=[1e-4,1e4],xtitle="!4R!x!lgas!n [M!lsolar!n/pc!u2!n]",ytitle="!4R!x!lSFR!n [M!lsolar!n/kpc!u2!n/yr]",xstyle=1,ystyle=1

END
