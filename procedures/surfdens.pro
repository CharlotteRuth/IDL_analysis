PRO surfdens, tipsyfile, deltar=deltar,OVER=over,VERBOSE=verbose,_EXTRA=_extra
; plot annularized surface density of a disk


timeunit=4.0435875d+10 ;to convert to years
kpcunit=1.d5
massunit=13.6d16

rtipsy,tipsyfile,h,g,d,s

rgas=SQRT(g.x^2+g.y^2)*kpcunit

IF (KEYWORD_SET(deltar) EQ 0) THEN deltar=.5 ;in kpc

rmin=0.0
rmax=20.0 ;max radius in kpc
nbin=FIX(((rmax-rmin)/deltar)+1)

rarray=FINDGEN(nbin+1)*deltar+rmin
rmid=FINDGEN(nbin)*deltar+rmin+deltar/2.0
sigmagas=FINDGEN(nbin) ;in Msol/pc^2

FOR i=0,nbin-1 DO BEGIN
    areakpc=!PI*rarray[i+1]^2-!PI*rarray[i]^2
    areapc=areakpc*1.e6

    indgas=WHERE(rgas GT rarray[i] AND rgas LT rarray[i+1],nindgas)

    IF (nindgas GT 0) THEN sigmagas[i]=TOTAL(g[indgas].mass*massunit)/areapc ELSE sigmagas[i]=0
    IF (keyword_set(verbose)) THEN print, rarray[i],areakpc,sigmagas[i]

ENDFOR

IF (keyword_set(over)) THEN oplot,rmid,sigmagas,_EXTRA=_extra ELSE $
  plot, rmid,sigmagas,/ylog,yrange=[1e-1,MAX(sigmagas)*2.],ytitle="!4R!x!lgas!n [M!lsolar!n/pc!u2!n]",xtitle="R [kpc]",xstyle=1,ystyle=1,_EXTRA=_extra
;  plot, sigmagas,sigmasfr,/xlog,/ylog,psym=1,xrange=[1e-1,5e5],yrange=[1e-4,1e4],xtitle="!4R!x!lgas!n [M!lsolar!n/pc!u2!n]",ytitle="!4R!x!lSFR!n [M!lsolar!n/kpc!u2!n/yr]",xstyle=1,ystyle=1

END
