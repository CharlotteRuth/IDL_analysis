PRO cosmoschmidt, tipsyfile, deltat, deltar=deltar,OVER=over,VERBOSE=verbose,COLOR=color
; try to find best-fitting schmidt-law for the galaxies.


timeunit=4.04313d+10 ;to convert to years
kpcunit=2.85714d4
massunit=3.171d15
dInitStarMass=9.53844e-11

rtipsy,tipsyfile,h,g,d,s
;r=read_rform(tipsyfile)
massformfile=tipsyfile+'.massform'
;RDFLOAT, massformfile, massform, skipline=1
;massform=fltarr(n_elements(s))+dInitStarMass*massunit
massform=dInitStarMass*massunit

print,h.time
rform=SQRT(s.x^2+s.y^2)*kpcunit
tform=s.tform*timeunit
print,max(tform)
rgas=SQRT(g.x^2+g.y^2)*kpcunit*h.time
g.mass=g.mass*massunit
;tcurrent=h.time*timeunit
tcurrent=1.346773e+10-lookback((1./h.time) - 1.)

IF (KEYWORD_SET(deltat) EQ 0) THEN deltat=100.e6 ; in yrs
IF (KEYWORD_SET(deltar) EQ 0) THEN deltar=.2 ;in kpc

tclip=tcurrent-deltat ; only use starformation events since tclip
print,tclip

rmin=0.0
rmax=5.0 ;max radius in kpc
nbin=FIX(((rmax-rmin)/deltar)+1)

rarray=FINDGEN(nbin+1)*deltar+rmin
rmid=FINDGEN(nbin)*deltar+rmin+deltar/2.0
sigmasfr=FINDGEN(nbin) ;in Msol/yr/kpc^2
sigmagas=FINDGEN(nbin) ;in Msol/pc^2

FOR i=0,nbin-1 DO BEGIN
    areakpc=!PI*rarray[i+1]^2-!PI*rarray[i]^2
    areapc=areakpc*1.e6

    indform=WHERE(rform GT rarray[i] AND rform LT rarray[i+1] AND tform GT tclip,nindform)
    indgas=WHERE(rgas GT rarray[i] AND rgas LT rarray[i+1],nindgas)

    IF (nindform GT 0) THEN sigmasfr[i]=n_elements(indform)*(massform)/(areakpc*deltat) ELSE sigmasfr[i]=0
    IF (nindgas GT 0) THEN sigmagas[i]=TOTAL(g[indgas].mass)/areapc ELSE sigmagas[i]=0
    IF (keyword_set(verbose)) THEN print, rarray[i],areakpc,sigmasfr[i],sigmagas[i]

ENDFOR

print,sigmasfr

IF (keyword_set(over)) THEN oplot,sigmagas,sigmasfr,psym=4,color=color ELSE $
  plot, sigmagas,sigmasfr,/xlog,/ylog,psym=1,xrange=[MIN(sigmagas)/2,MAX(sigmagas)*2.],yrange=[MIN(sigmasfr)/2,MAX(sigmasfr)*2],xtitle="!4R!x!lgas!n [M!lsolar!n/pc!u2!n]",ytitle="!4R!x!lSFR!n [M!lsolar!n/kpc!u2!n/yr]",xstyle=1,ystyle=1
;  plot, sigmagas,sigmasfr,/xlog,/ylog,psym=1,xrange=[1e-1,5e5],yrange=[1e-4,1e4],xtitle="!4R!x!lgas!n [M!lsolar!n/pc!u2!n]",ytitle="!4R!x!lSFR!n [M!lsolar!n/kpc!u2!n/yr]",xstyle=1,ystyle=1

xsigma=FINDGEN(10000)/10.
ysigma=2.5e-4*xsigma^1.4
oplot,xsigma,ysigma
;STOP
END
