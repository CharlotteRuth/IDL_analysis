PRO schmidtlaw, tipsyfile, deltat, deltar=deltar,OVER=over,VERBOSE=verbose,COLOR=color,OLD=old,_EXTRA=_extra
; try to find best-fitting schmidt-law for the galaxies.


timeunit=4.0435875d+10 ;to convert to years
kpcunit=1.d5
massunit=13.6d16

rtipsy,tipsyfile,h,g,d,s
if(keyword_set(old)) then r=read_old_rform(tipsyfile) else $
r=read_rform(tipsyfile)
massformfile=tipsyfile+'.massform'
RDFLOAT, massformfile, massform, skipline=1
massform=massform[(h.ngas+h.ndark):*]*massunit

rform=SQRT(r.x^2+r.y^2)*kpcunit
tform=s.tform*timeunit
rgas=SQRT(g.x^2+g.y^2)*kpcunit
g.mass=g.mass*massunit
tcurrent=h.time*timeunit

IF (KEYWORD_SET(deltat) EQ 0) THEN deltat=100.e6 ; in yrs
IF (KEYWORD_SET(deltar) EQ 0) THEN deltar=.5 ;in kpc

tclip=tcurrent-deltat ; only use starformation events since tclip

rmin=0.0
rmax=15.0 ;max radius in kpc
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

    IF (nindform GT 0) THEN sigmasfr[i]=TOTAL(massform[indform])/(areakpc*deltat) ELSE sigmasfr[i]=0
    IF (nindgas GT 0) THEN sigmagas[i]=TOTAL(g[indgas].mass)/areapc ELSE sigmagas[i]=0
    IF (keyword_set(verbose)) THEN print, rarray[i],areakpc,sigmasfr[i],sigmagas[i]

ENDFOR

IF (keyword_set(over)) THEN oplot,sigmagas,sigmasfr,psym=4,color=color,_EXTRA=_extra ELSE $
  plot, sigmagas,sigmasfr,/xlog,/ylog,psym=1,xrange=[0.8,MAX(sigmagas)*2.],yrange=[1e-6,MAX(sigmasfr)*2],xtitle=textoidl(' \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl(' \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]",xstyle=1,ystyle=1,_EXTRA=_extra
;  plot, sigmagas,sigmasfr,/xlog,/ylog,psym=1,xrange=[1e-1,5e5],yrange=[1e-4,1e4],xtitle="!4R!x!lgas!n [M!lsol!n pc!u-2!n]",ytitle="!4R!x!lSFR!n [M!lsolar!nkpc!u-2!n yr!u-1!n]",xstyle=1,ystyle=1

xsigma=FINDGEN(10000)/10.
ysigma=2.5e-4*xsigma^1.4
oplot,xsigma,ysigma
END
