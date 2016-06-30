function schmidtlaw, files, lengthunit, massunit, deltat, deltar=deltar,OVER=over,VERBOSE=verbose,COLOR=color,OLD=old,_EXTRA=_extra
; try to find best-fitting schmidt-law for the galaxies.


timeunit=4.122d+10 ;to convert to years
;kpcunit=50000.
kpcunit=lengthunit
;massunit=1.69875d16
readcol, files, filename, format='a'
massformfile=filename+'.massform'

sigmasfr=FINDGEN(n_elements(filename)) ;in Msol/yr/kpc^2
sigmasfr2=FINDGEN(n_elements(filename)) ;in Msol/yr
sigmagas=FINDGEN(n_elements(filename)) ;in Msol/pc^2
FOR j=0L,n_elements(filename)-1 do begin
  RDFLOAT, massformfile[j], massform, skipline=1
  rtipsy,filename[j]+'.bin',h,g,d,s
  rform=SQRT(s.x^2+s.y^2)*kpcunit
  tform=s.tform*timeunit
  rgas=SQRT(g.x^2+g.y^2)*kpcunit
  g.mass=g.mass*massunit
  tcurrent=13.73d9

IF (KEYWORD_SET(deltat) EQ 0) THEN deltat=100.e6 ; in yrs
;IF (KEYWORD_SET(deltar) EQ 0) THEN deltar=.5 ;in kpc

tclip=tcurrent-deltat ; only use starformation events since tclip

rmin=0.0
rmax = [max(g.x)*kpcunit, 15.]
rmax = min(rmax)
;rmax=10.0 ;max radius in kpc
;nbin=FIX(((rmax-rmin)/deltar)+1)

;rarray=FINDGEN(nbin+1)*deltar+rmin
;rmid=FINDGEN(nbin)*deltar+rmin+deltar/2.0

;FOR i=0,nbin-1 DO BEGIN
    areakpc=!PI*rmax^2
    areapc=areakpc*1.e6

    indform=WHERE(rform GT rmin AND rform LT rmax AND tform GT tclip,nindform)
    indgas=WHERE(rgas GT rmin AND rgas LT rmax,nindgas)

    IF (nindform GT 0) THEN sigmasfr[j]=TOTAL(massform[indform])/(areakpc*deltat) ELSE sigmasfr[j]=0
    IF (nindform GT 0) THEN sigmasfr2[j]=TOTAL(massform[indform])/(deltat) ELSE sigmasfr[j]=0
    IF (nindgas GT 0) THEN sigmagas[j]=TOTAL(g[indgas].mass)/areapc ELSE sigmagas[j]=0
;ENDFOR
ENDFOR

IF (keyword_set(over)) THEN oplot,sigmagas,sigmasfr,psym=4,color=color,_EXTRA=_extra ELSE $
  plot, sigmagas,sigmasfr,/xlog,/ylog,psym=1,xrange=[0.01,MAX(sigmagas)*2.],yrange=[1e-6,MAX(sigmasfr)*2],xtitle=textoidl(' \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl(' \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]",xstyle=1,ystyle=1,_EXTRA=_extra, pos=[.15,.15,.95,.95]

xsigma=FINDGEN(10000)/10.
ysigma=2.5e-4*xsigma^1.4
oplot,xsigma,ysigma

return, sigmasfr2

END
