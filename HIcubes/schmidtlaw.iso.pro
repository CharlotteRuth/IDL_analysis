PRO schmidtlaw, filename, binary=binary, deltat=deltat, deltar=deltar,OVER=over,VERBOSE=verbose,COLOR=color,OLD=old,_EXTRA=_extra
; try to find best-fitting schmidt-law for the galaxies.
; for a non-cosmological simulation (isolated galaxy)

set_plot, 'ps'
device, filename='kstest.ps', /color

loadct, 13

kpcunit=1.0
timeunit=1.e9 ;in yr
massunit=2.3262e5
massformfile=filename+'.massform'

FOR j=0L,n_elements(filename)-1 do begin
  IF keyword_set(binary) then begin
  readarr, massformfile[j], h, massform
  massform = massform[h.ngas+h.ndark:h.n-1]
  ENDIF ELSE rdfloat, massformfile[j], massform, skipline=1
  rtipsy,filename[j],h,g,d,s
  ;timeunit=13.7d9/max(s.tform) ;to convert to years
  rform=SQRT(s.x^2+s.y^2+s.z^2)*kpcunit
  tform=s.tform*timeunit
  rgas=SQRT(g.x^2+g.y^2+g.z^2)*kpcunit
  g.mass=g.mass*massunit
  ;z = (1/h.time)-1.
  ;tcurrent=13.7e9-wmap3_lookback(z)
  tcurrent=max(tform)

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

    IF (nindform GT 0) THEN sigmasfr[i]=TOTAL(massunit*massform[indform])/(areakpc*deltat) ELSE sigmasfr[i]=0
    IF (nindgas GT 0) THEN sigmagas[i]=TOTAL(g[indgas].mass)/areapc ELSE sigmagas[i]=0
    IF (keyword_set(verbose)) THEN print, rarray[i],areakpc,sigmasfr[i],sigmagas[i]

ENDFOR

print, sigmasfr
print, sigmagas

xtitle=textoidl(' \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"
ytitle=textoidl(' \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]"
;xtitle=textoidl(' \Sigma')+"!lgas!n [M"+textoidl('_o')+" pc!u-2!n]"
;ytitle=textoidl(' \Sigma')+"!lSFR!n [M"+textoidl('_o')+"kpc!u-2!n yr!u-1!n]"
;xrange=[0.8,MAX(sigmagas)*2.]

plotsym, 0, 1.5, /fill
if j ne 2 then psym=j+1 else psym=8
IF j eq 0 THEN begin
  plot, sigmagas,sigmasfr,/xlog,/ylog,psym=6,xrange=[0.8,20.], yrange=[1e-6,2.*MAX(sigmasfr)*2],xtitle=xtitle, ytitle=ytitle, xstyle=1,ystyle=1,_EXTRA=_extra, pos=[.15,.15,.95,.95], title=filename[j]
  oplot, sigmagas,sigmasfr, linestyle=0, color=240
ENDIF
IF j ne 0 then begin
  oplot,sigmagas,sigmasfr,psym=psym,color=color,_EXTRA=_extra 
  oplot, sigmagas,sigmasfr, linestyle=0, color=240/j
ENDIF
;  plot, sigmagas,sigmasfr,/xlog,/ylog,psym=1,xrange=[1e-1,5e5],yrange=[1e-4,1e4],xtitle="!4R!x!lgas!n [M!lsol!n pc!u-2!n]",ytitle="!4R!x!lSFR!n [M!lsolar!nkpc!u-2!n yr!u-1!n]",xstyle=1,ystyle=1

xsigma=FINDGEN(10000)/10.
ysigma=2.5e-4*xsigma^1.4
oplot,xsigma,ysigma

ENDFOR
device, /close

END
