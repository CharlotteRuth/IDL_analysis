function starprof,file,nbin=nbin,kpcunit=kpcunit,massunit=m,overplot=overplot,yz=yz,xz=xz,logbin=logbin,linbin=linbin,_extra=_extra

  rtipsy,file,h,g,d,s
;Use minimum potential as center
  xcen=s[where(s.phi EQ min(s.phi))].x
  ycen=s[where(s.phi EQ min(s.phi))].y
  zcen=s[where(s.phi EQ min(s.phi))].z
  s.x=s.x-xcen
  s.y=s.y-ycen
  s.z=s.z-zcen

  if (keyword_set(xz)) then Rsqr = (s.x^2. + s.z^2.) $
  else if (keyword_set(yz)) then Rsqr = (s.y^2. + s.z^2.) $
  else Rsqr = (s.x^2. + s.y^2.)
  R = sqrt(Rsqr)

  if (keyword_set(nbin) eq 0) then nbin=40
  if (keyword_set(kpcunit) eq 0) then kpcunit = 1.
  if (keyword_set(massunit) eq 0) then massunit = 2.325e5
; Can't have log scale if Rmin = 0
  Rmin = min(R[where(R GT 0)])
  Rmax = max(R)
; To create nbins, we need r's on both outer borders
  nrs = nbin+1
  Rarray = fltarr(nrs)
  midRarray = fltarr(nbin)
  sigma = fltarr(nbin)
  if (keyword_set(logbin) ) then begin
    logdeltar = (alog10(Rmax)-alog10(Rmin))/nrs
    binsize = 10.^(findgen(nrs)*logdeltar + alog10(Rmin)) - $ 
      10.^((findgen(nrs)-1)*logdeltar+alog10(Rmin))
    Rarray = 10.^(findgen(nrs)*logdeltar+alog10(Rmin))
  endif else begin
    if (keyword_set(linbin)) then begin
      deltar = (Rmax-Rmin)/nrs
      Rarray = findgen(nrs)*deltar+Rmin
    endif else begin
; Same number of stars per bin
; here's how many stars, make sure it's a float so we get out to the edge
      n = float(h.nstar)/nrs
      sri = sort(R)
      sortR = R[sri]
      for j=1,nrs do Rarray[j-1] = sortR[j*n-1]
    endelse
  endelse
  FOR i=0,nbin-1 DO BEGIN
    areakpc=!PI*(Rarray[i+1])^2-!PI*Rarray[i]^2
    areapc=areakpc*1.e6

    indform=WHERE(R GT Rarray[i] AND R LT Rarray[i+1],nindform)

    midRarray[i] = (Rarray[i]+Rarray[i+1])/2.
    IF (nindform GT 0) THEN sigma[i]=TOTAL(s[indform].mass)/(areakpc) $
     ELSE sigma[i]=0
  ENDFOR

  return,{r:midRarray,sigma:sigma}
END

