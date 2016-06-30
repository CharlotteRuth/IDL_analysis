function beta,s,nbin=nbin,kpcunit=kpcunit,massunit=m,yz=yz,xz=xz,logbin=logbin,linbin=linbin,_extra=_extra

  if (keyword_set(nbin) eq 0) then nbin=40
  if (keyword_set(kpcunit) eq 0) then kpcunit = 1.
  if (keyword_set(massunit) eq 0) then massunit = 2.325e5

; Calculate radial and tangential velocities
R = (s.x*s.x + s.y*s.y + s.z*s.z)
vrsq = (s.vx*s.x + s.vy*s.y +s.vz*s.z)/R
Larr =crossp([s[i].x,s[i].y,s[i].z],[s[i].vx,s[i].vy,s[i].vz])
;vthsq=transpose(Larr)#Larr
vthsq=Larr[2]*Larr[2]/R

; Can't have log scale if Rmin = 0
  Rmin = min(R[where(R GT 0)])
  Rmax = max(R)
; To create nbins, we need r's on both outer borders
  nrs = nbin+1
  Rarray = fltarr(nrs)
  midRarray = fltarr(nbin)
  anisotropy = fltarr(nbin)
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
      n = float(n_elements(mass))/nrs
      sortR = R[sort(R)]
      for j=1,nrs do Rarray[j-1] = sortR[j*n-1]
    endelse
  endelse

  FOR i=0,nbin-1 DO BEGIN
    indform=WHERE(R GT Rarray[i] AND R LT Rarray[i+1],nindform)

    midRarray[i] = (Rarray[i]+Rarray[i+1])/2.
    IF (nindform GT 0) THEN $
       anisotropy[i]=1-mean(vthsq[indform])/(mean(vrsq[indform])) $
     ELSE sigma[i]=0
  ENDFOR

  return,{r:midRarray,b:anisotropy}
END

