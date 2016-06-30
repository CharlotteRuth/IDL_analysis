pro gasprof,file,nbin=nbin,all=all,kpcunit=kpcunit,massunit=m,overplot=overplot,showfit=showfit,_extra=_extra

  rtipsy,file,h,g,d,s
  if (keyword_set(all)) then begin
    if (h.nstar GT 0) then rsqr = [(g.x^2. + g.y^2. + g.z^2.), $
            (d.x^2. + d.y^2. + d.z^2.), $
            (s.x^2. + s.y^2. + s.z^2.)] $
    ELSE rsqr = [(g.x^2. + g.y^2. + g.z^2.), $
            (d.x^2. + d.y^2. + d.z^2.)]
    r = sqrt(rsqr)
  endif else begin
    rsqr = (g.x^2. + g.y^2. + g.z^2.)
    r = sqrt(rsqr)
  endelse

  if (keyword_set(nbin) eq 0) then nbin=40
  if (keyword_set(kpcunit) eq 0) then kpcunit = 1.d0
  if (keyword_set(massunit) eq 0) then massunit = 2.325d5
  mconv = double(massunit * 1.20d57)
  vconv = (kpcunit * 3.086d21)^3.
  Rmin = min(R)
  Rmax = max(R)
  rho = fltarr(nbin)
  ;deltar = (Rmax-Rmin)/nbin
  logdeltar = (alog10(Rmax)-alog10(Rmin))/nbin
  binsize = 10.^(findgen(nbin)*logdeltar + alog10(Rmin)) - $
    10.^((findgen(nbin)-1)*logdeltar+alog10(Rmin))
  ;Rarray = findgen(nbin)*deltar+Rmin
  Rarray = 10.^(findgen(nbin)*logdeltar+alog10(Rmin))
  ngd = h.ngas+h.ndark
  FOR i=0,nbin-2 DO BEGIN
    volume=vconv*4./3.*!PI*((Rarray[i+1])^3-Rarray[i]^3)

    bini=WHERE(R GT Rarray[i] AND R LT Rarray[i+1],ninbin)

   IF (ninbin GT 0) then begin
    IF (keyword_set(all)) then begin
      gi = where(bini LT h.ngas,ngi)
      if (ngi GT 0) THEN binmass = total(g[bini[gi]].mass) ELSE binmass = 0
      di = where(bini GT h.ngas AND bini LT ngd,ndi)
      if (ndi GT 0) THEN binmass = binmass+total(d[bini[di] - h.ngas].mass)
      si = where(bini LT h.n AND bini GT ngd,nsi)
      if (nsi GT 0) THEN binmass = binmass+total(s[bini[si]-ngd].mass)
    ENDIF ELSE binmass = total(g[bini].mass)
    rho[i]=TOTAL(mconv*binmass)/(volume)
   ENDIF ELSE rho[i]=0
  ENDFOR

  minrho = min(rho[where(rho GT 0)])
 if (keyword_set(overplot)) then oplot,Rarray,rho,_extra=_extra $
 else plot,Rarray,rho,/ylog,yrange=[minrho,1.1*max(rho)],ystyle=1, $
     xrange = [Rmin,Rmax], /xlog, $
     ;xtit = 'r [kpc]',ytit='Density [M'+sunsymbol()+'/kpc!u3!n]',_extra=_extra
     xtit = 'r [kpc]',ytit='Density [cm!u-3!n]',_extra=_extra

 if (keyword_set(showfit)) then begin
   fit = fitprof(rho,Rarray)
 endif
end
