;
;
;
;   return a mapping of energy to angular momentum of a circular orbit
;   for a given particle distribution
;
;   We assume that the system has already been centered and aligned
;   and return an array of Energy values and Jc to be mapped against
;   each other 



pro getlcirc, rotfile,s, E, Jc, nbins=nbins,rmax=rmax,zmax=zmax,dist_units,v_units
; *****************************
; generate a potential profile
; *****************************

if not keyword_set(nbins) then nbins = 200.
if not keyword_set(rmax) then rmax = 70.
if not keyword_set(zmax) then zmax = 1


ind = where(abs(s.z) lt zmax)
stemp = s[ind]
r = sqrt(stemp.x^2+stemp.y^2)/dist_units
rbins = equalnbins(stemp,r,2000.)
nbins = n_elements(rbins)-1
rcentbins = fltarr(nbins)
pot = fltarr(nbins)
E = fltarr(nbins)
Jc = fltarr(nbins)


for i = 0L, nbins-1 do begin

    ind = where(r gt rbins[i] and r lt rbins[i+1], n)

    rcentbins[i] = (rbins[i]+rbins[i+1])/2.

    if n gt 0 then $
      pot[i] = total(stemp[ind].phi/v_units/v_units)/n $
    else pot[i] = 0.

endfor

; *****************************************************************
; what is the angular momentum of a circular orbit at each radius?
; Need to use the rotation curve
; *****************************************************************

rdfloat,rotfile, r_rot, vc_total

Vc = spline(r_rot, vc_total, rcentbins)
Jc = Vc*rcentbins

E = pot + 0.5*Jc^2./(rcentbins^2.)
Jc= Jc*v_units*dist_units
E = E*v_units*v_units
end


