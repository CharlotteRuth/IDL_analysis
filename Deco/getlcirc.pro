
;
;
;
;   return a mapping of energy to angular momentum of a circular orbit
;   for a given particle distribution
;
;   We assume that the system has already been centered and aligned
;   and return an array of Energy values and Jc to be mapped against
;   each other 
;


pro getlcirc, filename, h, g, d, s, E, Jc, $
              rotfile=rotfile, nbins=nbins,rmax=rmax,zmax=zmax, cosmo=cosmo, $
              nperbin = nperbin, log = log, equaln = equaln


; *****************************
; generate a potential profile
; *****************************

if not keyword_set(nbins) then nbins = 200.
if not keyword_set(rmax) then rmax = 20.
if not keyword_set(zmax) then zmax = 0.1
if not keyword_set(nperbin) then nperbin = 1000

rg = sqrt(g.x^2+g.y^2)
rd = sqrt(d.x^2+d.y^2)
rs = sqrt(s.x^2+s.y^2)

indg = where(abs(g.z) lt zmax and rg lt rmax, ng)
indd = where(abs(d.z) lt zmax and rd lt rmax, nd)
inds = where(abs(s.z) lt zmax and rs lt rmax, ns)

if keyword_set(log) then $
  rbins = makenlog(min(r), rmax, nbins+1) $
else if keyword_set(equaln) then begin
    rbins = equalnbins(ng+nd+ns,[rg[indg],rd[indd],rs[inds]],1000.)
    nbins = n_elements(rbins)-1
endif else $
  rbins = arrscl(findgen(nbins+1), min([rg[indg],rd[indd],rs[inds]]), rmax)

rcentbins = fltarr(nbins)
pot = fltarr(nbins)
E = fltarr(nbins)
Jc = fltarr(nbins)


for i = 0L, nbins-1 do begin

    subg = where(rg[indg] gt rbins[i] and rg[indg] lt rbins[i+1], ng)
    subd = where(rd[indd] gt rbins[i] and rd[indd] lt rbins[i+1], nd)
    subs = where(rs[inds] gt rbins[i] and rs[inds] lt rbins[i+1], ns)

    rcentbins[i] = (rbins[i]+rbins[i+1])/2.

    if ng gt 0 then pot[i] = total(g[indg[subg]].phi)
    if nd gt 0 then pot[i] = pot[i]+total(d[indd[subd]].phi)
    if ns gt 0 then pot[i] = pot[i]+total(s[inds[subs]].phi)
    if(ng eq 0 and nd eq 0 and ns eq 0) then pot[i] = 0
    pot[i] = pot[i]/(ng+nd+ns)
    
endfor

; *****************************************************************
; what is the angular momentum of a circular orbit at each radius?
; Need to use the rotation curve
; *****************************************************************

if not keyword_set(rotfile) then $
  rdfloat, filename+'.rotcur', r_rot, vc_total,/silent $
else $
  rdfloat, rotfile, r_rot, vc_total,/silent

if keyword_set(cosmo) then begin
    readcol, 'units', units,/silent
    r_rot = r_rot*units[0]*h.time
    vc_total = vc_total*units[2]*h.time
endif

Vc = spline(r_rot, vc_total, rcentbins)
Jc = Vc*rcentbins

E = pot + 0.5*Jc^2./(rcentbins^2.)


end


