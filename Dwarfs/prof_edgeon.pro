;
;
;
;
;
; make an edge on profile - based on peter's code
;
;
;


function prof_edgeon, s, nbins = nbins, rmax = rmax, orient = orient 


if(keyword_set(nbins) EQ 0) then nbins=100
if(keyword_set(rmax) eq 0) then rmax = 15.

bin=rmax/nbins

s = s[where(sqrt(s.x^2+s.y^2) lt rmax)]

if(keyword_set(orient) eq 0) then orient = 'x'

if orient eq 'x' then $
  im=better_hist2d(s.x,s.z,s.mass, binsize1=bin, binsize2=bin, $
                   obin1=obin1, obin2=obin2) $
else $
  im=better_hist2d(s.y,s.z,s.mass, binsize1=bin, binsize2=bin, $
                   obin1=obin1, obin2=obin2) 

if(keyword_set(clip) EQ 0) then clip=max(im)/5 ;max to clip at

;makex, im,x,y

;x=x(*,0)
;y=transpose(y(0,*))

;convert x  and y to kpc
;x=x*(obin1(1)-obin1(0))
;y=y*(obin2(1)-obin2(0))

rprof=total(im,2) ;edge-on radial profile

prof = {rho: rprof, x: obin1}

return, prof

end
