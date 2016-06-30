; Written by Peter Yoachim 5/2/2006
; Copyright PY 2006
function diskfit,s,NBINS=nbins,clip=CLIP,NOPLOT=noplot,TITLE=title,_extra=_EXTRA,rmax = rmax, radial = radial
if(keyword_set(nbins) EQ 0) then nbins=100
if not (keyword_set(rmax)) then rmax = max(abs(sqrt(s.x*s.x + s.y*s.y)))
s = s[where(sqrt(s.x*s.x + s.y*s.y) LT rmax)]
bin=rmax/nbins ;(max(s.x)-min(s.x))/nbins
if keyword_set(radial) then r = sqrt(s.x*s.x + s.y*s.y) else r = s.y

;addpath, '/astro/net/grads-1/stinson/trace'

;im=hist_2d(s.x*s.mass,s.z*s.mass,bin1=bin, bin2=bin) ;so this line is
                                ;          clearly wrong, really need
                                ;a loop to calc the stellar mass in
                                ;each bin


im=better_hist2d(r,s.z,s.mass, binsize1=bin, binsize2=bin, $
                 obin1=obin1, obin2=obin2)
if(keyword_set(clip) EQ 0) then clip=max(im)/5 ;max to clip at

makex, im,x,y
;center the galaxy in an easy way
maxp=where(im eq max(im))
im=shift(im, -1.*x(maxp), -1.*y(maxp))
;now we're centered on (0,0)

x=x(*,0)
y=transpose(y(0,*))

;convert x  and y to kpc
x=x*(obin1(1)-obin1(0))
y=y*(obin2(1)-obin2(0))
if(keyword_set(noplot) EQ 0) THEN BEGIN
loadct, 39
multiplot,/reset
window,1,xsize = 400,ysize = 400
contour,alog10(im),x,y,/fill,nlevels = 240,xrange = [-1*rmax,rmax],yrange = [-1*rmax,rmax],xstyle = 1,ystyle = 1
ENDIF

rprof=total(im,2) ;edge-on radial profile

good=where(rprof lt clip and rprof gt 0 and x ne 0) ;take out the
                                ;nasty bulge, I think Besselk blows up
                                ;at 0, so take that out too
 
guess=[max(rprof(good))*2., max(x(good))/10.,0] ;make random guess based
                                ;on the size of the image.

;assume Poisson errors for lack of anything better to do
rprof_params=mpfitfun('edge_exp', x(good),rprof(good), sqrt(rprof(good)), $
             guess, perror=perror, bestnorm=bestnorm, /quiet)
print, 'edge-on scale length = ', trim(rprof_params(1)), ' kpc'

if(keyword_set(noplot) EQ 0) THEN BEGIN
window,0
!p.multi=[0,1,2]
plot, x,rprof, yrange=[min(rprof(good)),max(rprof)], /ylog, $
tit=title,xtitle='R (kpc)', ytitle='Stellar mass',_extra=_EXTRA
oplot, x(good), edge_exp(x(good), rprof_params), color=250
oplot,[-1*rmax,rmax],[clip,clip],linestyle = 1
oplot,[rprof_params(1),rprof_params(1)],[1,clip*6],linestyle = 2
oplot,-1*[rprof_params(1),rprof_params(1)],[1,clip*6],linestyle = 2
ENDIF

;now to fit the scale height
good=where(x gt 4.*rprof_params(1) and x lt 8.*rprof_params(1))
; measure the scale height between 2 and 4 scale lengths

z_prof=total(im(good,*),1) ;edge-on vertical profile

good=where(z_prof ne 0)  ;eliminate where there are no particles

guess=[max(z_prof), max(y(good))/30.,0] ;guess params

;Poisson errors on the bins
z_prof_params=mpfitfun('oned_exp',y(good), z_prof(good), $
                       sqrt(z_prof(good)), guess, perror=perror, $
                       bestnorm=bestnorm, /quiet) 

print, 'exponential scale height=', trim(z_prof_params(1)), ' kpc'
 
if(keyword_set(noplot) EQ 0) THEN BEGIN
plot, y, z_prof, /ylog, $
yrange=[min(z_prof(good)), max(z_prof(good))], $
xrange=[-4.*z_prof_params(1), 4.*z_prof_params(1)],title=title, $
xtitle='z (kpc)', ytitle='Stellar mass',_extra=_EXTRA 
oplot, y(good), oned_exp(y(good), z_prof_params), color=250
oplot,[z_prof_params(1),z_prof_params(1)],[min(z_prof(good)), max(z_prof(good))],linestyle = 2
oplot,-1*[z_prof_params(1),z_prof_params(1)],[min(z_prof(good)), max(z_prof(good))],linestyle = 2
multiplot,/reset
stop
ENDIF


;do a sech^2 profile too while I'm at it
z_prof_params2=mpfitfun('sech_one',y(good), z_prof(good), $
                       sqrt(z_prof(good)), guess, perror=perror, $
                       bestnorm=bestnorm, /quiet) 

if(keyword_set(noplot) EQ 0) THEN BEGIN
oplot, y(good), sech_one(y(good), z_prof_params2), color=50
ENDIF

;print, 'sech^2 z_0 height=', trim(z_prof_params2(1)), ' kpc'
;print, 'Beware, scale height looks to flare ~20% at larger R'

!p.multi=0
return,{hl:rprof_params[1],hz:z_prof_params[1],sechz:z_prof_params2[1]}
end
