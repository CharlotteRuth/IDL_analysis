function meanrotcur,file,nbins=nbins,_extra=_extra

if (keyword_set(nbins) eq 0) then nbins=40
if (keyword_set(nsigs) eq 0) then nsigs=4.

rtipsy,file,h,g,d,s
;***********************************************
;FIND THE CENTRE OF MASS AND REPOSITION
;***********************************************
print,'Finding center of mass and recentering'
cx = 0.0
cy = 0.0
cz = 0.0
mt = 0.0
Rxy = sqrt(s.x*s.x +s.y*s.y)
limit = max(Rxy)

for i=1,4 do begin
rad = where(sqrt(s.x*s.x+s.y*s.x+s.z*s.z)   lt limit)
cx = total(s[rad].x*s[rad].mass)
cy = total(s[rad].y*s[rad].mass)
cz = total(s[rad].z*s[rad].mass)
mt = total(s[rad].mass)

cx = cx/mt
cy = cy/mt
cz = cz/mt

s.x =  s.x - cx
s.y =  s.y - cy
s.z =  s.z - cz
g.x = g.x - cx
g.y = g.y - cy
g.z = g.z - cz
d.x = d.x - cx
d.y = d.y - cy
d.z = d.z - cz


if  i eq 3 then limit=limit*0.5
if  i eq 6  then limit=limit*0.5
if  i eq 9  then limit=limit*0.5
if  i eq 12  then limit=limit*0.5
if  i eq 15  then limit=limit*0.5

endfor
print,'Done finding center of mass and recentering'

Rxy = sqrt(s.x*s.x +s.y*s.y)
Rxyzg = sqrt(g.x*g.x +g.y*g.y+g.z*g.z)
Rxyzd = sqrt(d.x*d.x +d.y*d.y+d.z*d.z)
Rxyz = sqrt(s.x*s.x +s.y*s.y+s.z*s.z)
V=(s.x*s.vy-s.y*s.vx)/Rxy
U=(s.x*s.vx+s.y*s.vy)/Rxy
W=(s.vz)
T =sqrt(U*U+W*W)

m=moment(Rxy)
Rxyrange=nsigs*sqrt(m[1])
;h=histogram(Rxy,nbins=nbins,max=Rxyrange,reverse_indices=r,locations=loc)
h=histogram(s.y,nbins=nbins,min=-Rxyrange,max=Rxyrange,reverse_indices=r,locations=loc)
meanV=fltarr(nbins)
maxV=fltarr(nbins)
errV=fltarr(nbins)
for i=0,nbins-1 do begin
  mom=moment(s[R[R[I] : R[i+1]-1]].vx)
  ;mom=moment(V[R[R[I] : R[i+1]-1]])
  meanV[i] = mom[0]
  errV[i]= sqrt(mom[1])
  if(mom[0] ge 0) then maxV[i]=meanV[i]+errV[i] $
  else if(mom[0] lt 0) then maxV[i]=mom[0]-errV[i]
endfor


print,'Calculating vc by totaling masses.'
posloc=loc[where(loc ge 0)]
nposbins = n_elements(posloc)
nnegbins = nbins-nposbins
vc=fltarr(nposbins)
mtot=0.
grord = sort(Rxyzg)
Rgord = Rxyzg[grord]
gord = g[grord]
drord = sort(Rxyzd)
Rdord = Rxyzd[drord]
dord = d[drord]
srord = sort(Rxyz)
Rsord = Rxyz[srord]
sord = s[srord]
jg=0L
jd=0L
js=0L
for i=0,nposbins-1 do begin
  rad=posloc[i]
  while(Rgord[jg] lt rad) do begin
    mtot=mtot+gord[jg].mass
    jg=jg+1
  endwhile
  while(Rdord[jd] lt rad) do begin
    mtot=mtot+dord[jd].mass
    jd=jd+1
  endwhile
  while(Rsord[js] lt rad) do begin
    mtot=mtot+sord[js].mass
    js=js+1
  endwhile
  vc[i]=sqrt(mtot/rad)
endfor


;plot,loc,meanv,/nodata,xtit="R [kpc]",ytit="Velocity [km/s]",_extra=_extra
ploterror,loc,meanV,errV,yrange=[1.5*min(maxV),1.5*max(maxV)],xtit="R [kpc]",ytit="Velocity [km/s]",line=1,_extra=_extra
oplot,loc,maxV
oplot,loc[nnegbins:*],-vc,line=2
;legend,['maxv','meanv','vc'],line=[0,1,2],/right

return,{maxV:max(maxV),r:abs(loc[where(maxV EQ max(maxV))])}
END

