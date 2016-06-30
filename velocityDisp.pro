pro meanrotcur,mstr=mstr,nbins=nbins,MASS=mass,_extra=_extra,RES = res

if (keyword_set(nbins) eq 0) then nbins=40
if (keyword_set(mass) eq 0) then mass = 1
if (keyword_set(res) eq 0) then res = '?'

Res=['5E1R','1E2R','1E3R','1E4R','1E5R']

for rct=0,4 do begin
    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
    file=Res[rct]+'/'+mstr+'/o'+mstr+'.00300'
    print,file
    rtipsy,file,h,g,d,s
if ((Size(s,/N_Elements) GT 1) OR (s[0].Mass NE 0)) then begin
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
;print,max(meanv),' ',min(meanv)
for i=1,4 do begin
rad = where(sqrt(s.x*s.x+s.y*s.x+s.z*s.z) lt limit)
if (total(rad) GT -1)then begin
    cx = total(s[rad].x*s[rad].mass)
    cy = total(s[rad].y*s[rad].mass)
    cz = total(s[rad].z*s[rad].mass)
    mt = total(s[rad].mass)
endif


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
Rxyz = sqrt(s.x*s.x +s.y*s.y+s.z*s.z)
V=(s.x*s.vy-s.y*s.vx)/Rxy
U=(s.x*s.vx+s.y*s.vy)/Rxy
W=(s.vz)
T =sqrt(U*U+W*W)

;print,Rxy
;print,total(Rxy)
if (Size(Rxy,/N_ELEMENTS) GT 1) then m=moment(Rxy) else m=[0,0]
Rxyrange=4.0*sqrt(m[1])
;h=histogram(Rxy,nbins=nbins,max=Rxyrange,reverse_indices=r)
h=histogram(s.y,nbins=nbins,min=-Rxyrange,max=Rxyrange,reverse_indices=r)
meanV=fltarr(nbins)
errV=fltarr(nbins)
loc=fltarr(nbins)
for i=0,nbins-1 do begin
  mom=moment(s[R[R[I] : R[i+1]-1]].vx)
  ;mom=moment(V[R[R[I] : R[i+1]-1]])
  loc[i]=-Rxyrange + i*(2.*Rxyrange/nbins)
  meanV[i] = mom[0]
  errV[i]= sqrt(mom[1])
endfor
print,max(meanv),' ',min(meanv)
miny=min(meanv)-errV(where(min(meanv)))
maxy=max(meanv)+errV(where(max(meanv)))
print,maxy,' ',miny

;set_plot,'ps'
filen=mstr+'VelocityDisp.eps'
print,filen
;device,filename=filen

if (rct EQ 0) then plot,loc,meanv,xtit="R [kpc]",ytit="Velocity [km/s]",_extra=_extra,title='Velocity Dispersion -- Dark Matter Mass: 10^'+trim(mass)+sunsymbol() + ' Resolution: '+res,psym=7,yrange=[miny,maxy]

errplot,loc,meanV-errV,meanV+errV,psym=rct+1

;ploterr,loc,meanV,errV,xtit="R [kpc]",ytit="Velocity
;[km/s]",tit='Velocity Dispersion -- Dark Matter Mass:
;                                    10^'+trim(mass)+sunsymbol()+
;                                    'Resolution: '+trim
;device,/close
endif
endfor
END

