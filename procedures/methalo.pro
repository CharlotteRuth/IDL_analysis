pro methalo, filename,h=h,g=g,d=d,s=s,psym=psym,_extra=_extra
;*********************
;Define some constants
;********************

ms_units=2.325e5
v_units=1
dist_units=1
t_units=1
XSOLfe=0.117E-2
XSOLO=0.96E-2
XSOLH=0.706
;*************************************
; READ THE FILES:
;**************************************
oxfile = filename+'.OxMassFrac'
fefile = filename+'.FeMassFrac'
HIfile = filename+'.HI'
HeIfile = filename+'.HeI'
HeIIfile = filename+'.HeII'

ox = rbarray(oxfile,/swap)
fe = rbarray(fefile,/swap)
HI = rbarray(HIfile,/swap)
HeI = rbarray(HeIfile,/swap)
HeII = rbarray(HeIIfile,/swap)
if(keyword_set(h) EQ 0 AND keyword_set(s) EQ 0) then rtipsy, filename, h,g,d,s

;*******************************************
;INDEX STARS FROM PARTICULAR HALO
;*******************************************
ind = where(fe GT -1.0)
indstars = ind(where(ind ge h.ngas+h.ndark))

oxMW = ox[indstars]
feMW = fe[indstars]
HIMW = HI[indstars]
HeIMW = HeI[indstars]
HeIIMW = HeII[indstars]
stars = s[indstars-h.ngas-h.ndark] 

;***********************************************************************
;MAKE NEW STRUCTURE CONTAINING ONLY THOSE STARS, WITH ALL THEIR INFO
;***********************************************************************
num = n_elements(stars)
MWstars = replicate({x:0.0,y:0.0,z:0.0,vx:0.0,vy:0.0,vz:0.0,mass:0.0,metals:0.0,fe:0.0,ox:0.0,tform:0.0,HI:0.0,HeI:0.0,HeII:0.0,eps:0.0,phi:0.0},num)
;rad = replicate({r:0.0},num)



MWstars.x = stars.x
MWstars.y = stars.y
MWstars.z = stars.z
MWstars.vx = stars.vx
MWstars.vy = stars.vy
MWstars.vz = stars.vz
MWstars.mass = stars.mass
MWstars.metals = stars.metals
MWstars.ox = oxMW
MWstars.fe = feMW
MWstars.HI = HIMW
MWstars.HeI = HeIMW
MWstars.HeII = HeIIMW
MWstars.tform = stars.tform
MWstars.eps = stars.eps
MWstars.phi = stars.phi

;***********************************************
;FIND THE CENTRE OF MASS AND REPOSITION
;***********************************************
cx = 0.0
cy = 0.0
cz = 0.0
mt = 0.0
limit = 500

for i=1,2 do begin
rad = where(sqrt(MWstars.x*MWstars.x+MWstars.y*MWstars.x+MWstars.z*MWstars.z)   lt limit)
cx = total(MWstars[rad].x*MWstars[rad].mass)
cy = total(MWstars[rad].y*MWstars[rad].mass) 
cz = total(MWstars[rad].z*MWstars[rad].mass) 
mt = total(MWstars[rad].mass) 

cx = cx/mt
cy = cy/mt
cz = cz/mt

MWstars.x =  MWstars.x - cx
MWstars.y =  MWstars.y - cy
MWstars.z =  MWstars.z - cz
g.x = g.x - cx
g.y = g.y - cy
g.z = g.z - cz


if  i eq 3 then limit=limit*0.5
if  i eq 6  then limit=limit*0.5
if  i eq 9  then limit=limit*0.5
if  i eq 12  then limit=limit*0.5
if  i eq 15  then limit=limit*0.5

endfor
 
;***********************************************
;eliminate proper motion of the halo/galaxy
;***********************************************
propx=mean(MWstars.vx)
propy=mean(MWstars.vy)
propz=mean(MWstars.vz)

MWstars.vx=MWstars.vx-propx
MWstars.vy=MWstars.vy-propy
MWstars.vz=MWstars.vz-propz
g.vx=g.vx-propx
g.vy=g.vy-propy
g.vz=g.vz-propz


;******************************************************
; find rotation curve and translate in XY plane
;*******************************************************
J=0.
Jxs=0.
Jys=0.
Jzs=0.
Js=0.
mt=0.
limit=200
rad = where(sqrt(MWstars.x*MWstars.x+MWstars.y*MWstars.y+MWstars.z*MWstars.z)   lt limit)
mt = total(MWstars[rad].mass)
Jxs = total(MWstars[rad].mass*(MWstars[rad].y*MWstars[rad].vz-MWstars[rad].z*MWstars[rad].vy))
Jys = total(MWstars[rad].mass*(MWstars[rad].z*MWstars[rad].vx-MWstars[rad].x*MWstars[rad].vz))
Jzs = total(MWstars[rad].mass*(MWstars[rad].x*MWstars[rad].vy-MWstars[rad].y*MWstars[rad].vx))

;****; if you want to use gas instead of stars
;rad = where(sqrt(g.x*g.x+g.y*g.y+g.z*g.z)   lt limit)
;mt = total(g[rad].mass)
;Jxs = total(g[rad].mass*(g[rad].y*g[rad].vz-g[rad].z*g[rad].vy))
;Jys = total(g[rad].mass*(g[rad].z*g[rad].vx-g[rad].x*g[rad].vz))
;Jzs = total(g[rad].mass*(g[rad].x*g[rad].vy-g[rad].y*g[rad].vx))
;*****************************

Js=sqrt(Jxs*Jxs+Jys*Jys+Jzs*Jzs)

if Js gt 0. then begin
           jjx=Jxs/Js
           jjy=Jys/Js
           jjz=Jzs/Js
           costh=jjz
           sinth=sqrt(1.0-jjz*jjz)
if  sinth gt 0.0 then begin
              sinph=jjy/sinth
              cosph=jjx/sinth
endif
endif 
if Js le 0.  then begin
           cosph = 1.0
           sinph = 0.0
endif

        ax=costh*cosph
        bx=costh*sinph
        cx=-sinth
        ay=-sinph
        by=cosph
        cy=0.0
        az=sinth*cosph
        bz=sinth*sinph
        cz=costh


        print, ax,bx,cx
        print, ay,by,cy
        print, az,bz,cz
; /**** translate star particles ****/

           txs=MWstars.x
           tys=MWstars.y
           tzs=MWstars.z
           MWstars.x=dist_units*(ax*txs+bx*tys+cx*tzs)
           MWstars.y=dist_units*(ay*txs+by*tys+cy*tzs)
           MWstars.z=dist_units*(az*txs+bz*tys+cz*tzs)

           txs=MWstars.vx
           tys=MWstars.vy
           tzs=MWstars.vz
           MWstars.vx=v_units*(ax*txs+bx*tys+cx*tzs)
           MWstars.vy=v_units*(ay*txs+by*tys+cy*tzs)
           MWstars.vz=v_units*(az*txs+bz*tys+cz*tzs)
; /**** translate gas particles ****/

           txs=g.x
           tys=g.y
           tzs=g.z
           g.x=dist_units*(ax*txs+bx*tys+cx*tzs)
           g.y=dist_units*(ay*txs+by*tys+cy*tzs)
           g.z=dist_units*(az*txs+bz*tys+cz*tzs)

           txs=g.vx
           tys=g.vy
           tzs=g.vz
           g.vx=v_units*(ax*txs+bx*tys+cx*tzs)
           g.vy=v_units*(ay*txs+by*tys+cy*tzs)
           g.vz=v_units*(az*txs+bz*tys+cz*tzs)
;*************************************************************
;change other units and define anundance ratios
;************************************************************
helium=MWstars.HeI+MWstars.HeII
hydrogen=(1.-MWstars.metals-helium)
mets=where(MWstars.metals gt 0.,comp=metalfree)
ofe=fltarr(n_elements(MWstars))
feh=ofe
feh[mets] = alog10(MWstars[mets].fe/hydrogen[mets])-alog10(XSOLFe/XSOLH) 
ofe[mets]= alog10(MWstars[mets].ox/MWstars[mets].fe)-alog10(XSOLO/XSOLFe)
feh[metalfree]=min(feh)
ofe[metalfree]=min(ofe)
MWstars.mass=ms_units*MWstars.mass
MWstars.tform=t_units*MWstars.tform


minfeh=-3.0
nbins=(max(feh)-minfeh)/0.1
notcrap=where(ofe GT -0.5)
ncofe=ofe[notcrap]
h=histogram(feh[notcrap],binsize=0.1,min=minfeh,reverse_indices=r,locations=loc)
meanV=fltarr(nbins)
errV=fltarr(nbins)
for i=0,nbins-1 do begin
  mom=moment(ncofe[R[R[I] : R[i+1]-1]])
  meanV[i] = mom[0]
  errV[i]= sqrt(mom[1])
endfor

oploterror,loc,meanV,errV,psym=psym
END

