pro alignHI,stars,dark,gas,limit, hi=hi,spinaxes = spinaxes, aligndark = aligndark,alignstars = alignstars

; edited to align with weights based on hi fraction.

;***************************************************************
; align the particles so angular momentun vector parallel to z axis
;***************************************************************

J=0.
Jxs=0.
Jys=0.
Jzs=0.
Js=0.
mt=0.

IF keyword_set(aligndark) THEN BEGIN
   gas_temp = gas
   gas = dark
ENDIF

IF NOT KEYWORD_SET(hi) THEN hi = fltarr(N_ELEMENTS(gas)) + 1.0

;*****************************
;if you want to use young stars
;******************************
r = where(sqrt(stars.x*stars.x+stars.y*stars.y+stars.z*stars.z) lt limit and stars.tform gt 0.2)
IF keyword_set(alignstars) AND n_elements(r) GT 1 THEN BEGIN
   mt = total(stars[r].mass)
   Jxs = total(stars[r].mass*(stars[r].y*stars[r].vz-stars[r].z*stars[r].vy))
   Jys = total(stars[r].mass*(stars[r].z*stars[r].vx-stars[r].x*stars[r].vz))
   Jzs = total(stars[r].mass*(stars[r].x*stars[r].vy-stars[r].y*stars[r].vx))
ENDIF ELSE BEGIN
;****; if you want to use gas instead of stars
   r = where(sqrt(gas.x*gas.x+gas.y*gas.y+gas.z*gas.z)   lt limit)
;if r[0] eq -1 then begin
;    print,'Using gas to find center of mass..'
;    foo = min(gas.phi,cmindexgas)
;    cxgas = gas[cmindexgas].x
;    cygas = gas[cmindexgas].y
;    czgas = gas[cmindexgas].z

;    stars.x = stars.x - cxgas
;    stars.y = stars.y - cygas
;    stars.z = stars.z - czgas

;    gas.x = gas.x - cxgas
;    gas.y = gas.y - cygas
;    gas.z = gas.z - czgas

;    dark.x = dark.x - cxgas
;    dark.y = dark.y - cygas
;    dark.z = dark.z - czgas

;    r = where(sqrt(gas.x*gas.x+gas.y*gas.y+gas.z*gas.z)   lt limit)
;endif

   mt = total(gas[r].mass*hi[r])
   Jxs = total(gas[r].mass*hi[r]*double((gas[r].y*gas[r].vz-gas[r].z*gas[r].vy)))
   Jys = total(gas[r].mass*hi[r]*double((gas[r].z*gas[r].vx-gas[r].x*gas[r].vz)))
   Jzs = total(gas[r].mass*hi[r]*double((gas[r].x*gas[r].vy-gas[r].y*gas[r].vx)))
ENDELSE
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

c = jjy/SQRT(jjy*jjy + jjz*jjz)
b = SQRT(1 - c*c)
print,'Spin Axes: <',jjx,',',jjy,',',jjz,'>'
print,'Aligned Axes: <0,',b,',',c,'>' 
spinaxes = [jjx,jjy,jjz]

ax=costh*cosph
bx=costh*sinph
cx=-sinth
ay=-sinph
by=cosph
cy=0.0
az=sinth*cosph
bz=sinth*sinph
cz=costh

IF keyword_set(aligndark) THEN gas = gas_temp

print, ax,bx,cx
print, ay,by,cy
print, az,bz,cz
; /**** translate star particles (and change units) ****/

txs=stars.x
tys=stars.y
tzs=stars.z
stars.x=(ax*txs+bx*tys+cx*tzs)
stars.y=(ay*txs+by*tys+cy*tzs)
stars.z=(az*txs+bz*tys+cz*tzs)

txs=stars.vx
tys=stars.vy
tzs=stars.vz
stars.vx=(ax*txs+bx*tys+cx*tzs)
stars.vy=(ay*txs+by*tys+cy*tzs)
stars.vz=(az*txs+bz*tys+cz*tzs)
; /**** translate gas particles  (and change units) ****/

txs=gas.x
tys=gas.y
tzs=gas.z
gas.x=(ax*txs+bx*tys+cx*tzs)
gas.y=(ay*txs+by*tys+cy*tzs)
gas.z=(az*txs+bz*tys+cz*tzs)

txs=gas.vx
tys=gas.vy
tzs=gas.vz
gas.vx=(ax*txs+bx*tys+cx*tzs)
gas.vy=(ay*txs+by*tys+cy*tzs)
gas.vz=(az*txs+bz*tys+cz*tzs)
; translate dark matter particles  (and change units)
txs=dark.x
tys=dark.y
tzs=dark.z
dark.x=(ax*txs+bx*tys+cx*tzs)
dark.y=(ay*txs+by*tys+cy*tzs)
dark.z=(az*txs+bz*tys+cz*tzs)

txs=dark.vx
tys=dark.vy
tzs=dark.vz
dark.vx=(ax*txs+bx*tys+cx*tzs)
dark.vy=(ay*txs+by*tys+cy*tzs)
dark.vz=(az*txs+bz*tys+cz*tzs)

RETURN
END
