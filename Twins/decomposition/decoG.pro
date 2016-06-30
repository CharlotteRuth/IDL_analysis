pro decoG1, filename, group, h=h,g=g,d=d,s=s
;*********************
;Define some constants
;********************
;XSOLfe=0.117E-2
;XSOLO=0.96E-2
;XSOLH=0.706

;*************************************
; READ THE FILES:
;**************************************
;oxfile = filename+'.OxMassFrac'
;fefile = filename+'.FeMassFrac'
;HIfile = filename+'.HI'
;HeIfile = filename+'.HeI'
;HeIIfile = filename+'.HeII'
;Tagfile = filename+'.iord'
;Taggfile = filename+'.igasorder'
grpfile = filename+'.amiga.grp'
gtpfile = filename+'.amiga.gtp'

;ox_in = read_ascii_array(oxfile)
;fe_in = read_ascii_array(fefile)
;HI_in = read_ascii_array(HIfile)
;HeI_in = read_ascii_array(HeIfile)
;HeII_in = read_ascii_array(HeIIfile)
;tag_in = read_ascii_array(Tagfile)
;tagg_in =  read_ascii_array(Taggfile)
grp = read_ascii_array(grpfile)

if(keyword_set(h) EQ 0 AND keyword_set(s) EQ 0) then rtipsy, filename, h,g,d,s
;*******************************************
;INDEX GAS,DM & STARS FROM PARTICULAR HALO
;*******************************************
ind = where(grp eq group,comp=indcomp)
inds = ind(where(ind ge h.ngas+h.ndark)) 
indg = ind(where(ind lt h.ngas))
indd = ind(where(ind ge h.ngas and ind lt h.ngas+h.ndark))
indsats= indcomp(where(indcomp ge h.ngas+h.ndark))

; ox = ox_in[inds]
; fe = fe_in[inds]
; HI = HI_in[inds]
; HeI = HeI_in[inds]
; HeII = HeII_in[inds]
; tag = tag_in[inds]
; tagg = tagg_in[inds]

index=lonarr(h.nstar)
for i=0L,h.nstar-1 do index[i] = i
indstars=index[inds-h.ngas-h.ndark]

stars = s[inds-h.ngas-h.ndark] 
gas = g[indg]
dark = d[indd-h.ngas] 

;**************************************************************************
;work out how many of each type, and put them into num (equiv to tipsy h)
;***********************************************************************
num=h
num.nstar=n_elements(stars)
num.ngas=n_elements(gas)
num.ndark=n_elements(dark)
num.n=num.ngas+num.ndark+num.nstar

;***************************************************************************
; get CENTRE OF MASS from amiga (first entry in gtp file): REPOSITION
;**************************************************************************
rtipsy, gtpfile, h1,g1,d1,s1

cx= s1[0].x
cy= s1[0].y
cz= s1[0].z
stars.x =  stars.x - cx
stars.y =  stars.y - cy
stars.z =  stars.z - cz
gas.x = gas.x - cx
gas.y = gas.y - cy
gas.z = gas.z - cz
dark.x = dark.x - cx
dark.y = dark.y - cy
dark.z = dark.z - cz
 
;***********************************************
;eliminate proper motion of the halo/galaxy
;***********************************************
propx=mean(stars.vx)
propy=mean(stars.vy)
propz=mean(stars.vz)

stars.vx=stars.vx-propx
stars.vy=stars.vy-propy
stars.vz=stars.vz-propz
gas.vx=gas.vx-propx
gas.vy=gas.vy-propy
gas.vz=gas.vz-propz
dark.vx=dark.vx-propx
dark.vy=dark.vy-propy
dark.vz=dark.vz-propz

;**************************************************************
; find rotation curve and translate in XY plane using align.pro
;*************************************************************
dist_units=1e5*h.time
limit=5./dist_units            ;use stars within this limit (kpc/dist_units)
;align,gas,dark,stars,limit
align,stars,dark,gas,limit

;*************************************************************
;change units and define anundance ratios
;************************************************************
unitsG,stars,gas,dark,h.time
; helium=HeI+HeII
; hydrogen=(1.-stars.metals-helium)
; mets=where(stars.metals gt 0.,comp=metalfree)
; ofe=fltarr(n_elements(stars))
; feh=fltarr(n_elements(stars))
; feh[mets] = alog10(fe[mets]/hydrogen[mets])-alog10(XSOLFe/XSOLH) 
; ofe[mets] = alog10(ox[mets]/fe[mets])-alog10(XSOLO/XSOLFe)
; feh[metalfree]=min(feh)
; ofe[metalfree]=max(ofe)
;********************************************************
; find J/Jz
;********************************************************
loadct,39
Jx = (stars.y*stars.vz-stars.z*stars.vy)
Jy = (stars.z*stars.vx-stars.x*stars.vz)
Jz = (stars.x*stars.vy-stars.y*stars.vx)
J=sqrt(Jx*Jx+Jy*Jy+Jz*Jz)

;for gas as well

Jxg = (gas.y*gas.vz-gas.z*gas.vy)
Jyg = (gas.z*gas.vx-gas.x*gas.vz)
Jzg = (gas.x*gas.vy-gas.y*gas.vx)
Jg = sqrt(Jxg*Jxg+Jyg*Jyg+Jzg*Jzg)


;*****************************************************
; cylindrical  velocities 
;*****************************************************

Rxy = sqrt(stars.x*stars.x+stars.y*stars.y) 
height = sqrt(stars.z*stars.z)
           V=(stars.x*stars.vy-stars.y*stars.vx)/Rxy 
           U=(stars.x*stars.vx+stars.y*stars.vy)/Rxy 
           W=(stars.vz)
           T =sqrt(U*U+W*W)

;***********************************************************
;find kinetic energy, potential energy, Jc of the stars
;**********************************************************

Rxyz = sqrt(stars.x*stars.x+stars.y*stars.y+stars.z*stars.z) 
v2=(stars.vx*stars.vx+stars.vy*stars.vy+stars.vz*stars.vz)
ke=0.5*v2
Energy=(ke+stars.phi)
Jc=sqrt(v2)*Rxyz

; for gas as well
Rg = sqrt(gas.x*gas.x+gas.y*gas.y+gas.z*gas.z) 
v2g=(gas.vx*gas.vx+gas.vy*gas.vy+gas.vz*gas.vz)
keg=0.5*v2g
Energyg=(ke+gas.phi)
Jcg=sqrt(v2g)*Rg

diskgas = where (Rg lt 20 and Jzg/Jcg gt .8)

;********************************************************************
;assuming amiga selects only bound stars, renormalise the potential
;*******************************************************************
stars.phi=stars.phi-max(Energy)
Energy=(ke+stars.phi)

JzJc=Jz/Jc
JzJ=Jz/J
JJc=J/Jc
JzJJcJc=J*Jz/Jc/Jc
JzJzJJc=Jz*Jz/Jc/J

;********************************************************************
; find JzJc cutoff so the sphere does not rotate, and define disk stars
;*******************************************************************
cut=.8
Vsphere = 10.
while (Vsphere gt 0.) do begin 
sphere=where(JzJc lt cut,comp=disk)
Vsphere=mean(V[sphere])
cut = cut -0.0002
endwhile

sphere=where(JzJc lt cut or  abs(stars.z) ge 6. or Rxyz gt 35,comp=rest1)
disk = where(JzJc gt 0.8 and abs(stars.z) lt 6.and Rxyz lt 35.,comp=rest2)
rest=intersect(rest1,rest2)

;****************************************************************
; check rotation/morphology of young spheroid: disk contaminants?
;***************************************************************
youngsphere=where(stars[sphere].tform gt 0.80*(max(stars.tform)))
print,'rotation of Young Spheroid ='
print,mean(V[sphere[youngsphere]])

if mean(V[sphere[youngsphere]]) gt 50. then begin
print,'looks like spheroid has disk contamination'
endif

;*************************************************
;make an energy cut to get halo/bulge
;*************************************
;Ecut=-3.2e5 ; cut for step 444     (=-4e5 for 512)
Ecut=-4e5
halo=where(Energy[sphere] gt Ecut,comp=bulge)   ; note now Rxyz[sphere[halo]]

; young=where(stars[disk].tform gt 0.8*max(stars.tform))
; mid=where(stars[disk].tform gt 4)   ; since thick disk era
; old=where(stars[disk].tform lt 4)

;plotGals,stars,disk,sphere,halo,bulge,rest,U,V,W,Rxy,Rxyz,mid,old,young,JzJc

;************************************************
;write files for compoenents
;*****************************************************

index=lonarr(h.n)
index[0:h.ngas+h.ndark-1] = 0

index[h.ngas+h.ndark+indstars[disk]] = 1
index[indg[diskgas]] = 1
index[h.ngas+h.ndark+indstars[sphere[halo]]] = 2
index[h.ngas+h.ndark+indstars[sphere[bulge]]] = 3
index[h.ngas+h.ndark+indstars[rest]] = 4
index[indsats] = 5

openw,1,filename+'.cmp'
printf,1,h.n
for i=0L,n_elements(index)-1 do printf,1,index[i]

close,1



stop


end

