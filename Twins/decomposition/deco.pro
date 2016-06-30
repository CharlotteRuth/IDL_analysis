;v_units = 2418.46 ;h603, h516
;dist_units = 50000 ;h603
;dist_units = 25000 ;h516
; ecut = -3e4 ;h603
;deco,'h603.cosmo50cmb.2304g2bwK.00512','h603.cosmo50cmb.2304g2bwK.00512.1',1,50000,2418.46,ecut= -3e4 
;deco,'h516.cosmo25cmb.2304g2bwK.00512','h516.cosmo25cmb.2304g2bwK.00512.1',1,25000,2418.46 
;.r ~/code/idl4tipsy/align.pro
;.r ~/code/Twins/decomposition/deco.pro
pro deco, filename,rotfile,halo,dist_units,v_units,Ecut = ECUT
;*************************************
; READ THE FILES:
;**************************************
gtpfile = filename+'.amiga.gtp'
grpfile = filename+'.amiga.grp'
grp = read_ascii_array(grpfile)
rotfile=rotfile+'.rot'

if(keyword_set(h) EQ 0 AND keyword_set(s) EQ 0) then rtipsy, filename, h,g,d,s
;*******************************************
;INDEX GAS,DM & STARS FROM PARTICULAR HALO
;*******************************************
ind = where(grp eq halo,comp=indcomp)
inds = ind(where(ind ge h.ngas+h.ndark)) 
indg = ind(where(ind lt h.ngas))
indd = ind(where(ind ge h.ngas and ind lt h.ngas+h.ndark))
indsats= indcomp(where(indcomp ge h.ngas+h.ndark))


index=lonarr(h.nstar)
for i=0L,h.nstar-1 do index[i] = i
indstars=index[inds-h.ngas-h.ndark]

stars = s[inds-h.ngas-h.ndark] 
gas = g[indg]
dark = d[indd-h.ngas] 

;***************************************************************************
; get CENTRE OF MASS from amiga (first entry in gtp file): REPOSITION
;**************************************************************************
rtipsy, gtpfile, h1,g1,d1,s1

cx= s1[halo-1].x
cy= s1[halo-1].y
cz= s1[halo-1].z
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
propx=s1[halo-1].vx
propy=s1[halo-1].vy
propz=s1[halo-1].vz

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
; translate into XY plane using align.pro
;*************************************************************
limit=5./dist_units            ;use gas within this limit (kpc/dist_units)
align,gas,dark,stars,limit

;*************************************************************
;change units 
;************************************************************
units,gas,dark,stars,h.time

;********************************************************
; find J/Jz
;********************************************************
loadct,39
Jx = (stars.y*stars.vz-stars.z*stars.vy)
Jy = (stars.z*stars.vx-stars.x*stars.vz)
Jz = (stars.x*stars.vy-stars.y*stars.vx)
J=sqrt(Jx*Jx+Jy*Jy+Jz*Jz)


;***********************************************************
;find kinetic, potential & total Energy of the stars
;**********************************************************

Rxyz = sqrt(stars.x*stars.x+stars.y*stars.y+stars.z*stars.z) 
v2=(stars.vx*stars.vx+stars.vy*stars.vy+stars.vz*stars.vz)
ke=0.5*v2
Etot=(ke+stars.phi)
;********************************************************************
;assuming amiga selects only bound stars, renormalise the potential
;*******************************************************************
stars.phi=stars.phi-max(Etot)
Etot=(ke+stars.phi)


;find the max radius of the rotation curve
lines=file_lines(rotfile)
readcol,rotfile,radius,skipline=lines-2
radius=radius[0]*dist_units
nstars = n_elements(stars.x)
JzJzmax=    replicate(-2.,nstars)
region = where(Rxyz lt radius) ; set to max radius of rotation curve

;*********************************************************************
; relate Energy and max Jz at given radius in the region or rot curve:
;**********************************************************************
getlcirc,rotfile,stars[region],E,Jc,rmax=radius,zmax=radius,dist_units,v_units
Energy = Etot[region]
sort=sort(Energy)
Energy = Energy[sort]
Jzmax = spline(E,Jc,Energy)
sort=sort(sort)
Jzmax = Jzmax[sort]
Energy = Energy[sort]
keep = where(Jzmax lt max(Jc) and Energy lt max(E)); get rid of weird outputs
Jzmax=Jzmax[keep]
Energy=Energy[keep]
region=region[keep]
JzJzmax[region]=Jz[region]/Jzmax


;********************************************************************
; find JzJc cutoff so the sphere does not rotate, and define disk stars
;*******************************************************************
Rxy = sqrt(stars.x*stars.x+stars.y*stars.y) 
           V=(stars.x*stars.vy-stars.y*stars.vx)/Rxy 
allstars = where(stars.x le max(stars.x))
outer = setdifference(allstars,region)
cut=.7
Vsphere = 10
while (Vsphere gt 0.) do begin 
test=where(JzJzmax[region] lt cut)
Vsphere=mean(V[region[test]])
cut = cut -0.0002
endwhile

sphere1=region[where(JzJzmax[region] lt cut,comp=rest1)]
sphere2= region[where(JzJzmax[region] gt cut and Jz[region]/J[region] lt 1./sqrt(3))]
sphere=[sphere1,sphere2]
sphere=[sphere,outer]; the unbound stuff is part of the spheroid
disk=region[where(JzJzmax[region] gt .8 and Jz[region]/J[region] gt 1./sqrt(3))]
rest=setdifference(allstars,[sphere,disk])

;***********************************************************
;make an energy cut to get halo/bulge & thick/pseudo
;**************************************************************
IF NOT(KEYWORD_SET(Ecut)) THEN BEGIN
    set_plot,'x'
    plot,Etot,Jz,psym=3,xtitle = "Etot",ytitle = 'Jz'         ; check the range for Ecut
    print,'Set Ecut to the energy at which the distribution begins to flare: '
    stop
ENDIF

 ;was -1.5e5
halo=where(Etot[sphere] gt Ecut,comp=bulge)   ; note now Rxyz[sphere[halo]]
;raise ecut to get less halo
thick = where(Etot[rest] gt Ecut,comp=PsBulge)

disktot = [disk,rest[thick]]
bulgetot= [sphere[bulge],rest[Psbulge]]

IF NOT KEYWORD_SET(noplot) THEN BEGIN
    set_plot,'x'
    loadct,39
    hist=hist1D(rxy[sphere],stars[sphere].mass,obin=xvals,nbins = 40)
    area=3.1459*((xvals+0.5)^2-(xvals-0.5)^2)
    plot,xvals,alog10(hist/area),xrange=[0,10],yrange=[4,12],title = filename + ' Ecut = '+strtrim(ecut,2),xtitle = 'Radius [kpc]'
    hist=hist1D(rxy[sphere[halo]],stars[sphere[halo]].mass,obin=xvals,nbins = 40)
    oplot,xvals,alog10(hist/area),color=250
    hist=hist1D(rxy[sphere[bulge]],stars[sphere[bulge]].mass,obin=xvals,nbins = 40)
    oplot,xvals,alog10(hist/area),color=150
    legend,['Sphereoid','Bulge','Halo'],color = [0,150,240],linestyle=[0,0,0]
    stop

    set_plot,'x'
    loadct,39
    hist=hist1D(rxy,stars.mass,obin=xvals)
    area=3.1459*((xvals+0.5)^2-(xvals-0.5)^2)
    plot,xvals,alog10(hist/area),xrange=[0,10],yrange=[4,12],title = filename,xtitle = 'Radius [kpc]'
    hist=hist1D(rxy[sphere],stars[sphere].mass,obin=xvals)
    oplot,xvals,alog10(hist/area),color=250
    hist=hist1D(rxy[disk],stars[disk].mass,obin=xvals)
    oplot,xvals,alog10(hist/area),color=50
    legend,['Sphereoid','Disk'],color = [250,50],linestyle=[0,0]
    stop

    plot,stars[sphere].x,stars[sphere].z,psym=3,xrange=[-40,40],yrange=[-40,40],title = filename + ' Ecut = '+strtrim(ecut,2),xtitle = 'Radius [kpc]'
    oplot,stars[sphere[halo]].x,stars[sphere[halo]].z,psym=3,color=250 
    oplot,stars[disk].x,stars[disk].z,psym=3,color=50
    oplot,stars[sphere[bulge]].x,stars[sphere[bulge]].z,psym=3,color=150
    oplot,stars[rest[thick]].x,stars[rest[thick]].z,psym=3,color=100
    oplot,stars[rest[psbulge]].x,stars[rest[psbulge]].z,psym=3,color=200

    set_plot,'ps'
    device,filename = rotfile + '_profile.ps',/color,bits_per_pixel=8
; plot the density profile and check it out:
    loadct,39
    hist=hist1D(rxy[sphere],stars[sphere].mass,obin=xvals)
    area=3.1459*((xvals+0.5)^2-(xvals-0.5)^2)
    plot,xvals,alog10(hist/area),xrange=[0,10],yrange=[4,12],title = filename + ' Ecut = '+strtrim(ecut,2)
    hist=hist1D(rxy[sphere[halo]],stars[sphere[halo]].mass,obin=xvals)
    oplot,xvals,alog10(hist/area),color=250
    hist=hist1D(rxy[sphere[bulge]],stars[sphere[bulge]].mass,obin=xvals)
    oplot,xvals,alog10(hist/area),color=150
    legend,['Bulge','Halo'],color = [150,240],linestyle=[0,0]
    device,/close


    device,filename = rotfile + '_galaxy.ps',/color,bits_per_pixel=8
    plot,stars[sphere].x,stars[sphere].z,psym=3,xrange=[-40,40],yrange=[-40,40],title = filename + ' Ecut = '+strtrim(ecut,2)
    oplot,stars[sphere[halo]].x,stars[sphere[halo]].z,psym=3,color=250 
    oplot,stars[disk].x,stars[disk].z,psym=3,color=50
    oplot,stars[sphere[bulge]].x,stars[sphere[bulge]].z,psym=3,color=150
    oplot,stars[rest[thick]].x,stars[rest[thick]].z,psym=3,color=100
    oplot,stars[rest[psbulge]].x,stars[rest[psbulge]].z,psym=3,color=200
    device,/close
ENDIF


;************************************************
;write iord values for compoenents
;*****************************************************

index=lonarr(h.n)
index[0:h.n-1] = 0
index[h.ngas+h.ndark+indstars[disk]] = 1
index[h.ngas+h.ndark+indstars[sphere[halo]]] = 2
index[h.ngas+h.ndark+indstars[sphere[bulge]]] = 3
index[h.ngas+h.ndark+indstars[rest[thick]]] = 4
index[h.ngas+h.ndark+indstars[rest[PsBulge]]] = 5
index[indsats] = 6

file = rotfile
openw,1,file+'.cmp'
printf,1,h.n
for i=0L,n_elements(index)-1 do printf,1,index[i]
close,1

stop

end



