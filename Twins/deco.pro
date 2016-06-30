pro deco, filename,rotfile,halo,dist_units,ms_units,Ecut,NOCHECK=nocheck

;PROCEDURE TO DECOMPOESE A GALAXY INTO BULGE, DISK, THICK DISK, PSEUDO
;BULGE AND HALO
;ORIGINALLY WRITTEN BY CHRIS BROOK
;MOST RECENTLY UPDATED BY LAUREN POPE <lapope@u.washington.edu> 7/27/2010

if (N_PARAMS() eq 0) then begin
    print, 'This is the procedure to kinematically decompose a galaxy simulation according to angular momentum.'
    print, 'Format iput as follows'
    print, 'deco, tipsyfile, rotation curve file, halo number (from amiga), dKpcUnit, dMsolUnit(from the param file), Energy Cut'
endif

print,'Deco.pro will output 5 files: range.ps, density_profile.ps,'
print,'decomposition_face.ps, decomposition_edge.ps and filename+.masses'


;*************************************
; READ THE FILES:
;**************************************
print, 'Reading in the Amiga files.....'
gtpfile = filename+'.amiga.gtp'
grpfile = filename+'.amiga.grp'
grp = read_lon_array(grpfile)
rotfile=rotfile                 ;this is the rotation file created using Tipsy
print,'Reading in the tipsy file. (This can take a few moments for large numbers of particles).....'
rtipsy, filename, h,g,d,s       ;read in the file 
dist_units=h.time*dist_units 

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
; Get CENTRE OF MASS from amiga (first entry in gtp file): REPOSITION
;**************************************************************************
print,'Getting the center of mass from the Amiga file and repositioning......'
rtipsy, gtpfile, h1,g1,d1,s1
cx= s1[halo-1].x
cy= s1[halo-1].y
cz= s1[halo-1].z

Print,'cx (Amiga): ',cx
print,'cv (Amiga): ',cv
print,'cz (Amiga): ',cz

gastemp = gas
gastemp.x = gas.x - cx
gastemp.y = gas.y - cy
limit=1./dist_units        ;use gas within this limit (kpc/dist_units)
radius = sqrt(gastemp.x*gastemp.x + gastemp.y*gastemp.y)
dens = MAX(g[radius lt limit].dens,subscript_min = subscript_min)
cx = gas[subscript_min].x
cy = gas[subscript_min].y
cz = gas[subscript_min].z

print,''
print,'cx (Rho Max): ',cx
print,'cv (Rho Max): ',cv
print,'cz (Rho Max): ',cz

Stars.x =  stars.x - cx
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


for p=0,1 do begin
print,'Entering Loop for Eliminating Proper Motion......'
    if p eq 0 then begin
        propx=s1[halo-1].vx
        propy=s1[halo-1].vy
        propz=s1[halo-1].vz
    endif
    
    if p eq 1 then begin
        propx=total(stars[disk].vx*stars[disk].mass)/total(stars[disk].mass)
        propy=total(stars[disk].vy*stars[disk].mass)/total(stars[disk].mass)
        propz=total(stars[disk].vz*stars[disk].mass)/total(stars[disk].mass)
    endif
    
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
    limit=5./dist_units    ;use gas within this limit (kpc/dist_units)
    align,gas,dark,stars,limit

    
;********************************************************
; find J/Jz
;********************************************************
;loadct,39

;***********************************************************
;find kinetic, potential & total Energy of the stars
;**********************************************************
    print,'Finding kinetic, potential, and total energy of the stars.....'
    Rxyz = sqrt(stars.x*stars.x+stars.y*stars.y+stars.z*stars.z) 
    v2=(stars.vx*stars.vx+stars.vy*stars.vy+stars.vz*stars.vz)
    ke=0.5*v2
    Etot=(ke+stars.phi)
    print, 'max Etot = ',max(Etot)

;********************************************************************
;assuming amiga selects only bound stars, renormalise the potential
;*******************************************************************
    print,'Renormalizing Potential for bound stars only....'
    stars.phi=stars.phi-max(Etot)
    gas.phi = gas.phi - max(Etot)
    dark.phi = dark.phi - max(Etot)
    Etot=(ke+stars.phi)
    
    Jz = (stars.x*stars.vy-stars.y*stars.vx)
    Jx = (stars.y*stars.vz-stars.z*stars.vy)
    Jy = (stars.z*stars.vx-stars.x*stars.vz)
    J=sqrt(Jx*Jx+Jy*Jy+Jz*Jz)
    
;find the max radius of the rotation curve
    print,'Finding the region included within the rotation curve.....'
    lines=file_lines(rotfile)
    readcol,rotfile,radius,skipline=lines-2,/silent
    radius=radius[0]
    nstars = n_elements(stars.x)
    JzJzmax=    replicate(-2.,nstars)
    region = where(Rxyz lt radius) ; set to max radius of rotation curve

;*********************************************************************
; relate Energy and max Jz at given radius in the region of rot curve:
;**********************************************************************
    print,'Relating energy and max J_z within the region of the rotation curve.....'
    getlcirc,'dummy',h,gas,dark,stars,E,Jc,rmax=radius,zmax=0.2/dist_units,rotfile=rotfile
    Energy = Etot[region]
    sort=sort(Energy)
    Energy = Energy[sort]
    Jzmax = spline(E,Jc,Energy)
    sort=sort(sort)
    Jzmax = Jzmax[sort]
    Energy = Energy[sort]
    keep = where(Jzmax lt max(Jc) and Energy lt max(E)) ; get rid of weird outputs
    Jzmax=Jzmax[keep]
    Energy=Energy[keep]
    region=region[keep]
    JzJzmax[region]=Jz[region]/Jzmax
    
    disk=region[where(JzJzmax[region] gt .8 and JzJzmax[region] lt 1.1)]
   
endfor


;********************************************************************
; find JzJc cutoff so the sphere does not rotate, and define disk stars
;*******************************************************************
print,'Finding J_z & J_C cutoff so the sphere does not rotate, also defining disk stars.....'
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

print,'Defining the spheroid.....'
sphere=region[where(JzJzmax[region] lt cut,comp=rest1)]
sphere=[sphere,outer]; the unbound stuff is part of the spheroid

rest=setdifference(allstars,[sphere,disk])

;***********************************************************
;make an energy cut to get halo/bulge & thick/pseudo
;**************************************************************


;Plot the energy range here, just in case it crashes. Make sure Ecut falls within the range or else bad things will happen (ie: a crash).

print,'Plotting the energy range, in case of crash make sure Ecut is within this range.....'

set_plot,'ps'
device,filename='range.ps',/color

plot,Etot,Jz,psym=3,title='Energy',xtitle='Etot',ytitle='J_z'; check the range for Ecut

device,/close
set_plot,'x'

IF NOT  keyword_set(NOCHECK) THEN BEGIN
    check = ''
    read,'Would you like to stop and check the energy cut? (y/n)    ',check

    IF ((check EQ 'y') OR (check EQ 'Y') OR (check EQ 'yes') OR (check EQ 'Yes') OR (check EQ 'YES')) THEN BEGIN
        print,'Plotting the Energy Distribution and Density Profile to check the Energy Cut'
        plot,Etot,Jz,psym=3,title='Energy',xtitle='Etot',ytitle='J_z' ; check the range for Ecut
        oplot,[Ecut,Ecut],[min(Jz),max(Jz)],color=80
        print,'Ecut = ',Ecut
        Change = ''
        read,'Would you like to change your guess for the Energy cut? (y/n)     ',change
        IF ((change EQ 'y') OR (change EQ 'yes')) THEN BEGIN
            read,'Enter new Ecut     ',Ecut
            oplot,[Ecut,Ecut],[min(Jz),max(Jz)],color=240
            print,"Continuing with new Ecut"
        ENDIF ELSE BEGIN 
            print,'Continuing with original Ecut =  ',Ecut
        print, 'Dont forget to look at the Density Profile plot to see how good your cut is!!!'

        ENDELSE

    ENDIF ELSE BEGIN 
        print,'Continuing with original Ecut =  ',Ecut
        print, 'Dont forget to look at the Density Profile plot to see how good your cut is!!!'

    ENDELSE
ENDIF

halo=where(Etot[sphere] gt Ecut,comp=bulge)   ; note now Rxyz[sphere[halo]]
thick = where(Etot[rest] gt Ecut,comp=PsBulge)

disktot = [disk,rest[thick]]
bulgetot= [sphere[bulge],rest[Psbulge]]
 
;************************************************
;write iord values for components
;*****************************************************
;stop

print,'Assigning an Index to each component:'
print,'Disk = 1'
print,'Halo = 2'
print,'Bulge = 3'
print,'Thick Disk = 4'
print,'Pseudo Bulge = 5'
print,'Satellites = 6'

index=lonarr(h.n)
index[0:h.n-1] = 0
index[h.ngas+h.ndark+indstars[disk]] = 1
index[h.ngas+h.ndark+indstars[sphere[halo]]] = 2
index[h.ngas+h.ndark+indstars[sphere[bulge]]] = 3
index[h.ngas+h.ndark+indstars[rest[thick]]] = 4
index[h.ngas+h.ndark+indstars[rest[PsBulge]]] = 5
index[indsats] = 6

openw,1,filename+'.cmp'
printf,1,h.n
for i=0L,n_elements(index)-1 do printf,1,index[i]

close,1

;stop


;*****************************************************
; Plot the density profile and check it out
;*****************************************************

set_plot,'ps'
device,/color,bits_per_pixel=8,filename='density_profile.ps',/inch,ysize=5,xsize=5,/times

loadct,39
interval=1./dist_units
hist=hist1D(rxy[sphere],stars[sphere].mass,bins=interval,obin=xvals)
area=3.1459*((xvals+0.5*interval)^2-(xvals-0.5*interval)^2)
plot,dist_units*xvals,alog10(ms_units*hist/(area*dist_units^2)),xrange=[0,10.],yrange=[3,10],xtit='Radius (kpc)',ytit='Log Density (Msol/Kpc^2)',title='Density Profile of Stars',/XSTYLE,/YSTYLE
hist=hist1D(rxy[sphere[halo]],stars[sphere[halo]].mass,obin=xvals,bins=interval)
oplot,dist_units*xvals,alog10(ms_units*hist/(area*dist_units^2)),color=250
hist=hist1D(rxy[sphere[bulge]],stars[sphere[bulge]].mass,obin=xvals,bins=interval)
oplot,dist_units*xvals,alog10(ms_units*hist/(area*dist_units^2)),color=150

legend,['All stars in the Sphere','Halo','Bulge'],textcolor=[0,250,150],linestyle=[0,0,0],colors=[0,240,150],/right



device,/close
set_plot,'x'


;*****************************************************
; Plot the components in a square so they look right.
;*****************************************************

set_plot,'ps'
device,/color,bits_per_pixel=8,filename='decomposition_edge.ps',/inch,ysize=5,xsize=5,/times

plot,dist_units*stars[sphere[halo]].x,dist_units*stars[sphere[halo]].z,psym=3,xrange=[-15.,15.],yrange=[-15.,15.],xtit='kpc (x direction)',ytit='kpc (z direction)',title='Edge on View',/NODATA,/XSTYLE,/YSTYLE
oplot,dist_units*stars[sphere[halo]].x,dist_units*stars[sphere[halo]].z,psym=3,color=250 
oplot,dist_units*stars[sphere[bulge]].x,dist_units*stars[sphere[bulge]].z,psym=3,color=30
oplot,dist_units*stars[rest[thick]].x,dist_units*stars[rest[thick]].z,psym=3,color=100
oplot,dist_units*stars[disk].x,dist_units*stars[disk].z,psym=3,color=50
oplot,dist_units*stars[rest[psbulge]].x,dist_units*stars[rest[psbulge]].z,psym=3,color=200

legend,['Bulge','Psbulge','Disk','Thickdisk','Halo'],textcolors=[30,200,50,100,250],position=[-18.5,12.5]

device,/close
set_plot,'x'


set_plot,'ps'
device,/color,bits_per_pixel=8,filename='decomposition_face.ps';,/inch,ysize=10,xsize=10,/times

plot,dist_units*stars[sphere[halo]].x,dist_units*stars[sphere[halo]].y,psym=3,xrange=[-15.,15.],yrange=[-15.,15.],xtit='kpc (x direction)',ytit='kpc (y direction)',title='Face on View',/NODATA,/XSTYLE,/YSTYLE
oplot,dist_units*stars[sphere[halo]].x,dist_units*stars[sphere[halo]].y,psym=3,color=250 
oplot,dist_units*stars[sphere[bulge]].x,dist_units*stars[sphere[bulge]].y,psym=3,color=30
oplot,dist_units*stars[rest[thick]].x,dist_units*stars[rest[thick]].y,psym=3,color=100
oplot,dist_units*stars[disk].x,dist_units*stars[disk].y,psym=3,color=50
oplot,dist_units*stars[rest[psbulge]].x,dist_units*stars[rest[psbulge]].y,psym=3,color=200

legend,['Bulge','Psbulge','Disk','Thickdisk','Halo'],textcolors=[30,200,50,100,250],position=[-18.5,12.5]

device,/close
set_plot,'x'


print,'Computing masses for the stars in each component.....'

BD=starmass(filename,stars,inds,ms_units,grp,index,s)

;stop
end
