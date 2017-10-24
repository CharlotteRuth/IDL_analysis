; Finds the r200, the circular velocities at various radii, the
; stellar mass within that radii, and the
; baryonic mass included within that radii (stars + gas where T < 40,000K)
;radiiANDv,'Data_files/halo_listX2X2.dat','X2X2'

pro radiiANDv, list, outputfile
msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
filebase = "/astro/net/scratch2/christensen/Database/figures/"

readcol, list, format = 'A60', files
velocity = fltarr(n_elements(files),3)
radii = fltarr(n_elements(files),3)
stellarmass = fltarr(n_elements(files),3)
baryonicmass = fltarr(n_elements(files),3)
flag = fltarr(n_elements(files)) ;flags whether a halo ever decresses in density to 200 X crit
nbins = 500
maxdistance = 1000
bind = maxdistance/dKpcUnit/nbins
bins = (findgen(nbins)+1)*bind  ;*(maxd-mind)/nbins+mind
tmass = fltarr(nbins)
density = fltarr(nbins)
set_plot,'x'
loadct,39

close,4
close,5
close,6
close,7
openw,4,'Data_files/'+outputfile+'r200.dat'
openw,5,'Data_files/'+outputfile+'rcirc1.dat'
openw,6,'Data_files/'+outputfile+'rcirc2.dat'
openw,7,'Data_files/'+outputfile+'randv.dat'
printf,7,format='(13A13)','HaloID ','R200 ','V200 ','StellarM200 ','BaryonM200 ','R1 ','Vcirc1 ','StellarM1 ','BaryonM1 ','R2 ','Vcirc2 ','StellarM2 ','BaryonM2 '

for ct = 0, n_elements(files)-1 do begin
    print,files[ct]
    rtipsy, files[ct], h, g, d, s

    distancesg = SQRT(g.x^2.0 + g.y^2.0 + g.z^2.0)
    distancesd = SQRT(d.x^2.0 + d.y^2.0 + d.z^2.0)
    massg = g.mass
    massd = d.mass
    IF (h.nstar GT 0) THEN BEGIN
        distancess = SQRT(s.x^2.0 + s.y^2.0 + s.z^2.0)
        distances = [distancesg,distancess,distancesd]
        masss = s.mass
        mass = [massg,masss,massd]
        maxd = MAX(distances)
        mind = MIN(distances)

        FOR i = 0, nbins-1 DO BEGIN
            ind = where(distances LE bins[i])
            if (ind[0] ne -1)then tmass[i] = TOTAL(mass[ind])else tmass[i] = 0
        ENDFOR
        density = tmass/(4.0/3.0*!PI*bins^3.0)

        IF((Where(density lt 200))[0] ne -1)  then begin 
            cross = MIN(WHERE(density lt 200)) ;Finds the first point density becomes less than 200
            IF (cross ne 0) THEN BEGIN

                temp = MIN(ABS(density[0:cross]  - 200.0),i200)
                radii[ct,0] = SPLINE(-1.*density[i200-2:i200+2],bins[i200-2:i200+2],-200.0)*dKpcUnit 
                IF((WHERE(distances*dkpcUnit le radii[ct,0]))[0] ne -1) THEN totalmass = TOTAL(mass[where(distances*dkpcUnit LE radii[ct,0])]) ELSE BEGIN
                    totalmass = 0
                    stop
                    ENDELSE
                velocity[ct,0] = SQRT(totalmass*dMsolUnit*msol*grav/radii[ct,0]/kpc)
                IF(h.nstar GT 0 AND (WHERE(distancess*dkpcUnit le radii[ct,0]))[0] ne -1) THEN  stellarmass[ct,0] = TOTAL(s[WHERE(distancess*dkpcUnit le radii[ct,0])].mass)*dMsolUnit ELSE  stellarmass[ct,0] = 0
                IF((WHERE(distancesg*dKpcUnit le radii[ct,0] AND g.tempg LE 40000.))[0] ne -1) THEN baryonicmass[ct,0] = stellarmass[ct,0] + TOTAL(g[WHERE(distancesg*dkpcUnit le radii[ct,0]) AND g.tempg LE 40000.].mass) * dMsolUnit ELSE baryonicmass[ct,0] = stellarmass[ct,0]

                radii[ct,1] = 20./220.*velocity[ct,0]/hubble
                IF((WHERE(distances*dkpcUnit le radii[ct,1]))[0] ne -1) THEN totalmass = TOTAL(mass[where(distances*dkpcUnit LE radii[ct,1])]) ELSE totalmass = 0
                velocity[ct,1] = SQRT(totalmass*dMsolUnit*msol*grav/radii[ct,1]/kpc)
                IF((WHERE(distancess*dkpcUnit le radii[ct,1]))[0] ne -1) THEN stellarmass[ct,1] = TOTAL(s[WHERE(distancess*dKpcUnit le radii[ct,1])].mass)*dMsolUnit ELSE stellarmass[ct,1] = 0
                IF((WHERE(distancesg*dKpcUnit le radii[ct,1] AND g.tempg LE 40000.))[0] ne -1) THEN baryonicmass[ct,1] = stellarmass[ct,1] + TOTAL(g[WHERE(distancesg*dKpcUnit le radii[ct,1] AND g.tempg LE 40000.)].mass) * dMsolUnit ELSE baryonicmass[ct,1] = stellarmass[ct,1]

                radii[ct,2] = 10./220.*velocity[ct,0]/hubble
                IF((WHERE(distances*dkpcUnit le radii[ct,2]))[0] ne -1) THEN totalmass = TOTAL(mass[where(distances*dkpcUnit LE radii[ct,2])]) ELSE totalmass = 0
                velocity[ct,2] = SQRT(totalmass*dMsolUnit*msol*grav/radii[ct,2]/kpc)
                IF((WHERE(distancess*dkpcUnit le radii[ct,2]))[0] ne -1) THEN stellarmass[ct,2] = TOTAL(s[WHERE(distancess*dKpcUnit le radii[ct,2])].mass)*dMsolUnit ELSE stellarmass[ct,2] = 0
                IF((WHERE(distancesg*dKpcUnit le radii[ct,2] AND g.tempg LE 40000.))[0] ne -1) THEN baryonicmass[ct,2] = stellarmass[ct,2] + TOTAL(g[WHERE(distancesg*dKpcUnit le radii[ct,2] AND g.tempg LE 40000.)].mass) * dMsolUnit ELSE baryonicmass[ct,2] = stellarmass[ct,2]

                plot,bins*dkpcUnit,density,/ylog,yrange=[100,1e9],title = 'Baryonic Mass Profile',xtitle='Radius (kpc)',ytitle = 'Total Mass Density (* p_crit)',/xlog
                oplot,[MIN(bins)*dkpcUnit,MAX(BINS)*dkpcUnit],[200,200],linestyle = 1
                halo_n = (STRSPLIT(files[ct],'.',/extract))[4]
   ;             IF (halo_n lt 50) THEN stop

                velocitycurve = SQRT(tmass*dMsolUnit*msol*grav/bins/dkpcUnit/kpc)

                set_plot,'x'    
;                halo_n = (STRSPLIT(files[ct],'.',/extract))[4]
                plot,bins*dkpcUnit,velocitycurve,xtitle = 'Radius',ytitle = 'Velocity (km/s)',title = 'Rotation Curve '+ outputfile+' '+halo_n,xrange = [0,150]
                oplot,[radii[ct,2],radii[ct,2]],[0,MAX(velocitycurve)],linestyle = 1
                oplot,[radii[ct,1],radii[ct,1]],[0,MAX(velocitycurve)]
                oplot,[radii[ct,0],radii[ct,0]],[0,MAX(velocitycurve)],linestyle = 2

                set_plot,'ps'
                print,filebase + 'rot_curves/'+outputfile+'_'+halo_n+'_RC.eps'
                device, filename=filebase + outputfile+"_"+halo_n+'_RC.eps'
                plot,bins*dkpcUnit,velocitycurve,xtitle = 'Radius',ytitle = 'Velocity (km/s)',title = 'Rotation Curve '+ outputfile+ ' '+halo_n,xrange = [0,150]
                oplot,[radii[ct,2],radii[ct,2]],[0,MAX(velocitycurve)],linestyle = 1
                oplot,[radii[ct,1],radii[ct,1]],[0,MAX(velocitycurve)]
                oplot,[radii[ct,0],radii[ct,0]],[0,MAX(velocitycurve)],linestyle = 2
                device,/close
                set_plot,'x'


;    IF ((Where (density le 200.0))[0] eq -1) THEN velocity[ct,0] =
;    1.0
                print,'Halo: ',halo_n
                print,'V200:           ',velocity[ct,0],' Vcirc:       ',velocity[ct,1],' Vcirc2:       ',velocity[ct,2]
                print,'R200:           ',radii[ct,0],' Rcirc:       ',radii[ct,1],' Rcirc2:       ',radii[ct,2]
                print,'StellarMass200: ',stellarmass[ct,0],' StellarMass: ',stellarmass[ct,1],' StellarMass2: ',stellarmass[ct,2]
                print,'BaryonMass200:  ',baryonicmass[ct,0],' BaryonMass:  ',baryonicmass[ct,1],' BaryonMass2:  ',baryonicmass[ct,2]
 ;              IF (halo_n lt 50) THEN stop
                if (radii[ct,0] le 0) then stop
                printf,4,halo_n," ",radii[ct,0]
                printf,5,halo_n," ",radii[ct,1]
                printf,6,halo_n," ",radii[ct,2]
                halo_n = halo_n + " "
                print,FORMAT='(A,2F,2E,2F,2E,2F,2E)',halo_n,radii[ct,0],velocity[ct,0],stellarmass[ct,0],baryonicmass[ct,0],radii[ct,1],velocity[ct,1],stellarmass[ct,1],baryonicmass[ct,1],radii[ct,2],velocity[ct,2],stellarmass[ct,2],baryonicmass[ct,2]
                print,halo_n
                printf,7,FORMAT='(A,2F,2E,2F,2E,2F,2E)',halo_n,radii[ct,0],velocity[ct,0],stellarmass[ct,0],baryonicmass[ct,0],radii[ct,1],velocity[ct,1],stellarmass[ct,1],baryonicmass[ct,1],radii[ct,2],velocity[ct,2],stellarmass[ct,2],baryonicmass[ct,2]
                flag[ct] = 1
            ENDIF ELSE BEGIN
                print,"Donut Profile"
                flag[ct] = 0 ;Flags this halo so that it is not plotted
            ENDELSE
        ENDIF ELSE BEGIN
            print,"Contaminated Subhalo"
            flag[ct] = 0    ;Flags this halo so that it is not plotted
        ENDELSE
    ENDIF ELSE BEGIN
        print,"No Stars"
        flag[ct] = 0        ;Flags this halo so that it is not plotted
    ENDELSE        
;    stop
endfor
close,4
close,5
close,6
close,7



!Y.style = 1
!X.style = 1
magfit = findgen((24.0-15.0)/0.1)*0.1 -24.0
slope = -0.128412
velocityfit = slope*(magfit + 15.5) + 1.5
;Fit given in Eke, Navarro and Steimetz 2001

;Estimates the luminosity of the galaxy based on a mass to light ratio
;of 4 Msol/Lsol
set_plot,'x'
plot,velocityfit,magfit,title = 'TF -- M/L = 1.5',ytitle = 'MI - 5Log(h)', xtitle = 'Log(Vcirc)',xrange=[0.1,2.5],yrange = [-5,-24]
oplot,ALOG10(velocity[where(flag eq 1),0]),4.08-2.5*ALOG10(stellarmass[where(flag eq 1),0]/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,ALOG10(velocity[where(flag eq 1),1]),4.08-2.5*ALOG10(stellarmass[where(flag eq 1),1]/1.5) - 5*ALOG10(hubble),psym = 4,color = 100
oplot,ALOG10(velocity[where(flag eq 1),2]),4.08-2.5*ALOG10(stellarmass[where(flag eq 1),2]/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['R200','20/220*v200/h','10/220*v200/h'],color=[50,100,240],psym = [4,4,4]
set_plot, 'ps'
device, filename='../figures/stellar_mlTF'+outputfile+'.eps',/color,bits_per_pixel=8
plot,velocityfit,magfit,title = 'TF -- M/L  = 1.5',ytitle = 'MI - 5Log(h)', xtitle = 'Log(Vcirc)',xrange=[1.5,2.8],yrange = [-15,-24]
oplot,ALOG10(velocity[where(flag eq 1),0]),4.08-2.5*ALOG10(stellarmass[where(flag eq 1),0]/1.5) - 5*ALOG10(hubble),psym = 4,color = 50
oplot,ALOG10(velocity[where(flag eq 1),1]),4.08-2.5*ALOG10(stellarmass[where(flag eq 1),1]/1.5) - 5*ALOG10(hubble),psym = 4,color = 100
oplot,ALOG10(velocity[where(flag eq 1),2]),4.08-2.5*ALOG10(stellarmass[where(flag eq 1),2]/1.5) - 5*ALOG10(hubble),psym = 4,color = 240
legend,['R200','20/220*v200/h','10/220*v200/h'],color=[50,100,240],psym = [4,4,4]
device,/close
stop

set_plot,'x'
plot,[10,400],[7e5,1e12],title = 'Baryonic TF',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[10,400],yrange = [1e5,1e12]
oplot,velocity[where(flag eq 1),0],baryonicmass[where(flag eq 1),0],psym = 4,color = 50
oplot,velocity[where(flag eq 1),1],baryonicmass[where(flag eq 1),1],psym = 4,color = 100
oplot,velocity[where(flag eq 1),2],baryonicmass[where(flag eq 1),2],psym = 4,color = 240
legend,['R200','20/220*v200/h','10/220*v200/h'],color=[50,100,240],psym = [4,4,4]
set_plot, 'ps'
device, filename='../figures/baryonicTF'+outputfile+'.eps',/color,bits_per_pixel=8
plot,[10,400],[7e5,1e12],title = 'Baryonic TF',xtitle = 'Log(Vcirc)', ytitle = 'Baryonic Mass',/ylog,/xlog,xrange=[10,400],yrange = [1e5,1e12]
oplot,velocity[where(flag eq 1),0],baryonicmass[where(flag eq 1),0],psym = 4,color = 50
oplot,velocity[where(flag eq 1),1],baryonicmass[where(flag eq 1),1],psym = 4,color = 100
oplot,velocity[where(flag eq 1),2],baryonicmass[where(flag eq 1),2],psym = 4,color = 240
legend,['R200','20/220*v200/h','10/220*v200/h'],color=[50,100,240],psym = [4,4,4]
device,/close
stop

end
