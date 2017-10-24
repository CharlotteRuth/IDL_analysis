function r200,file,output_graph = output_graph, nbins = nbins, maxdistance = maxdistance

;NAME:
;   r200.pro
;
;PURPOSE:
;   This function will find r200 for a tipsy file, and print the circular velocity
;   and mass inclosed.
;   If the density does not dip below 200 X the critical density within
;   1000kpc (or some other determined maxdistance), -1 is returned
;
;CALLING SEQUENCE:
;   IDL>r = r200("myfile.std")
;
;REQUIRED INPUTS:
;   "myfile.std" = Name of a TIPSY snapshot file
;
;OPTIONAL INPUTS:
;
;   output_graph = If you would like the rotation curve and the
;   density profile saved to a file, the is the root of that file
;
;   nbins = is the number of bins the density is divided into to find
;   r200
;
;   maxdistance = the maximum distance in kpc that r200 is searched to
;
;OUTPUTS:
;   r200, a graph of the density profile, a graph of the rotation
;   curve
;
;EXTERNAL ROUTINES:
;   rtipsy

;EXTERNAL DATA:
;
;COMMENTS:
;
;REVISION HISTORY:
;


msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
IF (keyword_set(nbins)) THEN nbins = nbins ELSE nbins = 100
IF (keyword_set(maxdistance)) THEN maxdistance = maxdistance ELSE maxdistance = 1000

bind = maxdistance/dKpcUnit/nbins
bins = (findgen(nbins)+1)*bind
tmass = fltarr(nbins)
density = fltarr(nbins)

rtipsy, file, h, g, d, s
distancesg = SQRT(g.x^2.0 + g.y^2.0 + g.z^2.0)
distancess = SQRT(s.x^2.0 + s.y^2.0 + s.z^2.0)
distancesd = SQRT(d.x^2.0 + d.y^2.0 + d.z^2.0)
distances = [distancesg,distancess,distancesd]
massg = g.mass
masss = s.mass
massd = d.mass
mass = [massg,masss,massd]
maxd = MAX(distances)
mind = MIN(distances)

FOR i = 0, nbins-1 DO BEGIN
    ind = where(distances LE bins[i])
    if (ind[0] ne -1)then tmass[i] = TOTAL(mass[ind])else tmass[i] = 0
ENDFOR

density = tmass/(4.0/3.0*!PI*bins^3.0)
temp = MIN(ABS(density  - 200.0),i200)
IF(i200 LT 2) THEN i200 = 2
IF(i200 GT nbins -2) THEN i200 = nbins - 2
IF((Where(density lt 200))[0] ne -1)  then begin 
    radius = SPLINE(-1.*density[i200-2:i200+2],bins[i200-2:i200+2],-200.0)*dKpcUnit 
    totalmass = TOTAL(mass[where(distances*dkpcUnit LE radius)])
    v200 = SQRT(totalmass*dMsolUnit*msol*grav/radius/kpc)

    print,file
    print,"R200: ",radius,' kpc'
    print,"M200: ",totalmass*dMsolUnit,' solar masses'
    print,"V200: ",v200,' km/s'

    set_plot,'x'
    plot,bins*dkpcUnit,density,/ylog,yrange=[100,1e9],title = 'Baryonic Mass Profile',xtitle='Radius (kpc)',ytitle = 'Total Mass Density (* p_crit)',xrange = [1,maxdistance],/xlog
    oplot,[MIN(bins)*50000,MAX(BINS)*50000],[200,200],linestyle = 1

    IF (keyword_set(output_graph)) THEN BEGIN
        set_plot,'ps'
        device,filename = output_graph+'_dens.eps'
        plot,bins*dkpcUnit,density,/ylog,yrange=[100,1e9],title = 'Baryonic Mass Profile',xtitle='Radius (kpc)',ytitle = 'Total Mass Density (* p_crit)',xrange = [1,maxdistance],/xlog
        oplot,[MIN(bins)*50000,MAX(BINS)*50000],[200,200],linestyle = 1 
        device,/close
    ENDIF
    print,'Type ".c" to continue'
    stop

    velocitycurve = SQRT(tmass*dMsolUnit*msol*grav/bins/dkpcUnit/kpc)
    set_plot,'x'
    plot,bins*dKpcUnit,velocitycurve,xtitle = 'Radius',ytitle = 'Velocity (km/s)',title = 'Rotation Curve '+file,xrange = [0,150] ,yrange = [0,MAX(velocitycurve)*1.1]
    oplot,[radius,radius],[0,MAX(velocitycurve)],linestyle = 1
    oplot,[20./220.*v200/hubble,20./220.*v200/hubble],[0,MAX(velocitycurve)]
    oplot,[10./220.*v200/hubble,10./220.*v200/hubble],[0,MAX(velocitycurve)],linestyle = 2
    
    IF (keyword_set(output_graph)) THEN BEGIN
        set_plot,'ps'
        device,filename = output_graph+'_RC.pro'
        plot,bins*dKpcUnit,velocitycurve,xtitle = 'Radius',ytitle = 'Velocity (km/s)',title = 'Rotation Curve '+file,xrange = [0,150] ,yrange = [0,MAX(velocitycurve)*1.1]
        oplot,[radius,radius],[0,MAX(velocitycurve)],linestyle = 1
        oplot,[20./220.*v200/hubble,20./220.*v200/hubble],[0,MAX(velocitycurve)]
        oplot,[10./220.*v200/hubble,10./220.*v200/hubble],[0,MAX(velocitycurve)],linestyle = 2
        set_plot,'x'
    ENDIF

ENDIF ELSE BEGIN
    print,"Density does not go below 200 X crit"
    radius = -1
ENDELSE
RETURN,radius
END
