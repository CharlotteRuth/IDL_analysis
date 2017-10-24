; Finds the r200, the circular velocities at various radii, the
; stellar mass within that radii, and the
; baryonic mass included within that radii (stars + gas where T < 40,000K)
;comp_RC,'halo_match_list.txt','matches'

pro comp_RC, list, outputfile
msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
filebasis = '/astro/net/scratch2/fabio/Xruns/'

readcol, list, format = 'I,A,A,A,A', num, files1, files2, files3, files4
maxdistance = 150
nbins = 1000
bind = maxdistance/dKpcUnit/nbins
tmass = fltarr(nbins)
density = fltarr(nbins)
loadct,39

for ct = 6, n_elements(num)-1 do begin
    print,files1[ct]
    rtipsy, filebasis+files1[ct], h, g, d, s
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
    bins = (findgen(nbins)+1)*bind;*(maxd-mind)/nbins+mind

    FOR i = 0, nbins-1 DO BEGIN
        ind = where(distances LE bins[i])
        if (ind[0] ne -1)then tmass[i] = TOTAL(mass[ind])else tmass[i] = 0
    ENDFOR
    density = tmass/(4.0/3.0*!PI*(bins+bind)^3.0)
    temp = MIN(ABS(density  - 200.0),i200)
;    plot,bins*dkpcUnit,density,/ylog,yrange=[100,1e9],title = 'Baryonic Mass Profile',xtitle='Radius (kpc)',ytitle = 'Total Mass Density (* p_crit)',/xlog
;    oplot,[MIN(bins)*50000,MAX(BINS)*50000],[200,200],linestyle = 1
;stop
    velocitycurve = SQRT(tmass*dMsolUnit*msol*grav/bins/dkpcUnit/kpc)
    print,"about to plot"
    set_plot,'x'
    set_plot,'ps';
    device, filename=outputfile + STRTRIM(num[ct],2) + '.eps',/color,bits_per_pixel=8
;    halo_n = (STRSPLIT(files[ct],'.',/extract))[4]
    plot,bins*dKpcUnit,velocitycurve,xtitle = 'Radius',ytitle = 'Velocity (km/s)',title = 'Rotation Curve '+STRTRIM(num[ct],2),xrange = [0,20],yrange = [120,170]

    IF((Where(density lt 200))[0] ne -1)  then begin 
        r200 = SPLINE(-1.*density[i200-2:i200+2],bins[i200-2:i200+2],-200.0)*dKpcUnit 
        totalmass = TOTAL(mass[where(distances*dkpcUnit LE r200)])
        v200 = SQRT(totalmass*dMsolUnit*msol*grav/r200/kpc)
        oplot,[r200,r200],[0,MAX(velocitycurve)],linestyle = 1
        oplot,[20./220.*v200/hubble,20./220.*v200/hubble],[0,MAX(velocitycurve)]
        oplot,[10./220.*v200/hubble,10./220.*v200/hubble],[0,MAX(velocitycurve)],linestyle = 2
        
    ENDIF
    print,files2[ct]
    if (files2[ct] ne "no") then begin
        rtipsy, filebasis+files2[ct], h, g, d, s
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
;        bins = findgen(nbins)*(maxd-mind)/nbins/mind

        FOR i = 0, nbins-1 DO BEGIN
            ind = where(distances LE bins[i])
            if (ind[0] ne -1)then tmass[i] = TOTAL(mass[ind])else tmass[i] = 0
        ENDFOR 
        density = tmass/(4.0/3.0*!PI*bins^3.0)  
        
        velocitycurve = SQRT(tmass*dMsolUnit*msol*grav/bins/dkpcUnit/kpc)
        oplot,bins*dkpcUnit,velocitycurve,color = 50
    ENDIF
        
    print,files3[ct]
    if (files3[ct] ne "no") then begin
        rtipsy, filebasis+files3[ct], h, g, d, s
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
;        bins = findgen(nbins)*(maxd-mind)/nbins+mind
        
        FOR i = 0, nbins-1 DO BEGIN
            ind = where(distances LE bins[i])
            if (ind[0] ne -1) then tmass[i] = TOTAL(mass[ind])else tmass[i] = 0
        ENDFOR 

        
;    plot,bins*dkpcUnit,density,/ylog,yrange=[100,1e9],title = 'Baryonic Mass Profile',xtitle='Radius (kpc)',ytitle = 'Total Mass Density (* p_crit)',/xlog
;    oplot,[MIN(bins)*50000,MAX(BINS)*50000],[200,200],linestyle = 1
;stop
        velocitycurve = SQRT(tmass*dMsolUnit*msol*grav/bins/dkpcUnit/kpc)
        oplot,bins*dkpcUnit,velocitycurve,color = 150
    ENDIF

    print,files4[ct]
    if (files4[ct] ne "no") then begin
        rtipsy, filebasis+files4[ct], h, g, d, s
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
 ;       bins = findgen(nbins)*(maxd-mind)/nbins+mind
        FOR i = 0, nbins-1 DO BEGIN
            ind = where(distances LE bins[i])
            if (ind[0] ne -1) then tmass[i] = TOTAL(mass[ind])else tmass[i] = 0
        ENDFOR        

        density = tmass/(4.0/3.0*!PI*bins^3.0)  
        
;    plot,bins*dkpcUnit,density,/ylog,yrange=[100,1e9],title = 'Baryonic Mass Profile',xtitle='Radius (kpc)',ytitle = 'Total Mass Density (* p_crit)',/xlog
;    oplot,[MIN(bins)*50000,MAX(BINS)*50000],[200,200],linestyle = 1
;stop
        velocitycurve = SQRT(tmass*dMsolUnit*msol*grav/bins/dkpcUnit/kpc)
        oplot,bins*dkpcUnit,velocitycurve,color = 240
    ENDIF
    
    legend,['X2X2','X2X2nf','X3X3','X5X5'],color=[0,50,150,240],/right,linestyle=[0,0,0,0]
;    stop
    device,/close
endfor

end
