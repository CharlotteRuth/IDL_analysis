pro rotcurve_hi, files, dMsolUnit, dKpcUnit, outfile = outfile, keys = keys, color = color,thicks = thicks, linestyle = linestyle,yrange = yrange,vfinal = vfinal,type = type,label = label,ctables = ctables,verbose = verbose,nbins = nbins,maxdistance = maxdistance
; type: gas, dark, star,all

;*********************************** Units ****************************************
hubble = 73.0
cm_per_kpc = 3.08568021d21
km_per_kpc = 3.08567758d16
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
sec_per_year = 31556926
grav = 6.67e-8

n = n_elements(files)
IF NOT keyword_set(nbins) THEN nbins = 100.0
IF NOT keyword_set(maxdistance) THEN  maxdistance = 8.0
bind = 1.0*(maxdistance)/nbins
vc_x = (findgen(nbins)+1)*bind 
tmass = fltarr(nbins)   
tmassd= tmass
tmassg = tmass
tmasss = tmass
formatplot,outplot = outfile

IF keyword_set(outfile) THEN  BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 400
    ysize = 266
ENDELSE
IF keyword_set(color) THEN  BEGIN
    loadct,39
    IF NOT keyword_set(ctables) THEN  ctables = [0,0,0]
    IF color[0] eq 1 THEN   colors = (findgen(n) + 1)*240/n else colors = color
    IF NOT keyword_set(thicks) THEN  thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyle) THEN  linestyle = fltarr(n) ;REVERSE(findgen(n)*2)
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN  ctables = [0,0,0]
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN  thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyle) THEN  linestyle = reverse(findgen(n)*2)   
ENDELSE
IF NOT keyword_set(yrange) THEN  yrange = [0,250]

vfinal = fltarr(n)

FOR ifile = 0, n - 1 DO BEGIN
    loadct,ctables[ifile]
    IF keyword_set(outfile) THEN device,filename=files[ifile] + '.HIrotcurve.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,0

    rtipsy,files[ifile],h,g,d,s
    readarr,files[ifile]+'.HI',h,HI,part = 'gas',/ascii

    a = h.time
    IF (n_elements(dKpcUnit) GT 1) THEN  kpc_per_syslength = dKpcUnit[ifile]*a ELSE kpc_per_syslength = (dKpcUnit)[0]*a
    IF (n_elements(dMsolUnit) GT 1) THEN  msol_per_sysmass = dMsolUnit[ifile] ELSE msol_per_sysmass = dMsolUnit[0]
    timeunit = sqrt((cm_per_kpc*kpc_per_syslength)^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year

    distancesg = sqrt(g.x^2.0 + g.y^2.0 + g.z^2.0)*kpc_per_syslength
    IF n_elements(s) GT 0 THEN  distancess = sqrt(s.x^2.0 + s.y^2.0 + s.z^2.0)*kpc_per_syslength
    distancesd = sqrt(d.x^2.0 + d.y^2.0 + d.z^2.0)*kpc_per_syslength
    IF n_elements(s) GT 0 THEN  distances = [distancesg,distancess,distancesd] ELSE distances = [distancesg,distancesd]    

    massg = g.mass*msol_per_sysmass
    IF n_elements(s) GT 0 THEN  masss = s.mass*msol_per_sysmass
    massd = d.mass*msol_per_sysmass
    IF n_elements(s) GT 0 THEN  mass = [massg,masss,massd] ELSE mass = [massg,massd]

    vr = sqrt(g.vx*g.vx + g.vy*g.vy)/timeunit/sec_per_year*kpc_per_syslength*km_per_kpc
 
    FOR i = 0, nbins-1 DO BEGIN
        ind = where(distances LE vc_x[i])
        IF (ind[0] NE -1)THEN  tmass[i] = total(mass[ind])ELSE tmass[i] = 0

        ind = where(distancesd LE vc_x[i])
        IF (ind[0] NE -1)THEN  tmassd[i] = total(massd[ind])ELSE tmassd[i] = 0

        ind = where(distancesg LE vc_x[i])
        IF (ind[0] NE -1)THEN  tmassg[i] = total(massg[ind])ELSE tmassg[i] = 0

        ind = where(distancess LE vc_x[i])
        IF (ind[0] NE -1)THEN  tmasss[i] = total(masss[ind])ELSE tmasss[i] = 0
    ENDFOR
    vc_y = sqrt(tmass*gm_per_msol*grav/vc_x/cm_per_kpc)/1e5
    vc_yd = sqrt(tmassd*gm_per_msol*grav/vc_x/cm_per_kpc)/1e5
    vc_yg = sqrt(tmassg*gm_per_msol*grav/vc_x/cm_per_kpc)/1e5
    vc_ys = sqrt(tmasss*gm_per_msol*grav/vc_x/cm_per_kpc)/1e5

    temp = weighted_histogram(distancesg,weight = massg*(HI),/cum,nbins = nbins,locations = rbins)
    temp = weighted_histogram(distancess,weight = masss,/cum,nbins = nbins,locations = rbins)
    IF keyword_set(verbose) THEN BEGIN
        plot,rbins,temp
        stop
    ENDIF
    temp2 = min(abs(temp/max(temp) - 0.85),indmaxdisk) ;0.9
    disk90 = rbins[indmaxdisk]
    temp2 = min(abs(temp/max(temp) - 0.75),indmaxdisk) ;0.8
    disk80 = rbins[indmaxdisk]
    eps = min(d.eps)*kpc_per_syslength
    ind_outerdisk = where(distancesg GE disk80 AND distancesg LE disk90)
    vfinal[ifile] = total(vr[ind_outerdisk]*HI[ind_outerdisk]*g[ind_outerdisk].mass)/total(HI[ind_outerdisk]*g[ind_outerdisk].mass)

 ;   IF ifile EQ 0 THEN  BEGIN
        plot,vc_x,vc_y,thick = thicks[iFile],linestyle = linestyle[iFile],xtitle = 'Radius [kpc]',ytitle = textoidl('V_{circ} [km/s]'),yrange = yrange,xrange = [0,disk90 + 5],title = label,/nodata
        eps = min(d.eps)*kpc_per_syslength
        oplot,[2.0*eps,2.0*eps],[0,300],linestyle = 1            
 ;   ENDIF
    oplot,distancesg,vr,psym = 3,color = 100
    oplot,vc_x,vc_y,thick = thicks[iFile],linestyle = linestyle[iFile],color = colors[iFile]
    oplot,[disk90,disk90],yrange,linestyle = 2
    oplot,[disk80,disk80],yrange,linestyle = 2
    oplot,[0,100],[vfinal[ifile],vfinal[ifile]],linestyle = 2
     print,'Average Outer Disk Velocity: ',vfinal[ifile]

    IF keyword_set(verbose) THEN BEGIN
         stop
    ENDIF

ENDFOR
IF keyword_set(keys) THEN  legend,keys,color = colors,linestyle = linestyle,thick = thicks,ctables = ctables,/right,/bottom,box = 0
IF keyword_set(outfile) THEN  device,/close else stop
IF keyword_set(multiframe) THEN  multiplot,/reset
END
