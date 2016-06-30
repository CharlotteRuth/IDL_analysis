FUNCTION atan_rotc,x,P
RETURN,P[0]*(2/!PI)*atan(X/P[1])
END

FUNCTION neg_atan_rotc,x,P
RETURN,P[1] - P[0]*(2/!PI)*atan(X/P[1])
END

pro vcirc_h603
halo = 1
msol_per_sysmass       = 1.84793e16
kpc_per_syslength      = 50000.

files = ['/astro/net/mondo-2/nbody/fabio/e11Gals/h603.cosmo50cmb.2304g2bwK.BUG/h603.cosmo50cmb.2304g2bwK.00328/h603.cosmo50cmb.2304g2bwK.00328.1.std',$
         '/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.00324/h603.cosmo50cmb.3072gs1MbwK.00324.1.std', $
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.00324/h603.cosmo50cmb.3072g14HBWK.00324.halo.1.std']

keys = ['Original','Metals','Metals + H2']

outfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/rotcurve.eps'

vcirc,files,msol_per_sysmass,kpc_per_syslength,keys = keys,/color,outfile = outfile;,/color
end

pro vrot, files, halos, dMsolUnit, dKpcUnit, outfile = outfile, keys = keys, color = color,thicks = thicks, linestyle = linestyle,yrange = yrange,maxdistance = maxdistance,vreturn = vreturn,type = type,multiframe = multiframe, label = label,ctables = ctables,verbose = verbose,nbins = nbins
; type: gas, dark, star,all

;*********************************** Units ****************************************
hubble = 73.0
cm_per_kpc = 3.08568021d21
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
vrotbin = fltarr(nbins)
vrotbin_std = fltarr(nbins)
formatplot,outplot = outfile

IF keyword_set(outfile) THEN  BEGIN
    fgcolor = 0 
    bgcolor = 255
   IF keyword_set(multiframe) THEN  BEGIN
        xsize = 10*n
        ysize = 12
        mxTitSize = 1.5
        mxTitOffset = 2
    ENDIF ELSE BEGIN
        xsize = 18
        ysize = 12
    ENDELSE
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    IF keyword_set(multiframe) THEN  BEGIN
        xsize = 400*n
        ysize = 475
        mxTitSize = 1.5
        mxTitOffset = 1
    ENDIF ELSE BEGIN
        xsize = 400
        ysize = 266
    ENDELSE
ENDELSE
IF keyword_set(color) THEN  BEGIN
    loadct,39
;    IF NOT keyword_set(ctables) THEN  ctables = [39,39,39]
    IF NOT keyword_set(ctables) THEN  ctables = [0,0,0]
    IF color[0] eq 1 THEN   colors = (findgen(n) + 1)*240/n else colors = color
    IF NOT keyword_set(thicks) THEN  thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyle) THEN  linestyle = fltarr(n) ;REVERSE(findgen(n)*2)
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN  ctables = fltarr(n)
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN  thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyle) THEN  linestyle = reverse(findgen(n)*2)   
ENDELSE
IF NOT keyword_set(yrange) THEN  yrange = [0,250]

vfinal = fltarr(n)
vfinalrot = fltarr(n)
vflat = fltarr(n)
vreturn = fltarr(n)

IF keyword_set(outfile) THEN      device,filename=outfile,/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,0
FOR ifile = 2, n - 1 DO BEGIN
    rtipsy,files[ifile] + '.halo.' + halos[ifile],h,g,d,s
    readarr,files[ifile] + '.halo.' + halos[ifile]+'.HI',h,HI,part = 'gas',/ascii
    IF (file_test(files[ifile]+".H2")) THEN  readarr,files[ifile] + '.halo.' + halos[ifile]+".H2",h,H2,/ascii ELSE H2 = fltarr(n_elements(HI))
    H2 = H2*2
    loadct,ctables[ifile]
    a = h.time
    IF (n_elements(dKpcUnit) GT 1) THEN  kpc_per_syslength = dKpcUnit[ifile]*a ELSE kpc_per_syslength = (dKpcUnit)[0]*a
    IF (n_elements(dMsolUnit) GT 1) THEN  msol_per_sysmass = dMsolUnit[ifile] ELSE msol_per_sysmass = dMsolUnit[0]
    h = 0.73
    vunit = 100.0* h * (kpc_per_syslength / 1000.0)/2.894405
    timeunit = sqrt((cm_per_kpc*dKpcUnit[ifile])^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year
    g.x = g.x*kpc_per_syslength
    g.y = g.y*kpc_per_syslength
    g.z = g.z*kpc_per_syslength
    g.mass = g.mass*msol_per_sysmass
    s.x = s.x*kpc_per_syslength
    s.y = s.y*kpc_per_syslength
    s.z = s.z*kpc_per_syslength
    s.mass = s.mass*msol_per_sysmass
    d.x = d.x*kpc_per_syslength
    d.y = d.y*kpc_per_syslength
    d.z = d.z*kpc_per_syslength
    d.mass = d.mass*msol_per_sysmass

    distancesg = sqrt(g.x^2.0 + g.y^2.0 + g.z^2.0)
    IF n_elements(s) GT 0 THEN  distancess = sqrt(s.x^2.0 + s.y^2.0 + s.z^2.0)
    distancesd = sqrt(d.x^2.0 + d.y^2.0 + d.z^2.0)
    IF n_elements(s) GT 0 THEN  distances = [distancesg,distancess,distancesd] 

;--------------------- 90% HI radius -------------
    HIhist = weighted_histogram(distancesg[where(distancesg LE 30)],weight = g[where(distancesg LE 30)].mass*HI[where(distancesg LE 30)],/cum,locations = xHIhist,nbins = 100)
    HIhist = HIhist/max(HIhist)
    r90HI = spline(HIhist,xHIhist,0.9)
    r80HI = spline(HIhist,xHIhist,0.8)

;---------------------- Rotational Velocity ---------------------------------
    g.vx = g.vx*vunit
    g.vy = g.vy*vunit
    g.vz = g.vz*vunit
    diskg = g[where((HI + H2) GT 0.1 AND abs(g.z) LT 3)]
    vel = [[diskg.vx],[diskg.vy]]
    tan = [[-1.0*diskg.y/sqrt(diskg.y*diskg.y + diskg.x*diskg.x)],[diskg.x/sqrt(diskg.y*diskg.y + diskg.x*diskg.x)]]
    tanvel = diskg.vx*tan[*,0] + diskg.vy*tan[*,1]
    diskgr = sqrt(diskg.x^2 +diskg.y^2)
   
;----------------------- Circular Velocity -----------------------------
    massg = g.mass
    IF n_elements(s) GT 0 THEN  masss = s.mass
    massd = d.mass
    IF n_elements(s) GT 0 THEN  mass = [massg,masss,massd] 

    FOR i = 0, nbins-1 DO BEGIN
        ind = where(distances LE vc_x[i])
        IF (ind[0] NE -1)THEN  tmass[i] = total(mass[ind])ELSE tmass[i] = 0

        ind = where(distancesd LE vc_x[i])
        IF (ind[0] NE -1)THEN  tmassd[i] = total(massd[ind])ELSE tmassd[i] = 0

        ind = where(distancesg LE vc_x[i])
        IF (ind[0] NE -1)THEN  tmassg[i] = total(massg[ind])ELSE tmassg[i] = 0

        ind = where(distancess LE vc_x[i])
        IF (ind[0] NE -1)THEN  tmasss[i] = total(masss[ind])ELSE tmasss[i] = 0
        
        ind = where(diskgr LE vc_x[i] + bind/2 AND diskgr GT vc_x[i] - bind/2)
        vrotbin[i] = mean(tanvel[ind])
        vrotbin_std[i] = stdev(tanvel[ind])
    ENDFOR
    vc_y = sqrt(tmass*gm_per_msol*grav/vc_x/cm_per_kpc)/1e5
    vc_yd = sqrt(tmassd*gm_per_msol*grav/vc_x/cm_per_kpc)/1e5
    vc_yg = sqrt(tmassg*gm_per_msol*grav/vc_x/cm_per_kpc)/1e5
    vc_ys = sqrt(tmasss*gm_per_msol*grav/vc_x/cm_per_kpc)/1e5
;    vfinal[ifile] = vc_y[nbins - 1]

;------------------------------ Plot ------------------------------
    loadct,39
    window,0
    plot,diskgr,tanvel,psym = 3,xrange = [0,maxdistance],xtitle = 'Radius [kpc]',ytitle = textoidl('Velocity [km/s]')
    oplot,vc_x,vrotbin,color = 200
    oploterror,vc_x,vrotbin,vrotbin_std,errcolor = 200,color = 200,/nohat
    oplot,vc_x,vc_y,color = 60

;----------------------------------- Prep fit ---------------------------------
    fit0 = [vc_y[nbins - 1],vc_x[nbins - 1]] ;100 km/s is tangental velocity, 50 kpc is the diameter it reaches it at
;outer edge of disk
;    temp = weighted_histogram(distancesg,weight = massg*(HI + H2),/cum,nbins = nbins,locations = rbins)
;    temp2 = min(abs(temp/max(temp) - 0.9),indmaxdisk)
;    maxdisk = rbins[indmaxdisk]
;minimum of disk
    eps = min(d.eps)*kpc_per_syslength
    mindisk = 2.0*eps
;determine if profile is peaked
    temp3 = max(vc_y,indmindisk)
    IF vc_x[indmindisk] LT 2 THEN BEGIN
        mindisk = vc_x[indmindisk] 
        print,'peaked'
    ENDIF

;-------------- Fit arctan function to rotational velocity
    indr = where(diskgr GT mindisk  AND diskgr LT r90HI)
    fitrot = mpfitfun('atan_rotc',diskgr[indr],tanvel[indr],fltarr(n_elements(indr)) + 1,fit0,/quiet)
    vfinalrot[ifile] = fitrot[0]
 
; if the curve is strongly peaked
    IF vrotbin[nbins - 1] LE max(vrotbin[0:where(vc_x EQ 5)] - 10) THEN  BEGIN
        fitrot = mpfitfun('atan_rotc',diskgr[indr],vrotbin[nbins - 1] + vrotbin[nbins - 1] - tanvel[indr],fltarr(n_elements(indr)) + 1,fit0,/quiet)
        vfinalrot[ifile] = vrotbin[nbins - 1] + vrotbin[nbins - 1] - fitrot[0]
    ENDIF

;   IF fitrot[0] LE max(vc_y - 30)
    IF vrotbin[nbins - 1] LE max(vrotbin[0:where(vc_x EQ 5)] - 10)  THEN oplot,vc_x,vrotbin[nbins - 1] + vrotbin[nbins - 1] - atan_rotc(vc_x,fitrot),linestyle = 2,color = 200 $
    ELSE oplot,vc_x,atan_rotc(vc_x,fitrot),linestyle = 2,color = 200

;------------------------------ Fit Arctan function to vcirc ------------------
;    indr = where(vc_x GT 2.0*eps)
    indr = where(vc_x GT mindisk  AND vc_x LT r90HI)
    fit = mpfitfun('atan_rotc',vc_x[indr],vc_y[indr],fltarr(n_elements(indr)) + 1,fit0,/quiet)
    vfinal[ifile] = fit[0]

    IF vc_y[nbins - 1] LE max(vc_y - 10) THEN  BEGIN
        fit = mpfitfun('atan_rotc',vc_x[indr],vc_y[nbins - 1] + vc_y[nbins - 1] - vc_y[indr],fltarr(n_elements(indr)) + 1,fit0,/quiet)
        vfinal[ifile] = vc_y[nbins - 1] + vc_y[nbins - 1] - fit[0]
    ENDIF

    IF fit[0] LE max(vc_y - 10) THEN  oplot,vc_x,vc_y[nbins - 1] + vc_y[nbins - 1] - atan_rotc(vc_x,fit),thick = thicks[iFile],linestyle = 2,color = 60 $
    ELSE oplot,vc_x,atan_rotc(vc_x,fit),linestyle = 2,color = 60

    IF keyword_set(verbose) THEN  BEGIN
        print,'Circular Velocity: ',vfinal[ifile]
        print,'Rotational Velocity: ',vfinalrot[ifile]
        print,'Final Velocity: ',vc_y[nbins - 1]
    ENDIF
    oplot,[2.0*eps,2.0*eps],[-300,300],linestyle = 2            
    oplot,[r90HI,r90HI],[-300,300],linestyle = 2
    oplot,[r80HI,r80HI],[-300,300],linestyle = 2

;--------------------- Get radius containing 80% of i-band flux ------------------
    r = eighty_flux(files[ifile],halo_str = halos[ifile]);,/verbose)
    temp = min(abs(vc_x - r),rind)
    oplot,[r,r],[-300,300],linestyle = 2 
    IF vc_y[nbins - 1] LE max(vc_y - 10) THEN print,'V Circ fit: ',vc_y[nbins - 1] + vc_y[nbins - 1] - atan_rotc(r,fit),vc_y[rind] ELSE print,'V Circ (80%): ',atan_rotc(r,fit),vc_y[rind]
    IF vc_y[nbins - 1] LE max(vrotbin[0:where(vc_x EQ 5)] - 10) THEN print,'V Rot fit: ',vrotbin[nbins - 1] + vrotbin[nbins - 1] - atan_rotc(r,fitrot),vrotbin[rind] ELSE print,'V Rot (80%):  ',atan_rotc(r,fitrot),vrotbin[rind]

    vreturn[iFile] = atan_rotc(r,fitrot)

    indr = where(diskgr GE r80HI AND diskgr LT r90HI)
    vflat[iFile] = mean(tanvel[indr])
    print,'V Flat: ',vflat[iFile]
    oplot,[r80HI,r90HI],[vflat[iFile],vflat[iFile]],color = 245,linestyle = 2
    stop
ENDFOR
IF keyword_set(keys) THEN  legend,keys,color = colors,linestyle = linestyle,thick = thicks,ctables = ctables,/right,/bottom,box = 0
IF keyword_set(outfile) THEN  device,/close; else stop
IF keyword_set(multiframe) THEN  multiplot,/reset
END
