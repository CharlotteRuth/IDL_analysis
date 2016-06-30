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

pro vcirc, files, dMsolUnit, dKpcUnit, outfile = outfile, keys = keys, color = color,thicks = thicks, linestyle = linestyle,yrange = yrange,maxdistance = maxdistance,vfinal = vfinal,type = type,multiframe = multiframe, label = label,ctables = ctables,verbose = verbose,nbins = nbins
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
    IF NOT keyword_set(ctables) THEN  ctables = [0,0,0]
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN  thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyle) THEN  linestyle = reverse(findgen(n)*2)   
ENDELSE
IF NOT keyword_set(yrange) THEN  yrange = [0,250]

vfinal = fltarr(n)

IF keyword_set(outfile) THEN      device,filename=outfile,/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,0
FOR ifile = 0, n - 1 DO BEGIN
    rtipsy,files[ifile],h,g,d,s
    IF (file_test(files[ifile]+".H2")) THEN readarr,files[ifile]+'.HI',h,HI,part = 'gas',/ascii ELSE HI = fltarr(n_elements(g)) + 1
    IF (file_test(files[ifile]+".H2")) THEN  readarr,files[ifile]+".H2",h,H2,/ascii ELSE H2 = fltarr(n_elements(HI))
    H2 = H2*2
    loadct,ctables[ifile]
    a = h.time
    IF (n_elements(dKpcUnit) GT 1) THEN  kpc_per_syslength = dKpcUnit[ifile]*a ELSE kpc_per_syslength = (dKpcUnit)[0]*a
    IF (n_elements(dMsolUnit) GT 1) THEN  msol_per_sysmass = dMsolUnit[ifile] ELSE msol_per_sysmass = dMsolUnit[0]
    timeunit = sqrt((cm_per_kpc*kpc_per_syslength)^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year

    distancesg = sqrt(g.x^2.0 + g.y^2.0 + g.z^2.0)*kpc_per_syslength
    IF n_elements(s) GT 0 THEN  distancess = sqrt(s.x^2.0 + s.y^2.0 + s.z^2.0)*kpc_per_syslength
    distancesd = sqrt(d.x^2.0 + d.y^2.0 + d.z^2.0)*kpc_per_syslength

    IF NOT keyword_set(type) THEN  type = 'all'
    IF type NE 'all' THEN  BEGIN
        IF type EQ 'dark' THEN  distances = distancesd
        IF type EQ 'gas'  THEN  distances = distancesg
        IF type EQ 'star' THEN  distances = distancess
    ENDIF ELSE $
      IF n_elements(s) GT 0 THEN  distances = [distancesg,distancess,distancesd] ELSE distances = [distancesg,distancesd]    

    massg = g.mass*msol_per_sysmass
    IF n_elements(s) GT 0 THEN  masss = s.mass*msol_per_sysmass
    massd = d.mass*msol_per_sysmass
    IF keyword_set(type) AND type NE 'all' THEN  BEGIN
        IF type EQ 'dark' THEN  mass = massd
        IF type EQ 'gas'  THEN  mass = massg
        IF type EQ 'star' THEN  mass = masss
    ENDIF ELSE $
    IF n_elements(s) GT 0 THEN  mass = [massg,masss,massd] ELSE mass = [massg,massd]

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
;    vfinal[ifile] = vc_y[nbins - 1]

    fit0 = [100,50] ;100 km/s is tangental velocity, 50 kpc is the diameter it reaches it at
    temp = weighted_histogram(distancesg,weight = massg*(HI + H2),/cum,nbins = nbins,locations = rbins)
    temp2 = min(abs(temp/max(temp) - 0.9),indmaxdisk)
    maxdisk = rbins[indmaxdisk]

    eps = MIN(d.eps)*kpc_per_syslength
    mindisk = 2.0*eps
    temp3 = max(vc_y,indmindisk)
    IF vc_x[indmindisk] LT 2 THEN BEGIN
        mindisk = vc_x[indmindisk] 
        print,'peaked'
        stop
    ENDIF

;    indr = where(vc_x GT 2.0*eps)
    indr = where(vc_x GT mindisk  AND vc_x LT maxdisk)
    fit = mpfitfun('atan_rotc',vc_x[indr],vc_y[indr],fltarr(n_elements(indr)) + 1,fit0)
    vfinal[ifile] = fit[0]

    IF vc_y[nbins - 1] LE max(vc_y - 30) THEN  BEGIN
        fit = mpfitfun('atan_rotc',vc_x[indr],vc_y[nbins - 1] + vc_y[nbins - 1] - vc_y[indr],fltarr(n_elements(indr)) + 1,fit0)
        vfinal[ifile] = vc_y[nbins - 1] + vc_y[nbins - 1] - fit[0]
    ENDIF

    IF keyword_set(multiframe) THEN  BEGIN
        multiplot,[n,1],/square,mxtitle = 'Radius [kpc]',mxTitSize = mxTitSize, mxTitOffset = mxTitOffset
        plot,vc_x,vc_y,thick = thicks[iFile],linestyle = linestyle[iFile],xtitle = 'Radius [kpc]',ytitle = textoidl('V_{circ} [km/s]'),yrange = yrange[0],xrange = [0,maxdistance]
    ENDIF ELSE BEGIN
        IF ifile EQ 0 THEN  BEGIN
            plot,vc_x,vc_y,thick = thicks[iFile],linestyle = linestyle[iFile],xtitle = 'Radius [kpc]',ytitle = textoidl('V_{circ} [km/s]'),yrange = yrange,xrange = [0,maxdistance],title = label,/nodata
            eps = min(d.eps)*kpc_per_syslength
            oplot,[2.0*eps,2.0*eps],[0,300],linestyle = 1            
        ENDIF
        oplot,vc_x,vc_y,thick = thicks[iFile],linestyle = linestyle[iFile],color = colors[iFile]
        IF keyword_set(verbose) THEN  BEGIN
            IF fit[0] LE MAX(vc_y - 30) THEN  oplot,vc_x,vc_y[nbins - 1] + vc_y[nbins - 1] - atan_rotc(vc_x,fit),thick = thicks[iFile],linestyle = 1,color = colors[iFile] $
                                        ELSE oplot,vc_x,atan_rotc(vc_x,fit),thick = thicks[iFile],linestyle = 1,color = colors[iFile]
            oplot,[maxdisk,maxdisk],yrange,linestyle = 2
        ENDIF
    ENDELSE
    IF keyword_set(type) THEN  $
      IF type eq 'all' THEN  BEGIN
;            oplot,vc_x,vc_yd,thick = thicks[iFile],linestyle = 2,color = colors[iFile]
;        oplot,vc_x,vc_yg,thick = thicks[iFile],linestyle = 2,color = colors[iFile]
;        oplot,vc_x,vc_ys,thick = thicks[iFile],linestyle = 3,color = colors[iFile]
    ENDIF
    IF keyword_set(multiframe)AND i NE n - 1 THEN  multiplot

    IF keyword_set(verbose) THEN BEGIN
        print,vfinal[ifile], vc_y[nbins - 1]
        stop
    ENDIF

ENDFOR
IF keyword_set(keys) THEN  legend,keys,color = colors,linestyle = linestyle,thick = thicks,ctables = ctables,/right,/bottom,box = 0
IF keyword_set(outfile) THEN  device,/close else stop
IF keyword_set(multiframe) THEN  multiplot,/reset
END
