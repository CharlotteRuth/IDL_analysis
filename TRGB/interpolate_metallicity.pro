; For an input stellar photometry structure, interpolate on isochrones to 
; return the metallicity of each star (assuming zero photometric errors!).
;
; stars that fall outside the region traced by isochrones default to Z=-1
;
; assumes input stars are a structure returned from READDOLPHOT
; isochrone file is assumed to contain ONLY RGB stars.
; requires fixed version of HIST_ND.PRO from: turtle.as.arizona.edu/idl/hist_nd.pro
; if requested return a vector corresponding to the inferred initial stellar mass
;
; Example:
;
;  @compilecodes.pro
;  file='/local/jd/Angst/Phot/KDG73_ap1_fs2_frc1.bi2.fits'
;  isofile = '/net/faculty-1/jd/Angst/Isochrones/Girardi07/rgb_isochrones_12Gyr.dat'
;  stars = readdolphot(file)
;  z = interpolate_metallicity(stars,isochrones=isofile,chatty=1)

FUNCTION interpolate_metallicity, stars, isochrones=isochrones, m_TRGB=m_TRGB, massvec=massvec, filt=filt, chatty=chatty, noplot=noplot, overplotiso=overplotiso

IF (NOT(keyword_set(filt))) THEN filt='bi2'
IF (NOT(keyword_set(isochrones))) THEN BEGIN
    isochrones='/net/faculty-1/jd/Angst/Isochrones/Girardi07/rgb_isochrones_12Gyr.dat'
ENDIF
IF (NOT(keyword_set(m_TRGB))) THEN m_TRGB=24.0
IF (NOT(keyword_set(chatty))) THEN chatty=0

; set the distance modulus to shift isochrones by

distmod = m_TRGB - (-4.)

; read in the isochrone file -- NOTE: Columns do _not_ correspond to Girardi default!
;                               (Processed Girardi default by removing narrow band columns,
;                               keeping only lines between RGBb & RGBt and adding column
;                               of the metallicity (using emacs rectangle-replace))
; note: readcol maxes out at 25 vectors!
READCOL, isochrones, z, age, M_init, M_act, logL_over_Lo, logTe, logG, mbol, $
         F435W, F475W, F550M, F555W, F606W, F625W, F658N, F660N, F775W, F814W, F850LP,    $
         C_over_O, M_hec, period, pmode, logMdot, int_IMF,           $
         FORMAT='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,X,D,D,D,I,D,D', silent=NOT(chatty)

; set default filters to use for interpolation
IF ((filt EQ 'bi2') OR (filt EQ 'bi')) THEN BEGIN
    x = F475W - F814W
    y = F814W + distmod
ENDIF

IF ((filt EQ 'vi2') OR (filt EQ 'vi')) THEN BEGIN
    x = F606W - F814W
    y = F814W + distmod
ENDIF

; triangulate isochrones in color-magnitude space
TRIANGULATE, x, y, objtriangles, boundaryvec, CONNECTIVITY=c

; make a default grid onto which to interpolate metallicities

xmin = min(x)
xmax = max(x)
dx = 0.005
;nx = (xmax - xmin)/dx + 1.0
nx = (xmax - xmin)/dx 
nx = fix(nx)
xgrid = xmin + findgen(nx)*dx

ymin = min(y)
ymax = max(y)
dy = 0.025
;ny = (ymax - ymin)/dy + 1.0
ny = (ymax - ymin)/dy
ny = fix(ny)
ygrid = ymin + findgen(ny)*dy

; make regularized grid of metallicity and initial stellar mass

zmap = TRIGRID(x, y, z, objtriangles, XOUT=xgrid, YOUT=ygrid, MISSING=-1)

massmap = TRIGRID(x, y, M_init, objtriangles, XOUT=xgrid, YOUT=ygrid, MISSING=-1)

; bin stars into the same regularized grid, so that correct metallicity can be assigned

starsx = stars.mag1_acs - stars.mag2_acs
starsy = stars.mag2_acs
xminhist = min(xgrid) - dx/2.
yminhist = min(ygrid) - dy/2.
xmaxhist = max(xgrid) + dx
ymaxhist = max(ygrid) + dy
hist=hist_nd([1.#starsx, 1.#starsy],[dx,dy],min=[xminhist,yminhist],max=[xmaxhist,ymaxhist],reverse_indices=ri)
nxhist = (xmaxhist - xminhist)/dx
nyhist = (ymaxhist - yminhist)/dy
nxhist = fix(nxhist+0.1)
nyhist = fix(nyhist+0.1)
print,nx,nxhist,ny,nyhist

; loop through reverse indices and assign each star to the proper metallicity bin

zvec = -1.0 + 0.*stars.x    ; initialize metallicity vector
mvec = -1.0 + 0.*stars.x    ; initialize metallicity vector
FOR j=0L,nyhist-2 DO BEGIN
    FOR i=0L,nxhist-2 DO BEGIN
        ind = i + (nxhist+1)*j
        IF (ri[ind] NE ri[ind+1]) THEN BEGIN
            starind = ri[ri[ind]:ri[ind+1]-1]
            zvec[starind]=zmap[i,j]
            mvec[starind]=massmap[i,j]
        ENDIF
    ENDFOR
ENDFOR

; plot data and check if all is sensible

loadct,39

IF (keyword_set(massvec) AND NOT(keyword_set(noplot))) THEN !P.MULTI=[0,2,2]

; plot metallicities
IF (NOT(keyword_set(noplot))) THEN BEGIN
    plot, starsx, starsy, psym=3, /ynozero, xrange=[-1.0,5], yrange=[28,20], /xstyle, /ystyle
ENDIF
uniqz = z[uniq(z)]  ; get unique isochrone metallicity values
nz = n_elements(uniqz)
startcolor = 25.
colorvec = startcolor + (255.-startcolor) * findgen(nz) / (float(nz-1))
FOR i=0,nz-2 DO BEGIN
    foo = where((zvec GT uniqz[i]) AND (zvec LE uniqz[i+1]),npts)
    IF ((keyword_set(chatty)) AND (npts GT 0)) THEN print,uniqz[i],uniqz[i+1],npts
    IF ((npts GT 1) AND (NOT(keyword_set(noplot)))) THEN oplot, starsx[foo], starsy[foo], color=fix(colorvec[i]), psym=2, symsize=0.5
ENDFOR
;
IF (keyword_set(overplotiso)) THEN BEGIN
    FOR i=0,nz-1 DO BEGIN
        goo = where(z EQ uniqz[i],npts)
        IF ((npts GT 1) AND (NOT(keyword_set(noplot)))) THEN oplot, x[goo], y[goo], thick=2
    ENDFOR
ENDIF   

medz = median(zvec[where((stars.mag2_acs - distmod) LT -3)])
print,'Median Metallicity: ',medz

; plot masses and return mass values if requested

IF (keyword_set(massvec)) THEN BEGIN
    IF (NOT(keyword_set(noplot))) THEN BEGIN
        plot, starsx, starsy, psym=3, /ynozero, xrange=[-1.0,5], yrange=[28,20], /xstyle, /ystyle
    ENDIF
    nm = 50.
    mmin = min(M_init)
    mmax = max(M_init)
    print,'Mass Range: ',mmin, mmax
    mstep = mmin + (mmax-mmin)*findgen(nm)/(float(nm-1))
    startcolor = 25.
    colorvec = startcolor + (255.-startcolor) * findgen(nm) / (float(nm-1))
    FOR i=0,nm-2 DO BEGIN
        foo = where((mvec GT mstep[i]) AND (mvec LE mstep[i+1]),npts)
        IF (keyword_set(chatty)) THEN print,mstep[i],mstep[i+1],npts
        IF ((npts GT 1) AND (NOT(keyword_set(noplot)))) THEN oplot, starsx[foo], starsy[foo], color=fix(colorvec[i]), psym=2, symsize=0.5
    ENDFOR

    IF (NOT(keyword_set(noplot))) THEN !P.MULTI=[0,1,1]

    massvec=mvec    ; set up return vector

ENDIF

RETURN,zvec

END
