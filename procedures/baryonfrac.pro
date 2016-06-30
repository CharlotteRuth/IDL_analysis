PRO baryonfrac, files, kcorrect = kcorrect, imf = imf, outplot = outplot, color = color, pysm = psym, key = key, xrange = xrange, yrange = yrange, symsize = symsize, title = title, obscolor = obscolor, obssym = obssym, ctables = ctables

IF NOT keyword_set(imf) THEN imf = 1
IF keyword_set(outplot) THEN fgcolor = 0 ELSE fgcolor = 255

IF keyword_set(color) THEN BEGIN
    IF NOT keyword_set(bw) THEN loadct,39 ELSE loadct,0
    IF color[0] EQ 1 THEN  color  = (findgen(n) + 1)*240/n ELSE color = color
    IF NOT keyword_set(psym) THEN psym = fltarr(n) + 4
    IF NOT keyword_set(obscolor) THEN obscolor = 0 ;fgcolor
    IF NOT keyword_set(obssym) THEN obssym = 4
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(obscolor) THEN obscolor = 100
    color = (findgen(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(n_elements(broadband)) + 5
    IF NOT keyword_set(psym) THEN  psym = (findgen(n)+2)*2
    IF NOT keyword_set(obssym) THEN obssym = 4
ENDELSE
IF NOT keyword_set(symsize) THEN symsize = 1.5
IF n_elements(symsize) eq 1 THEN symsize = fltarr(n_elements(files)) + symsize 
IF NOT keyword_set(xrange) THEN xrange = [9.5,13]
IF NOT keyword_set(yrange) THEN yrange = [1e-4,1]


FOR j = 0, n_elements(files) -1 DO BEGIN
    file = files[j]
    print,file
 
    fmt = 'I,X,X,X,X,F,X,F,F,F,X,X,X,X,X,X,X,X,X,A,A'
    readcol,file + '.amiga.stat',F=fmt,haloid,vir,gas,star,dark,cont,sat,/silent
    IF ((where(cont eq 'clean'))[0] eq -1) THEN BEGIN
        fmt = 'I,X,X,X,X,F,X,F,F,F,X,X,X,X,X,X,X,X,X,X,A,A'
        readcol,file + '.amiga.stat',F=fmt,haloid,vir,gas,star,dark,cont,sat,/silent
    ENDIF
    IF ((where(sat eq 'no'))[0] eq -1) THEN BEGIN
        fmt = 'I,X,X,X,X,F,X,F,F,F,F,X,F,X,X,X,X,X,X,A'
        readcol,file + '.amiga.stat',F=fmt,haloid,vir,gas,star,dark,cont,/silent
    ENDIF
    IF keyword_set(sat) THEN BEGIN
        central = where(sat EQ 'no')
        print,n_elements(sat),n_elements(central)
        haloid = haloid[central]
        vir = vir[central]
        gas = gas[central]
        star = star[central]
        dark = dark[central]
        cont = cont[central]
        sat = sat[central]
    ENDIF
    IF (keyword_set(kcorrect)) THEN BEGIN
        IF IMF eq 1 THEN mags = mrdfits(file + '.amiga_r200.halos.star99_K_ab.Mv.fits',1) ELSE mags = mrdfits(file + '.amiga_vir.halos.star99_MS_ab.Mv.fits',1)
        stop
        haloid2 = haloid[hstar]
        vir = vir[hstar]
        dark = dark[hstar]
        star = star[hstar]
        mags = mags[hstar]

        mags_bes = transpose([[mags.u],[mags.b],[mags.v],[mags.r],[mags.i]])
        mags_errs = fltarr(5,n_elements(mags.u))+0.02
        mags_ivar=1./mags_errs^2
        dmod = 19.4576       ;distance modulus for z = 0 in this model
        z = fltarr(n_elements(mags.u))
        mgy =(10.D)^(-(0.4D)*(mags_bes + dmod))
        mgy_ivar = mags_ivar/(0.4*alog(10.)*mgy)^2.
        kcorrect, mgy, mgy_ivar, z, kcorrect, mass = star, mtol = mtol_k, absmag = absmag_k, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']
    ENDIF ELSE BEGIN
;        hstar = where(cont eq 'clean')
;        hstar = where(cont eq 'clean' AND velmax gt veldisp/2.0)
;        hstar = where(velmax gt veldisp/2.0)
;        IF ((where(sat eq 'no'))[0] ne -1) then hstar = where(sat eq 'no' AND gas ne 0) else hstar = where(gas ne 0) 
;        IF (size(haloids))[0] eq 2 THEN match,haloids[*,j],haloid,temp,hstar ELSE match,haloids,haloid,temp,hstar
        IF (size(haloids))[0] gt 1 THEN temphaloids = haloids[*,j] ELSE temphaloids = haloids
        hstar = fltarr(n_elements(temphaloids))
        FOR i = 0, n_elements(temphaloids) - 1 DO BEGIN
            k = where(temphaloids[i] eq haloid)
            if k ne -1 THEN hstar[i] = k ELSE hstar[i] = 0
        ENDFOR
        stop
        haloid = haloid[hstar]
        star = star[hstar]
        vir = vir[hstar]
        dark = dark[hstar]
        gas = gas[hstar]
    ENDELSE

END
