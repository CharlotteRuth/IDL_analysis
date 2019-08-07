; Relies on the output from mzr.pro

PRO galMS,filenames, halos = halos, outfile = outfile, key = key, psym  = psym, color = color, symsizes = symsizes, ctables = ctables, thicks = thicks, simStellarMass = simStellarMass, noObsStellarMass = noObsStellarMass, xrange = xrange, yrange = yrange, readfile = readfile,multiz = multiz,redshift = redshift,specific = specific

IF NOT keyword_set(multiz) THEN multiz = 1
n = fix(n_elements(filenames)/multiz)
zsolar = 0.0130215
IF keyword_set(color) THEN BEGIN
    white = 13
    black = 0
    IF NOT keyword_set(obsct) THEN obsct = 0;39
    loadct,obsct
    distinct_colors,n_colors = 12
ENDIF ELSE BEGIN
    black = 0
    white = 255
    IF NOT keyword_set(obsct) THEN obsct = 0 
    loadct,obsct
;    IF NOT keyword_set(ctables) THEN ctables = fltarr(n)
ENDELSE
IF keyword_set(outfile) THEN BEGIN
    l_charsize = 1.2 ;1 ;0.75
    fgcolor = black
    bgcolor = white
    xsize = 14;7
    ysize = 14;5
    formatplot,/outplot,thick = formatthick
ENDIF ELSE BEGIN
    l_charsize = !p.charsize ;1.0
    fgcolor = white
    bgcolor = black
    xsize = 500
    ysize = 500
    formatplot,thick = formatthick
ENDELSE
IF keyword_set(color) THEN BEGIN
;    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF color[0] EQ 1 THEN  colors = (fltarr(n) + 1)*fgcolor else colors = fltarr(n) + color
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(psym) THEN psym = fltarr(n) + 4 ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
    IF NOT keyword_set(obscolor) THEN obscolor = [2,6,11];fgcolor
    obssymsize = 1.25;2;.5
    obssym = 16
    obsthick = thicks[0];2
ENDIF ELSE BEGIN
    colors = (fltarr(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(psym) THEN psym = findgen(n) + 4;(findgen(n) + 2)*2   
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
    IF NOT keyword_set(obscolor) THEN obscolor = [150,150,150]
    obssymsize = 1.25;2;.5
    obssym = 16
    obsthick = thicks[0];2
ENDELSE
spawn,'hostname',hostname

IF (keyword_set(outfile)) THEN device,filename=outfile+'_' + strtrim(redshift,2)+'_galMS.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2,/encapsul ELSE window,0,xsize = xsize,ysize = ysize
IF keyword_set(specific) THEN ytitle = textoidl('log sSFR (Gyr^{-1})') ELSE ytitle = textoidl('log SFR/M_\odot (Gyr^{-1})')

;plot,[10,10],[-2,-1],xtitle = textoidl('log(M_{*} [M') + sunsymbol() + '])', ytitle = textoidl('log sSFR (Gyr^{-1})'),/nodata,xrange = xrange,yrange = yrange
FOR iz = 0, multiz - 1 DO BEGIN ;Counting through timesteps
plot,[10,10],[-2,-1],xtitle = textoidl('log(M_{*} [M') + sunsymbol() + '])', ytitle = ytitle,/nodata,xrange = xrange,yrange = yrange,title = 'z = ' + strtrim(redshift[iz],2)
gamma = 0
gamma1 = 0
gamma2 = 0
gamma3 = 0
alpha1 = 0
alpha2 = 0
alpha3 = 0
IF redshift[iz] GE 3 AND redshift[iz] LT 4 THEN BEGIN ;z = 3
    alpha = 1.02 ;Santini 2017, Table 2 ****
    beta = 1.37 + alpha*(-9.7)
ENDIF
IF redshift[iz] GE 2 AND redshift[iz] LT 3 THEN BEGIN ; z = 2.2
    alpha = 1.16 ;Santini 2017, Table 2; HST frontier Feilds (2 citations)
    beta = 1.22 + alpha*(-9.7)
    alpha1 = -0.311 + 1 ;Feng 2017, Table 2 CANDELS, 0 citations
    beta1 = -8.714 + (-0.311)*(-10)
    alpha2 = 3.44  ;Whitaker 2014, Table 1; 3D-HST, 218 citations ****
    beta2 = -19.99
    gamma2 = -0.13
    alpha3 = 0.46 ;Zahid+ 2012, Table 2, SDSS (82 citations)
    beta3 = 1.608 - 10*alpha3
;    alpha4 = ;Daddi+ 2007 (936 citations)
;    beta4 = 
ENDIF
IF redshift[iz] GE 0.5 AND redshift[iz] LT 1 THEN BEGIN ; z = 0.8
    alpha = -0.063 + 1 ;Feng 2017, Table 2, Table 2 CANDELS, 0 citations
    beta = -8.987 + (-0.063)*(-10) 
    alpha1 = 5.02 ;Whitaker 2014, Table 1; 3D-HST, 218 citations ****
    beta1 = -27.40
    gamma1 = -0.22
    alpha2 = 0.67; Zahid+ 2012, Table 2
    beta2 = 0.787 - 10*alpha2
    alpha3 = 0.9; Elbaz+ 2007, Eq 4 (849 citations)
    beta3 = alog10(7.2) - 10*alpha3
ENDIF
IF redshift[iz] EQ 0 THEN BEGIN ;Look up Elbaz+ 2007 ; z= 0
    alpha = -0.35 + 1 ;Salim 2007, Eq 11 ******
    beta = -0.35*(-10) - 9.82
    alpha1 = 0.71 ; Zahid+ 2012, Table 2
    beta1 = 0.317 - 10*alpha1
    alpha2 = 0.77; Elbaz+ 2007, Eq 5 (849 citations)
    beta2 = alog10(8.7) - 11*alpha2
ENDIF
;Show fit to star forming galaxies from Salim 2007 (Eq 11)
logmstar_arr = findgen(100)/100*6+6
;sSFR_arr = -0.35*(logmstar_arr - 10) - 9.83
SFR_arr = alpha*logmstar_arr + beta + gamma*logmstar_arr^2
sSFR_arr = (alpha - 1)*logmstar_arr + beta + gamma*logmstar_arr^2
IF keyword_set(specific) THEN oplot,logmstar_arr,sSFR_arr ELSE oplot,logmstar_arr,SFR_arr
IF alpha1 NE 0 THEN BEGIN
    SFR_arr1 = alpha1*logmstar_arr + beta1 + gamma1*logmstar_arr^2
    sSFR_arr1 = (alpha1 - 1)*logmstar_arr + beta1 + gamma1*logmstar_arr^2
    IF keyword_set(specific) THEN oplot,logmstar_arr,sSFR_arr1 ELSE oplot,logmstar_arr,SFR_arr1,linestyle = 1
ENDIF
IF alpha2 NE 0 THEN BEGIN
    SFR_arr2 = alpha2*logmstar_arr + beta2 + gamma2*logmstar_arr^2
    sSFR_arr2 = (alpha2 - 1)*logmstar_arr + beta2 + gamma2*logmstar_arr^2
    IF keyword_set(specific) THEN oplot,logmstar_arr,sSFR_arr2 ELSE oplot,logmstar_arr,SFR_arr2,linestyle = 2
ENDIF
IF alpha3 NE 0 THEN BEGIN
    SFR_arr3 = alpha3*logmstar_arr + beta3 + gamma3*logmstar_arr^2
    sSFR_arr3 = (alpha3 - 1)*logmstar_arr + beta3 + gamma3*logmstar_arr^2
    IF keyword_set(specific) THEN oplot,logmstar_arr,sSFR_arr3 ELSE oplot,logmstar_arr,SFR_arr3,linestyle = 3
ENDIF
FOR i = 0, n - 1 DO BEGIN ;Counting through simulations
;    loadct,ctables[i]
    IF keyword_set(onehalos) THEN metalfile = filenames[i + iz*n] + '.halo.' + strtrim(onehalos[i],2) ELSE  metalfile = filenames[i + iz*n]
    IF keyword_set(readfile) AND file_test(metalfile + '.metals.fits') THEN BEGIN 
        metals = mrdfits(metalfile + '.metals.fits',1,/silent) 
    ENDIF ELSE BEGIN
        IF keyword_set(onehalos) THEN IF keyword_set(noObsStellarMass) THEN metals = mzr(filenames[i + iz*n], onehalo = onehalos[i]) ELSE metals = mzr(filenames[i + iz*n],/obs, onehalo = onehalos[i]) ELSE $
          IF keyword_set(noObsStellarMass) THEN metals = mzr(filenames[i + iz*n]) ELSE metals = mzr(filenames[i + iz*n],/obs) 
        mwrfits,metals,metalfile + '.metals.fits',/create
    ENDELSE
    IF n_elements(tag_names(metals)) GT 34 THEN BEGIN
    IF keyword_set(halos) THEN BEGIN
       FOR j = 0, n_elements(halos[*,i,iz]) - 1 DO BEGIN
           IF halos[j,i,iz] NE 0 THEN BEGIN
               ind = where(metals.grp EQ halos[j,i,iz])
               IF keyword_set(specific) THEN weights = metals.smassO ELSE weights = fltarr(n_elements(metals.smassO)) + 1
               IF ind[0] NE -1 THEN BEGIN
                   IF NOT keyword_set(noObsStellarMass) THEN $
                       oplot,[alog10(metals[ind].smassO),alog10(metals[ind].smassO)],[alog10(metals[ind].sfr/weights[ind]),alog10(metals[ind].sfr/weights[ind])],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i]
                   IF keyword_set(simStellarMass) THEN $
                     oplot,[alog10(metals[ind].smass),alog10(metals[ind].smass)],[alog10(metals[ind].sfr/weights[ind]),alog10(metals[ind].sfr/weights[ind])],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i]
                   print,alog10(metals[ind].smass),alog10(metals[ind].sfr),alog10(metals[ind].sfr/metals[ind].smass)
;                   stop
               ENDIF
           ENDIF
       ENDFOR
   ENDIF ELSE BEGIN
       indsat = where(metals.sat EQ 1, complement = indnosat)
       IF (indsat[0] EQ -1) THEN indnosat = indgen(n_elements(metals))
       IF NOT keyword_set(noObsStellarMass) THEN BEGIN
           IF (indnosat[0] NE -1) THEN BEGIN
               IF n_elements(indnosat) EQ 1 THEN oplot,[alog10(metals[indnosat].smassO),alog10(metals[indnosat].smassO)],[metals[indnosat].sfr/weights[indnosat],metals[indnosat].sfr/weights[indnosat]],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] $
               ELSE oplot,alog10(metals[indnosat].smassO),metals[indnosat].sfr/weights[indnosat],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i]
               IF keyword_set(simStellarMass) THEN oplot,alog10(metals[indnosat].smass),metals[indnosat].sfr/weights[indnosat],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] $
               ELSE oplot,alog10(metals[indnosat].smassO),metals[indnosat].sfr/weights[indnosat],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i]
           ENDIF
       ENDIF
;       stop
   ENDELSE
ENDIF ELSE BEGIN
    print,'Missing SFR tag: ',filenames[i + iz*n]
    metals = mzr(filenames[i + iz*n],/obs) 
    mwrfits,metals,filenames[i + iz*n] + '.metals.fits',/create
ENDELSE
ENDFOR
stop
ENDFOR
END
