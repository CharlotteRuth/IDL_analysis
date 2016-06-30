;Observational Data from Geha06
PRO tully_fisher_obs, broadband, velocities, gmass, key = key,ctables = ctables, color = color, obscolor = obscolor,outfile = outfile, filternums = filternums, symbols = symbols, symsizes = symsizes, thicks = thicks_psym, kcorrect = kcorrect, halo = halo

!Y.STYLE = 1
!X.STYLE = 1
IF NOT keyword_set(halo) THEN halo = strtrim(fltarr(n) + '1',2)
IF keyword_set(outfile) THEN BEGIN
    set_plot,'ps' 
    nbins=100.0
    linestyles = [0,2]
    !P.THICK = 3.5
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25
;    !p.font=0 
    l_charsize = 0.75
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    set_plot,'x'
    nbins=100.0
    linestyles = [0,2]
    !P.THICK = 1.5
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    l_charsize = 1.0
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    fgcolor = 255
    bgcolor = 0
ENDELSE
;!p.multi = [0,4,1]
!p.multi = 0
n = n_elements(broadband)
IF keyword_set(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  colors = (findgen(n_elements(broadband)) + 1)*240/n_elements(broadband) else colors = color
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n_elements(broadband)) + 6;2
    IF NOT keyword_set(symbols) THEN symbols = fltarr(n_elements(broadband)) + 4
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n_elements(broadband)) + 2
    obscolor = fgcolor
    obssymsize = 2
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n_elements(broadband)) + fgcolor
;    colors = (findgen(n_elements(broadband)) + 1)*10.0 + 5.0;  fltarr(n_elements(broadband)) + 5
;    thicks = (findgen(n_elements(broadband)) + 1)*6/n_elements(broadband) - 1
    IF NOT keyword_set(THICKS) THEN thicks = fltarr(n_elements(broadband)) + 6;2
    IF NOT keyword_set(symbols) THEN symbols = (findgen(n_elements(broadband))+2)*2
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n_elements(broadband)) + 2
    obscolor = 150
    obssymsize = 2
ENDELSE
IF not keyword_set(filternums) THEN filternums = intarr(n) + 13
IF n_elements(filternums) eq 1 THEN filternums = intarr(n) + filternums

Imag = fltarr(n_elements(broadband))
smass = fltarr(n_elements(broadband))
IF NOT keyword_set(kcorrect) THEN BEGIN
    Bmag = fltarr(n_elements(broadband))
    isdssmag = fltarr(n_elements(broadband))
    gmag = fltarr(n_elements(broadband))
    rmag = fltarr(n_elements(broadband))
    r_solarMag = 4.76           ; Blanton et al. 2003, ApJ, 592, 819
;r_zpmag = 24.80
    MLa = 0.306                 ;Bell 03, page 296
    MLb = 1.097
ENDIF 

FOR i = 0, n_elements(broadband) - 1 DO BEGIN
    IF keyword_set(kcorrect) THEN BEGIN
         mags = mrdfits(broadband[i] + '.amiga_vir.halos.star99_K_ab.Mv.fits',1)
         ind = where(fix(mags.id) EQ halo[i])
         mags = mags[ind]
         mags_bes = transpose([[mags.u],[mags.b],[mags.v],[mags.r],[mags.i]])
         mags_errs = fltarr(5,n_elements(mags.u))+0.02
         mags_ivar=1./mags_errs^2
         dmod = 19.4576         ;distance modulus for z = 0 in this model
         z = fltarr(n_elements(mags.u))
         mgy =(10.D)^(-(0.4D)*(mags_bes + dmod))
         mgy_ivar = mags_ivar/(0.4*alog(10.)*mgy)^2.
         kcorrect, mgy, mgy_ivar, z, kcorrect, mass = star, mtol = mtol_k, absmag = absmag_k, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']
         smass[i] = star
         Imag[i] = mags.I
    ENDIF ELSE BEGIN
        lums = mrdfits(broadband[i], filternums[i], header)
        In = where(STRCMP(strtrim(lums.filter,2),'I_Cousins.res') eq 1)
        Imag[i] = lums[In].AB_mag1
        Bn = where(STRCMP(strtrim(lums.filter,2),'B_Johnson.res') eq 1)
        Bmag[i] = lums[Bn].AB_mag1
        isdssn = where(STRCMP(strtrim(lums.filter,2),'i_SDSS.res') eq 1)
        isdssmag[i] = lums[isdssn].AB_mag1
        gn = where(STRCMP(strtrim(lums.filter,2),'g_SDSS.res') eq 1)
        gmag[i] = lums[gn].AB_mag1
        rn = where(STRCMP(strtrim(lums.filter,2),'r_SDSS.res') eq 1)
        rmag[i] = lums[rn].AB_mag1  
;    L_r_band = lums[rn].L_lambda_eff1 ;Watts/m
;    r_band_width = lums[rn].ewidth_lambda ;m
;    L_r = L_r_band*r_band_width ;Watts
        g_r = gmag - rmag
        MLfit = MLa + MLb*g_r
        smass[i] = MLfit[i]*10^((rmag[i] - r_solarMag)/(-2.5))      
    ENDELSE
ENDFOR
bmass = smass + gmass

;IF (keyword_set(outfile)) THEN device,filename=outfile+'_tfgeha.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 10,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 712,ysize = 392;,xsize = 950,ysize = 500;392

IF(1) THEN BEGIN ;From Geha 2006 (see marla.pro)
    h = .70
    spawn,'hostname',hostname
    IF hostname EQ 'ozma' THEN datadir = '/home/christensen/Code/Datafiles/HIcubes/' $
    ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN datadir = '/nobackupp2/crchrist/MolecH/' $
    ELSE datadir = '/astro/users/christensen/code/Datafiles/HIcubes/'
    sfi = mrdfits(datadir + 'SFI_data.fits',1)
    q=where(sfi.W20_err/sfi.W20 le 0.01 and sfi.absmag_err le 0.2 and $
        sfi.sini ge 0.7,n)
    sfi = sfi[q]
    verh = mrdfits(datadir + 'verh_data.fits',1)
    q=where(verh.w20_err lt 10)
    verh=verh[q]
    matt = mrdfits(datadir + 'matthews_data.fits',1)
    q=where(matt.sini ge 0.6)
    matt=matt[q]
    mcg = mrdfits(datadir + 'mcg_data.fits',1)
    readcol, datadir + 'marla_data.txt', galnum,mw20, mw20err, mi, mierr,/silent
    readcol, datadir + 'geha06.dat', galnumAll, a1,a2,a3,d1,d2,d3,distance,mr,m_r,g_r,mu,reff, reff2,ba,Mstar, MHI, sigmaHI, mw202orig, mw202, mw20err2,/silent 
    match,LONG(galnum),LONG(galnumAll),temp,selectedGeha
    MHI = 10^MHI
    Mstar = 10^Mstar
    sfidata = replicate( $
                     create_struct('w20', 0.0, $
                                   'w20_err', 0.0,$
                                   'i', 0.0, $
                                   'i_err', 0.0), $
                     n_elements(sfi))

    verhdata = replicate($
                      create_struct('w20', 0.0, $
                                    'w20_err', 0.0,$
                                    'i', 0.0, $
                                    'i_err', 0.0), $
                      n_elements(verh))

    mattdata = replicate( $
                      create_struct('w20', 0.0, $
                                    'w20_err', 0.0,$
                                    'i', 0.0, $
                                    'i_err', 0.0), $
                    n_elements(matt))

    mcgdata = replicate( $
                     create_struct('w20', 0.0, $
                                   'w20_err', 0.0,$
                                   'i', 0.0, $
                                   'i_err', 0.0), $
                     n_elements(mcg))

    marladata = replicate( $
                           create_struct('w20', 0.0, $
                                         'w20_err', 0.0, $
                                         'i', 0.0, $
                                         'i_err', 0.0), $
                           n_elements(mw20))

    sfidata.w20 = sfi.w20/2.
    sfidata.w20_err = sfi.w20_err
    sfidata.i = sfi.absmagi

    verhdata.w20 = verh.w20/2.
    verhdata.w20_err = verh.w20_err
    verhdata.i = verh.absmag_bvri[3]

    mattdata.w20 = matt.w20/2.
    mattdata.w20_err = matt.w20_err
    mattdata.i = matt.absmagi

    mcgdata.w20 = mcg.w20
    mcgdata.w20_err = mcg.w20_err
    mcgdata.i = mcg.absmagi

    marladata.w20 = mw20 - 6
    marladata.w20_err = mw20err
    marladata.i = mi
    marladata.i_err = mierr

    colobs = 100
ENDIF ELSE BEGIN 
;From Giovanelli, Haynes et al 1997 (To high of TF) 
;Form: y = M - 5 log h, where M is the I band mag and h = H_0/(100 km s^-1 Mpc^-1)
;      x = log W - 2.5, where W is twice the maximum rotational velocity, from single-dish 21cm line data
; y(x) = c0 + c1 x + c2 x^2 (eq 20)
    h = 73.0/100.0
    logW = findgen(100)/100.0 + 1.9
    x = logW - 2.5
    
    if 1 THEN BEGIN             ;Bivariate case, Table 4
        a = -21.00              ; in, alphaLW = -0.5, a_bi   
        b = -7.76               ; b_bi
        sigma = 0.34
        y = a + b*x 
        ymax = a + b*x + sigma
        ymin = a + b*x - sigma
    ENDIF ELSE BEGIN            ;Quadratic fit, Table 5
        c0 = -21.03              
        c1 = -6.84
        c2 = 1.93
        sigma = 0.33
        y = c0 + c1*x + c2 * x^2
        ymax = c0 + c1*x + c2 * x^2 + sigma
        ymin = c0 + c1*x + c2 * x^2 - sigma
    ENDELSE
ENDELSE

IF (keyword_set(outfile)) THEN device,filename=outfile+'_stfgeha.eps',/color,bits_per_pixel= 8,/times,xsize = 12,ysize = 20,xoffset =  2,yoffset =  2  ELSE $
  window,0,xsize = 392,ysize = 500 ;xsize = 712,ysize = 392;,xsize = 950,ysize = 500;392

;------------------------------------ Stellar Mass vs Velocity ----------------------------------------
plot,[velocities[0],velocities[0]],[bmass[0],bmass[0]],psym = symbols[0],xrange = [9,105],yrange = [6.9,9.2],/xlog,/nodata, ystyle=1, xstyle=1, $
      xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
      ytitle='log Stellar Mass [M' + sunsymbol() + textoidl('h_{70}^{-2}]'),symsize = symsizes[0];/ylog,yrange = [8e6,1.5e9
;coeff = [-0.14,0.43]
;y = [5e6, 5e9]
;x = 10^(coeff[0] + coeff[1]*alog10(y/1e8))*50;   ( alog10(x/50.0) - coeff[0])/coeff[1]
x = [20.5,97]
y = [2e7,1e9]
coeff = line(alog10(x), alog10(y))
x = [8, 500]
y = coeff[0] + coeff[1] * alog10(x)
oplot,x,y
hierror = 0.06
oploterror, mw202[selectedGeha] - 6, alog10(Mstar[selectedGeha]), mw20err2[selectedGeha], hierror[selectedGeha],psym=symcat(16),errcolor=obscolor, /nohat  , color=obscolor;,symsize = obssymsize

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL
IF keyword_set(COLOR) THEN BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO oplot,[velocities[i],velocities[i]],alog10([smass[i],smass[i]]),color = colors[i],psym = 4,SYMSIZE = symsizes[i], thick = thicks[i]
    IF keyword_set(key) THEN legend, key, color = colors, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),charsize = l_charsize,psym = 4*(fltarr([n_elements(broadband)]) + 1),symsize = symsizes
ENDIF ELSE BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO oplot,[velocities[i],velocities[i]],alog10([smass[i],smass[i]]),psym = symbols[i],SYMSIZE = symsizes[i], thick = thicks[i]
    IF keyword_set(key) THEN legend, key, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),psym = symbols, charsize = l_charsize,symsize = symsizes
ENDELSE

IF (keyword_set(outfile)) THEN BEGIN
   device,/close
   device,filename=outfile+'_gtfgeha.eps',/color,bits_per_pixel= 8,/times,xsize = 12,ysize = 20,xoffset =  2,yoffset =  2  
ENDIF ELSE BEGIN
    stop
    window,0,xsize = 392,ysize = 500;,xsize = 712,ysize = 392;,xsize = 950,ysize = 500 ;392
ENDELSE

;------------------------------------ Gass Mass vs Velocity ----------------------------------------
plot,[velocities[0],velocities[0]],[bmass[0],bmass[0]],psym = symbols[0],xrange = [9,105],yrange = [6.9,9.2],/xlog,/nodata, ystyle=1, xstyle=1, $
      xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
      ytitle='Log HI Mass [M' + sunsymbol() + textoidl('h_{70}^{-2}]'),symsize = symsizes[0];/ylog,yrange = [8e6,1.5e9]
x = [10,100]
y = [alog10(2e7),alog10(1e9)]
coeff = line(alog10(x), y)
x = [8, 500]
y = coeff[0] + coeff[1] * alog10(x)
oplot,x,y
;coeff = [0.37,-0.09]
;y = [5e6, 5e9]
;x = 10^(coeff[0] + coeff[1]*alog10(y/1e8))*50;   ( alog10(x/50.0) - coeff[0])/coeff[1]
;oplot,x,y
hierror_coeff = [2.05340,0.656304]
hierror = 10^(hierror_coeff[0] + hierror_coeff[1]*alog10(Mstar))
;oploterror, mw202, MHI*1.4, mw20err2, hierror,psym=3, color=100,/lobar,/nohat
;oploterror, mw202[selectedGeha] - 6, MHI[selectedGeha]*1.4, mw20err2[selectedGeha], hierror[selectedGeha],psym=3, errcolor=obscolor,/nohat  
oploterror, mw202[selectedGeha] - 6, alog10(MHI[selectedGeha]), mw20err2[selectedGeha], sigmaHI[selectedGeha],psym=symcat(16),errcolor=obscolor,/nohat, color=obscolor;, symsize = obssymsize
; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL
IF keyword_set(COLOR) THEN BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO oplot,[velocities[i],velocities[i]],alog10([gmass[i],gmass[i]]),psym = symbols[i],SYMSIZE = symsizes[i],color = colors[i], thick = thicks[i]
    IF keyword_set(key) THEN legend, key, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),charsize = l_charsize,psym = symbols,symsize = symsizes, color = colors
ENDIF ELSE BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO oplot,[velocities[i],velocities[i]],alog10([gmass[i],gmass[i]]),psym = symbols[i],SYMSIZE = symsizes[i], thick = thicks[i]
    IF keyword_set(key) THEN legend, key, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),charsize = l_charsize,psym = symbols,symsize = symsizes
ENDELSE

IF (keyword_set(outfile)) THEN BEGIN
   device,/close
   device,filename=outfile+'_btfgeha.eps',/color,bits_per_pixel= 8,/times,xsize = 12,ysize = 20,xoffset =  2,yoffset =  2  
ENDIF ELSE BEGIN
    stop
    window,0,xsize = 392,ysize = 500;,xsize = 712,ysize = 392;,xsize = 950,ysize = 500 ;392
ENDELSE

;------------------------------------ Baryonic TF ----------------------------------------
plot,[velocities[0],velocities[0]],[bmass[0],bmass[0]],psym = symbols[0],xrange = [9,105],yrange = [6.9,9.2],/xlog,/nodata, ystyle=1, xstyle=1, $
      xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
      ytitle='Log Baryonic Mass [M' + sunsymbol() + textoidl('h_{70}^{-2}]'), symsize = symsizes[0];,yrange = [alog10(8e6),alog10(1.5e9)
x = [9,88]
y = [alog10(5e7),alog10(1e9)]
coeff = line(alog10(x), y)
x = [8, 500]
y = coeff[0] + coeff[1] * alog10(x)
oplot,x,y
;coeff = [-0.19,0.27]
;y = [5e6, 5e9]
;x = 10^(coeff[0] + coeff[1]*alog10(y/1e8))*50;   ( alog10(x/50.0) - coeff[0])/coeff[1]
;oplot,x,y
hierror_coeff = [-0.698970,1.00000]
hierror = 10^(hierror_coeff[0] + hierror_coeff[1]*alog10(Mstar))
;oploterror, mw202[selectedGeha] - 6, Mstar[selectedGeha] + MHI[selectedGeha]*1.4, mw20err2[selectedGeha], hierror[selectedGeha] + hierror[selectedGeha]*10*sigmaHI,psym=3, errcolor=obscolor,/nohat  
oploterror, mw202[selectedGeha] - 6, alog10(Mstar[selectedGeha] + MHI[selectedGeha]*1.4), mw20err2[selectedGeha],sigmaHI[selectedGeha],psym=symcat(16), errcolor=obscolor,/nohat, color=obscolor;, symsize = obssymsize;
IF keyword_set(COLOR) THEN BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO oplot,[velocities[i],velocities[i]],alog10([bmass[i],bmass[i]]),psym = symbols[i],SYMSIZE = symsizes[i],color = colors[i], thick = thicks[i]
;    legend key, color = colors, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),charsize = l_charsize,psym = 4*(fltarr([n_elements(broadband)]) + 1)
ENDIF ELSE BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO oplot,[velocities[i],velocities[i]],alog10([bmass[i],bmass[i]]),psym = symbols[i],SYMSIZE = symsizes[i], thick = thicks[i]
;    legend, key, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),psym = symbols, charsize = l_charsize
ENDELSE
IF (keyword_set(outfile)) THEN BEGIN
   device,/close
   device,filename=outfile+'_tfgeha.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  
ENDIF ELSE BEGIN
    stop
    window,0,xsize = 712,ysize = 392;,xsize = 950,ysize = 500 ;392
ENDELSE

;------------------------------------ Baryonic TF AND I-band Mag vs Velocity ----------------------------------------
!p.multi = [0,2,1]
IF 1 THEN BEGIN  ;Geha
    plot, marladata.w20, marladata.i, /xlog, yrange=[-13.8,-17.1], $
      xrange=[9, 100], /nodata, ystyle=1, xstyle=1, $
      xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
      ytitle='M' + textoidl('_I') + ' - 5' + $
      textoidl(' log_{10} h_{70}')
    oplot,marladata.w20, marladata.i, psym=symcat(16) , color=obscolor;,symsize = obssymsize
;    oplot, mcgdata.w20, mcgdata.i, psym=2, color=obscolor 
;    oplot, mattdata.w20, mattdata.i, psym=2, color=obscolor 
;    oplot, sfidata.w20, sfidata.i, psym=2, color=obscolor 
;    oplot, verhdata.w20, verhdata.i, psym=2, color=obscolor 
    oploterror, marladata.w20, marladata.i, marladata.w20_err, marladata.i_err,psym=3, errcolor=obscolor,/nohat ;   ,symsize = obssymsize    
    x = [30,300]
    y = [-14.5, -24]
    coeff = line(alog10(x), y)
    x = [10, 500]
    y = coeff[0] + coeff[1] * alog10(x)
    oplot,x,y
ENDIF ELSE BEGIN                ;Giovanelli
    plot,logW - alog10(2),y + 5*alog10(h),xtitle = textoidl('log_{10}W'),ytitle = 'M - 5 log h',xrange = [1.6,2.6],yrange = [-17.5,-25] ;yrange = [-15,-24],xrange = [1.9,2.9]
    oplot,logW + 5*alog10(h),ymax,linestyle = 1
    oplot,logW + 5*alog10(h),ymin,linestyle = 1
ENDELSE
; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL
IF keyword_set(COLOR) THEN BEGIN 
    FOR i = 0, n_elements(broadband) - 1 DO BEGIN
        oplot,[velocities[i],velocities[i]],[Imag[i] - 5.0*alog10(h),Imag[i] - 5.0*alog10(h)], $
          psym = symbols[i],SYMSIZE = symsizes[i],color = colors[i], thick = thicks[i]
    ENDFOR
;    IF keyword_set(key) THEN legend, key, color = colors, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),charsize = l_charsize,psym = 4*(fltarr([n_elements(broadband)]) + 1)
ENDIF ELSE BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO BEGIN
        oplot,[velocities[i],velocities[i]],[Imag[i] - 5.0*alog10(h),Imag[i] - 5.0*alog10(h)], $
          psym = symbols[i],SYMSIZE = symsizes[i], thick = thicks[i]
    ENDFOR
;    IF keyword_set(key) THEN legend, key, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),psym = symbols, charsize = l_charsize
ENDELSE

plot,[velocities[0],velocities[0]],[bmass[0],bmass[0]],psym = symbols[0],xrange = [9,105],yrange = [7.7,9.2],/xlog,/nodata, ystyle=1, xstyle=1, $
      xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
      ytitle='log Baryonic Mass [M' + sunsymbol() + textoidl(' h_{70}^{-2}]');,yrange = [alog10(8e6),alog10(1.5e)9
x = [9,88]
y = [alog10(5e7),alog10(1e9)]
coeff = line(alog10(x), y)
x = [8, 500]
y = coeff[0] + coeff[1] * alog10(x)
oplot,x,y
;coeff = [-0.19,0.27]
;y = [5e6, 5e9]
;x = 10^(coeff[0] + coeff[1]*alog10(y/1e8))*50;   ( alog10(x/50.0) - coeff[0])/coeff[1]
;oplot,x,y
hierror_coeff = [-0.698970,1.00000]
hierror = 10^(hierror_coeff[0] + hierror_coeff[1]*alog10(Mstar))
;oploterror, mw202[selectedGeha] - 6, Mstar[selectedGeha] + MHI[selectedGeha]*1.4, mw20err2[selectedGeha], hierror[selectedGeha] + hierror[selectedGeha]*10*sigmaHI,psym=3, errcolor=obscolor,/nohat  
oploterror, mw202[selectedGeha] - 6, alog10(Mstar[selectedGeha] + MHI[selectedGeha]*1.4), mw20err2[selectedGeha],sigmaHI[selectedGeha],psym=symcat(16), errcolor=obscolor,/nohat, color=obscolor; ,symsize = obssymsize
IF keyword_set(color) THEN BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO oplot,[velocities[i],velocities[i]],alog10([bmass[i],bmass[i]]),psym = symbols[i],SYMSIZE = symsizes[i],color = colors[i], thick = thicks[i]
    IF keyword_set(key) THEN legend, [key,"Geha et al. 2006"], /top, /left, linestyle = [fltarr([n_elements(broadband)]),1],psym  = [symbols,16],charsize = l_charsize, color = [colors,obscolor],thick = [thicks,1];, SYMSIZE = [symsizes,1]
ENDIF ELSE BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO oplot,[velocities[i],velocities[i]],alog10([bmass[i],bmass[i]]),psym = symbols[i],SYMSIZE = symsizes[i], thick = thicks[i]
    IF keyword_set(key) THEN legend, [key,"Geha et al. 2006"], /top, /left, linestyle = [fltarr([n_elements(broadband)]),1],psym = [symbols,16], charsize = l_charsize,  color = [colors,obscolor],thick = [thicks,1];, SYMSIZE = [symsizes,1],
ENDELSE

IF (keyword_set(outfile)) THEN BEGIN
   device,/close
   device,filename=outfile+'_tfgehaMassiv.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  
ENDIF ELSE BEGIN
    stop
    window,0,xsize = 392,ysize = 500 ;392
ENDELSE

!p.multi = 0
IF 1 THEN BEGIN  ;Geha
    plot, marladata.w20, marladata.i, /xlog, yrange=[-13.5,-24.1], $
      xrange=[9, 500], /nodata, ystyle=1, xstyle=1, $
      xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
      ytitle='M' + textoidl('_I') + ' - 5' + $
      textoidl('log h_{70}');, symsize = obssymsize
    oplot,marladata.w20, marladata.i, psym=symcat(16) , color=obscolor;, symsize = obssymsize ;open circle
    oplot, mcgdata.w20, mcgdata.i, psym=symcat(14), color=obscolor;, symsize = obssymsize ;filled diamond
    oplot, mattdata.w20, mattdata.i, psym=symcat(15), color=obscolor;, symsize = obssymsize ;filled square
    oplot, sfidata.w20, sfidata.i, psym=symcat(17), color=obscolor;, symsize = obssymsize ;filled triangle
    oplot, verhdata.w20, verhdata.i, psym=symcat(46), color=obscolor;, symsize = obssymsize ;filled star
    oploterror, marladata.w20, marladata.i, marladata.w20_err, marladata.i_err,psym=3, errcolor=obscolor,/nohat;, symsize = obssymsize             
    x = [30,300]
    y = [-14.5, -24]
    coeff = line(alog10(x), y)
    x = [10, 500]
    y = coeff[0] + coeff[1] * alog10(x)
    oplot,x,y
ENDIF ELSE BEGIN                ;Giovanelli
    plot,logW - alog10(2),y + 5*alog10(h),xtitle = textoidl('log_{10}W'),ytitle = 'M - 5 log h',xrange = [1.6,2.6],yrange = [-17.5,-25] ;yrange = [-15,-24],xrange = [1.9,2.9]
    oplot,logW + 5*alog10(h),ymax,linestyle = 1
    oplot,logW + 5*alog10(h),ymin,linestyle = 1
ENDELSE
; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL
IF keyword_set(COLOR) THEN BEGIN 
    FOR i = 0, n_elements(broadband) - 1 DO BEGIN
        oplot,[velocities[i],velocities[i]],[Imag[i] - 5.0*alog10(h),Imag[i] - 5.0*alog10(h)], $
          psym = symbols[i],SYMSIZE = symsizes[i],color = colors[i], thick = thicks[i]
    ENDFOR
    IF keyword_set(key) THEN legend, key, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),psym = symbols, charsize = l_charsize, SYMSIZE = symsizes, color = colors
ENDIF ELSE BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO BEGIN
        oplot,[velocities[i],velocities[i]],[Imag[i] - 5.0*alog10(h),Imag[i] - 5.0*alog10(h)], $
          psym = symbols[i],SYMSIZE = symsizes[i], thick = thicks[i]
    ENDFOR
    IF keyword_set(key) THEN legend, key, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),psym = symbols, charsize = l_charsize,SYMSIZE = symsizes
ENDELSE
IF (keyword_set(outfile)) THEN device,/close ELSE stop

END
