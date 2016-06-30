;Observational Data from Geha06
PRO tully_fisher_obs, broadband, velocities, gmass, key = key,ctables = ctables, color = color, obscolor = obscolor,obssymsize = obssymsize, obspsym = obspsym, outfile = outfile, filternums = filternums, symbols = symbols, symsizes = symsizes, thicks = thicks, do_kcorrect = do_kcorrect, halo = halo,obsct = obsct,formatthick = formatthick,pizagno = pizagno

formatplot,outplot = outfile,thick = formatthick
;!Y.STYLE = 1
;!X.STYLE = 1
IF NOT keyword_set(halo) THEN halo = strtrim(fltarr(n) + '1',2)
IF keyword_set(outfile) THEN BEGIN
    l_charsize = !p.charsize;1.0;0.75
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    l_charsize = !p.charsize;1.0
    fgcolor = 255
    bgcolor = 0
ENDELSE

nfiles = n_elements(broadband)
IF keyword_set(color) THEN BEGIN
    IF NOT keyword_set(obsct) THEN obsct = 0;39
    IF color[0] EQ 1 THEN  colors = fltarr(n_elements(broadband)) + fgcolor else colors = fltarr(n_elements(broadband)) + color;olors = (findgen(nfiles) + 1)*240/nfiles
    IF NOT keyword_set(thicks) THEN thicks = fltarr(nfiles) + 6;2
    IF NOT keyword_set(symbols) THEN symbols = fltarr(nfiles) + 4
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(nfiles) + 2
    IF NOT keyword_set(obscolor) THEN obscolor = fgcolor
    IF NOT keyword_set(obssymsize) THEN obssymsize = 1.25
    IF NOT keyword_set(obspsym) THEN obspsym = 2
    IF NOT keyword_set(ctables) THEN ctables = fltarr(nfiles) + 39
ENDIF ELSE BEGIN
    IF NOT keyword_set(obsct) THEN obsct = 0    
    colors = fltarr(nfiles) + fgcolor
;    colors = (findgen(nfiles) + 1)*10.0 + 5.0;  fltarr(nfiles) + 5
;    thicks = (findgen(nfiles) + 1)*6/nfiles - 1
    IF NOT keyword_set(THICKS) THEN thicks = fltarr(nfiles) + 6;2
    IF NOT keyword_set(symbols) THEN symbols = (findgen(nfiles)+2)*2
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(nfiles) + 2
    IF NOT keyword_set(obscolor) THEN obscolor = 150
    IF NOT keyword_set(obssymsize) THEN obssymsize = 1.25
    IF NOT keyword_set(obspsym) THEN obspsym = 2
    IF NOT keyword_set(ctables) THEN ctables = fltarr(nfiles)
ENDELSE
IF NOT keyword_set(filternums) THEN filternums = intarr(nfiles) + 13
IF n_elements(filternums) EQ 1 THEN filternums = intarr(nfiles) + filternums

Imag = fltarr(nfiles)
smass = fltarr(nfiles)
isdssmag = fltarr(nfiles)
IF NOT keyword_set(do_kcorrect) THEN BEGIN
    Bmag = fltarr(nfiles)
    gmag = fltarr(nfiles)
    rmag = fltarr(nfiles)
    r_solarMag = 4.76           ; Blanton et al. 2003, ApJ, 592, 819
;r_zpmag = 24.80
    MLa = 0.306                 ;Bell 03, page 296
    MLb = 1.097
ENDIF 

FOR i = 0, nfiles - 1 DO BEGIN
    IF keyword_set(do_kcorrect) OR NOT file_test(broadband[i] + '.' + strtrim(halo[i],2) + '/broadband.fits') THEN BEGIN
        print,'kcorrect: ',broadband[i],halo[i]
        mags = mrdfits(broadband[i] + '.amiga_r200.halos.star99_K_ab.Mv.fits',1)
        ind = where(fix(mags.id) EQ halo[i])
        mags = mags[ind]
        mags_bes = transpose([[mags.u],[mags.b],[mags.v],[mags.r],[mags.i]])
        mags_errs = fltarr(5,n_elements(mags.u))+0.02
        mags_ivar=1./mags_errs^2
        dmod = 19.4576          ;distance modulus for z = 0 in this model
        z = fltarr(n_elements(mags.u))
        mgy =(10.D)^(-(0.4D)*(mags_bes + dmod))
        mgy_ivar = mags_ivar/(0.4*alog(10.)*mgy)^2.
        kcorrect, mgy, mgy_ivar, z, kcorrect, mass = star, mtol = mtol_k, absmag = absmag_k, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']
        smass[i] = star
        Imag[i] = mags.I
        isdssmag[i] = mags.iprime
    ENDIF ELSE BEGIN
        print,'Sunrise',broadband[i],halo[i]
        lums = mrdfits(broadband[i] + '.' + strtrim(halo[i],2) + '/broadband.fits', filternums[i], header)
        In = where(strcmp(strtrim(lums.filter,2),'I_Cousins.res') eq 1)
        Imag[i] = lums[In].AB_mag1
        Bn = where(strcmp(strtrim(lums.filter,2),'B_Johnson.res') eq 1)
        Bmag[i] = lums[Bn].AB_mag1
        isdssn = where(strcmp(strtrim(lums.filter,2),'i_SDSS.res') eq 1)
        isdssmag[i] = lums[isdssn].AB_mag1
        gn = where(strcmp(strtrim(lums.filter,2),'g_SDSS.res') eq 1)
        gmag[i] = lums[gn].AB_mag1
        rn = where(strcmp(strtrim(lums.filter,2),'r_SDSS.res') eq 1)
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
    IF NOT keyword_set(pizagno) THEN BEGIN
        sfi = mrdfits(datadir + 'SFI_data.fits',1)
        q=where(sfi.W20_err/sfi.W20 le 0.01 and sfi.absmag_err le 0.2 and $
                sfi.sini ge 0.7,n)
        sfi = sfi[q]
        verh = mrdfits(datadir + 'verh_data.fits',1) ;Verheijen 2001
        q=where(verh.w20_err lt 10)
        verh=verh[q]
        matt = mrdfits(datadir + 'matthews_data.fits',1) ;Matthews, van Driel, Gallagher 1998
        q=where(matt.sini ge 0.6)
        matt=matt[q]
        mcg = mrdfits(datadir + 'mcg_data.fits',1) ;McGaugh 2000
    ENDIF ELSE BEGIN
        readcol,datadir + 'PizagnoData.txt',galname,mi_p,mierr_p,vcirc80_p,vcirc80err_p,format = '(A19,D,D,D,D)'
    ENDELSE
    readcol, datadir + 'marla_data.txt', galnum,mw20, mw20err, mi, mierr,/silent
    readcol, datadir + 'geha06.dat', galnumAll, a1,a2,a3,d1,d2,d3,distance,mr,m_r,g_r,mu,reff, reff2,ba,Mstar, MHI, sigmaHI, mw202orig, mw202, mw20err2,/silent
;    readcol, datadir + 'geha06_color.txt',objID,ra,dec,camcol,field,type,unknown,type2,type_str,m_r2,m_i2,format='(F,F,F,F,F,F,F,F,A,F,F)'
    match,LONG(galnum),LONG(galnumAll),temp,selectedGeha
;    match2,LONG(galnum),LONG(galnumAll),ind1,ind2
    mi = mi + 5*alog10(0.7)
    r_i = m_r[selectedGeha] - mi[temp]
    mi_j = m_r[selectedGeha] - 1.2444*r_i - 0.3820 ;convert using Lupton 2005
    MHI = 10^MHI
    Mstar = 10^Mstar
    IF NOT keyword_set(pizagno) THEN BEGIN
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
    ENDIF ELSE BEGIN
        marladata = replicate( $
                    create_struct('w20', 0.0, $
                                  'w20_err', 0.0, $
                                  'i', 0.0, $
                                  'i_err', 0.0), $
                    n_elements(mi_j))
        
        pizagnodata = replicate( $
                      create_struct('w20', 0.0, $
                                    'w20_err', 0.0, $
                                    'i', 0.0, $
                                    'i_err', 0.0), $
                      n_elements(mi_p))
    ENDELSE
    IF NOT keyword_set(pizagno) THEN BEGIN
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
    ENDIF ELSE BEGIN
        pizagnodata.w20 = vcirc80_p
        pizagnodata.w20_err = vcirc80err_p
        pizagnodata.i = mi_p
        pizagnodata.i_err = mierr_p

        marladata.w20 = mw20[temp]; - 6
        marladata.w20_err = mw20err[temp]
        marladata.i = mi_j
        marladata.i_err = mierr[temp]
ENDELSE
    colobs = 100
ENDIF ELSE BEGIN 
;From Giovanelli, Haynes et al 1997 (Too high of TF) 
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

;------------------------------------ Stellar Mass vs Velocity ----------------------------------------
!p.multi = 0
loadct,obsct
IF (keyword_set(outfile)) THEN device,filename=outfile+'_stfgeha.eps',/color,bits_per_pixel= 8,/times,xsize = 10,ysize = 12,xoffset =  2,yoffset =  2  ELSE $
  window,0,xsize = 392,ysize = 500 ;xsize = 712,ysize = 392;,xsize = 950,ysize = 500;392
plot,[velocities[0],velocities[0]],[bmass[0],bmass[0]],psym = symcat(symbols[0]),symsize = symsizes[0],thick = thicks[0],$
     xrange = [9,105],yrange = [6.9,9.2],/xlog,ystyle=1,xstyle=1, $
     xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
     ytitle='log Stellar Mass [M' + sunsymbol() + textoidl('h_{70}^{-2}]'),/nodata ;/ylog,yrange = [8e6,1.5e9
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
oplot, mw202[selectedGeha] - 6, alog10(Mstar[selectedGeha]),psym = symcat(obspsym),color = obscolor,symsize = obssymsize
oploterror, mw202[selectedGeha] - 6, alog10(Mstar[selectedGeha]), mw20err2[selectedGeha], hierror[selectedGeha],psym=3,errcolor=obscolor,/nohat,color=obscolor
FOR i = 0, nfiles - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[velocities[i],velocities[i]],alog10([smass[i],smass[i]]),color = colors[i],psym = symcat(symbols[i]),symsize = symsizes[i],thick = thicks[i]
ENDFOR
IF keyword_set(key) THEN legend,key,color = colors,thick = thicks,psym = symbols,symsize = symsizes,ctables = ctables,linestyle = fltarr(nfiles),/top,/left;,charsize = l_charsize

IF (keyword_set(outfile)) THEN BEGIN
   device,/close
   device,filename=outfile+'_gtfgeha.eps',/color,bits_per_pixel= 8,/times,xsize = 10,ysize = 12,xoffset =  2,yoffset =  2  
ENDIF ELSE BEGIN
    stop
    window,0,xsize = 392,ysize = 500;,xsize = 712,ysize = 392;,xsize = 950,ysize = 500 ;392
ENDELSE

;------------------------------------ Gass Mass vs Velocity ----------------------------------------
loadct,obsct
plot,[velocities[0],velocities[0]],[bmass[0],bmass[0]],psym = symcat(symbols[0]),symsize = symsizes[0],thick = thicks[0],$
     xrange = [9,105],yrange = [6.9,9.2],/xlog, ystyle=1, xstyle=1, $
     xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
     ytitle='Log HI Mass [M' + sunsymbol() + textoidl('h_{70}^{-2}]'),/nodata ;/ylog,yrange = [8e6,1.5e9]
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
oplot,mw202[selectedGeha] - 6, alog10(MHI[selectedGeha]),psym = symcat(obspsym),color = obscolor,symsize = obssymsize
oploterror, mw202[selectedGeha] - 6, alog10(MHI[selectedGeha]), mw20err2[selectedGeha], sigmaHI[selectedGeha],psym=3,errcolor=obscolor,/nohat,color=obscolor
FOR i = 0, nfiles - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[velocities[i],velocities[i]],alog10([gmass[i],gmass[i]]),psym = symcat(symbols[i]),symsize = symsizes[i],color = colors[i], thick = thicks[i]
ENDFOR
IF keyword_set(key) THEN legend,key,color = colors,thick = thicks,psym = symbols,symsize = symsizes,ctables = ctables,linestyle = fltarr(nfiles),/top,/left;,charsize = l_charsize


IF (keyword_set(outfile)) THEN BEGIN
   device,/close
   device,filename=outfile+'_btfgeha.eps',/color,bits_per_pixel= 8,/times,xsize = 10,ysize = 12,xoffset =  2,yoffset =  2  
ENDIF ELSE BEGIN
    stop
    window,0,xsize = 392,ysize = 500;,xsize = 712,ysize = 392;,xsize = 950,ysize = 500 ;392
ENDELSE

;------------------------------------ Baryonic TF ----------------------------------------
loadct,obsct
plot,[velocities[0],velocities[0]],[bmass[0],bmass[0]],psym = symcat(symbols[0]),symsize = symsizes[0],thick = thicks[0],$
     xrange = [9,105],yrange = [6.9,9.2],/xlog,ystyle=1,xstyle=1, $
     xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
     ytitle='Log Baryonic Mass [M' + sunsymbol() + textoidl('h_{70}^{-2}]'),/nodata ;,yrange = [alog10(8e6),alog10(1.5e9)
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
oplot,mw202[selectedGeha] - 6, alog10(Mstar[selectedGeha] + MHI[selectedGeha]*1.4),psym = symcat(obspsym),color = obscolor,symsize = obssymsize
oploterror,mw202[selectedGeha] - 6, alog10(Mstar[selectedGeha] + MHI[selectedGeha]*1.4), mw20err2[selectedGeha],sigmaHI[selectedGeha],psym=3, errcolor=obscolor,/nohat,color=obscolor, symsize = obssymsize;

FOR i = 0, nfiles - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[velocities[i],velocities[i]],alog10([bmass[i],bmass[i]]),psym = symcat(symbols[i]),symsize = symsizes[i],color = colors[i], thick = thicks[i]
ENDFOR
IF keyword_set(key) THEN legend, key,color = colors,thick = thicks,psym = symbols,symsize = symsizes,ctables = ctables,linestyle = fltarr(nfiles),/top,/left;,charsize = l_charsize

IF (keyword_set(outfile)) THEN BEGIN
   device,/close
   device,filename=outfile+'_tfgeha.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  
ENDIF ELSE BEGIN
    stop
    window,0,xsize = 712,ysize = 392;,xsize = 950,ysize = 500 ;392
ENDELSE

;------------------------------------ Baryonic TF AND I-band Mag vs Velocity ----------------------------------------
loadct,obsct
!p.multi = [0,2,1]
IF 1 THEN BEGIN  ;Geha
    plot, marladata.w20, marladata.i, $
          xrange=[9, 100],yrange=[-13.8,-17.1],/xlog,ystyle=1, xstyle=1, $
          xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
          ytitle='M' + textoidl('_I') + ' - 5' + textoidl(' log_{10} h_{70}'),/nodata
    oplot,marladata.w20, marladata.i,psym=symcat(obspsym),color=obscolor, symsize = obssymsize ;16
;    oplot, mcgdata.w20, mcgdata.i, psym=2, color=obscolor 
;    oplot, mattdata.w20, mattdata.i, psym=2, color=obscolor 
;    oplot, sfidata.w20, sfidata.i, psym=2, color=obscolor 
;    oplot, verhdata.w20, verhdata.i, psym=2, color=obscolor 
    oploterror, marladata.w20, marladata.i, marladata.w20_err, marladata.i_err,psym=3, errcolor=obscolor,/nohat,symsize = obssymsize    
    x = [30,300]
    y = [-14.5, -24]
    coeff = line(alog10(x), y)
    x = [10, 500]
    y = coeff[0] + coeff[1] * alog10(x)
    oplot,x,y
ENDIF ELSE BEGIN                ;Giovanelli
    plot,logW - alog10(2),y + 5*alog10(h),$
         xrange = [1.6,2.6],yrange = [-17.5,-25],$
         xtitle = textoidl('log_{10}W'),$
         ytitle = 'M - 5 log h',/nodata ;yrange = [-15,-24],xrange = [1.9,2.9]
    oplot,logW - alog10(2),y + 5*alog10(h),psym=symcat(obspsym), color=obscolor, symsize = obssymsize
    oplot,logW + 5*alog10(h),ymax,linestyle = 1
    oplot,logW + 5*alog10(h),ymin,linestyle = 1
ENDELSE
FOR i = 0, nfiles - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[velocities[i],velocities[i]],[Imag[i] - 5.0*alog10(h),Imag[i] - 5.0*alog10(h)], $
          psym = symcat(symbols[i]),symsize = symsizes[i],color = colors[i], thick = thicks[i]
ENDFOR
;    IF keyword_set(key) THEN legend, key, color = colors, thick = thicks, /top, /left, linestyle = fltarr([nfiles]),charsize = l_charsize,psym = 4*(fltarr(nfiles) + 1)
loadct,obsct
plot,[velocities[0],velocities[0]],[bmass[0],bmass[0]],psym = symbols[0],symsize = symsizes[0],thick = thicks[0],$
     xrange = [9,105],yrange = [7.7,9.2],/xlog,ystyle=1,xstyle=1, $
     xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
     ytitle='log Baryonic Mass [M' + sunsymbol() + textoidl(' h_{70}^{-2}]'),/nodata ;,yrange = [alog10(8e6),alog10(1.5e)9
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
oplot,mw202[selectedGeha] - 6, alog10(Mstar[selectedGeha] + MHI[selectedGeha]*1.4),psym=symcat(obspsym),color=obscolor, symsize = obssymsize
oploterror, mw202[selectedGeha] - 6, alog10(Mstar[selectedGeha] + MHI[selectedGeha]*1.4), mw20err2[selectedGeha],sigmaHI[selectedGeha],psym=3, errcolor=obscolor,/nohat, color=obscolor,symsize = obssymsize
IF keyword_set(color) THEN BEGIN
    FOR i = 0, nfiles - 1 DO BEGIN
            loadct,ctables[i]
            oplot,[velocities[i],velocities[i]],alog10([bmass[i],bmass[i]]),psym = symcat(symbols[i]),symsize = symsizes[i],color = colors[i], thick = thicks[i]
        ENDFOR
    IF keyword_set(key) THEN legend, [key,"Geha et al. 2006"], /top, /left, linestyle = [fltarr(nfiles),1],psym  = [symbols,16], color = [colors,obscolor],thick = [thicks,1], symsize = [symsizes,1];charsize = l_charsize
ENDIF ELSE BEGIN
    FOR i = 0, nfiles - 1 DO BEGIN
        loadct,ctables[i]
        oplot,[velocities[i],velocities[i]],alog10([bmass[i],bmass[i]]),psym = symcat(symbols[i]),symsize = symsizes[i], thick = thicks[i]
    ENDFOR
    IF keyword_set(key) THEN legend, [key,"Geha et al. 2006"], /top, /left, linestyle = [fltarr(nfiles),1],psym = [symbols,16],  color = [colors,obscolor],thick = [thicks,1], symsize = [symsizes,1];, charsize = l_charsize
ENDELSE

IF (keyword_set(outfile)) THEN BEGIN
   device,/close
   device,filename=outfile+'_tfgehaMassiv.eps',/color,bits_per_pixel= 8,/times,xsize = 14.4,ysize = 17.25,xoffset =  2,yoffset =  2  
ENDIF ELSE BEGIN
    stop
    window,0,xsize = 392,ysize = 500 ;392
ENDELSE

;-------------------------------------------------------------------------------------------------
!p.multi = 0
loadct,obsct
print,obsct
IF 1 THEN BEGIN  ;Geha
    IF NOT keyword_set(pizagno) THEN BEGIN
        xrange = [10,400]
        yrange = [-11.5,-24]
        plot, marladata.w20, marladata.i, $
              xrange=xrange,yrange=yrange,/xlog,ystyle=1, xstyle=1, $
              xtitle='W' + textoidl('_{HI}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
              ytitle='M' + textoidl('_I') + ' - 5' + textoidl('log_{10}h_{70}'),/nodata, symsize = obssymsize
        oploterror, marladata.w20, marladata.i, marladata.w20_err, marladata.i_err,psym=symcat(19), errcolor=obscolor,color=obscolor,/nohat, symsize = obssymsize  ;Filled square                  
        oplot, mcgdata.w20, mcgdata.i,    psym=symcat(14),color=obscolor,symsize = obssymsize  ;filled diamond,14 ;McGaugh et al., 2000
        oplot, mattdata.w20, mattdata.i,  psym=symcat(17),color=obscolor,symsize = obssymsize  ;filled upward triangle, 17 ;Matthews, van Driel, Gallagher 1998
        oplot, sfidata.w20, sfidata.i,    psym=symcat(18),color=obscolor,symsize = obssymsize ;filled downward triangle,18 ;Haynes et al., 1999
        oplot, verhdata.w20, verhdata.i,  psym=symcat(46),color=obscolor,symsize = obssymsize ;filled star, 46 ;Verheijen 2001
        FOR j  = 0, nfiles - 1 DO BEGIN
            loadct,ctables[j]
            oplot,[velocities[j],velocities[j]],[Imag[j] - 5.0*alog10(h),Imag[j] - 5.0*alog10(h)], $
                  psym = symcat(symbols[j]),symsize = symsizes[j],color = colors[j], thick = thicks[j]
        ENDFOR
    ENDIF ELSE BEGIN
        xrange = [10,400]
        yrange = [-11.5,-24]
        plot, marladata.w20, marladata.i, $
              xrange=xrange,yrange=yrange,/xlog,ystyle=1, xstyle=1, $
              xtitle=textoidl('Velocity [km s^{-1}]'), $
              ytitle=textoidl('M_i'),/nodata
        oploterror, marladata.w20, marladata.i, marladata.w20_err, marladata.i_err,psym=symcat(19), errcolor=obscolor,color=obscolor,/nohat, symsize = obssymsize             
        oploterror,pizagnodata.w20, pizagnodata.i, pizagnodata.w20_err, pizagnodata.i_err,psym=symcat(16), errcolor=obscolor,color = obscolor,/nohat, symsize = obssymsize
        FOR j  = 0, nfiles - 1 DO BEGIN
            loadct,ctables[j]
            oplot,[velocities[j],velocities[j]],[isdssmag[j],isdssmag[j]], $
                  psym = symcat(symbols[j]),symsize = symsizes[j],color = colors[j], thick = thicks[j]
        ENDFOR
    ENDELSE
    x = [30,300]
    y = [-14.5, -24]
    coeff = line(alog10(x), y)
    x = [10, 500]
    y = coeff[0] + coeff[1] * alog10(x)
;    oplot,x,y
ENDIF ELSE BEGIN                ;Giovanelli
    plot,logW - alog10(2),y + 5*alog10(h),$
         xrange = [1.6,2.6],yrange = [-13,-25],$
         xtitle = textoidl('Velocity [km s^{-1}]'),$ ;textoidl('log_{10}W'),$
         ytitle = textoidl('M - 5 log_{10}h_{70}'),/nodata ;yrange = [-15,-24],xrange = [1.9,2.9]
    oplot,logW - alog10(2),y + 5*alog10(h),psym=symcat(obspsym), color=obscolor, symsize = obssymsize
    oplot,logW + 5*alog10(h),ymax,linestyle = 1
    oplot,logW + 5*alog10(h),ymin,linestyle = 1
    FOR j  = 0, nfiles - 1 DO BEGIN
        loadct,ctables[j]
;    oplot,[velocities[j],velocities[j]],[Imag[j] - 5.0*alog10(h),Imag[j] - 5.0*alog10(h)], $
        oplot,[velocities[j],velocities[j]],[isdssmag[j],isdssmag[j]], $
              psym = symcat(symbols[j]),symsize = symsizes[j],color = colors[j], thick = thicks[j]
    ENDFOR
ENDELSE
uniqind0 = (uniq(symbols,sort(symbols)))[0]
IF n_elements((uniq(symbols,sort(symbols)))) GT 1 THEN uniqind1 = (uniq(symbols,sort(symbols)))[1]
IF NOT keyword_set(pizagno) THEN BEGIN
    obskey = ['Geha et al., 2006','McGaugh et al., 2000','Matthews et al., 1998','Haynes et al., 1999','Verheijen 2001']
    obscolorl = obscolor + fltarr(n_elements(obskey))
    obssyml = [19,14,17,18,46]
    obssymsizel = obssymsize + fltarr(n_elements(obskey))
    obsctl = obsct + fltarr(n_elements(obskey))
    obssymthick = 1 + fltarr(n_elements(obskey))
ENDIF ELSE BEGIN
    obskey = ['Geha et al., 2006','Pizagno et al., 2007']
    obssyml = [19,16]
    obscolorl = obscolor + fltarr(n_elements(obskey))
    obssymsizel = obssymsize + fltarr(n_elements(obskey))
    obsctl = obsct + fltarr(n_elements(obskey))
    obssymthick = 1 + fltarr(n_elements(obskey))
ENDELSE
IF keyword_set(key) THEN legend, key,color = colors, thick = thicks,psym = symbols, symsize = symsizes, ctables = ctables, linestyle = fltarr(nfiles),/top, /left ELSE $
  IF n_elements((uniq(symbols,sort(symbols)))) GT 1 $
  THEN legend,[obskey,'Med-res Sims','High-res Sims'],color = [obscolorl,colors[uniqind0],colors[uniqind1]],psym=[obssyml,symbols[uniqind0],symbols[uniqind1]],thick = [obssymthick,thicks[uniqind0],thicks[uniqind1]],symsize = [obssymsizel,symsizes[uniqind0],symsizes[uniqind1]],ctables = [obsctl,ctables[uniqind0],ctables[uniqind1]],/left,/top,box = 0,charsize = l_charsize $
  ELSE legend,[obskey,'Sims'],color = [obscolorl,colors[0]],psym=[obssyml,symbols[0]],thick = [obssymthick,thicks[0]],symsize = [obssymsizel,symsizes[0]],ctables = [obsctl,ctables[0]],/left,/top,box = 0;,charsize = l_charsize

IF (keyword_set(outfile)) THEN device,/close ELSE stop

END
