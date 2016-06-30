
;Observational Data from McGaugh 08
PRO tully_fisher_obs_btf, broadband, velocities, gmass, smass_true = smass_true,velocities_true = velocities_true,filternums = filternums,key = key, color = color,symbols = symbols,symsizes = symsizes,thicks = thicks,obscolor = obscolor,ctables = ctables,formatthick = formatthick,ctobs = ctobs, outfile = outfile,kcorrect = kcorrect,halo = halo,obssym = obssym, obsthick = obssymthick, obssize = obssize,obsct = obsct
formatplot,outplot = outfile,thick = formatthick
IF keyword_set(outfile) THEN BEGIN
    nbins=100.0
    linestyles = [0,2]
    l_charsize = !p.charsize;1.0
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    nbins=100.0
    linestyles = [0,2]
    l_charsize = !p.charsize;0.75
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF keyword_set(color) THEN BEGIN
    IF NOT keyword_set(obscolor) THEN obscolor = fgcolor
    IF NOT keyword_set(obsct) THEN obsct = 0 ;39
    IF NOT keyword_set(ctables) then ctables = fltarr(n_elements(broadband)) + 39
    IF color[0] eq 1 then  colors = fltarr(n_elements(broadband)) + fgcolor ELSE colors = fltarr(n_elements(broadband)) + color ;colors = (findgen(n_elements(broadband)) + 1)*240/n_elements(broadband)
    IF NOT keyword_set(symbols) THEN symbols = fltarr(n_elements(broadband)) + 4
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n_elements(broadband)) + 2;1.5
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n_elements(broadband)) + 5
ENDIF ELSE BEGIN  
    IF NOT keyword_set(obscolor) THEN obscolor = 150
    IF NOT keyword_set(obsct) THEN obsct = 0
    IF NOT keyword_set(ctables) then ctables = fltarr(n_elements(broadband))
    colors = fgcolor + fltarr(n_elements(broadband))
;fgcolor + (findgen(n_elements(broadband)) + 1)*10.0 + 5.0 ELSE colors = fgcolor - ((findgen(n_elements(broadband)) + 1)*10.0 + 5.0);  fltarr(n_elements(broadband)) + 5
    IF NOT keyword_set(symbols) THEN symbols = (findgen(n_elements(broadband))+2)*2
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n_elements(broadband)) + 2;1.5
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n_elements(broadband)) + 5
ENDELSE
IF NOT keyword_set(obssym) THEN obssym = 16 ;2
IF NOT keyword_set(obssize) THEN obssize = 1.25 ;1
IF NOT keyword_set(obsthick) THEN obsthick = 2
n = n_elements(broadband)

IF NOT keyword_set(filternums) THEN filternums = intarr(n) + 13
IF n_elements(filternums) eq 1 THEN filternums = intarr(n) + filternums

L_B_sol = 4.67e25 ;Watts, from B+M, p 53, Table 2.1
MLa = -0.942
MLb = 1.69
smass = fltarr(n_elements(broadband))
FOR i = 0, n_elements(broadband) - 1 DO BEGIN
    IF (keyword_set(kcorrect)) OR NOT file_test(broadband[i] + '.' + strtrim(halo[i],2) + '/broadband.fits') THEN BEGIN
        print,'Kcorrect',broadband[i],halo[i]
        IF 0 THEN BEGIN
;file_test(broadband[i] + '.' + strtrim(halo[i],2) + '/broadband.fits') THEN BEGIN
            lums = mrdfits(broadband[i] + '.' + strtrim(halo[i],2) + '/broadband.fits', filternums[i], header)
            Un = where(strcmp(strtrim(lums.filter,2),'U_Johnson.res') eq 1)
            Bn = where(strcmp(strtrim(lums.filter,2),'B_Johnson.res') eq 1)
            Vn = where(strcmp(strtrim(lums.filter,2),'V_Johnson.res') eq 1)
            Rn = where(strcmp(strtrim(lums.filter,2),'R_Cousins.res') eq 1)
            In = where(strcmp(strtrim(lums.filter,2),'I_Cousins.res') eq 1)
            mags_bes = transpose([[lums[Un].AB_mag1],[lums[Bn].AB_mag1],[lums[Vn].AB_mag1],[lums[Rn].AB_mag1],[lums[In].AB_mag1]])
        ENDIF ELSE BEGIN
            mags = mrdfits(broadband[i] + '.amiga_r200.halos.star99_K_ab.Mv.fits',1)
            ind = where(fix(mags.id) EQ halo[i])
            mags = mags[ind]
            mags_bes = transpose([[mags.u],[mags.b],[mags.v],[mags.r],[mags.i]])
        ENDELSE
        mags_errs = fltarr(5,n_elements(mags.u))+0.02
        mags_ivar=1./mags_errs^2
        dmod = 19.4576          ;distance modulus for z = 0 in this model
        z = fltarr(n_elements(mags.u))
        mgy =(10.D)^(-(0.4D)*(mags_bes + dmod))
        mgy_ivar = mags_ivar/(0.4*alog(10.)*mgy)^2.
        kcorrect, mgy, mgy_ivar, z, kcor, mass = star, mtol = mtol_k, absmag = absmag_k, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']
        smass[i] = star
    ENDIF ELSE BEGIN
        print,'Sunrise',broadband[i],halo[i]
        lums = mrdfits(broadband[i] + '.' + strtrim(halo[i],2) + '/broadband.fits', filternums[i], header)
        Bn = where(STRCMP(strtrim(lums.filter,2),'B_Johnson.res') eq 1)
        Vn = where(STRCMP(strtrim(lums.filter,2),'V_Johnson.res') eq 1)
        Bmag = lums[Bn].AB_mag1
        Vmag = lums[Vn].AB_mag1
        Bmag = Bmag + 0.09                    ; to Vega
        Vmag = Vmag + 0.03                    ; to Vega
        L_B_band = lums[Bn].L_lambda_eff1     ;Watts/m
        B_band_width = lums[Bn].ewidth_lambda ;m
        L_B = L_B_band*B_band_width           ;Watts
        B_V = Bmag - Vmag
        MLfit = 10^(MLa + MLb*B_V)
        smass[i] = MLfit*L_B/L_B_sol 
    ENDELSE
ENDFOR
bmass = smass + gmass
do_true = keyword_set(smass_true) OR keyword_set(velocities_true)
IF keyword_set(smass_true) THEN bmass_true = smass_true + gmass ELSE BEGIN 
    smass_true = smass
    bmass_true = bmass
ENDELSE
IF NOT keyword_set(velocities_true) THEN velocities_true = velocities

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN datadir = '/home/christensen/Code/Datafiles/HIcubes/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN datadir = '/nobackupp2/crchrist/MolecH/' $
ELSE datadir = '/astro/users/christensen/code/Datafiles/HIcubes/'
datafile = datadir + 'BTF.dat'
readcol,datafile,name,V_f,M_star,M_gas,mu_0,R_d,B_V,gamma_max,gamma_pop,gamma_acc,FORMAT = '(A9D)'
M_star = M_star*1e10
M_gas = M_gas*1e10
M_disk = M_star + M_gas

;------------------------------------------ Gas Fraction ---------------
!p.multi = [0,2,1]
IF (keyword_set(outfile)) THEN device,filename=outfile+'_gfrac_btf.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 712,ysize = 392

;Figure 5, top left
loadct,obsct
yrange = [0,1]
xrange = [1e8,1e12]
xrange = [8e6,1e12]
plot, M_disk,M_star/M_disk,psym = 1,xtitle = 'Disk Mass [Msol]',ytitle = 'Stellar Fraction',/xlog,xrange= xrange,yrange = yrange,/nodata
oplot,M_disk,M_star/M_disk,psym = symcat(obssym),color = obscolor,thick = obsthick,symsize = obssize
FOR i = 0, n_elements(broadband) - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[bmass[i],bmass[i]],[smass[i]/bmass[i],smass[i]/bmass[i]],color = colors[i],psym = symcat(symbols[i]),thick = thicks[i],symsize = symsizes[i]
    IF keyword_set(do_true) THEN oplot,[bmass_true[i],bmass_true[i]],[smass_true[i]/bmass_true[i],smass_true[i]/bmass_true[i]],color = colors[i],psym = symcat(9),thick = obsthick,symsize = symsizes[i]
ENDFOR
IF keyword_set(key) THEN legend,['McGaugh 2005',key],color = [obscolor,colors],psym=[obssym,symbols],thick = [obsthick,thicks],symsize = [obssize,symsizes],ctables = [obsct,obsctctables],/right,/bottom,box = 0,charsize = l_charsize

;--------------------------------------------- Baryonic TF ------------
loadct,obsct
yrange = [3e8,4e11]
yrange = [8e6,4e11]
xrange = [50,300]
xrange = [20,300]
;xrange = [10,300]
plot,V_f,M_star + M_gas,psym = symcat(obssym),/xlog,/ylog,xrange = xrange,yrange = yrange,xtitle = textoidl('Velocity [km s^{-1}]'),ytitle = textoidl('Baryonic Mass [M')+sunsymbol()+']',/nodata
oplot,V_f,M_star + M_gas,psym = symcat(obssym),color = obscolor,thick = obsthick,symsize = obssize
FOR i = 0, n_elements(broadband) - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[velocities[i],velocities[i]],[bmass[i],bmass[i]],psym = symcat(symbols[i],thick = thicks[i]),symsize = symsizes[i],color = colors[i]
    IF keyword_set(do_true) THEN  oplot,[velocities_true[i],velocities_true[i]],[bmass_true[i],bmass_true[i]],color = colors[i],psym = symcat(9),thick = thicks[i],symsize = symsizes[i]/2.0  
ENDFOR
IF (keyword_set(outfile)) THEN device,/close ELSE stop

;----------------------------------------------- Baryonic TF -------------
!p.multi = 0
IF (keyword_set(outfile)) THEN device,filename=outfile+'_btf.eps',/color,bits_per_pixel= 8,/times,xsize = 14.4,ysize = 17.25,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 392,ysize = 500
loadct,obsct
yrange = [3e8,4e11]
yrange = [5e6,5e11]
xrange = [50,300]
xrange = [20,300]
;xrange = [10,300]
plot,V_f,M_star + M_gas,psym = symcat(obssym),/xlog,/ylog,xrange = xrange,yrange = yrange,xtitle = textoidl('Velocity [km s^{-1}]'),ytitle = textoidl('Baryonic Mass [M')+sunsymbol()+']',/nodata
oplot,V_f,M_star + M_gas,psym = symcat(obssym),color = obscolor,thick = obsthick,symsize = obssize
FOR i = 0, n_elements(broadband) - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[velocities[i],velocities[i]],[bmass[i],bmass[i]],psym = symcat(symbols[i],thick = thicks[i]),symsize = symsizes[i],color = colors[i]
    IF keyword_set(do_true) THEN  oplot,[velocities_true[i],velocities_true[i]],[bmass_true[i],bmass_true[i]],color = colors[i],psym = symcat(9),thick = thicks,symsize = symsizes[i]/2.0
ENDFOR
uniqind0 = (uniq(symbols,sort(symbols)))[0]
IF n_elements((uniq(symbols,sort(symbols)))) GT 1 THEN uniqind1 = (uniq(symbols,sort(symbols)))[1]
IF keyword_set(key) THEN legend,['McGaugh 2005',key],color = [obscolor,colors],psym=[obssym,symbols],thick = [obsthick,thicks],symsize = [obssize,symsizes],ctables = [obsct,obsctctables],/right,/bottom,box = 0,charsize = l_charsize ELSE $
  IF n_elements((uniq(symbols,sort(symbols)))) GT 1 $
  THEN legend,['McGaugh 2005','Med-res Sims','High-res Sims'],color = [obscolor,colors[uniqind0],colors[uniqind1]],psym=[obssym,symbols[uniqind0],symbols[uniqind1]],thick = [obsthick,thicks[uniqind0],thicks[uniqind1]],symsize = [obssize,symsizes[uniqind0],symsizes[uniqind1]],ctables = [obsct,ctables[uniqind0],ctables[uniqind1]],/right,/bottom,box = 0,charsize = l_charsize $
  ELSE legend,['McGaugh 2005','Simulations'],color = [obscolor,colors[0]],psym=[obssym,symbols[0]],thick = [obsthick,thicks[0]],symsize = [obssize,symsizes[0]],ctables = [obsct,ctables[0]],/top,/left,box = 0,charsize = l_charsize
IF (keyword_set(outfile)) THEN device,/close ELSE stop

;------------------------------------------------ Stellar TF -----------------
!p.multi = 0
IF (keyword_set(outfile)) THEN device,filename=outfile+'_stf.eps',/color,bits_per_pixel= 8,/times,xsize = 15,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 392,ysize = 500
loadct,obsct
yrange = [3e8,4e11]
yrange = [8e6,4e11]
xrange = [50,300]
xrange = [20,300]
;xrange = [10,300]
plot,V_f,M_star,psym = symcat(obssym),/xlog,/ylog,xrange = xrange,yrange = yrange,xtitle = textoidl('Velocity [km s^{-1}]'),ytitle = textoidl('Stellar Mass [M')+sunsymbol()+']',/nodata
oplot,V_f,M_star ,psym = symcat(obssym),color = obscolor,thick = obsthick,symsize = obssize
FOR i = 0, n_elements(broadband) - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[velocities[i],velocities[i]],[smass[i],smass[i]],psym = symcat(symbols[i],thick = thicks[i]),symsize = symsizes[i],color = colors[i]
    IF keyword_set(do_true) THEN  oplot,[velocities_true[i],velocities_true[i]],[smass_true[i],smass_true[i]],color = colors[i],psym = symcat(9),thick = thicks[0],symsize = symsizes[i]/2.0
ENDFOR
uniqind0 = (uniq(symbols,sort(symbols)))[0]
IF n_elements((uniq(symbols,sort(symbols)))) GT 1 THEN uniqind1 = (uniq(symbols,sort(symbols)))[1]
;IF keyword_set(key) THEN legend,['McGaugh',key],color = [obscolor,colors],psym=[obssym,symbols],thick = [obsthick,thicks],symsize = [obssize,symsizes],ctables = [obsct,obsctctables],/right,/bottom,box = 0,charsize = l_charsize ELSE $
;  IF n_elements((uniq(symbols,sort(symbols)))) GT 1 $
;  THEN legend,['McGaugh 2005','Low-res Sims','High-res Sims'],color = [obscolor,colors[uniqind0],colors[uniqind1]],psym=[obssym,symbols[uniqind0],symbols[uniqind1]],thick = [obsthick,thicks[uniqind0],thicks[uniqind1]],symsize = [obssize,symsizes[uniqind0],symsizes[uniqind1]],ctables = [obsct,ctables[uniqind0],ctables[uniqind1]],/right,/bottom,box = 0,charsize = l_charsize $
;  ELSE legend,['McGaugh 2005','Simulations'],color = [obscolor,colors[0]],psym=[obssym,symbols[0]],thick = [obsthick,thicks[0]],symsize = [obssize,symsizes[0]],ctables = [obsct,ctables[0]],/top,/left,box = 0,charsize = l_charsize
IF (keyword_set(outfile)) THEN device,/close ELSE stop

;-------------------------------------------- Baryonic and Stellar TF
IF (keyword_set(outfile)) THEN device,filename=outfile+'_bstf.eps',/color,bits_per_pixel= 8,/times,xsize = 15,ysize = 12,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 500,ysize = 392
loadct,obsct
yrange = [3e8,4e11]
yrange = [8e6,4e11]
xrange = [50,300]
xrange = [20,300]
;xrange = [10,300]
multiplot,[2,1],mytitle = 'Mass [M' + sunsymbol() + ']',mxtitle =  textoidl('Velocity [km s^{-1}]'),mxtitsize = 1.5,mytitsize = 1.5
plot,V_f,M_star + M_gas,psym = symcat(obssym),/xlog,/ylog,xrange = xrange,yrange = yrange,/nodata
oplot,V_f,M_star + M_gas,psym = symcat(obssym),color = obscolor,thick = obsthick,symsize = obssize
FOR i = 0, n_elements(broadband) - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[velocities[i],velocities[i]],[bmass[i],bmass[i]],psym = symcat(symbols[i],thick = thicks[i]),symsize = symsizes[i],color = colors[i]
    IF keyword_set(do_true) THEN  oplot,[velocities_true[i],velocities_true[i]],[bmass_true[i],bmass_true[i]],color = colors[i],psym = symcat(9),thick = 5,symsize = symsizes[i]/2.0
ENDFOR
uniqind0 = (uniq(symbols,sort(symbols)))[0]
IF n_elements((uniq(symbols,sort(symbols)))) GT 1 THEN uniqind1 = (uniq(symbols,sort(symbols)))[1]

legend,['Baryonic Mass'],box = 0
multiplot
loadct,obsct
plot,V_f,M_star,psym = symcat(obssym),/xlog,/ylog,xrange = xrange,yrange = yrange,/nodata
oplot,V_f,M_star ,psym = symcat(obssym),color = obscolor,thick = obsthick,symsize = obssize
FOR i = 0, n_elements(broadband) - 1 DO BEGIN
    loadct,ctables[i]
    oplot,[velocities[i],velocities[i]],[smass[i],smass[i]],psym = symcat(symbols[i],thick = thicks[i]),symsize = symsizes[i],color = colors[i]
    IF keyword_set(do_true) THEN  oplot,[velocities_true[i],velocities_true[i]],[smass_true[i],smass_true[i]],color = colors[i],psym = symcat(9),thick = 5,symsize = symsizes[i]/2.0
ENDFOR
uniqind0 = (uniq(symbols,sort(symbols)))[0]
IF n_elements((uniq(symbols,sort(symbols)))) GT 1 THEN uniqind1 = (uniq(symbols,sort(symbols)))[1]
legend,['Stellar Mass'],box = 0
IF keyword_set(key) THEN legend,['McGaugh 2005',key],color = [obscolor,colors],psym=[obssym,symbols],thick = [obsthick,thicks],symsize = [obssize,symsizes],ctables = [obsct,obsctctables],/right,/bottom,box = 0,charsize = l_charsize ELSE $
  IF n_elements((uniq(symbols,sort(symbols)))) GT 1 $
  THEN legend,['McGaugh 2005','Med-res Sims','High-res Sims'],color = [obscolor,colors[uniqind0],colors[uniqind1]],psym=[obssym,symbols[uniqind0],symbols[uniqind1]],thick = [obsthick,thicks[uniqind0],thicks[uniqind1]],symsize = [obssize,symsizes[uniqind0],symsizes[uniqind1]],ctables = [obsct,ctables[uniqind0],ctables[uniqind1]],/right,/bottom,box = 0,charsize = l_charsize $
  ELSE IF keyword_set(smass_true) $
  THEN legend,['McGaugh 2005','Sims., Obs. Stellar Mass','Sims., True Stellar Mass'],color = [obscolor,colors[0],colors[0]],psym=[obssym,symbols[0],9],thick = [obsthick,thicks[0],5],symsize = [obssize,symsizes[0],symsizes[0]/2],ctables = [obsct,ctables[0],ctables[0]],/bottom,/right,box = 0,charsize = l_charsize*0.8 $
  ELSE legend,['McGaugh 2005','Simulations'],color = [obscolor,colors[0]],psym=[obssym,symbols[0]],thick = [obsthick,thicks[0]],symsize = [obssize,symsizes[0]],ctables = [obsct,ctables[0]],/bottom,/right,box = 0,charsize = l_charsize
multiplot,/reset

IF (keyword_set(outfile)) THEN device,/close ELSE stop
END
