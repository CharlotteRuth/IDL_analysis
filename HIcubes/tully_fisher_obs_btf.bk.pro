
;Observational Data from McGaugh 08
PRO tully_fisher_obs_btf, broadband, velocities, gmass, smass_true = smass_true,velocities_true = velocities_true,filternums = filternums,key = key, color = color,symbols = symbols,symsizes = symsizes,thicks = thicks,obscolor = obscolor,ctables = ctables,formatthick = formatthick,ctobs = ctobs, outfile = outfile,kcorrect = kcorrect,halo = halo
formatplot,outplot = outfile,thick = formatthick
IF keyword_set(outfile) THEN BEGIN
    nbins=100.0
    linestyles = [0,2]
    l_charsize = 0.75
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    nbins=100.0
    linestyles = [0,2]
    fgcolor = 255
    bgcolor = 0
ENDELSE
!p.multi = [0,2,1]
IF keyword_set(color) THEN BEGIN
    IF NOT keyword_set(obscolor) THEN obscolor = fgcolor
    IF NOT keyword_set(obsct) THEN obsct = 39
    obssym = 2
    obssymsize = 1
    IF NOT keyword_set(ctables) then ctables = fltarr(n_elements(broadband)) + 39
    if color[0] eq 1 then  colors = (findgen(n_elements(broadband)) + 1)*240/n_elements(broadband) else colors = color
    IF NOT keyword_set(symbols) THEN symbols = fltarr(n_elements(broadband)) + 4
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n_elements(broadband)) + 1.5
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n_elements(broadband)) + 5
ENDIF ELSE BEGIN  
    IF NOT keyword_set(obscolor) THEN obscolor = 150
    IF NOT keyword_set(obsct) THEN obsct = 0
    obssym = 2
    obssymsize = 1
    IF NOT keyword_set(ctables) then ctables = fltarr(n_elements(broadband))
;    IF keyword_set(outfile)  THEN 
    colors = fgcolor + fltarr(n_elements(broadband))
;fgcolor + (findgen(n_elements(broadband)) + 1)*10.0 + 5.0 ELSE colors = fgcolor - ((findgen(n_elements(broadband)) + 1)*10.0 + 5.0);  fltarr(n_elements(broadband)) + 5
    IF NOT keyword_set(symbols) THEN symbols = (findgen(n_elements(broadband))+2)*2
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n_elements(broadband)) + 1.5
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n_elements(broadband)) + 5
ENDELSE
obssymthick = 2
n = n_elements(broadband)
IF not keyword_set(filternums) THEN filternums = intarr(n) + 13
IF n_elements(filternums) eq 1 THEN filternums = intarr(n) + filternums

nk = n_elements(key)
lb_keys = strarr(nk)
lb_color = intarr(nk) + bgcolor
lb_linestyle = intarr(nk) + 1
lb_symbols = intarr(nk) + 3
lb_symsizes = intarr(nk + 1) + 1
lb_symthick = intarr(nk + 1) + 1

L_B_sol = 4.67e25 ;Watts, from B+M, p 53, Table 2.1
MLa = -0.942
MLb = 1.69
smass = fltarr(n_elements(broadband))
FOR i = 0, n_elements(broadband) - 1 DO BEGIN
    IF keyword_set(kcorrect) THEN BEGIN
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
    ENDIF ELSE BEGIN
        lums = mrdfits(broadband[i], filternums[i], header)
        Bn = where(STRCMP(strtrim(lums.filter,2),'B_Johnson.res') eq 1)
        Vn = where(STRCMP(strtrim(lums.filter,2),'V_Johnson.res') eq 1)
        Bmag = lums[Bn].AB_mag1
        Vmag = lums[Vn].AB_mag1
        L_B_band = lums[Bn].L_lambda_eff1 ;Watts/m
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
;datafile = '/Users/christensen/Code/IDL/MolecH/BTF.dat'
readcol,datafile,name,V_f,M_star,M_gas,mu_0,R_d,B_V,gamma_max,gamma_pop,gamma_acc,FORMAT = '(A9D)'
M_star = M_star*1e10
M_gas = M_gas*1e10
M_disk = M_star + M_gas

IF (keyword_set(outfile)) THEN device,filename=outfile+'_btf.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 712,ysize = 392
loadct,obsct
plot, M_disk,M_star/M_disk,psym = 1,xtitle = 'Disk Mass [Msol]',ytitle = 'Stellar Fraction',/xlog,xrange= [1e8,1e12],yrange = [-0.2,1.2],/nodata
oplot,M_disk,M_star/M_disk,psym = obssym,color = obscolor,thick = obssymthick
IF keyword_set(key) THEN,legend,['McGaugh et al. 2008',key],color = [obscolor,colors
IF keyword_set(key) THEN BEGIN
    l_keys = lb_keys
    l_keys[0] = 'McGaugh et al. 2008'
    l_color = lb_color
    l_color[0] = obscolor
    l_symbols = lb_symbols
    l_symbols[0]  = obssym
    l_symsizes = lb_symsizes
    l_symsizes[0]  = obssymsize
    l_symthick = lb_symthick
    l_symthick[0] = obssymthick 
    legend,l_keys,color = l_color,symsize = l_symsizes,psym = l_symbols,thick = l_symthick,/right,/bottom,box = 0
ENDIF

;Figure 5, top left
IF keyword_set(COLOR) THEN BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO BEGIN
        loadct,ctables[i]
        oplot,[bmass[i],bmass[i]],[smass[i]/bmass[i],smass[i]/bmass[i]],color = colors[i],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i]
        IF keyword_set(do_true) THEN oplot,[bmass_true[i],bmass_true[i]],[smass_true[i]/bmass_true[i],smass_true[i]/bmass_true[i]],color = colors[i],psym = symcat(symbols[i],thick = obssymthick),SYMSIZE = symsizes[i]
    ENDFOR
    IF keyword_set(key) THEN legend, key, color = colors, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),charsize = l_charsize,psym = symbols,symsize = symsizes,ctables = ctables
ENDIF ELSE BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO BEGIN
        loadct,ctables[i]
        oplot,[bmass[i],bmass[i]],[smass[i]/bmass[i],smass[i]/bmass[i]],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i],color = colors[i]
        IF keyword_set(do_true) THEN oplot,[bmass_true[i],bmass_true[i]],[smass_true[i]/bmass_true[i],smass_true[i]/bmass_true[i]],psym = symcat(symbols[i]),SYMSIZE = symsizes[i],color = colors[i] 
    ENDFOR
    IF keyword_set(key) THEN legend, key, thick = thicks, /top, /left, linestyle = fltarr([n_elements(broadband)]),psym = symbols, charsize = l_charsize,symsize = symsizes,ctables = ctables
ENDELSE

;IF keyword_set(color) THEN loadct,39 ELSE loadct,0
loadct,obsct
plot,V_f,M_star + M_gas,psym = obssym,/xlog,/ylog,xrange = [30,500],yrange = [1e8,4e11],xtitle = textoidl('V_{f} [km/s]'),ytitle = textoidl('Baryonic Mass [M')+sunsymbol()+']',/nodata
oplot,V_f,M_star + M_gas,psym = obssym,color = obscolor,thick = obssymthick
IF keyword_set(COLOR) THEN BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO BEGIN
        loadct,ctables[i]
        oplot,[velocities[i],velocities[i]],[bmass[i],bmass[i]],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i],color = colors[i]
      IF keyword_set(do_true) THEN  oplot,[velocities_true[i],velocities_true[i]],[bmass_true[i],bmass_true[i]],color = colors[i],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i]/2.0  
    ENDFOR
ENDIF ELSE BEGIN
  FOR i = 0, n_elements(broadband) - 1 DO BEGIN
      loadct,ctables[i]
      oplot,[velocities[i],velocities[i]],[bmass[i],bmass[i]],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i],color = colors[i]
      IF keyword_set(do_true) THEN oplot,[velocities_true[i],velocities_true[i]],[bmass_true[i],bmass_true[i]],color = colors[i],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i]/2.0 
  ENDFOR
ENDELSE
IF (keyword_set(outfile)) THEN device,/close ELSE stop

;-----------------------------------------------------------------------------------
!p.multi = 0
IF (keyword_set(outfile)) THEN device,filename=outfile+'_btf2.eps',/color,bits_per_pixel= 8,/times,xsize = 10,ysize = 12,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 712,ysize = 392
loadct,obsct
plot,V_f,M_star + M_gas,psym = obssym,/xlog,/ylog,xrange = [30,500],yrange = [1e8,4e11],xtitle = textoidl('V_{f} [km/s]'),ytitle = textoidl('Baryonic Mass [M')+sunsymbol()+']',/nodata
oplot,V_f,M_star + M_gas,psym = obssym,color = obscolor,thick = obssymthick
IF keyword_set(key) THEN BEGIN
    l_keys = lb_keys
    l_keys[0] = 'McGaugh et al. 2008'
    l_color = lb_color
    l_color[0] = obscolor
    l_symbols = lb_symbols
    l_symbols[0]  = obssym
    l_symsizes = lb_symsizes
    l_symsizes[0]  = obssymsize  
    l_symthick = lb_symthick
    l_symthick[0] = thick
    legend,l_keys,color = l_color,symsize = l_symsizes,psym = l_symbols,thick = l_symthick,/right,/bottom,box = 0
ENDIF
IF keyword_set(COLOR) THEN BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO BEGIN
        loadct,ctables[i]
        oplot,[velocities[i],velocities[i]],[bmass[i],bmass[i]],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i],color = colors[i]
      IF keyword_set(do_true) THEN  oplot,[velocities_true[i],velocities_true[i]],[bmass_true[i],bmass_true[i]],color = colors[i],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i]/2.0
      IF keyword_set(key) and i lt n_elements(key) THEN BEGIN
          l_keys = lb_keys
          l_keys[i + 1] = key[i]
          l_color = lb_color
          l_color[i + 1] = colors[i]
          l_symbols = lb_symbols
          l_symbols[i + 1]  = symbols[i]
          l_symsizes = lb_symsizes
          l_symsizes[i + 1]  = symsizes[i]   
          l_symthick = lb_symthick
          l_symthick[i + 1] = thick
          legend,l_keys,color = l_color,psym = l_symbols,thick = l_symthick,/right,/bottom,box = 0
      ENDIF
    ENDFOR
ENDIF ELSE BEGIN
    FOR i = 0, n_elements(broadband) - 1 DO BEGIN
        loadct,ctables[i]
        oplot,[velocities[i],velocities[i]],[bmass[i],bmass[i]],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i],color = colors[i]
        IF keyword_set(do_true) THEN oplot,[velocities_true[i],velocities_true[i]],[bmass_true[i],bmass_true[i]],color = colors[i],psym = symcat(symbols[i],thick = thick),SYMSIZE = symsizes[i]/2.0
        IF keyword_set(key) and i lt n_elements(key) THEN BEGIN
            l_keys = lb_keys
            l_keys[i + 1] = key[i]
            l_color = lb_color
            l_color[i + 1] = colors[i]
            l_symbols = lb_symbols
            l_symbols[i + 1]  = symbols[i]
            l_symsizes = lb_symsizes
            l_symsizes[i + 1]  = symsizes[i] 
            l_symthick = lb_symthick
            l_symthick[i + 1] = thick
            legend,l_keys,color = l_color,psym = l_symbols,thick = thick,/right,/bottom,box = 0
        ENDIF
    ENDFOR
ENDELSE
IF (keyword_set(outfile)) THEN device,/close ELSE stop
END
