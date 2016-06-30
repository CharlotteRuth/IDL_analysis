;prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/'
;prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/'
;files = prefix +['h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.out/h516.cosmo25cmb.3072g14HBWK.halo.1','h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.out/h516.cosmo25cmb.3072g1MBWK.00444.halo.1']
;coolrate,files,[25000.,25000.],[2.310e15,2.310e15],keys = ['DH2','DnoH2'],outplot = '~/Plots/h516.cosmo25cmb.3072g.paper'

pro coolrate,files,distunits,massunits,outplot = outplot,color = color,keys = keys
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
ZSOLAR    = 0.0130215

formatplot,outplot = outplot
loadct,0
n = N_ELEMENTS(files)
IF KEYWORD_SET(outplot) THEN BEGIN
    fgcolor = 0
    grey = 100 ;20
    device,filename = outplot + '_coolcurve.eps',/color,bits_per_pixel= 8,/times,ysize=18,xsize=18,xoffset =  2,yoffset =  2
ENDIF ELSE BEGIN
    fgcolor = 255
    grey = 200    
    window,0,xsize = 800, ysize = 800
ENDELSE
IF keyword_set(color) THEN BEGIN
    loadct,39
    grey = 240
    fgcolor = 50
ENDIF
color = [grey,fgcolor]
;X = [-0.5, 0, 0.5, 0, -0.5]  
;Y = [0, 0.5, 0, -0.5, 0]  
;USERSYM, X, Y,/FILL

multiplot,[2,1],/doxaxis,/doyaxis,mxtitle = 'T [K]',mytitle = textoidl('Cooling Rate [erg cm^3/s]'),mxtitsize = 1.5,mytitsize = 1.5,mxtitoffset = 2, mytitoffset = 3
FOR i = 0, n-1 DO BEGIN
   rtipsy,files[i],h,g,d,s
   read_tipsy_arr,files[i] + '.eCool',h,cooling,part = 'gas'
;   readarr,files[i] + '.HI',h,HI,part = 'gas',/ascii
;   readarr,files[i] + '.H2',h,H2,part = 'gas',/ascii
;   readarr,files[i] + '.OxMassFrac',h,ox,part = 'gas',/ascii
;   atomicgas_mass = (HI + 2.0*H2)/MAX(HI + 2.0*H2)*g.mass*massunits[i]
;   h_mass = (HI + 2.0*H2)*g.mass*massunits[i]
;   nh = (h_mass/1.6726d-24)*(2.0d33)
;   n_ox = (ox*atomicgas_mass)*2.0d33/2.66d-23 
;   ox_abund = n_ox/nh
;   cooling_ox = 2.89e-20*10^(1.07 + 0.43*alog10(ox_abund))*ox_abund
   
;units of ergs per gm per s
   cooling = cooling/5.9753790e+23 ;ergs per s
   dens_convert =  massunits[i] * gm_per_msol * 5.9753790e+23/distunits[i]^3/cm_per_kpc^3/h.time^3
;    if i eq 0 then histogramp,alog10(g.zmetal/zsolar),nbins = 100,min = -4,max = 1 else histogramp,alog10(g.zmetal/zsolar),/overplot,color = color[i],nbins = 100,min = -4,max = 1
   indlow = where(alog10(g.zmetal/zsolar) gt -2.5 AND alog10(g.zmetal/zsolar) lt -1.5)
   indhigh = where(alog10(g.zmetal/zsolar) gt -1.5 AND alog10(g.zmetal/zsolar) lt -0.5)
   
   if i eq 0 then plot, g.tempg,cooling/g.dens/dens_convert,psym = symcat(16),xrange = [10,6e7],yrange = [1e-30,1e-20],symsize = 0.2,/nodata,/ylog,/xlog;,xtitle = 'T [K]',ytitle = textoidl('Cooling Rate [erg cm^3/s]')
   oplot,g.tempg,cooling/g.dens/dens_convert,psym = symcat(16),color = color[i],symsize = 0.2
   
;   if i eq 0 then plot, g[indhigh].tempg,cooling[indhigh]/g[indhigh].dens/dens_convert,psym = symcat(16),/ylog,/xlog,xtitle = 'T [K]',ytitle = textoidl('Cooling Rate [erg cm^3/s]'),xrange = [10,6e7],yrange = [1e-30,1e-20],symsize = 0.2,/nodata
;    oplot,g[indhigh].tempg,cooling[indhigh]/g[indhigh].dens/dens_convert,psym = symcat(16),color = color[i],symsize = 0.2
   
;    if i eq 0 then plot, g[indhigh].dens,g[indhigh].tempg*dens_convert,psym= symcat(16),/ylog,/xlog
;    oplot,g[indhigh].dens,g[indhigh].tempg*dens_convert,psym = symcat(16),color =    color[i]
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = color,psym = [8,8],/right,/bottom,symsize = [0.6, 0.6]
multiplot,/doxaxis
FOR i = 0, 0 DO BEGIN
    rtipsy,files[i],h,g,d,s
    read_tipsy_arr,files[i] + '.eCool',h,cooling,part = 'gas'
    read_tipsy_arr,files[i] + '.eCoolH2',h,coolingH2,part = 'gas'
    readarr,files[i] + '.HI',h,HI,part = 'gas',/ascii
    readarr,files[i] + '.H2',h,H2,part = 'gas',/ascii
    readarr,files[i] + '.OxMassFrac',h,ox,part = 'gas',/ascii
    atomicgas_mass = (HI + 2.0*H2)/MAX(HI + 2.0*H2)*g.mass*massunits[i]
    h_mass = (HI + 2.0*H2)*g.mass*massunits[i]
    nh = (h_mass/1.6726d-24)*(2.0d33)
    n_ox = (ox*atomicgas_mass)*2.0d33/2.66d-23 
    ox_abund = n_ox/nh
    cooling_ox = 2.89e-20*10^(1.07 + 0.43*alog10(ox_abund))*ox_abund

;units of ergs per gm per s
    cooling = cooling/5.9753790e+23 ;ergs per s
    coolingH2 = coolingH2/5.9753790e+23 ;ergs per s
    dens_convert =  massunits[i] * gm_per_msol * 5.9753790e+23/distunits[i]^3/cm_per_kpc^3/h.time^3
    indlow = where(alog10(g.zmetal/zsolar) gt -2.5 AND alog10(g.zmetal/zsolar) lt -1.5)
    indhigh = where(alog10(g.zmetal/zsolar) gt -1.5 AND alog10(g.zmetal/zsolar) lt -0.5)

    if i eq 0 then plot, g.tempg,cooling/g.dens/dens_convert,psym = symcat(16),/ylog,/xlog,xrange = [10,6e7],yrange = [1e-30,1e-20],symsize = 0.2,/nodata;,xtitle = 'T [K]',ytitle = textoidl('Cooling Rate [erg cm^3/s]')
    indcold = where(g.tempg LE 8e3)
 ;   cooling_oxcold = cooling_ox*0
 ;   cooling_oxcold[indcold] = cooling_ox[indcold]
    oplot,g.tempg,cooling/g.dens/dens_convert,psym = symcat(16),symsize = 0.2,color = color[i]
    oplot,g[indcold].tempg,cooling_ox[indcold]/g[indcold].dens/dens_convert,psym = symcat(16),symsize = 0.2,color = color[1]
;    oplot,g[indcold].tempg,coolingH2[indcold]/g[indcold].dens/dens_convert,psym = symcat(16),symsize = 0.2,color = grey
ENDFOR
IF KEYWORD_SET(keys) THEN legend,[keys[0],'Estimated CII cooling'],color = color,psym = [8,8],/right,/bottom,symsize = [0.6, 0.6]
multiplot,/reset
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
IF 0 THEN BEGIN
multiplot,[1,1]
IF keyword_set(outplot) THEN device,filename = outplot + '_coolcomp.eps',/color,bits_per_pixel= 8,/times,ysize=18,xsize=18,xoffset =  2,yoffset =  2 ELSE window,1,xsize = 600,ysize = 600
plot,cooling_ox[indcold]/g[indcold].dens/dens_convert,coolingH2[indcold]/g[indcold].dens/dens_convert,/xlog,/ylog,xrange = [1e-29,1e-23],yrange = [1e-29,1e-23],xtitle = 'Estimated CII Cooling [erg cm^3/s]',ytitle = textoidl('H_2') + ' [erg cm^3/s]',psym = symcat(16),symsize = 0.2
oplot,[1e-30,1e-20],[1e-30,1e-20],color = grey
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
multiplot,/reset
ENDIF
end
