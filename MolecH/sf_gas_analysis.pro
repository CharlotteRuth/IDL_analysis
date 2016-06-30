pro sf_gas_analysis_master
base = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
slfiles = base + ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.starlog', $
           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g15HBWK/h516.cosmo25cmb.1536g15HBWK.starlog', $
           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK_noJeans.starlog', $
           'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/h516.cosmo25cmb.1536g1MBWK.starlog']
unitfiles = base + ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/h516.cosmo25cmb.1536g14HBWK.param', $
             'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g15HBWK/h516.cosmo25cmb.1536g15HBWK.param', $
             'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK_noJeans.param', $
             'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/h516.cosmo25cmb.1536g1MBWK.param']
mh = [1,1,0,0]
densunit = DOUBLE(mh)
FOR i = 0, N_ELEMENTS(unitfiles) - 1 DO BEGIN
    units = tipsyunits(unitfiles[i])
    densunit[i] = units.rhounit
ENDFOR

sf_gas_analysis,slfiles,densunit,mh,/color

end

pro sf_gas_analysis,slfiles,densunit,mh,color = color
IF KEYWORD_SET(outfile) THEN BEGIN
    set_plot,'ps' 
    nbins=100.0
    linestyles = [0,2]
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25
;    !p.font=0 
    cb_charsize = 0.75
    l_charsize = 0.75
ENDIF ELSE BEGIN
    set_plot,'x'
    nbins=100.0
    linestyles = [0,2]
     !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    cb_charsize = 1.0
    l_charsize = 1.0
ENDELSE
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
   if color[0] eq 1 then  colors = (findgen(N_ELEMENTS(slfiles)) + 1)*240/N_ELEMENTS(slfiles) else colors = color
    thicks = fltarr(N_ELEMENTS(slfiles)) + 2
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(N_ELEMENTS(slfiles)) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(files)) + 5
    thicks = (findgen(N_ELEMENTS(slfiles)) + 1)*6/N_ELEMENTS(slfiles) - 1
ENDELSE
FOR i = 0, N_ELEMENTS(slfiles) - 1 DO BEGIN
    data = rstarlog(slfiles[i],molecularh = mh[i])
    if i eq 0 then plot,data.rhoform*densunit[i],data.tempform,psym = 3,/xlog,/ylog,xrange = [0.1,1e4],yrange = [10,1e4]
    oplot,data.rhoform*densunit[i],data.tempform,psym = 3,color = colors[i]
    stop
ENDFOR
end
