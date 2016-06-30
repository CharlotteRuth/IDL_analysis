PRO rhoSF_master
;file0 = 'steps_00406/h516.cosmo25cmb.2304g3MBWK_00406.00024.dir/h516.cosmo25cmb.2304g3MBWK_00406.00024'
;starlog0 = 'h516.cosmo25cmb.2304g3MBWK_00406.starlog'
;file0 = 'steps/h516.cosmo25cmb.2304g3MBWK.00456.dir/h516.cosmo25cmb.2304g3MBWK.00456'
;starlog0 = 'h516.cosmo25cmb.2304g3MBWK.starlog'
;dPhysDenMin = 100
;dCStar = 0.1
;dESN = 1e51

;file1 = 'h516.cosmo25cmb.2304g2bwK/h516.cosmo25cmb.2304g2bwK.00512/h516.cosmo25cmb.2304g2bwK.00512'
;starlog1 = 'h516.cosmo25cmb.2304g2bwK/h516.cosmo25cmb.2304g2bwK.00512/h516.cosmo25cmb.2304g2bwK.starlog'
;dPhysDenMin = .1
;dCStar = 0.05
;dESN = 0.4e51

;file2 = 'h516.cosmo25cmb.2304g4bwdK/h516.cosmo25cmb.2304g4bwdK.00512/h516.cosmo25cmb.2304g4bwdK.00512'
;starlog2 = 'h516.cosmo25cmb.2304g4bwdK/h516.cosmo25cmb.2304g4bwdK.starlog'
;dPhysDenMin = 100.
;dCStar = 0.1
;dESN = .4e51

;file3 = 'h516.cosmo25cmb.1536g1MBWK/h516.cosmo25cmb.1536g1MBWK.00512.dir/h516.cosmo25cmb.1536g1MBWK.00512'
;starlog3 = 'h516.cosmo25cmb.1536g1MBWK/h516.cosmo25cmb.1536g1MBWK.starlog'
;dPhysDenMin = 100.
;dCStar = 0.1
;dESN = .8e51

;file4 = 'h516.cosmo25cmb.3072g1MBWK.00512/h516.cosmo25cmb.3072g1MBWK.00512'
;file4 ='h516.cosmo25cmb.3072g1MBWK.00512/../h516.cosmo25cmb.3072g1MBWK.00108/h516.cosmo25cmb.3072g1MBWK.00108'
;starlog4 = 'h516.cosmo25cmb.3072g1MBWK.00512/../h516.cosmo25cmb.3072g1MBWK.starlog'
;dPhysDenMin = 100.
;dCStar = 0.1
;dESN = 1.0e51


dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g3MBWK_SFrho'
cd,dir
pfile = 'h516.cosmo25cmb.2304g3MBWK.param'
units = tipsyunits(pfile)
base = ['h516.cosmo25cmb.2304g3MBWK','h516.cosmo25cmb.2304g2MBWK']
stepname = ['steps','steps.2']
color = [50,240]
key = ['2304g3MBWK Rho SF','2304g2MBWK']
step = ['00024','00036','00048','00072','00084','00096','00108','00120','00132','00156','00180','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00456']
nsteps = N_ELEMENTS(step) - 1

FOR i = nsteps, nsteps DO BEGIN
    files = stepname + '/' + base + '.' + step[i] + '.dir/' + base + '.' + step[i] + '.halo.1'
    rhoSF,files,base,units,color = color, key = key, title = step[i],outfile = 'h516.2304_RhoSF_Rho1.5SF_'+step[i]
    rotCurve,files,units.massunit,units.lengthunit,keys = keys,color = color,outfile = 'h516.2304_RhoSF_Rho1.5SF_'+step[i]
ENDFOR

END

PRO rhoSF,files,base,units,color = color,key = key,title = title, outfile = outfile
loadct,39
IF KEYWORD_SET(outfile) THEN BEGIN
    set_plot,'ps' 
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.0 ;0.75
    !y.charsize=1.0 ;0.75
    l_charsize = 0.75
;    !p.font=0 
ENDIF ELSE BEGIN
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5 
     l_charsize = 1.0
ENDELSE

IF NOT keyword_set(color) THEN color = (findgen(N_ELEMENTS(files)) + 1)*240/N_ELEMENTS(files)
IF NOT keyword_set(key) THEN key = 'file ' + STRTRIM(findgen(N_ELEMENTS(files)))

rtipsy,files[0],h0,g0,d0,s0
starlog0 = base[0] + '.starlog'
rtipsy,files[1],h1,g1,d1,s1
starlog1 = base[1] + '.starlog'

IF KEYWORD_SET(outfile) THEN resSchmidtLaw,['.','.'],files,[units.lengthunit,units.lengthunit],[units.massunit,units.massunit],key,color = color, outplot = outfile+'_resSK.eps' ELSE BEGIN
    window,0
    resSchmidtLaw,['.','.'],files,[units.lengthunit,units.lengthunit],[units.massunit,units.massunit],key,color = color
ENDELSE

;------------------------------------------------------------------------------------
IF KEYWORD_SET(outfile) THEN device,filename=outfile+'_histRhoSF.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 ELSE window,3,xsize = 712,ysize = 392

stardata0 = rstarlog(starlog0)
stardata1 = rstarlog(starlog1)
histrho0 = histogram(alog10(stardata0.rhoform*units.RHOUNIT),locations = x,min = 0.75,max = 3.75,nbins = 50)/1.0/N_ELEMENTS(stardata0)
histrho1 = histogram(alog10(stardata1.rhoform*units.RHOUNIT),locations = x,min = 0.75,max = 3.75,nbins = 50)/1.0/N_ELEMENTS(stardata0)

rtipsy,files[0],h0,g0,d0,s0
rtipsy,files[1],h1,g1,d1,s1
histrhog0 = histogram(alog10(g0.dens*units.RHOUNIT/h0.time/h0.time/h0.time),locations = x,min = 0.75,max = 3.75,nbins = 75)/1.0/N_ELEMENTS(g0)
histrhog1 = histogram(alog10(g1.dens*units.RHOUNIT/h1.time/h1.time/h1.time),locations = x,min = 0.75,max = 3.75,nbins = 75)/1.0/N_ELEMENTS(g1)

maxSF = max(s0.tform*units.TIMEUNIT)
minSF = maxSF - 5e8;min(stardata0.timeform) ;0.07
ind = where(stardata0.timeform*units.TIMEUNIT ge minSF AND stardata0.timeform*units.TIMEUNIT le maxSF)
stardata0new = stardata0[ind]
ind = where(stardata1.timeform*units.TIMEUNIT ge minSF AND stardata1.timeform*units.TIMEUNIT le maxSF)
stardata1new = stardata1[ind]

histrho0new = histogram(alog10(stardata0new.rhoform*units.RHOUNIT),locations = x,min = 0.75,max = 3.75,nbins = 50)/1.0/N_ELEMENTS(stardata0new)
histrho1new = histogram(alog10(stardata1new.rhoform*units.RHOUNIT),locations = x,min = 0.75,max = 3.75,nbins = 50)/1.0/N_ELEMENTS(stardata1new)

plot,x,histrho0,/ylog,xrange = [1.0,4.0],yrange = [1e-6,1],psym = 10,xstyle = 1,xtitle = 'SF density [amu/cc]',ytitle = 'N_stars/TOTAL(stars)',title = title
oplot,x,histrho0,psym = 10,color = color[0]
oplot,x,histrho1,psym = 10,color = color[1]
oplot,x,histrhog0,psym = 10,linestyle = 2,color = color[0]
oplot,x,histrhog1,psym = 10,linestyle = 2,color = color[1]
oplot,x,histrho0new,psym = 10,linestyle = 1, color = color[0]
oplot,x,histrho1new,psym = 10,linestyle = 1, color = color[1]
legend,key,color = color,linestyle = fltarr(N_ELEMENTS(color)),charsize = l_charsize
legend,['Star Forming Gas','Recent Star Forming Gas','Current Gas'],linestyle = [0,1,2],/right,charsize = l_charsize

;-----------------------------------------------------------------------------------
IF KEYWORD_SET(outfile) THEN BEGIN
    device,/close
    device,filename=outfile+'_sfh_haloall.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 
ENDIF ELSE window,2,xsize = 712,ysize = 392

sfr,s0,massunit = units.massunit,timeunit = units.timeunit,binsize = 1e8,yrange = [0,0.2],xrange = [0,15],title = title;,xrange = [minSF*units.timeunit/1e9,maxSF*units.timeunit/1e9],xstyle = 1
sfr,s0,/overplot,color = color[0],massunit = units.massunit,timeunit = units.timeunit,binsize = 1e8
sfr,s1,/overplot,color = color[1],massunit = units.massunit,timeunit = units.timeunit,binsize = 1e8

sdat0 = stardata0
sdat1 = stardata1
sfr,sdat0,massunit = units.massunit,timeunit = units.timeunit,binsize = 1e8,yrange = [0,0.2],/starlog,linestyle = 1,/overplot
sfr,sdat0,/overplot,color = color[0],massunit = units.massunit,timeunit = units.timeunit,binsize = 1e8,linestyle = 1
sfr,sdat1,/overplot,color = color[1],massunit = units.massunit,timeunit = units.timeunit,binsize = 1e8,/starlog,linestyle = 1
legend,[key[0] + ' ' + STRMID(STRTRIM(TOTAL(s0.mass)*units.massunit,2),0,4)+STRMID(STRTRIM(TOTAL(s0.mass)*units.massunit,2),3,/REVERSE_OFFSET), $
        key[1] + ' ' + STRMID(STRTRIM(TOTAL(s1.mass)*units.massunit,2),0,4)+STRMID(STRTRIM(TOTAL(s1.mass)*units.massunit,2),3,/REVERSE_OFFSET)]$
  ,color = color,thick = thicks,/top,/right,linestyle = fltarr([N_ELEMENTS(color)]),charsize = l_charsize 
;oplot,[maxSF,maxSF]/1e9,[0,1]
IF KEYWORD_SET(outfile) THEN  device,/close ELSE stop
close,/all
END
