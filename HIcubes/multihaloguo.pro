PRO multihaloguo, files, outplot = outplot, color = color, psym = psym, key = key,xrange = xrange, yrange = yrange, imf = imf, kcorrect = kcorrect, symsize = symsize,title = title,halomasses = halomasses,haloids = haloids,overplot = overplot,obscolor = obscolor,bw = bw,moster = moster, things = things,tags = tags,offset = offset,mainhalos = mainhalos,obssym = obssym,ctables = ctables
;Kroupa: IMF = 1, MS: IMF = 2 
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23

IF NOT keyword_set(IMF) THEN IMF = 1
IF NOT keyword_set(offset) THEN offset = [0.1,-0.2]

IF NOT keyword_set(overplot) THEN BEGIN
;    formatplot,outplot = outplot
    IF keyword_set(outplot) THEN BEGIN
        IF NOT strmatch(outplot,'*eps') THEN $
          IF keyword_set(kcorrect) THEN outplot = outplot + '_MstarMhalo_Obs.eps' ELSE outplot = outplot + '_MstarMhalo.eps'
        formatplot,outplot = outplot
        nbins=100.0
;   device, filename=outbase+'LG.'+step[0]+'.guo.ps',/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch
        device, filename=outplot,/COLOR,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch
    ENDIF
ENDIF

IF keyword_set(outplot) THEN fgcolor = 0 ELSE fgcolor = 255

n = n_elements(files)
mainhalos = fltarr(n,2)
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
IF n_elements(symsize) EQ 1 THEN symsize = fltarr(n_elements(files)) + symsize 
IF NOT keyword_set(xrange) THEN xrange = [9.5,13]
IF NOT keyword_set(yrange) THEN yrange = [6,11]

x = findgen(700)/100 + 9
y = alog10(0.1257*((10^x/10^(11.36))^(-0.9147) + (10^x/10^(11.36))^(0.2485))^(-2.574)*10^x) ;guo 09, (3)
;plot,x,y,xstyle = 1, ystyle = 1,xrange = xrange,yrange = yrange,xtitle = textoidl('log(M_{halo}[M')+sunsymbol()+'])',ytitle = textoidl('log(M_{*}[M')+sunsymbol()+'])',title = title,xmargin = xmargin,ymargin = ymargin

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Code/Datafiles/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN prefix = '/home1/crchrist/Datafiles/' $
ELSE prefix = '/astro/users/christensen/code/Datafiles/'

IF NOT keyword_set(overplot) THEN plot,[0,0],[0,0],xstyle = 1, ystyle = 1,xrange = xrange,yrange = yrange,xtitle = textoidl('log(M_{halo}[M')+sunsymbol()+'])',ytitle = textoidl('log(M_{*}[M')+sunsymbol()+'])',title = title,xmargin = xmargin,ymargin = ymargin,/nodata
IF keyword_set(moster) THEN BEGIN
    readcol,prefix + 'HIcubes/'+moster,mhalo,mstar_mhalo,mstar,logmstar
    oplot,mhalo,logmstar,color = obscolor
ENDIF
readcol,prefix + 'HIcubes/ThingsDwarfs.txt',name,D,inc,z0,MB,Vmax,Rmax,Mdyn_THINGS,Mhalo_THINGS,Mstar_THINGS,Mgas_THINGS,format = '(A,F10)'
IF keyword_set(things) THEN oplot,alog10(Mhalo_THINGS), alog10(Mstar_THINGS),color = obscolor,psym = symcat(obssym),symsize = symsize[0]
FOR j = 0, n_elements(files) -1 DO BEGIN
    IF NOT keyword_set(haloids) THEN haloids_arr = indgen(1000) ELSE haloids_arr = haloids[*,j]
    file = files[j]
;    print,file
    stat = read_stat_struc_amiga(file + '.amiga.stat')
    IF keyword_set(sat) THEN BEGIN
        central = where(stat.sat EQ 'no')
 ;       print,n_elements(stat),n_elements(central)
        stat = stat[central]
    ENDIF
    match,stat.group,haloids_arr,ind1,ind2
    haloids_arr = haloids_arr[ind2]
    stat = stat[ind1]

    nhalos = n_elements(stat)
    guo_mass = fltarr(3,nhalos)    
    mainhalos[j,*] = [stat[0].m_dark,stat[0].m_star]
    FOR i = 0, nhalos - 1 DO BEGIN
        guo_mass[2,i] = [stat[i].m_star]
    ENDFOR

    IF (keyword_set(kcorrect)) THEN BEGIN
        IF IMF eq 1 THEN mags = mrdfits(file + '.amiga_r200.halos.star99_K_ab.Mv.fits',1,/silent) ELSE mags = MRDFITS(file + '.amiga_vir.halos.star99_MS_ab.Mv.fits',1,/silent)
;        hstar = where(mags.u ne 0 AND cont eq 'clean')
        match,stat.group,mags.id,ind1,ind2
        mags = mags[ind2]
        stat = stat[ind1]

        mags_bes = transpose([[mags.u],[mags.b],[mags.v],[mags.r],[mags.i]])
        mags_errs = fltarr(5,n_elements(mags.u))+0.02
        mags_ivar=1./mags_errs^2
        dmod = 19.4576          ;distance modulus for z = 0 in this model
        z = fltarr(n_elements(mags.u))
        mgy =(10.D)^(-(0.4D)*(mags_bes + dmod))
        mgy_ivar = mags_ivar/(0.4*alog(10.)*mgy)^2.
        kcorrect, mgy, mgy_ivar, z, kcorrect, mass = star, mtol = mtol_k, absmag = absmag_k, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']
        stat.m_star = star
    ENDIF

    FOR i = 0, nhalos - 1 DO BEGIN
        guo_mass[0,i] = stat[i].m_dark
        guo_mass[1,i] = stat[i].m_star
    ENDFOR
    IF keyword_set(ctables) THEN loadct,ctables[j]
    FOR k = 0, nhalos - 1 DO BEGIN
        IF alog10(stat[k].m_star) GE 0 THEN BEGIN
            oplot,[alog10(guo_mass[0,k]),alog10(guo_mass[0,k])],[alog10(guo_mass[1,k]),alog10(guo_mass[1,k])],psym = symcat(psym[j]),color = color[j],symsize = symsize[j] 
            oplot,[alog10(guo_mass[0,k]),alog10(guo_mass[0,k])],[alog10(guo_mass[2,k]),alog10(guo_mass[2,k])],psym = 2,color = 250,symsize = symsize[j]/2
;           print,k,stat[k].group,alog10(stat[k].m_dark),alog10(stat[k].m_star)
        ENDIF
    ENDFOR
    IF keyword_set(tags) THEN xyouts,alog10(guo_mass[0,0]) + offset[0],alog10(guo_mass[1,0]) + offset[1],tags[j]
    IF j eq 0 THEN BEGIN
        halomasses = [[stat.m_dark],[stat.m_star]]
    ENDIF $
    ELSE BEGIN
        IF n_elements(halomasses[*,1]) LT nhalos THEN halomasses[*,1] = stat[0:n_elements(halomasses[*,1])-1].m_star/halomasses[*,1] $
        ELSE BEGIN
            tstar = fltarr(n_elements(halomasses[*,1]))
            tstar[0:nhalos -1 ] = stat.m_star
            halomasses[0:nhalos -1,1] = stat.m_star/halomasses[0:nhalos -1,1] 
        ENDELSE
    ENDELSE
    print,[[stat.group],[reform(alog10(guo_mass[0,0:nhalos - 1]))],[reform(alog10(guo_mass[1,0:nhalos - 1]))],[reform(alog10(guo_mass[2,0:nhalos - 1]))],[reform(guo_mass[1,0:nhalos - 1]/guo_mass[2,0:nhalos - 1])]]
ENDFOR
;IF keyword_set(key) THEN legend,['Moster','THINGS',key],psym = [-3,symcat(obssym),symcat(psym)],color = [fgcolor,obscolor,color],/bottom,/right
IF keyword_set(outplot) THEN device,/close ;ELSE stop
END
