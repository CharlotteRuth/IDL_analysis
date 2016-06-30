PRO baryonicfrac, files, haloids_arr, gmass, imf = imf, kcorrect = kcorrect, xrange = xrange,  color = color, ctables = ctables, psym = psym,symsizes = symsizes, thicks = thicks, key = key, outplot = outplot,yrange = yrange

f_bar = 0.16510

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Code/Datafiles/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN prefix = '/home1/crchrist/Datafiles/' $
ELSE prefix = '/astro/users/christensen/code/Datafiles/'

IF NOT keyword_set(IMF) THEN IMF = 1

formatplot,outplot = outplot
IF keyword_set(outplot) THEN BEGIN
    IF NOT strmatch(outplot,'*eps') THEN $
      IF keyword_set(kcorrect) THEN outplot = outplot + '_baryfrac_obs.eps' ELSE outplot = outplot + '_baryfrac.eps'
    device, filename=outplot,/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18;ysize=5,xsize=7,/inch
    fgcolor = 0
ENDIF ELSE BEGIN
    fgcolor = 255
ENDELSE

n = n_elements(files)
mainhalos = fltarr(n,2)
IF keyword_set(color) THEN BEGIN
    IF NOT keyword_set(bw) THEN loadct,39 ELSE loadct,0
    IF color[0] EQ 1 THEN  color  = (findgen(n) + 1)*240/n ELSE color = color
    bcolor = [254,60,140]
    IF NOT keyword_set(psym) THEN psym = fltarr(n) + 4
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 4
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(obscolor) THEN obscolor = 100
    color = (findgen(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(n_elements(broadband)) + 5
    bcolor = [fgcolor,fgcolor,fgcolor]
    IF NOT keyword_set(psym) THEN  psym = fltarr(n) + 4;(findgen(n)+2)*2
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 4
ENDELSE
IF NOT keyword_set(symsize) THEN symsize = 1.5
IF n_elements(symsize) EQ 1 THEN symsize = fltarr(n_elements(files)) + symsize 
IF NOT keyword_set(xrange) THEN xrange = [9.5,13]
IF NOT keyword_set(yrange) THEN yrange = [1e-2,1]
psym2 = psym
psym[where(psym EQ 6)] = 4
psym2[where(psym EQ 14)] = 15
psym3 = psym
psym3[where(psym2 EQ 6)] = 9
psym3[where(psym EQ 14)] = 16
last = 0
plot,[0,0],[0,0],xstyle = 1, ystyle = 1,xrange = xrange,yrange = yrange,xtitle = textoidl('log(M_{halo}[M')+sunsymbol()+'])',ytitle = textoidl('M/f_B M_{vir}'),title = title,xmargin = xmargin,ymargin = ymargin,/nodata,/ylog
FOR j = 0, n_elements(files) -1 DO BEGIN
    stat = read_stat_struc_amiga(files[j] + '.amiga.stat')
    haloids = haloids_arr[*,j]
    match,stat.group,haloids,ind1,ind2
    haloids = haloids[ind2]
    stat = stat[ind1]

    IF (keyword_set(kcorrect)) THEN BEGIN
        IF IMF eq 1 THEN mags = mrdfits(files[j] + '.amiga_r200.halos.star99_K_ab.Mv.fits',1) ELSE mags = MRDFITS(file + '.amiga_vir.halos.star99_MS_ab.Mv.fits',1)
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
        kcorrect, mgy, mgy_ivar, z, kcorrect, mass = smass, mtol = mtol_k, absmag = absmag_k, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']
    ENDIF ELSE smass = stat.m_star

    loadct,ctables[j]
    IF n_elements(stat) NE 1 THEN $
      oplot,alog10(stat.m_tot),(stat.m_gas + stat.m_star)/(stat.m_tot*f_bar),color = bcolor[0],psym = symcat(psym[j]),symsize = symsize[j],thick = thicks[j] ELSE $
        oplot,[alog10(stat.m_tot),alog10(stat.m_tot)],[(stat.m_gas + stat.m_star)/(stat.m_tot*f_bar),(stat.m_gas + stat.m_star)/(stat.m_tot*f_bar)],color = bcolor[0],psym = symcat(psym[j]),symsize = symsize[j],thick = thicks[j]

    IF n_elements(stat) NE 1 THEN $
      oplot,alog10(stat.m_tot),(smass + gmass[(last):(last + n_elements(stat) - 1)])/(stat.m_tot*f_bar),color = bcolor[1],psym = symcat(psym2[j]),symsize = symsize[j],thick = thicks[j] ELSE $
        oplot,[alog10(stat.m_tot),alog10(stat.m_tot)],[(smass + gmass[(last):(last + n_elements(stat) - 1)])/(stat.m_tot*f_bar),(smass + gmass[(last):(last + n_elements(stat) - 1)])/(stat.m_tot*f_bar)],color = bcolor[1],psym = symcat(psym2[j]),symsize = symsize[j],thick = thicks[j]

    IF n_elements(stat) NE 1 THEN $
      oplot,alog10(stat.m_tot),(smass)/(stat.m_tot*f_bar),color = bcolor[2],psym = symcat(psym3[j]),symsize = symsize[j],thick = thicks[j] ELSE $
        oplot,[alog10(stat.m_tot),alog10(stat.m_tot)],[(smass)/(stat.m_tot*f_bar),(smass)/(stat.m_tot*f_bar)],color = bcolor[2],psym = symcat(psym3[j]),symsize = symsize[j],thick = thicks[j]

    last = last + n_elements(stat.m_tot)
ENDFOR
legend,['All Bayrons','Disk Baryons','Stars'],psym = [4,6,9],color = bcolor,symsize = symsize[0]*(fltarr(3) + 1),thick = thicks[0]*(fltarr(3) + 1),linestyle = fltarr(3),/right,/bottom,box = 0
IF keyword_set(outplot) THEN device,/close ELSE stop
END
