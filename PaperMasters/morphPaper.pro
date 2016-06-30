;morphPaper,outplot = '~/Plots/Bulges/all'
;morphPaper,outplot = '~/Plots/Bulges/e11'

PRO morphPaper,outplot = outplot,color = color,formatthick = formatthick

formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF keyword_set(outplot) THEN bgcolor = 255 ELSE bgcolor = 0
loadct,39
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN prefix = '/nobackupp2/crchrist/MolecH/' $
ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
IF keyword_set(outplot) THEN BEGIN
    l_charsize = 1.0
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    l_charsize = 0.75
    fgcolor = 255
    bgcolor = 0
ENDELSE

lu_c50 = 50000
mu_c50 = 1.84793e16
tu_c50 = 3.87815e+10

find_vrot = 0
find_imag = 0
find_gmass = 0
plot_tully_fisher_obs = 0
size_mag_bulge = 0
size_mag_disk = 0
mu_size = 0
mu_mag = 0
mu_size_mag = 0
mu_size_mag2 = 1
plot_bulgeSFH = 0

x = 2
CASE x OF
   1: BEGIN
        dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
        file986 = 'h986.cosmo50cmb.3072g14HBWK'
        key986 = 'h986'
        dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
        file603 = 'h603.cosmo50cmb.3072g14HBWK'
        key603 = 'h603'
        dir277 = prefix + 'h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
        file277 = 'h277.cosmo50cmb.3072g14HMbwK'
        key277 = 'h277'
        dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
        file258 = 'h258.cosmo50cmb.3072g14HMbwK'
        key258 = 'h258'
        dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
        file285 = 'h285.cosmo50cmb.3072g14HMbwK'
        key285 = 'h285'
        dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
        file239 = 'h239.cosmo50cmb.3072g14HMbwK'
        key239 = 'h239'
        dirs  = [dir986,dir603,dir277,dir258,dir285,dir239]
        files = [file986,file603,file277,file258,file285,file239]
        steps = ['00512','00512','00512','00512','00512','00512']
        haloids = ['1','1','1','1','1','1']
        distunits = [fltarr(6) + lu_c50]
        massunits = [fltarr(6) + mu_c50]
        timeunits = [fltarr(6) + tu_c50]
        key  = [key986,key603,key277,key258,key285,key239] ;+ ', ' + haloid

        colors = fltarr(n_elements(files)) + 254
        ctables = fltarr(n_elements(files)) + 39
        psym = fltarr(n_elements(files)) + 7
        symsizes = fltarr(n_elements(files)) + 2
        thicks = fltarr(n_elements(files)) + 5

        obscolor1 = 60
        obscolor2 = 100
        obspsym = 14
        obssymsize = 1.5
        obssymsize2 = 2.5
        OBSSYMTHICK = 1
        obsct = 39

        camera = 14             ;FO with dust
        band = 7                ; 2MASS H
        filter = 13
        galfitname = 'galfit.H_2MASS'
;        galfitname = 'galfit_H_2MASS.res.feedme'

                  ;[h986,     h603,    h277,    h258,    h285,    h239]
 ;       vfinal = [100.913, 111.279, 213.539, 208.034, 203.057, 204.288] ;Made using vcirc.pro
        vfinal = [103.148, 115.102, 188.999, 181.909, 164.350, 165.470]
        vfinal = [100.913, 120.11,  201.204, 186.46,  189.241, 186.63] ;from fitting a tangent to v_circular
        vfinal = [109.015, 111.279, 205.083, 197.09,  169.28,  143.36] ;fitting a tangent to v_rotational
        vfinal = [101.391, 107.881, 206.796, 193.385, 194.573, 186.637] ;Rotational velocity at 80% i
        vfinal = [99.4566, 110.209, 203.130, 202.557, 189.019, 175.625] ;Fit to rotational velocity at 80% i, (Don't like fit for h239
        vfinal = [99.4566, 110.209, 203.130, 202.557, 189.019, 186.637] ;Fit to rotational velocity at 80% i, Rotational velocity @ 80% for h239
        vflat =  [106.895, 116.062, 193.159, 206.368, 183.045, 156.8];h239 is extrememly uncertain.  Between 156.8 and 195.377
        vHI    = [105.455, 118.880, 205.101, 209.299, 202.902, 211.828]
        imag   = [-19.7434, -20.2580, -21.7940, -21.9892, -21.9059, -22.1313]
        gmass = [ 3.45968e+09, 4.23651e+09, 5.76264e+09, 5.68696e+09, 8.47653e+09, 6.19468e+09]
                ;[h986,     h603,    h277,    h258,    h285,    h239
        outfiles2 = dirs + files + '.' + steps + '/' + files + '.' + steps
        outfiles = dirs + files + '.' + steps + '/' + files + '.' + steps + '.' + haloids + '/broadband.fits'
    END
   2: BEGIN
        dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
        file986 = 'h986.cosmo50cmb.3072g14HBWK'
        key986 = 'h986'
        dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
        file603 = 'h603.cosmo50cmb.3072g14HBWK'
        key603 = 'h603'
        dirs  = [dir986,dir603]
        files = [file986,file603]
        steps = ['00512','00512']
        haloids = ['1','1']
        distunits = [fltarr(n_elements(files)) + lu_c50]
        massunits = [fltarr(n_elements(files)) + mu_c50]
        timeunits = [fltarr(6) + tu_c50]
        key  = [key986,key603] ;+ ', ' + haloid

        colors = fltarr(n_elements(files)) + 254
        ctables = fltarr(n_elements(files)) + 39
        psym = fltarr(n_elements(files)) + 16
        symsizes = fltarr(n_elements(files)) + 2
        thicks = fltarr(n_elements(files)) + 5

        obscolor1 = 60
        obscolor2 = 100
        obspsym = 14
        obssymsize = 2
        obssymsize2 = 2.5
        OBSSYMTHICK = 1
        obsct = 39

        camera = 14             ;FO with dust
        band = 7                ; 2MASS H
        filter = 13
        galfitname = 'galfit.H_2MASS'
;        galfitname = 'galfit_H_2MASS.res.feedme'

                  ;[h986,     h603,    h277,    h258,    h285,    h239]
 ;       vfinal = [100.913, 111.279, 213.539, 208.034, 203.057, 204.288] ;Made using vcirc.pro
        vfinal = [103.148, 115.102]
        vfinal = [100.913, 120.11] ;from fitting a tangent to v_circular
        vfinal = [109.015, 111.279];fitting a tangent to v_rotational
        vfinal = [101.391, 107.881];Rotational velocity at 80% i
        vfinal = [99.4566, 110.209];Fit to rotational velocity at 80% i, (Don't like fit for h239
        vfinal = [99.4566, 110.209] ;Fit to rotational velocity at 80% i, Rotational velocity @ 80% for h239
        vflat =  [106.895, 116.062];h239 is extrememly uncertain.  Between 156.8 and 195.377
        vHI    = [105.455, 118.880]
        imag   = [-19.7434, -20.2580]
        smass_true = [ 4.5257366e+09, 7.8307607e+09]
        gmass = [ 3.45968e+09, 4.23651e+09]

        outfiles2 = dirs + files + '.' + steps + '/' + files + '.' + steps
        outfiles = dirs + files + '.' + steps + '/' + files + '.' + steps + '.' + haloids + '/broadband.fits'
    END
ENDCASE

nfiles = n_elements(files)

fits = {M:0.0,M_B:0.0,r_e:0.0,r_0:0.0,mu:0.0,sersic:0.0,M_D:0.0,r_h:0.0,B_D:0.0}
fits = replicate(fits,nfiles)
fits_galfit = {M:0.0,M_B:0.0,r_e:0.0,r_0:0.0,mu:0.0,sersic:0.0,M_D:0.0,r_h:0.0,B_D:0.0}
fits_galfit = replicate(fits_galfit,nfiles)

h603fits = {M:0.0,M_B:-20.86,r_e:2.92,r_0:2.92,mu:16.35,sersic:10^0.22,M_D:-21.17,r_h:3.53,B_D:0.0}
h986fits = {M:0.0,M_B:-20.81,r_e:3.34,r_0:3.34,mu:18.78,sersic:10^0.40,M_D:-20.68,r_h:3.73,B_D:0.0}
h277fits = {M:0.0,M_B:-22.39,r_e:2.77,r_0:2.77,mu:14.25,sersic:10^0.50,M_D:-23.22,r_h:3.49,B_D:0.0}
h258fits = {M:0.0,M_B:-22.99,r_e:2.74,r_0:2.74,mu:13.42,sersic:10^0.41,M_D:-23.31,r_h:3.64,B_D:0.0}
h285fits = {M:0.0,M_B:-22.90,r_e:2.85,r_0:2.85,mu:14.28,sersic:10^0.41,M_D:-23.14,r_h:3.56,B_D:0.0}
h239fits = {M:0.0,M_B:-22.67,r_e:2.71,r_0:2.71,mu:13.94,sersic:10^0.40,M_D:-23.62,r_h:3.53,B_D:0.0}
h603fits.M = -2.5*alog10(10^(-0.4*h603fits.M_B) + 10^(-0.4*h603fits.M_D))
h986fits.M = -2.5*alog10(10^(-0.4*h986fits.M_B) + 10^(-0.4*h986fits.M_D))
h277fits.M = -2.5*alog10(10^(-0.4*h277fits.M_B) + 10^(-0.4*h277fits.M_D))
h258fits.M = -2.5*alog10(10^(-0.4*h258fits.M_B) + 10^(-0.4*h258fits.M_D))
h285fits.M = -2.5*alog10(10^(-0.4*h285fits.M_B) + 10^(-0.4*h285fits.M_D))
h239fits.M = -2.5*alog10(10^(-0.4*h239fits.M_B) + 10^(-0.4*h239fits.M_D))
h603fits.B_D = 10^(-0.4*h603fits.M_B)/10^(-0.4*h603fits.M_D)
h986fits.B_D = 10^(-0.4*h986fits.M_B)/10^(-0.4*h986fits.M_D)
h277fits.B_D = 10^(-0.4*h277fits.M_B)/10^(-0.4*h277fits.M_D)
h258fits.B_D = 10^(-0.4*h258fits.M_B)/10^(-0.4*h258fits.M_D)
h285fits.B_D = 10^(-0.4*h285fits.M_B)/10^(-0.4*h285fits.M_D)
h239fits.B_D = 10^(-0.4*h239fits.M_B)/10^(-0.4*h239fits.M_D)
fits[0] = h986fits
fits[1] = h603fits
IF n_elements(fits) GT 2 THEN BEGIN
    fits[2] = h285fits
    fits[3] = h277fits
    fits[4] = h258fits
    fits[5] = h239fits
ENDIF

formatplot,outplot = outplot,thick = formatthick

;------------------------------------------- Vrot -------------------------------------
IF find_vrot THEN BEGIN
;    vmax = fltarr(nfiles)
;    FOR i = 0, nfiles - 1 DO BEGIN
;        data = read_stat_struc_amiga(outfiles[i] + '.amiga.stat')
;        ind = where(data.group eq haloids[i])
;        vmax[i] = data[ind].vc
;    ENDFOR

    vrot,dirs+files+'.'+steps+'/'+files+'.'+steps, haloids, massunits, distunits, maxdistance = 25, vreturn = vreturn, /verbose
    print,vreturn
ENDIF

;------------------------------------------- SDSS i mag -------------------------------------
IF find_imag THEN BEGIN
    imags = fltarr(nfiles)
    FOR i = 0, nfiles -1 DO imags[i] = sunrise_mag(dirs[i]+files[i]+'.'+steps[i]+'/'+files[i]+'.'+steps[i],halo_str = haloids[i],band = 13,filtername = "i_SDSS.res")
    print,imags
ENDIF

;------------------------------------------- Gas Mass -------------------------
IF find_gmass THEN gmass = tully_fisher_obs_gasmass(outfiles + '.halo.' + haloid, massunits)

; ------------------------------------------ Tully Fisher -----------------------------------------------------
IF plot_tully_fisher_obs THEN BEGIN
;    tully_fisher_pizagno,vfinal,imag,color = colors,symbols = psym,symsizes = symsizes,thicks = thicks,ctables = ctables,obscolor = obscolor1,obspsym = obspsym,obssymsize = obssymsize,outplot = outplot,keys = ['Simulations']

;    tully_fisher_obs,    outfiles,vHI,gmass,halo = haloids,symbols = psym,symsizes = symsizes,thicks = thicks,ctables = ctables,obscolor = obscolor1,color = colors,outfile = outplot
    tully_fisher_obs_btf,outfiles,vflat,gmass,halo = haloids,symbols = psym,symsizes = symsizes,thicks = thicks,ctables = ctables,obscolor = obscolor1,obssym = obspsym,obssize = obssymsize,color = colors,outfile = outplot,smass_true = smass_true;,key = key
ENDIF

;-------------------------------------------- Bulge Morphology ------------------
IF 0 THEN BEGIN
;    size_mag_bulge OR size_mag_disk OR mu_size OR mu_mag OR mu_size_mag THEN BEGIN
    formatplot
    FOR i = 0, nfiles - 1 DO BEGIN
        cd,dirs[i]+files[i]+'.'+steps[i]+'/'+files[i]+'.'+steps[i]+'.'+haloids[i]
        checkfit,'broadband.fits',camera,filter,band,key[i]+galfitname,/recenter,parameters = parameters
        fits_galfit[i] = parameters
        fits_galfit[i].M = fits_galfit[i].M - 1.39
        fits_galfit[i].mu = fits_galfit[i].mu - 1.39
        fits_galfit[i].M_D = fits_galfit[i].M_D - 1.39
        fits_galfit[i].M_B = fits_galfit[i].M_B - 1.39        
    ENDFOR
    formatplot,outplot = outplot,thick = formatthick
ENDIF

IF size_mag_bulge OR size_mag_disk OR mu_size OR mu_mag OR mu_size_mag OR mu_size_mag2 THEN BEGIN
    readcol,'/home/christensen/Code/Datafiles/DecompObs/Fisher.dat',name,logn,M_B,mu_B,r_e,r_h,mu_D,M_D,format='(A8,D,D,D,D,D,D,D)'
    spiral  = where(mu_D NE 999,complement = ellip)
    r_e_spiral = r_e[spiral]
    r_e_ellip = r_e[ellip]
    M_B_spiral = M_B[spiral]
    M_B_ellip = M_B[ellip]
    r_h_spiral = r_h[spiral]
    M_D_spiral = M_D[spiral]
    mu_spiral = mu_B[spiral]
    mu_ellip = mu_B[ellip]
    n_spiral = 10^logn[spiral]
    n_ellip = 10^logn[ellip]
    B_D = 10^(-0.4*M_B_spiral)/10^(-0.4*M_D_spiral)
    M_spiral = -2.5*alog10(10^(-0.4*M_B_spiral) + 10^(-0.4*M_D_spiral))
;file277,file258,file285,file239

;h239,h258,h277,h285
;    fits[2:5].r_e = [2.77,2.74,2.85,2.71];][2.71,2.74,2.77,2.85]
;    fits[2:5].M_B = [-22.39,-22.99,-22.90,-22.67]+1.39 ;[-22.67,-22.99,-22.39,-22.90]
;    fits[2:5].mu = [14.25,13.42,14.28,13.94]+1.39;[13.94,13.42,14.25,14.28]
ENDIF
;stop
IF size_mag_bulge THEN BEGIN
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/size_magnitude_spiral.txt',r_e_spiral,M_B_spiral
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/size_magnitude_elliptical.txt',r_e_ellip,M_B_ellip
    IF keyword_set(outplot) THEN device,filename=outplot + '_size_mag.eps', /color,xsize = 15, ysize = 15,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 600,ysize = 600
    plot,r_e_spiral,M_B_spiral,xtitle = textoidl('log(r_e/[pc])'),ytitle = textoidl('M_H(Bulge)'),/nodata,xrange = [1,5],yrange = [-17,-27.5]
    oplot,r_e_ellip, M_B_ellip, psym = symcat(obspsym),color = obscolor2,symsize =  obssymsize2
    oplot,r_e_spiral,M_B_spiral,psym = symcat(obspsym),color = obscolor1,symsize = obssymsize2
    oploterror,[1.5,1.5],[-23,-23],[0.25,0.25],[0.45,0.45],psym = 3
    oplot,fits.r_e,fits.M_B,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits_galfit.r_e,fits_galfit.M_B,psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    IF n_elements(fits) GT 2 THEN oplot,fits[2:nfiles-1].r_e,fits[2:nfiles-1].M_B - 1.39,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;   oplot,fits[0:1       ].r_e,fits[0:1       ].M_B - 1.39,psym = symcat(psym[0]), color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h603fits.r_e,h603fits.r_e],[h603fits.M_B,h603fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h986fits.r_e,h986fits.r_e],[h986fits.M_B,h986fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
    legend,['Observed Bulges in Spirals','Observed Ellipticals','Simulations'],color = [obscolor1,obscolor2,colors[0]],psym=[obspsym,obspsym,psym[0]],thick = [obssymthick,obssymthick,thicks[0]],symsize = [obssymsize,obssymsize,symsizes[0]],ctables = [obsct,obsct,ctables[0]],/top,/left,box = 0,charsize = l_charsize
    IF keyword_set(outplot) THEN device,/close ELSE stop
ENDIF

IF size_mag_disk THEN BEGIN
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/h_mdisk.txt',r_h_spiral,M_D_spiral
    IF keyword_set(outplot) THEN device,filename=outplot + '_h_mag.eps', /color,xsize = 15, ysize = 15,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 600,ysize = 600
    plot,r_h_spiral,M_D_spiral,xtitle = textoidl('log(h/[pc])'),ytitle = textoidl('M_H(Disk)'),/nodata,xrange = [2.2,4.5],yrange = [-17.5,-26]
    oplot,r_h_spiral,M_D_spiral,psym = symcat(obspsym),color = obscolor1,symsize = obssymsize2
    oploterror,[2.4,2.4],[-24,-24],[0.05,0.05],[0.25,0.25]
    oplot,fits.r_h,fits.M_D,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oploterror,[2.3,2.3],[],[],[]
;    oplot,fits_galfit.r_h,fits_galfit.M_D,psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    IF n_elements(fits) GT 2 THEN oplot,fits[2:nfiles-1].r_h,fits[2:nfiles-1].M_D - 1.39,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits[0:1       ].r_h,fits[0:1       ].M_D - 1.39,psym = symcat(psym[0]) ,color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h603fits.r_h,h603fits.r_h],[h603fits.M_D,h603fits.M_D],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h986fits.r_h,h986fits.r_h],[h986fits.M_D,h986fits.M_D],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
    oploterror,[h986fits.r_h,h986fits.r_h],[h986fits.M_D,h986fits.M_D],[h986fits.r_h - alog10(10^h986fits.r_h/2),h986fits.r_h - alog10(10^(h986fits.r_h)/2)],[0,0],errcolor = colors[0],/lobar;,thick = thicks[0]
    oploterror,[h986fits.r_h,h986fits.r_h],[h986fits.M_D,h986fits.M_D],[alog10(2*10^h986fits.r_h) - h986fits.r_h,alog10(2*10^(h986fits.r_h)) - h986fits.r_h],[0,0],errcolor = colors[0],/highbar;,thick = thicks[0]
    legend,['Observed Bulges in Spirals','Simulations'],color = [obscolor1,colors[0]],psym=[obspsym,psym[0]],thick = [obssymthick,thicks[0]],symsize = [obssymsize,symsizes[0]],ctables = [obsct,ctables[0]],/top,/left,box = 0,charsize = l_charsize
    IF keyword_set(outplot) THEN device,/close ELSE stop
ENDIF
tHisletter = "155B
mu = '!9' + String(thisLetter) + '!X'

IF mu_size THEN BEGIN
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/size_mu_spiral.txt',r_e_spiral,mu_spiral
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/size_mu_elliptical.txt',r_e_ellip,mu_ellip
    IF keyword_set(outplot) THEN device,filename=outplot + '_mu_size.eps', /color,xsize = 15, ysize = 15 ELSE window,0,xsize = 600,ysize = 600
    plot,mu_spiral,r_e_spiral,xtitle = cgGreek('mu') + textoidl(' (r_e) [H mag arcsec^{-2}]'),ytitle = textoidl('log(r_e/[pc])'),/nodata,xrange = [24,10],yrange = [1.7,5.2],symsize = 2
    oplot,mu_spiral,r_e_spiral,psym = symcat(obspsym),color = obscolor1, symsize = obssymsize2
    oplot,mu_ellip, r_e_ellip, psym = symcat(obspsym),color = obscolor2,symsize = obssymsize2
    oploterror,[22,22],[2.3,2.3],[0.4,0.4],[0.25,0.25]
    oplot,fits.mu,fits.r_e,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits_galfit.mu,fits_galfit.r_e,psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    IF n_elements(fits) GT 2 THEN oplot,fits[2:nfiles-1].mu - 1.39,fits[2:nfiles-1].r_e,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits[0:1       ].mu - 1.39,fits[0:1       ].r_e,psym = symcat(psym[0]), color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h603fits.mu,h603fits.mu],[h603fits.r_e,h603fits.r_e],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h986fits.mu,h986fits.mu],[h986fits.r_e,h986fits.r_e],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
    legend,['Observed Bulges in Spirals','Observed Ellipticals','Simulations'],color = [obscolor1,obscolor2,colors[0]],psym=[obspsym,obspsym,psym[0]],thick = [obssymthick,obssymthick,thicks[0]],symsize = [obssymsize,obssymsize,symsizes[0]],ctables = [obsct,obsct,ctables[0]],/top,/right,box = 0,charsize = l_charsize
    IF keyword_set(outplot) THEN device,/close ELSE stop
ENDIF

IF mu_mag THEN BEGIN
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/mu_magnitude_spiral.txt',mu_spiral,M_B_spiral
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/mu_magnitude_elliptical.txt',mu_ellip,M_B_ellip
    IF keyword_set(outplot) THEN device,filename=outplot + '_mu_mag.eps', /color,xsize = 15, ysize = 15,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 600,ysize = 600
    plot,mu_spiral,M_B_spiral,xtitle = cgGreek('mu') + textoidl(' (r_e) [H mag arcsec^{-2}]'),ytitle = textoidl('M_H(Bulge)'),/nodata,xrange = [24,10],yrange = [-17,-27.5] 
    oplot,mu_spiral,M_B_spiral,psym = symcat(obspsym),color = obscolor1, symsize = obssymsize2
    oplot,mu_ellip, M_B_ellip ,psym = symcat(obspsym),color = obscolor2,symsize = obssymsize2
    oploterror,[22,22],[-19,-19],[0.4,0.4],[0.45,0.45]
    oplot,fits.mu,fits.M_B,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits_galfit.mu,fits_galfit.M_B,psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    IF n_elements(fits) GT 2 THEN oplot,fits[2:nfiles-1].mu - 1.39,fits[2:nfiles-1].M_B - 1.39,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits[0:1       ].mu - 1.39,fits[0:1       ].M_B - 1.39,psym = symcat(psym[0]), color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h603fits.mu,h603fits.mu],[h603fits.M_B,h603fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h986fits.mu,h986fits.mu],[h986fits.M_B,h986fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
    legend,['Observed Bulges in Spirals','Observed Ellipticals','Simulations'],color = [obscolor1,obscolor2,colors[0]],psym=[obspsym,obspsym,psym[0]],thick = [obssymthick,obssymthick,thicks[0]],symsize = [obssymsize,obssymsize,symsizes[0]],ctables = [obsct,obsct,ctables[0]],/top,/right,box = 0,charsize = l_charsize
    IF keyword_set(outplot) THEN device,/close ELSE stop
ENDIF

IF mu_size_mag THEN BEGIN
    x_mu_size = findgen(100)/5+10
    y_mu_size = 0.29674*(x_mu_size - 11.7) + 1.4

;    readcol,'/home/christensen/Code/Datafiles/DecompObs/size_mu_spiral.txt',r_e_spiral,mu_spiral
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/size_mu_elliptical.txt',r_e_ellip,mu_ellip
    IF keyword_set(outplot) THEN device,filename=outplot + '_mu_mag.eps', /color,xsize = 12, ysize = 20,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 500,ysize = 920
    multiplot,[1,2],/square
    plot,mu_spiral,r_e_spiral,ytitle = textoidl('log(r_e/[pc])'),/nodata,xrange = [24,10],yrange = [1.7,5.2],symsize = 2;,xtitle = textoidl('\mu (r_e) [H mag arcsec^{-2}]')
    oplot,mu_spiral,r_e_spiral,psym = symcat(obspsym),color = obscolor1, symsize = obssymsize2
    oplot,mu_ellip, r_e_ellip, psym = symcat(obspsym),color = obscolor2,symsize = obssymsize2
    oploterror,[22,22],[2.3,2.3],[0.4,0.4],[0.25,0.25]
    oplot,fits.mu,fits.r_e,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits_galfit.mu,fits_galfit.r_e,psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    IF n_elements(fits) GT 2 THEN oplot,fits[2:nfiles-1].mu - 1.39,fits[2:nfiles-1].r_e,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits[0:1       ].mu - 1.39,fits[0:1       ].r_e,psym = symcat(psym[0]), color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h603fits.mu,h603fits.mu],[h603fits.r_e,h603fits.r_e],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h986fits.mu,h986fits.mu],[h986fits.r_e,h986fits.r_e],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
    oplot,x_mu_size,y_mu_size,linestyle = 2
    legend,['Observed Bulges in Spirals','Observed Ellipticals','Simulations'],color = [obscolor1,obscolor2,colors[0]],psym=[obspsym,obspsym,psym[0]],thick = [obssymthick,obssymthick,thicks[0]],symsize = [obssymsize,obssymsize,symsizes[0]],ctables = [obsct,obsct,ctables[0]],/top,/right,box = 0,charsize = l_charsize

    x_mu_mag = x_mu_size
    y_mu_mag = -0.592857*(x_mu_mag - 10) - 18.8
    multiplot
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/mu_magnitude_spiral.txt',mu_spiral,M_B_spiral
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/mu_magnitude_elliptical.txt',mu_ellip,M_B_ellip
    plot,mu_spiral,M_B_spiral,xtitle = cgGreek('mu') + textoidl(' (r_e) [H mag arcsec^{-2}]'),ytitle = textoidl('M_H(Bulge)'),/nodata,xrange = [24,10],yrange = [-17,-27.5] ;,xtitle = textoidl('\mu (r_e) [H mag arcsec^{-2}]')
    oplot,mu_spiral,M_B_spiral,psym = symcat(obspsym),color = obscolor1, symsize = obssymsize2
    oplot,mu_ellip, M_B_ellip ,psym = symcat(obspsym),color = obscolor2,symsize = obssymsize2
    oploterror,[22,22],[-19,-19],[0.4,0.4],[0.45,0.45]
    oplot,fits.mu,fits.M_B,psym = symcat (psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits_galfit.mu,fits_galfit.M_B,psym = symcat (1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    IF n_elements(fits) GT 2 THEN oplot,fits[2:nfiles-1].mu - 1.39,fits[2:nfiles-1].M_B - 1.39,psym = symcat (psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits[0:1       ].mu - 1.39,fits[0:1       ].M_B - 1.39,psym = symcat(psym[0]), color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h603fits.mu,h603fits.mu],[h603fits.M_B,h603fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h986fits.mu,h986fits.mu],[h986fits.M_B,h986fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
    oplot,x_mu_mag,y_mu_mag,linestyle = 2
    multiplot,/reset
    IF keyword_set(outplot) THEN device,/close ELSE stop
ENDIF


IF mu_size_mag2 THEN BEGIN
    x_mu_size = findgen(100)/5+10
    y_mu_size = 0.29674*(x_mu_size - 11.7) + 1.4

;    readcol,'/home/christensen/Code/Datafiles/DecompObs/size_mu_spiral.txt',r_e_spiral,mu_spiral
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/size_mu_elliptical.txt',r_e_ellip,mu_ellip
    IF keyword_set(outplot) THEN device,filename=outplot + '_mu_mag.eps', /color,xsize = 20.5, ysize = 20,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 920,ysize = 920
    multiplot,[2,2],/square
    plot,mu_spiral,r_e_spiral,ytitle = textoidl('log(r_e/[pc])'),/nodata,xrange = [24,10],yrange = [1.7,5.2],symsize = 2;,xtitle = textoidl('\mu (r_e) [H mag arcsec^{-2}]')
    oplot,mu_ellip, r_e_ellip, psym = symcat(obspsym),color = obscolor2,symsize = obssymsize2
    oplot,mu_spiral,r_e_spiral,psym = symcat(obspsym),color = obscolor1, symsize = obssymsize2
    oploterror,[22,22],[2.3,2.3],[0.4,0.4],[0.25,0.25]
    oplot,fits.mu,fits.r_e,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits_galfit.mu,fits_galfit.r_e,psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    IF n_elements(fits) GT 2 THEN oplot,fits[2:nfiles-1].mu - 1.39,fits[2:nfiles-1].r_e,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits[0:1       ].mu - 1.39,fits[0:1       ].r_e,psym = symcat(psym[0]), color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h603fits.mu,h603fits.mu],[h603fits.r_e,h603fits.r_e],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h986fits.mu,h986fits.mu],[h986fits.r_e,h986fits.r_e],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
    oplot,x_mu_size,y_mu_size,linestyle = 2
    legend,['Observed Bulges in Spirals','Observed Ellipticals','Simulations'],color = [obscolor1,obscolor2,colors[0]],psym=[obspsym,obspsym,psym[0]],thick = [obssymthick,obssymthick,thicks[0]],symsize = [obssymsize,obssymsize,symsizes[0]],ctables = [obsct,obsct,ctables[0]],/top,/right,box = 0,charsize = l_charsize

    x_mu_mag = x_mu_size
    y_mu_mag = -0.592857*(x_mu_mag - 10) - 18.8
    multiplot
    multiplot
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/mu_magnitude_spiral.txt',mu_spiral,M_B_spiral
;    readcol,'/home/christensen/Code/Datafiles/DecompObs/mu_magnitude_elliptical.txt',mu_ellip,M_B_ellip
    plot,mu_spiral,M_B_spiral,xtitle = cgGreek('mu') + textoidl(' (r_e) [H mag arcsec^{-2}]'),ytitle = textoidl('M_H(Bulge)'),/nodata,xrange = [24,10],yrange = [-17,-27.5] ;,xtitle = textoidl('\mu (r_e) [H mag arcsec^{-2}]')
    oplot,mu_ellip, M_B_ellip ,psym = symcat(obspsym),color = obscolor2,symsize = obssymsize2
    oplot,mu_spiral,M_B_spiral,psym = symcat(obspsym),color = obscolor1, symsize = obssymsize2
    oploterror,[22,22],[-19,-19],[0.4,0.4],[0.45,0.45]
    oplot,fits.mu,fits.M_B,psym = symcat (psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits_galfit.mu,fits_galfit.M_B,psym = symcat (1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    IF n_elements(fits) GT 2 THEN oplot,fits[2:nfiles-1].mu - 1.39,fits[2:nfiles-1].M_B - 1.39,psym = symcat (psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits[0:1       ].mu - 1.39,fits[0:1       ].M_B - 1.39,psym = symcat(psym[0]), color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h603fits.mu,h603fits.mu],[h603fits.M_B,h603fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h986fits.mu,h986fits.mu],[h986fits.M_B,h986fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
    oplot,x_mu_mag,y_mu_mag,linestyle = 2

    multiplot
;    IF keyword_set(outplot) THEN device,filename=outplot + '_size_mag.eps', /color,xsize = 15, ysize = 15,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 600,ysize = 600
    plot,r_e_spiral,M_B_spiral,xtitle = textoidl('log(r_e/[pc])'),/nodata,xrange = [1,5],yrange = [-17,-27.5];,ytitle = textoidl('M_H(Bulge)'),/nodata,xrange = [1,5],yrange = [-17,-27.5]
    oplot,r_e_ellip, M_B_ellip, psym = symcat(obspsym),color = obscolor2,symsize =  obssymsize2
    oplot,r_e_spiral,M_B_spiral,psym = symcat(obspsym),color = obscolor1,symsize = obssymsize2
    oploterror,[1.5,1.5],[-23,-23],[0.25,0.25],[0.45,0.45],psym = 3
    oplot,fits.r_e,fits.M_B,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,fits_galfit.r_e,fits_galfit.M_B,psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    IF n_elements(fits) GT 2 THEN oplot,fits[2:nfiles-1].r_e,fits[2:nfiles-1].M_B - 1.39,psym = symcat(psym[0]),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;   oplot,fits[0:1       ].r_e,fits[0:1       ].M_B - 1.39,psym = symcat(psym[0]), color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h603fits.r_e,h603fits.r_e],[h603fits.M_B,h603fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    oplot,[h986fits.r_e,h986fits.r_e],[h986fits.M_B,h986fits.M_B],psym = symcat(1),color = colors[0],symsize = symsizes[0],thick = thicks[0]
;    legend,['Observed Bulges in Spirals','Observed Ellipticals','Simulations'],color = [obscolor1,obscolor2,colors[0]],psym=[obspsym,obspsym,psym[0]],thick = [obssymthick,obssymthick,thicks[0]],symsize = [obssymsize,obssymsize,symsizes[0]],ctables = [obsct,obsct,ctables[0]],/top,/left,box = 0,charsize = l_charsize

    multiplot,/reset
    IF keyword_set(outplot) THEN device,/close ELSE stop
ENDIF


IF plot_bulgeSFH THEN BEGIN
    bulgeSFH,dirs+files+'.'+steps+'/'+files+'.'+steps + '.halo.' + haloids,massunits,lengthunits,timeunits,key,outplot = outplot
ENDIF

END
