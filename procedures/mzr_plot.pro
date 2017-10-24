
;filenames=['/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir/h516.cosmo25cmb.3072g14HBWK.00492','/home/christensen/Storage2/UW/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512','/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512','/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512']
;psym=[4,5,6,9] ;11,12,13x
;outfile = '~/mz_all.eps'

;filenames=['/astro/store/nbody2/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir/h516.cosmo25cmb.3072g14HBWK.00492','/astro/store/nbody2/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512','/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512','/astro/store/nbody2/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512','/astro/store/nbody2/abrooks/h277.cosmo50cmb.3072g14HMbwK/h277.cosmo50cmb.3072g14HMbwK.00512/h277.cosmo50cmb.3072g14HMbwK.00512','/astro/store/nbody3/abrooks/h285.cosmo50cmb.3072g14HMbwK/h285.cosmo50cmb.3072g14HMbwK.00512/h285.cosmo50cmb.3072g14HMbwK.00512']
;psym=[4,5,6,9,11,12,13]
;psym=[6,6,6,6,6, 6, 6]
;outfile = '~/plots/all'
;mzr_plot,filename,psym = psym, /readfile, outfile = outfile

PRO mzr_plot,filenames, halos = halos, outfile = outfile, key = key, psym  = psym, color = color, symsizes = symsizes, ctables = ctables, thicks = thicks, simStellarMass = simStellarMass, noObsStellarMass = noObsStellarMass, xrange = xrange, yrange = yrange, readfile = readfile, obscolor = obscolor, onehalos = onehalos, obsct = obsct, formatthick = formatthick,redshift = redshift,stellar = stellar
n = n_elements(filenames)
zsolar = 0.0130215

IF keyword_set(color) THEN BEGIN
    white = 13
    black = 0
    IF NOT keyword_set(obsct) THEN obsct = 0;39
    loadct,obsct
    distinct_colors,n_colors = 12
ENDIF ELSE BEGIN
    black = 0
    white = 255
    IF NOT keyword_set(obsct) THEN obsct = 0 
    loadct,obsct
;    IF NOT keyword_set(ctables) THEN ctables = fltarr(n)
ENDELSE
IF keyword_set(outfile) THEN BEGIN
    l_charsize = 1.2 ;1 ;0.75
    fgcolor = black
    bgcolor = white
    xsize = 18;7
    ysize = 12;5
    formatplot,/outplot,thick = formatthick
ENDIF ELSE BEGIN
    l_charsize = !p.charsize ;1.0
    fgcolor = white
    bgcolor = black
    xsize = 600
    ysize = 400
    formatplot,thick = formatthick
ENDELSE
IF keyword_set(color) THEN BEGIN
;    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF color[0] EQ 1 THEN  colors = (fltarr(n) + 1)*fgcolor else colors = fltarr(n) + color
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(psym) THEN psym = fltarr(n) + 4 ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
    IF NOT keyword_set(obscolor) THEN obscolor = [2,6,11];fgcolor
    obssymsize = 1.25;2;.5
    obssym = 16
    obsthick = thicks[0];2
ENDIF ELSE BEGIN
    colors = (fltarr(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(psym) THEN psym = findgen(n) + 4;(findgen(n) + 2)*2   
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
    IF NOT keyword_set(obscolor) THEN obscolor = [150,150,150]
    obssymsize = 1.25;2;.5
    obssym = 16
    obsthick = thicks[0];2
ENDELSE
IF NOT keyword_set(xrange) THEN xrange = [5.5,11.5]
IF NOT keyword_set(yrange) THEN $
  IF keyword_set(stellar) THEN yrange = [-2.8,0.5] ELSE yrange = [7,9.1]

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN datafile_base = '/home/christensen/Code/Datafiles/' ELSE $
  IF strcmp(hostname,'pfe',3) THEN datafile_base = '/home6/crchrist/Datafiles/' ELSE datafile_base = '/astro/users/christensen/code/Datafiles/'
readcol,datafile_base + 'Lee06.txt',gal,source,aorkey,f,ferr,mag,magerr,disref,magA,magAerr,logOH,logOHerr,OPHref,B,logMstar,format = '(A,A,I,F,F,F,F,I,F,F,F,F,F,F,F,F)'
readcol,datafile_base + 'Tremonti04.txt',LogMstarT,P_25,P_16,P_50,P_84,P_975

IF (keyword_set(outfile)) THEN BEGIN
   IF keyword_set(stellar) THEN device,filename=outfile+'_mz_stellar.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2,/encapsul $
     ELSE device,filename=outfile+'_' + strtrim(redshift,2)+'_mz.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2,/encapsul
ENDIF ELSE BEGIN
   window,0,xsize = xsize,ysize = ysize ;392
ENDELSE

LogMstar_arr = findgen(101)/100.0*(8.5 - 5.5) + 5.5

uniqind0 = (uniq(psym,sort(psym)))[0]
uniqind1 = (uniq(psym,sort(psym)))[1]

IF NOT keyword_set(redshift) THEN redshift = 0
IF keyword_set(stellar) $
THEN plot,logMstar,logOH,xtitle = textoidl('log(M_{*} [M') + sunsymbol() + '])', ytitle = textoidl('log(Z_*/Z')+sunsymbol()+')',/nodata,xrange = xrange,yrange = yrange $
ELSE plot,logMstar,logOH,xtitle = textoidl('log(M_{*} [M') + sunsymbol() + '])', ytitle = '12 + log(O/H)',/nodata,xrange = xrange,yrange = yrange
IF redshift EQ 0 THEN BEGIN
    IF NOT keyword_set(stellar) THEN BEGIN
;Adjust Lee data
        oploterror,logMstar,logOH,logMstar*0,logOHerr,psym = symcat(obssym), color = obscolor[0], symsize = obssymsize, thick = obsthick,errcolor=obscolor[0]
       a = 319.7570 ;-0.26 ;Used in Christensen 2016
        b = -107.13160 ; 1
        c = 12.208670 ; 0
        d = -0.4606539 ; 0
;        correction = a + b*LogMstarT + c*LogMstarT^2 + d*LogMstarT^3
        oplot,LogMstarT,a + b*P_25 + c*P_25^2 + d*P_25^3,color = obscolor[1], thick = obsthick, linestyle = 1
        oplot,LogMstarT,a + b*P_16 + c*P_16^2 + d*P_16^3,color = obscolor[1], thick = obsthick, linestyle = 2
        oplot,LogMstarT,a + b*P_50 + c*P_50^2 + d*P_50^3,color = obscolor[1], thick = obsthick, linestyle = 0
        oplot,LogMstarT,a + b*P_84 + c*P_84^2 + d*P_84^3,color = obscolor[1], thick = obsthick, linestyle = 2
        oplot,LogMstarT,a + b*P_975 + c*P_975^2 + d*P_975^3,color = obscolor[1], thick = obsthick, linestyle = 1
;oplot,LogMstar_arr, -1.492 + 1.847*LogMstar_arr - 0.08026*LogMstar_arr*LogMstar_arr - 0.26,color = obscolor, thick = obsthick, linestyle = 3

        LogMstarAM = findgen(32)/10+7.4
        mzrAM = 8.798 - alog10(1 + (10^8.901/10^LogMstarAM)^0.64) ;fit from Andrew & Martini 2013, see table 4
        oplot,LogMstarAM,mzrAM,color = obscolor[2], thick = obsthick, linestyle = 3
        IF keyword_set(key) THEN legend,[key,"Lee et al., 2006","Tremonti et al., 2004",'Andrew & Martini, 2013'],color = [color,obscolor[0],obscolor[1],obscolor[2]],linestyle = intarr(n_elements(key) + 3),psym = [psym,obssym,0,0],thick = [thicks,obsthick,obsthick,obsthick] $
        ELSE legend,['Lee et al., 2006','Tremonti et al., 2004','Andrew & Martini, 2013','Med-res Sims','High-res Sims'],color = [obscolor[0],obscolor[1],obscolor[2],colors[uniqind0],colors[uniqind1]],linestyle = [0,0,3,2,0],psym = [obssym,0,0,psym[uniqind0],psym[uniqind1]],thick = [obsthick,obsthick,obsthick,thicks[uniqind0],thicks[uniqind1]],box = 0,symsize = [obssymsize,obssymsize,obssymsize,symsizes[uniqind0],symsizes[uniqind1]],charsize = l_charsize,/bottom,/right,position = [0.934145,0.250005],/norm;,ctables = [obsct,obsct,obsct,ctables[uniqind0],ctables[uniqind1]] ;corners = corners;0.934145     0.503365     0.460109     0.271045
        legend,['z = 0'],box = 0,/left,/top
    ENDIF ELSE BEGIN
        ;Gallazzi et al, 2005
        readcol,datafile_base + 'Gallazzi2005.txt',LogMstarG,LogZ_50,LogZ_16,LogZ_84,Logt_50,Logt_16,Logt_84
        oplot,LogMstarG,LogZ_50,color = obscolor[0], thick = obsthick
        oplot,LogMstarG,LogZ_16,color = obscolor[0], thick = obsthick, linestyle = 1
        oplot,LogMstarG,LogZ_84,color = obscolor[0], thick = obsthick, linestyle = 1
        ;Kirby et al., 2013
        readcol,datafile_base + 'Kirby2013.txt',gal,nstars,LogL_V,LogL_V_err,LogM_star,LogM_star_err,LogFe,LogFe_err,sigma,sigma2,median,mad,iqr,skewness,skewness_err,Kurtosis,Kurtosis_err,format='(A,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)'
        oplot,LogM_star,LogFe,psym = symcat(obssym),color = obscolor[1], thick = obsthick
        oploterror,logM_star,LogFe,LogM_star_err,LogFe_err,psym = symcat(obssym), color = obscolor[1], symsize = obssymsize, thick = obsthick,errcolor=obscolor
        IF keyword_set(key) THEN legend,[key,"Kirby et al., 2013","Gallazzi et al., 2005"],color = [color,obscolor[1],obscolor[0]],linestyle = intarr(n_elements(key) + 2),psym = [psym,obssym,0],thick = [thicks,obsthick,obsthick] $
        ELSE legend,["Kirby et al., 2013","Gallazzi et al., 2005",'Med-res Sims','High-res Sims'],color = [obscolor[1],obscolor[0],colors[uniqind0],colors[uniqind1]],linestyle = [0,0,0,0],psym = [obssym,0,psym[uniqind0],psym[uniqind1]],thick = [obsthick,obsthick,thicks[uniqind0],thicks[uniqind1]],box = 0,symsize = [obssymsize,obssymsize,symsizes[uniqind0],symsizes[uniqind1]],charsize = l_charsize,/bottom,/right,position = [0.934145,0.250005],/norm;,,ctables = [obsct,obsct,ctables[uniqind0],ctables[uniqind1]] ;corners = corners;0.934145     0.503365     0.460109     0.271045

    ENDELSE
ENDIF

IF redshift EQ 0.8 AND NOT keyword_set(stellar) THEN BEGIN
;Also consider: z = 0.6\u20130.8 (Lamareille et al. 2009, labeled L09). The dot-dashed curve is for galaxies with z = 0.4\u20130.98 (Savaglio et al. 2005, labeled S05). The triple-dot-dashed curve is for galaxies with z = 0.475\u20130.9 (Cowie & Barger 2008, labeled C08). The diamonds are the median instead of the mean of the data from Cowie & Barger (2008). 
    a = 916.7484
    b = -309.54480
    c = 35.051680
    d = -1.3188
    
    LogMstarZ = findgen(16)/10+9.1
    mzrZ = 8.923 + 0.24*(LogMstarZ-10) - 0.06*(LogMstarZ-10)^2 ;fit from Zahid 2011, see eq 8, with adjustment
    oplot,LogMstarZ,a + b*mzrZ + c*mzrZ^2 + d*mzrZ^3,color = obscolor[0], thick = obsthick, linestyle = 3
    IF keyword_set(key) THEN legend,[key,"Zahid et al., 2011"],color = [color,obscolor[0]],linestyle = intarr(n_elements(key) + 1),psym = [psym,0,0],thick = [thicks,obsthick] $
    ELSE legend,['Zahid et al., 2011','Med-res Sims','High-res Sims'],color = [obscolor[0],colors[uniqind0],colors[uniqind1]],linestyle = [3,2,0],psym = [0,psym[uniqind0],psym[uniqind1]],thick = [obsthick,thicks[uniqind0],thicks[uniqind1]],box = 0,symsize = [obssymsize,symsizes[uniqind0],symsizes[uniqind1]],charsize = l_charsize,/bottom,/right,position = [0.934145,0.250005],/norm;,,ctables = [obsct,ctables[uniqind0],ctables[uniqind1]]
    legend,['z = 0.8'],box = 0,/left,/top
ENDIF

IF (redshift EQ 2.2 OR redshift EQ 2.3) AND NOT keyword_set(stellar) THEN BEGIN
;Sanders+ 2015, Table 1, z = 2.3
    LogMstarSand = [9.45,9.84,10.11,10.56]
    LogMstarSand_low = LogMstarSand - [9.15,9.68,9.99,10.29]
    LogMstarSand_high= [9.68,9.94,10.27,11.11] - LogMstarSand
    Z_N2_Sand = [8.18,8.30,8.44,8.52]
    Z_N2_Sand_low = [0.1,0.05,0.04,0.02]
    Z_N2_Sand_high =[0.07,0.04,0.04,0.02]
    Z_O3N2_Sand = [8.11,8.20,8.31,8.42]
    Z_O3N2_Sand_low = [0.06,0.03,0.03,0.02] 
    Z_O3N2_Sand_high = [0.04,0.02,0.02,0.02]
    oplot,LogMstarSand,Z_N2_Sand,psym = symcat(17), color = obscolor[0], symsize = obssymsize, thick = obsthick
    oploterror,LogMstarSand,Z_N2_Sand,LogMstarSand_low,Z_N2_Sand_low,psym = symcat(obssym), color = obscolor[0], symsize = obssymsize, thick = obsthick,errcolor=obscolor[0],/lobar,errthick = 1;,errstyle = 2
    oploterror,LogMstarSand,Z_N2_Sand,LogMstarSand_high,Z_N2_Sand_high,psym = symcat(obssym), color = obscolor[0], symsize = obssymsize, thick = obsthick,errcolor=obscolor[0],/hibar,errthick = 1;,errstyle = 2
;    oploterror,LogMstarSand,Z_O3N2_Sand,LogMstarSand_low,Z_N2_Sand_low,psym = symcat(obssym), color = obscolor, symsize = obssymsize, thick = obsthick,errcolor=obscolor,/lobar,errstyle = 1
;    oploterror,LogMstarSand,Z_O3N2_Sand,LogMstarSand_high,Z_N2_Sand_high,psym = symcat(obssym), color = obscolor, symsize = obssymsize, thick = obsthick,errcolor=obscolor,/hibar,errstyle = 1

;Steidel+ 2015, 
;" We have also shown (see also Newman et al. 2013) that there is a
;systematic offset of \u0394 sime 0.13\u2009dex between metallicities
;inferred from the PP04 N2 and O3N2 indices when they are applied to
;the same galaxies in the KBSS-MOSFIRE z sime 2.3 sample, in spite of
;the fact that these calibrations were established using direct Te
;abundances of the same local galaxy sample."
;Includes a conversion between O3N2 and N2, Eq 18
    LogMstarS = findgen(30)/10+8.4
    Z_N2_S = 8.41 + 0.2*(logMstarS - 10)
    Z_O3N2_S = 8.27 + 0.19*(logMstarS - 10)
    oplot,LogMstarS,Z_N2_S,color = obscolor[1], thick = obsthick, linestyle = 3
;    oplot,LogMstarS,Z_O3N2_S,color = obscolor[1], thick = obsthick, linestyle = 3

;Erb+ 2006
    LogMstarE = alog10([0.26,0.71,1.5,2.6,4.1,10.5]*1e10)
    LogMstarE_err = [0.15,0.17,0.3,0.4,0.6,5.4]
    Z_N2_E = [8.2,8.33,8.42,8.46,8.52,8.58]
    Z_N2_E_low = [0.07,0.05,0.05,0.05,0.04]
    Z_N2_E_high = [0.07,0.06,0.06,0.06,0.06]
    oplot,LogMstarE,Z_N2_E,psym = symcat(16), color = obscolor[2], symsize = obssymsize, thick = obsthick
    oploterror,LogMstarE,Z_N2_E,LogMstarE_err,Z_N2_E_low,psym = symcat(obssym), color = obscolor[2], symsize = obssymsize, thick = obsthick,errcolor=obscolor[2],/lobar,errthick = 1
    oploterror,LogMstarE,Z_N2_E,LogMstarE_err,Z_N2_E_high,psym = symcat(obssym), color = obscolor[2], symsize = obssymsize, thick = obsthick,errcolor=obscolor[2],/hibar,errthick = 1 

    IF keyword_set(key) THEN legend,[key,"Sanders et al., 2015",'Erb et al., 2006',"Steidel et al., 2015"],color = [color,obscolor[0],obscolor[2],obscolor[1]],linestyle = intarr(n_elements(key) + 3),psym = [psym,14,16,0],thick = [thicks,obsthick,obsthick,obsthick] $
    ELSE legend,["Sanders et al., 2015",'Erb et al., 2006',"Steidel et al., 2015",'Med-res Sims','High-res Sims'],color = [obscolor[0],obscolor[2],obscolor[1],colors[uniqind0],colors[uniqind1]],linestyle = [0,0,3,2,0],psym = [17,16,0,psym[uniqind0],psym[uniqind1]],thick = [obsthick,obsthick,obsthick,thicks[uniqind0],thicks[uniqind1]],box = 0,symsize = [obssymsize,obssymsize,obssymsize,symsizes[uniqind0],symsizes[uniqind1]],charsize = l_charsize,/bottom,/right,position = [0.934145,0.250005],/norm;,,ctables = [obsct,obsct,obsct,ctables[uniqind0],ctables[uniqind1]]
    legend,['z = 2.2'],box = 0,/left,/top
ENDIF

IF redshift EQ 3 AND NOT keyword_set(stellar) THEN BEGIN
;Mannucci+ 2009 (Probably over larger redshift range)
    LogMstarM = findgen(20)/10+9.0
    Z_M = -0.0864*(LogMstarM - 12.28)^2 + 8.69
    oplot,LogMstarM,Z_M,color = obscolor[0], thick = obsthick, linestyle = 3
    IF keyword_set(key) THEN legend,[key,"Mannucci et al., 2009"],color = [color,obscolor[0]],linestyle = intarr(n_elements(key) + 1),psym = [psym,0,0],thick = [thicks,obsthick] $
    ELSE legend,['Mannucci et al., 2009','Med-res Sims','High-res Sims'],color = [obscolor[0],colors[uniqind0],colors[uniqind1]],linestyle = [3,2,0],psym = [0,psym[uniqind0],psym[uniqind1]],thick = [obsthick,thicks[uniqind0],thicks[uniqind1]],box = 0,symsize = [obssymsize,symsizes[uniqind0],symsizes[uniqind1]],charsize = l_charsize,/bottom,/right,position = [0.934145,0.250005],/norm;,,ctables = [obsct,ctables[uniqind0],ctables[uniqind1]]
    legend,['z = 3'],box = 0,/left,/top
ENDIF

FOR i = 0, n - 1 DO BEGIN
;    loadct,ctables[i]
    IF keyword_set(onehalos) THEN metalfile = filenames[i] + '.halo.' + strtrim(onehalos[i],2) ELSE  metalfile = filenames[i]
    IF keyword_set(readfile) AND file_test(metalfile + '.metals.fits') THEN BEGIN 
        metals = mrdfits(metalfile + '.metals.fits',1,/silent) 
        zstellar = alog10((2.09*metals.ox_stellar + 1.06*metals.fe_stellar)/zsolar)
        zstellar_rband = alog10((2.09*metals.ox_stellar_rband + 1.06*metals.fe_stellar_rband)/zsolar)
    ENDIF ELSE BEGIN
        IF keyword_set(onehalos) THEN IF keyword_set(noObsStellarMass) THEN metals = mzr(filenames[i], onehalo = onehalos[i]) ELSE metals = mzr(filenames[i],/obs, onehalo = onehalos[i]) ELSE $
          IF keyword_set(noObsStellarMass) THEN metals = mzr(filenames[i]) ELSE metals = mzr(filenames[i],/obs) 
        mwrfits,metals,metalfile + '.metals.fits',/create
        zstellar = alog10((2.09*metals.ox_stellar + 1.06*metals.fe_stellar)/zsolar)
        zstellar_rband = alog10((2.09*metals.ox_stellar_rband + 1.06*metals.fe_stellar_rband)/zsolar)
    ENDELSE
    IF keyword_set(halos) THEN BEGIN
       FOR j = 0, n_elements(halos[*,i]) - 1 DO BEGIN
           IF halos[j,i] NE 0 THEN BEGIN
               ind = where(metals.grp EQ halos[j,i])
               IF ind[0] NE -1 THEN BEGIN
                   metallicity = metals[ind].ox_sfr
                   IF NOT finite(metals[ind].ox_sfr) THEN metallicity = metals[ind].ox_cold
                   IF (metals[ind].ncold LE 10) THEN metallicity = 0
                   IF NOT keyword_set(noObsStellarMass) THEN BEGIN
;                   IF metals[ind].sat EQ 0 THEN 
                       IF keyword_set(stellar) $
                       THEN oplot,[alog10(metals[ind].smassO),alog10(metals[ind].smassO)],[zstellar_rband[ind],zstellar_rband[ind]],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] $
                       ELSE oplot,[alog10(metals[ind].smassO),alog10(metals[ind].smassO)],[metallicity,metallicity],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] ;ELSE oplot,[alog10(metals[ind].smassO),alog10(metals[ind].smassO)],[metals[ind].ox_sfr,metals[ind].ox_sfr],psym = symcat(psym[i]) + 1, color = colors[i], symsize = symsizes[i], thick = thicks[i]
                   ENDIF
;                   IF keyword_set(simStellarMass) THEN oplot,[alog10(metals[ind].smass),alog10(metals[ind].smass)],[metals[ind].ox_sfr,metals[ind].ox_sfr],psym = symcat(psym[i]), color = colors[i], symsize = 0.5*symsizes[i], thick = thicks[i]              
                   IF keyword_set(simStellarMass) THEN $
                     IF keyword_set(stellar) $ 
                     THEN oplot,[alog10(metals[ind].smass),alog10(metals[ind].smass)],[zstellar[ind],zstellar[ind]],psym = symcat(psym[i]), color = colors[i], symsize = 0.5*symsizes[i], thick = thicks[i] $
                     ELSE oplot,[alog10(metals[ind].smass),alog10(metals[ind].smass)],[metallicity,metallicity],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i]*0.75, thick = thicks[i]
                   print,filenames[i],' ',strtrim(halos[j,i],2),' ',strtrim(alog10(metals[ind].smassO),2),' ',strtrim(alog10(metals[ind].smass),2),' ',strtrim(metals[ind].ox_sfr,2),' ',strtrim(zstellar[ind],2),strtrim(zstellar_rband[ind],2)
               ENDIF
           ENDIF
       ENDFOR
;       stop
   ENDIF ELSE BEGIN
       indsat = where(metals.sat EQ 1, complement = indnosat)
       IF (indsat[0] EQ -1) THEN indnosat = indgen(n_elements(metals))
       IF NOT keyword_set(noObsStellarMass) THEN BEGIN
           IF (indnosat[0] NE -1) THEN BEGIN
               IF n_elements(indnosat) EQ 1 THEN oplot,[alog10(metals[indnosat].smassO),alog10(metals[indnosat].smassO)],[metals[indnosat].ox_sfr,metals[indnosat].ox_sfr],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] $
               ELSE oplot,alog10(metals[indnosat].smassO),metals[indnosat].ox_sfr,psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i]
               ind_plotable = where(finite(metals[indnosat].ox_sfr))
               IF NOT ind_plotable[0] EQ -1 THEN BEGIN
                   IF n_elements(indnosat[ind_plotable]) EQ 1 THEN BEGIN
                       IF keyword_set(simStellarMass) THEN $
                         oplot,[alog10(metals[indnosat[ind_plotable]].smass),alog10(metals[indnosat[ind_plotable]].smass)],[metals[indnosat[ind_plotable]].ox_inner,metals[indnosat[ind_plotable]].ox_inner],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] $
                       ELSE $
                         oplot,[alog10(metals[indnosat[ind_plotable]].smassO),alog10(metals[indnosat[ind_plotable]].smassO)],[metals[indnosat[ind_plotable]].ox_inner,metals[indnosat[ind_plotable]].ox_inner],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i]
                   ENDIF ELSE BEGIN
                       IF keyword_set(simStellarMass) THEN $
                         oplot,alog10(metals[indnosat[ind_plotable]].smass),metals[indnosat[ind_plotable]].ox_inner,psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] $
                       ELSE $
                         oplot,alog10(metals[indnosat[ind_plotable]].smassO),metals[indnosat[ind_plotable]].ox_inner,psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i]
                   ENDELSE
                   ox_plot = metals[indnosat[ind_plotable]].ox_inner
               ENDIF ELSE ox_plot = metals[indnosat].ox_sfr
               IF keyword_set(onehalos) THEN print,filenames[i],' ',strtrim(onehalos[i],2),alog10(metals[indnosat].smassO),ox_plot ELSE print,filenames[i],alog10(metals[indnosat[ind_plotable]].smassO),ox_plot
           ENDIF
;           IF (indsat[0] NE -1) THEN $
;             IF n_elements(indnosat) EQ 1 THEN oplot,[alog10(metals[indsat].smassO),alog10(metals[indsat].smassO)],[metals[indsat].ox_sfr,metals[indsat].ox_sfr],psym = symcat(psym[i] + 1), color = colors[i], symsize = symsizes[i], thick = thicks[i] ELSE oplot,alog10(metals[indsat].smassO),metals[indsat].ox_sfr,psym = symcat(psym[i] + 1), color = colors[i], symsize = symsizes[i], thick = thicks[i]
       ENDIF
;       IF keyword_set(simStellarMass) THEN oplot,[alog10(metals[indsat].smass),alog10(metals[indsat].smass)],[metals[indsat].ox_sfr,metals[indsat].ox_sfr],psym = symcat(psym[i] + 1), color = colors[i], symsize = 0.5, thick = thicks[i]
;       oplot,alog10(metals.smassO),metals.ox_sfr,psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i]*0.75, thick = thicks[i]       
;       print,alog10(metals[indnosat].smassO),alog10(metals.smass),metals[indnosat].ox_sfr
   ENDELSE
   IF NOT keyword_set(readfile) OR NOT file_test(filenames[i] + '.metals.fits') THEN mwrfits,metals,filenames[i] + '.metals.fits',/create
ENDFOR

IF (keyword_set(outfile)) THEN device,/close ELSE stop
END


