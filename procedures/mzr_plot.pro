
;filenames=['/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir/h516.cosmo25cmb.3072g14HBWK.00492','/home/christensen/Storage2/UW/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512','/home/christensen/Storage2/UW/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512','/home/christensen/Storage2/UW/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512']
;psym=[4,5,6,9] ;11,12,13x
;outfile = '~/mz_all.eps'

;filenames=['/astro/store/nbody2/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir/h516.cosmo25cmb.3072g14HBWK.00492','/astro/store/nbody2/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512','/astro/store/nbody2/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512','/astro/store/nbody2/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512','/astro/store/nbody2/abrooks/h277.cosmo50cmb.3072g14HMbwK/h277.cosmo50cmb.3072g14HMbwK.00512/h277.cosmo50cmb.3072g14HMbwK.00512','/astro/store/nbody3/abrooks/h285.cosmo50cmb.3072g14HMbwK/h285.cosmo50cmb.3072g14HMbwK.00512/h285.cosmo50cmb.3072g14HMbwK.00512']
;psym=[4,5,6,9,11,12,13]
;psym=[6,6,6,6,6, 6, 6]
;outfile = '~/plots/all'
;mzr_plot,filename,psym = psym, /readfile, outfile = outfile

PRO mzr_plot,filenames, halos = halos, outfile = outfile, key = key, psym  = psym, color = color, symsizes = symsizes, ctables = ctables, thicks = thicks, simStellarMass = simStellarMass, noObsStellarMass = noObsStellarMass, xrange = xrange, yrange = yrange, readfile = readfile, obscolor = obscolor, onehalos = onehalos, obsct = obsct, formatthick = formatthick
n = n_elements(filenames)

IF keyword_set(outfile) THEN BEGIN
    l_charsize = 1.2 ;1 ;0.75
    fgcolor = 0
    bgcolor = 255
    xsize = 18;7
    ysize = 12;5
    formatplot,/outplot,thick = formatthick
ENDIF ELSE BEGIN
    l_charsize = !p.charsize ;1.0
    fgcolor = 255
    bgcolor = 0
    xsize = 600
    ysize = 400
    formatplot,thick = formatthick
ENDELSE
IF keyword_set(color) THEN BEGIN
    IF NOT keyword_set(obsct) THEN obsct = 0;39
    IF NOT keyword_set(ctables) THEN ctables = 39 + fltarr(n)
    IF color[0] EQ 1 THEN  colors = (findgen(n) + 1)*240/n else colors = fltarr(n) + color
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(psym) THEN psym = fltarr(n) + 4 ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
    IF NOT keyword_set(obscolor) THEN obscolor = fgcolor
    obssymsize = 1.25;2;.5
    obssym = 16
    obsthick = thicks[0];2
ENDIF ELSE BEGIN
    IF NOT keyword_set(obsct) THEN obsct = 0 
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n)
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(psym) THEN psym = findgen(n) + 4;(findgen(n) + 2)*2   
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
    IF NOT keyword_set(obscolor) THEN obscolor = 150
    obssymsize = 1.25;2;.5
    obssym = 16
    obsthick = thicks[0];2
ENDELSE
IF NOT keyword_set(xrange) THEN xrange = [5.5,11.5]
IF NOT keyword_set(yrange) THEN yrange = [7,9.1]

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN datafile_base = '/home/christensen/Code/Datafiles/' ELSE datafile_base = '/astro/users/christensen/code/Datafiles/'
readcol,datafile_base + 'Lee06.txt',gal,source,aorkey,f,ferr,mag,magerr,disref,magA,magAerr,logOH,logOHerr,OPHref,B,logMstar,format = '(A,A,I,F,F,F,F,I,F,F,F,F,F,F,F,F)'
readcol,datafile_base + 'Tremonti04.txt',LogMstarT,P_25,P_16,P_50,P_84,P_975

IF (keyword_set(outfile)) THEN BEGIN
   device,filename=outfile+'_mz.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2,/encapsul
ENDIF ELSE BEGIN
   window,0,xsize = xsize,ysize = ysize ;392
ENDELSE

LogMstar_arr = findgen(101)/100.0*(8.5 - 5.5) + 5.5
loadct,obsct
plot,logMstar,logOH,xtitle = textoidl('log(M_{star} [M') + sunsymbol() + '])', ytitle = '12 + log(O/H)',/nodata,xrange = xrange,yrange = yrange
oploterror,logMstar,logOH,logMstar*0,logOHerr,psym = symcat(obssym), color = obscolor, symsize = obssymsize, thick = obsthick,errcolor=obscolor
oplot,LogMstarT,P_25  - 0.26,color = obscolor, thick = obsthick, linestyle = 1
oplot,LogMstarT,P_16  - 0.26,color = obscolor, thick = obsthick, linestyle = 2
oplot,LogMstarT,P_50  - 0.26,color = obscolor, thick = obsthick, linestyle = 0
oplot,LogMstarT,P_84  - 0.26,color = obscolor, thick = obsthick, linestyle = 2
oplot,LogMstarT,P_975 - 0.26,color = obscolor, thick = obsthick, linestyle = 1
;oplot,LogMstar_arr, -1.492 + 1.847*LogMstar_arr - 0.08026*LogMstar_arr*LogMstar_arr - 0.26,color = obscolor, thick = obsthick, linestyle = 3

LogMstarAM = findgen(32)/10+7.4
mzrAM = 8.798 - alog10(1 + (10^8.901/10^LogMstarAM)^0.64) ;fit from Andrew & Martini 2013, see table 4
oplot,LogMstarAM,mzrAM,color = obscolor, thick = obsthick, linestyle = 3
FOR i = 0, n - 1 DO BEGIN
    loadct,ctables[i]
    IF keyword_set(onehalos) THEN metalfile = filenames[i] + '.halo.' + strtrim(onehalos[i],2) ELSE  metalfile = filenames[i]
    IF keyword_set(readfile) AND file_test(metalfile + '.metals.fits') THEN metals = mrdfits(metalfile + '.metals.fits',1) ELSE BEGIN
        IF keyword_set(onehalos) THEN IF keyword_set(noObsStellarMass) THEN metals = mzr(filenames[i], onehalo = onehalos[i]) ELSE metals = mzr(filenames[i],/obs, onehalo = onehalos[i]) ELSE $
          IF keyword_set(noObsStellarMass) THEN metals = mzr(filenames[i]) ELSE metals = mzr(filenames[i],/obs) 
        mwrfits,metals,metalfile + '.metals.fits',/create
    ENDELSE
    IF keyword_set(halos) THEN BEGIN
       FOR j = 0, n_elements(halos[*,i]) - 1 DO BEGIN
           IF halos[j,i] NE 0 THEN BEGIN
               ind = where(metals.grp EQ halos[j,i])
               IF ind[0] NE -1 THEN BEGIN
                   IF NOT keyword_set(noObsStellarMass) THEN BEGIN
;                   IF metals[ind].sat EQ 0 THEN 
                       oplot,[alog10(metals[ind].smassO),alog10(metals[ind].smassO)],[metals[ind].ox_sfr,metals[ind].ox_sfr],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] ;ELSE oplot,[alog10(metals[ind].smassO),alog10(metals[ind].smassO)],[metals[ind].ox_sfr,metals[ind].ox_sfr],psym = symcat(psym[i]) + 1, color = colors[i], symsize = symsizes[i], thick = thicks[i]
                   ENDIF
                   IF keyword_set(simStellarMass) THEN oplot,[alog10(metals[ind].smass),alog10(metals[ind].smass)],[metals[ind].ox_sfr,metals[ind].ox_sfr],psym = symcat(psym[i]), color = colors[i], symsize = 0.5*symsizes[i], thick = thicks[i]              
;               oplot,[alog10(metals[ind].smassO),alog10(metals[ind].smassO)],[metals[ind].ox_sfr,metals[ind].ox_sfr],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i]*0.75, thick = thicks[i]
                   print,halos[j,i],alog10(metals[ind].smassO),alog10(metals[ind].smass),metals[ind].ox_sfr
                   stop
               ENDIF
           ENDIF
       ENDFOR
   ENDIF ELSE BEGIN
       indsat = where(metals.sat EQ 1, complement = indnosat)
       IF (indsat EQ -1) THEN indnosat = indgen(n_elements(metals))
       IF NOT keyword_set(noObsStellarMass) THEN BEGIN
           IF (indnosat[0] NE -1) THEN BEGIN
               IF n_elements(indnosat) EQ 1 THEN oplot,[alog10(metals[indnosat].smassO),alog10(metals[indnosat].smassO)],[metals[indnosat].ox_sfr,metals[indnosat].ox_sfr],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] $
               ELSE oplot,alog10(metals[indnosat].smassO),metals[indnosat].ox_sfr,psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i]
               IF NOT finite(metals[indnosat].ox_sfr) THEN BEGIN
                   IF n_elements(indnosat) EQ 1 THEN oplot,[alog10(metals[indnosat].smassO),alog10(metals[indnosat].smassO)],[metals[indnosat].ox_inner,metals[indnosat].ox_inner],psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i] $
                   ELSE oplot,alog10(metals[indnosat].smassO),metals[indnosat].ox_inner,psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i], thick = thicks[i]      
               ENDIF
           ENDIF
           IF (indsat[0] NE -1) THEN $
             IF n_elements(indnosat) EQ 1 THEN oplot,[alog10(metals[indsat].smassO),alog10(metals[indsat].smassO)],[metals[indsat].ox_sfr,metals[indsat].ox_sfr],psym = symcat(psym[i] + 1), color = colors[i], symsize = symsizes[i], thick = thicks[i] ELSE oplot,alog10(metals[indsat].smassO),metals[indsat].ox_sfr,psym = symcat(psym[i] + 1), color = colors[i], symsize = symsizes[i], thick = thicks[i]
       ENDIF
       IF keyword_set(simStellarMass) THEN oplot,[alog10(metals[indsat].smass),alog10(metals[indsat].smass)],[metals[indsat].ox_sfr,metals[indsat].ox_sfr],psym = symcat(psym[i] + 1), color = colors[i], symsize = 0.5, thick = thicks[i]
;       oplot,alog10(metals.smassO),metals.ox_sfr,psym = symcat(psym[i]), color = colors[i], symsize = symsizes[i]*0.75, thick = thicks[i]       
       print,alog10(metals[indnosat].smassO),alog10(metals.smass),metals[indnosat].ox_sfr
   ENDELSE
   IF NOT keyword_set(readfile) OR NOT file_test(filenames[i] + '.metals.fits') THEN mwrfits,metals,filenames[i] + '.metals.fits',/create
;   stop
ENDFOR
uniqind0 = (uniq(psym,sort(psym)))[0]
uniqind1 = (uniq(psym,sort(psym)))[1]
IF keyword_set(key) THEN legend,[key,"Lee et al., 2006","Tremonti et al., 2004",'Andrew & Martini, 2013'],color = [color,obscolor,obscolor,obscolor],linestyle = intarr(n_elements(key) + 3),psym = [psym,obssym,0,0],thick = [thicks,obsthick,obsthick,obsthick] $
ELSE legend,['Lee et al., 2006','Tremonti et al., 2004','Andrew & Martini, 2013','Med-res Sims','High-res Sims'],color = [obscolor,obscolor,obscolor,colors[uniqind0],colors[uniqind1]],linestyle = [0,0,3,2,0],psym = [obssym,0,0,psym[uniqind0],psym[uniqind1]],thick = [obsthick,obsthick,obsthick,thicks[uniqind0],thicks[uniqind1]],ctables = [obsct,obsct,obsct,ctables[uniqind0],ctables[uniqind1]],box = 0,symsize = [obssymsize,obssymsize,obssymsize,symsizes[uniqind0],symsizes[uniqind1]],charsize = l_charsize,/bottom,/right,position = [0.934145,0.250005],/norm ;corners = corners;0.934145     0.503365     0.460109     0.271045
stop
IF (keyword_set(outfile)) THEN device,/close ELSE stop
END


