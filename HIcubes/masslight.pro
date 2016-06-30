;This program will call kcorrect to find the observationally
;determined stellar masses given the magnitude of the galaxy.
;It will then print a plot of the stellar mass to dark matter mass in
;each halo

;For it to run, you will need to have kcorrect installed on your system

;Call by referenceing a datafile.
;This datafile should contain one or more lines of text, the text
;containing the information for each simulation you would like
;included
;The format of the data file is:
;dir,sunrise,file,massunit,kpcunit,massform,color,psym,key
;whre dir is the directory of the simulation
;sunrise is the subdirectory containing broadband.fits (it dosen't matter what you put here if you are not using SR mags)
;file is the filename of the simulation
;massunit is the sim. mass unit in solar masses
;kpcunit is the sim. length unit in kpc
;massform is the initial mass of star particles in sim. units
;color is the color you would like this galaxy plotted with
;psym is the symbol to be used when plotting the galaxy
;key is what you would like to call the galaxy in the legend

;If you would like to save your plot as a .ps file, set outplot to the name of the file
;Use the sunrise flag if you would like to use sunrise generated magnitudes
;Use the sdss flag if you would like base you masses off of sdss magnitudes as opposed to Johnson/Bessel (only use with Sunrise)
;Use the calculate flag if you would like the program to use the get_halo_mags idl script to calculate the galaxy magnitudes

;Examples:
;masslight,'schmidtlaw_files_sdss.txt',/sunrise,/sdss
;masslight,'schmidtlaw_files_bessel.txt',/calculate
;masslight,datafile


; .r /astro/users/christensen/code/kcorrect/pro/fit/sdss_kcorrect.pro

pro masslight,datafile,outplot = outplot,sunrise = sunrise,sdss = sdss,calculate = calculate
loadct,39

readcol,'/astro/users/christensen/code/HIcubes/'+datafile,dir,sunrise,file,massunit,kpcunit,massform,color,psym,key,format = '(A60,A33,A33,F,F,F,F,F,D,A18)'
outbase = '/astro/net/nbody1/christensen/Schmidtlaw/data/'
mags = fltarr(5,N_ELEMENTS(dir))
mags_bes = fltarr(5,N_ELEMENTS(dir))
mags_errs = fltarr(5,N_ELEMENTS(dir))
guo_mass = fltarr(3,N_ELEMENTS(dir))
FMT = 'X,X,X,X,X,F,X,F,F,F,X,X,X,X,X,X,X,X,X,X,A'

IF (NOT KEYWORD_SET(sr)) AND KEYWORD_SET(calculate) THEN BEGIN
    FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
        get_halo_mags,massunit[i],kpcunit[i],iso_set = 0, init_mass = massform[i],mag_sys='ab',tipsy_file = dir[i]+'/'+file[i]
    ENDFOR
ENDIF

FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
    IF KEYWORD_SET(sr) THEN BEGIN
        IF KEYWORD_SET(sdss) THEN BEGIN
            filters = mrdfits(dir[i]+'/'+sunrise[i]+'/broadband.fits',13)
            mags[*,i] = [filters[2].AB_MAG_NONSCATTER0,filters[3].AB_MAG_NONSCATTER0,filters[4].AB_MAG_NONSCATTER0,filters[5].AB_MAG_NONSCATTER0,filters[6].AB_MAG_NONSCATTER0]
            mags_errs[*,i] = [filters[2].AB_MAG0,filters[3].AB_MAG0,filters[4].AB_MAG0,filters[5].AB_MAG0,filters[6].AB_MAG0]

            filtersbes = mrdfits(dir[i] + '/' + file[i]+ '.star99_ab.Mv.fits',1)
            mags_bes[*,i] = [filtersbes.U,filtersbes.B,filtersbes.V,filtersbes.R,filtersbes.I]
        ENDIF
    ENDIF ELSE BEGIN
        filters = mrdfits(dir[i] + '/' + file[i]+ '.star99_ab.Mv.fits',1)
        mags[*,i] = [filters.U,filters.B,filters.V,filters.R,filters.I]
    ENDELSE
    mags_errs[*,i] = [0.05, 0.02, 0.02, 0.02, 0.03] ;From the kcorrect website ;ABS(sdss_errs[*,i] - sdss_mags[*,i])
    IF i eq 0 THEN plot,[1,4,6,8,10],mags[*,i],ytitle = "Magnitudes",xtitle = "Filters",xrange = [0,11],yrange = [-12,-25],psym = psym[i] + 3,xtickname = ['','u','U','B','g','V','r','R','i','I','z',''],xticks = 11,xstyle = 1
    oplot,[1,4,6,8,10],mags[*,i],psym = psym[i]+3,color = color[i]
    oplot,[2,3,5,7,9],mags_bes[*,i],psym = psym[i]+1,color = color[i]
;    oplot,[3,5,7],[filters[22].AB_MAG_NONSCATTER0,filters[23].AB_MAG_NONSCATTER0,filters[24].AB_MAG_NONSCATTER0],psym = psym[i]+3,color = color[i]
;    oploterr,[1,2,3,4,5],mags[*,i],mags_errs[*,i],color = color[i]
    readcol,dir[i]+'/'+file[i] + '.amiga.stat',F=FMT,vir,gas,star,dark,contam,/silent 
    guo_mass[*,i] = [dark[0],star[0],0]
ENDFOR
legend,key,psym = psym + 3,color = color
z = fltarr(N_ELEMENTS(dir))
dmod = 19.4576 ;distance modulus for z = 0 in this model


;****************** Using kcorrect to find the mass from the magnitudes.  Mag/mgy are the magnitudes and the output we are interested in is mass
IF KEYWORD_SET(sdss) THEN BEGIN
    kcorrect = sdss_kcorrect(z,mag = mags + dmod, err = mags_errs,mass = mass, mtol = mtol,absmag = absmag,rmaggies = rmaggies)
    kcorrect = sdss2bessell(z,mag = mags + dmod, err = mags_errs,mass = massbes, absmag = absmag_bessel)
ENDIF ELSE BEGIN
    mgy=(10.D)^(-(0.4D)*(mags + dmod))
    mags_ivar=1./mags_errs^2
    mgy_ivar= mags_ivar/(0.4*alog(10.)*mgy)^2.
    kcorrect, mgy, mgy_ivar, z, kcorrect, mass = mass, mtol = mtol, absmag = absmag, filterlist = ['bessell_U.par','bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par']
ENDELSE

FOR i = 0, N_ELEMENTS(dir) - 1 DO oplot,[2,3,5,7,9],absmag_bessel[*,i],psym = psym[i]+2,color = color[i]

guo_mass[2,*] = mass

IF KEYWORD_SET(outplot) THEN BEGIN
   set_plot, 'ps'
   !P.THICK=1.5                 ;4
   !P.CHARTHICK=1.5             ;4
   !X.THICK=1.5                 ;4
   !Y.THICK=1.5                 ;4
   !p.charsize=1.2              ;1.8
   !p.font=0
   device, filename=outbase+'guo.ps',/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch
ENDIF ELSE window,1

x = findgen(700)/100 + 9
y = alog10(0.1257*((10^x/10^(11.36))^(-0.9147) + (10^x/10^(11.36))^(0.2485))^(-2.574)*10^x);guo 09, (3)
plot,x,y,xstyle = 1, ystyle = 1,xrange = [9.5,13],yrange = [6,11],xtitle = textoidl('log(M_{halo}[M')+sunsymbol()+'])',ytitle = textoidl('log(M_{*}[M')+sunsymbol()+'])'
FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
;    oplot,[alog10(guo_mass[0,i]),alog10(guo_mass[0,i])],[alog10(guo_mass[1,i]),alog10(guo_mass[1,i])],psym = psym[i]+3,color = color[i]
    oplot,[alog10(guo_mass[0,i]),alog10(guo_mass[0,i])],[alog10(guo_mass[2,i]),alog10(guo_mass[2,i])],psym = psym[i]+3,color = color[i]
;    oplot,[alog10(guo_mass[0,i]),alog10(guo_mass[0,i])],[alog10(massbes[i]),alog10(massbes[i])],psym = psym[i]+1,color = color[i]
ENDFOR
legend,key,psym = psym+3,color = color,/bottom,/right
IF KEYWORD_SET(outplot) THEN device,/close
stop
END
