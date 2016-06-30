;CC, 12/12/11
;This program reads in stat files to create an array (t,z,x_c,y_c,z_c)

FUNCTION read_center,files = files,halos = halos
ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)

;find files
IF NOT keyword_set(files) THEN spawn,'ls */*amiga.stat',files
center = fltarr(N_ELEMENTS(files),5)


FOR i = 0, N_ELEMENTS(files) - 1 DO BEGIN
    dotpos = strsplit(files[i],'.')
    base = strmid(files[i],0,dotpos[N_ELEMENTS(dotpos) - 2] - 1)
    rtipsy,base,h,g,d,s,/justhead
    a = h.time
    z = (1 - a)/a
    time = (ageUniverse - wmap3_lookback(z))/1e9 ;in Gyrs
    satsdata = read_stat_struc_AMIGA(files[i])
    IF keyword_set(halos) THEN j = where(satsdata.group EQ halos[i]) ELSE j = 0
    center[i,*] = [time,z,satsdata[j].xc,satsdata[j].yc,satsdata[j].zc]
ENDFOR

RETURN,center
END
