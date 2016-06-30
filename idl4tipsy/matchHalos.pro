FUNCTION matchHalos,files,minmass = minmass, dmmass = dmmass

halos = 0
array0 = read_stat_struc_amiga(files[0])
array0 = array0[where(array0.sat ne 'yes')]

array1 = read_stat_struc_amiga(files[1])
array1 = array1[where(array1.sat ne 'yes')]
center1all = TRANSPOSE([[array1.xc],   [array1.yc],   [array1.zc]])

IF NOT KEYWORD_SET(minmass) THEN minmass = 0.01 * array0[0].m_tot
temp = MIN(ABS(array0.m_tot - minmass),maxhalo)
matches = intarr(maxhalo,2)
mindistances = fltarr(maxhalo)
massratios = fltarr(maxhalo)
dmparticles = fltarr(maxhalo)
stars = fltarr(maxhalo)
FOR i = 0, maxhalo - 1 DO BEGIN
    IF KEYWORD_SET(dmmass) THEN BEGIN
        center0 = [array0[i].xc,array0[i].yc,array0[i].zc]
        distances = fltarr(N_ELEMENTS(center1all[0,*]))
        FOR j = 0, N_ELEMENTS(distances) - 1 DO distances[j] = sqrt((center1all[0,j] - center0[0])^2 + (center1all[1,j] - center0[1])^2 + (center1all[2,j] - center0[2])^2)        
        temp = min(abs(array0[i].m_dark - array1.m_dark),j)
        matches[i,*] = [array0[i].group,array1[j].group]
        mindistances[i] = [distances[j]]
        massratios[i] = array1[j].m_dark/array0[i].m_dark
        stars[i] = MIN([array1[j].m_star,array0[i].m_star])
        dmparticles[i] = MIN([array1[j].ndark,array0[i].ndark]) 
    ENDIF ELSE BEGIN
        center0 = [array0[i].xc,array0[i].yc,array0[i].zc]
        closemass = where(array1.m_dark ge array0[i].m_dark*0.33 AND array1.m_dark le array0[i].m_dark*3.0)
        IF closemass[0] ne -1 THEN BEGIN
            center1 = center1all[*,closemass]
            array1close = array1[closemass]           
            distances = fltarr(N_ELEMENTS(center1[0,*]))
            FOR j = 0, N_ELEMENTS(distances) - 1 DO distances[j] = sqrt((center1[0,j] - center0[0])^2 + (center1[1,j] - center0[1])^2 + (center1[2,j] - center0[2])^2)
            temp = min(distances,j)
            matches[i,*] = [array0[i].group,array1close[j].group]
            mindistances[i] = [distances[j]]
            massratios[i] = array1close[j].m_dark/array0[i].m_dark
            stars[i] = MIN([array1close[j].m_star,array0[i].m_star])
            dmparticles[i] = MIN([array1close[j].ndark,array0[i].ndark])
        ENDIF ELSE matches[i,*] = [array0[i].group,0] 
    ENDELSE
ENDFOR
FOR i = 0, N_ELEMENTS(matches[*,0]) - 1 DO BEGIN
    print,matches[i,0],matches[i,1],mindistances[i],massratios[i],stars[i],dmparticles[i]
ENDFOR
return,matches
END
