FUNCTION VERTICAL_PRO,particles, NBINS = NBINS, ZMAX = ZMAX
IF NOT KEYWORD_SET(nbins) THEN nbins = 500.
IF NOT KEYWORD_SET(zmax) THEN zmax = 20.
dz = zmax/nbins
kpc_unit = 1
particles.z = particles.z*kpc_unit

profile = fltarr(2,nbins)
profile[0,*] = (findgen(nbins)+1)*dz
between = where(ABS(particles.z) LE zmax)
IF (between[0] eq -1) THEN RETURN, profile
totalmass = total(particles[between].mass)
FOR i = 0, nbins-1 DO BEGIN
    ind = where(ABS(particles.z) LE (profile[0,i]))
    IF (ind[0] ne -1) THEN profile[1,i] = TOTAL(particles[ind].mass)/totalmass ELSE  profile[1,i] = 0
ENDFOR

RETURN,profile
END


;files = ['../ResMassTests/1E5R/10M/o10M_1.00300']
;comp = 's'
;nbins = '10'
;maxradii = '10'
PRO VERTICAL_SL, files, comp, NBINS = NBINS, MAXRADII = maxradii, OUTPLOT = OUTPLOT, TITLE = TITLE, XTITLE = XTITLE, YTITLE = YTITLE, LEGEND_NAME = LEGEND_NAME, YRANGE = YRANGE, COLORS = COLORS, PRINT = PRINT, LINESTYLE = LINESTYLE, XRANGE=XRANGE,thickness=thickness
;This will plot half mass and ninety-percent mass heights of stellar
;or gas disks
msol = 2.362e5
kpc_unit = 1
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
loadct,39
fudge = 0.0001

IF NOT KEYWORD_SET(nbins) THEN nbins = 10.
IF NOT KEYWORD_SET(maxradii) THEN maxradii = 10.
;IF (keyword_set(print)) THEN BEGIN 
;    set_plot,'ps'
;    print,"plot"
;    device,filename = outplot,/color,bits_per_pixel = 8
;ENDIF ELSE set_plot,'x'

radii = findgen(nbins + 1)*maxradii/nbins
sl90 = findgen(nbins)
sl50 = findgen(nbins)
nfiles = N_ELEMENTS(files)

FOR i = 0, nfiles-1 DO BEGIN
    print,files[i]
    rtipsy, files[i], h, g, d, s
    cofm, h, g, d, s, /pot
    particles = gal_align( h, g, d, s, rbox = 1, lunit = 1)
    h = particles.h
    g = particles.g
    d = particles.d
    s = particles.s
    FOR ibin = 0, nbins - 1 DO BEGIN
        print,radii[ibin],":",radii[ibin+1]
        between = where((SQRT(s.x^2 + s.y^2)*kpc_unit GT radii[ibin]) AND (SQRT(s.x^2 + s.y^2)*kpc_unit LT radii[ibin+1]))
        IF (between[0] eq -1 ) THEN BEGIN
            sl90[ibin] = -1 
            sl50[ibin] = -1
        ENDIF ELSE BEGIN
            IF (comp eq 's') THEN particles = s[between]
            IF (comp eq 'g') THEN particles = g[between]
            profile = vertical_pro(particles)
            IF (MAX(profile) EQ 0) THEN BEGIN
                sl90[ibin] = -1 
                sl50[ibin] = -1
            ENDIF ELSE BEGIN
                temp50 = MIN(ABS(profile[1,*]  - 0.5),i50)
                IF (i50 eq 0 OR i50 eq nbins - 1) THEN sl50[ibin] = profile[1,i50] ELSE BEGIN
                    WHILE(profile[1,i50+1]lt 0.5 AND i50+1 lt nbins -1) DO i50 = i50 + 1
                    WHILE(profile[1,i50-1]gt 0.5 AND i50-1 gt 0) DO i50 = i50 - 1
                    IF(profile[1,i50] EQ profile[1,i50 + 1]) THEN sl50[ibin] = profile[0,i50-1]*(0.5 - profile[1,i50-1])/(profile[1,i50] - profile[1,i50-1]) + profile[0,i50]*(profile[1,i50] - 0.5)/(profile[1,i50] - profile[1,i50-1]) ELSE BEGIN
                        IF(profile[1,i50] EQ profile[1,i50 - 1]) THEN sl50[ibin] = profile[0,i50]*(0.5 - profile[1,i50])/(profile[1,i50+1] - profile[1,i50]) + profile[0,i50+1]*(profile[1,i50+1] - 0.5)/(profile[1,i50+1] - profile[1,i50]) ELSE sl50[ibin] = SPLINE(profile[1,i50-1:i50+1],profile[0,i50-1:i50+1],0.5)*Kpc_Unit 
                    ENDELSE
                ENDELSE

                temp90 = MIN(ABS(profile[1,*]  - 0.9),i90)
                IF (i90 eq 0 OR i90 eq nbins - 1) THEN sl90[ibin] = profile[1,i90] ELSE BEGIN   
                    WHILE(profile[1,i90+1]lt 0.9 AND i90+1 lt nbins - 1) DO i90 = i90 + 1
                    WHILE(profile[1,i90-1]gt 0.9 AND i90-1 gt 0) DO i90 = i90 - 1
                    IF(profile[1,i90] EQ profile[1,i90 + 1]) THEN sl90[ibin] = profile[0,i90-1]*(0.9 - profile[1,i90-1])/(profile[1,i90] - profile[1,i90-1]) + profile[0,i90]*(profile[1,i90] - 0.9)/(profile[1,i90] - profile[1,i90-1]) ELSE BEGIN
                        IF(profile[1,i90] EQ profile[1,i90 - 1]) THEN sl90[ibin] = profile[0,i90]*(0.9 - profile[1,i90])/(profile[1,i90+1] - profile[1,i90]) + profile[0,i90+1]*(profile[1,i90+1] - 0.9)/(profile[1,i90+1] - profile[1,i90]) ELSE sl90[ibin] = SPLINE(profile[1,i90-1:i90+1],profile[0,i90-1:i90+1],0.9)*Kpc_Unit 
                    ENDELSE
                ENDELSE
            ENDELSE
        ENDELSE
;        print,sl50[ibin],sl90[ibin]
;        oplot,[sl50[ibin],sl50[ibin]],[0,1],linestyle = 2
;        oplot,[sl90[ibin],sl90[ibin]],[0,1],linestyle = 3       
;        stop
    ENDFOR
    IF (i eq 0) THEN plot,radii[where(sl90 ne -1)]+maxradii/nbins/2,sl90[where(sl90 ne -1)],title = title, xtitle = xtitle, ytitle = ytitle ,linestyle = linestyle[i],yrange = yrange,color = colors[i],thick=thickness[i],xrange = xrange,ytickinterval = 5;[0,11]
    oplot,radii[where(sl90 ne -1)]+maxradii/nbins/2,sl90[where(sl90 ne -1)],linestyle = linestyle[i], color = colors[i],thick=thickness[i]
;    oplot,radii[where(sl50 ne -1)]+maxradii/nbins/2,sl50[where(sl50 ne -1)],linestyle = 1, color = colors[i]
ENDFOR
IF (KEYWORD_SET(legend_name)) THEN legend,legend_name,color = colors,linestyle = linestyle,/top,/right,charsize = 1,thick=thickness,pspacing = 1.5*1.25
;IF (keyword_set(print)) THEN device,/close
END

PRO master_vertical_sl
cd,'/astro/net/scratch1/christensen/DwarfResearch/procedures'
;This master plotter will loop through interesting files

!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1
!Y.Charsize = 1
!X.style = 1
base = ['../ResMassTests/','../ResSpecTests/1E5R/']
imf = ['k']
type = ['s']
type_str = ['Stellar']
legend_name0 = [textoidl('10^6'),textoidl('10^5'),textoidl('10^4'),textoidl('10^3')]
legend_name1 = [textoidl('2.5\times10^{-4} R_{vir}'),textoidl('10^{-3} R_{vir} '),textoidl('10^{-2} R_{vir} '),textoidl('4\times10^{-2} R_{vir}')]
res_str0 = ['1E6R','1E5R','1E4R','1E3R']
res_str1 = ['50pc','200pc','2kpc','8kpc']

colors0 = [0,0,0,0]
linestyle0 = [0,0,2,1]
colors1 = [0,0,0,0]
linestyle1 = [0,0,2,1]
yrange0 = [0.0001,15]
yrange1 = [0,5]
mass_str = ['12','10']
maxradii = [10,6]
xrange=[0,11]
thickness = [5,1,5,1]

;set_plot,'x'
set_plot,'ps'
device,filename = '../results/s_Spec_scaleHeight_k.eps',/color,bits_per_pixel = 8,/times
!P.MULTI = [0]         
!P.MULTI = [0,2,2,0,0] 
!P.FONT = 0

FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
    FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
        FOR tct = 0, N_ELEMENTS(type) - 1 DO BEGIN
            FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN
                IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00300' ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
  ;          outplot = '../results/' + mass_str[mct] + 'M_' + type[tct] + 'Num_scaleHeight_' + imf[imfct] + '.eps'
 ;           outplot = '../results/' + mass_str[mct] + 'M_' + type[tct] + 'Spec_scaleHeight_' + imf[imfct] + '.eps'
 ;           print,outplot

;            title = type_str[tct]+' Scale Height, '+textoidl('10^{'+mass_str[mct]+'}') + ' M' +sunsymbol()+ ' Halo'
                xtitle = 'Radius [kpc]'
                ytitle = textoidl('z_{*,90} [kpc]')
;            print,"type: ",type[tct],", mass: ",mass_str[mct]
;            print,files
                multiplot
                IF (mct eq 0) THEN BEGIN 
                    IF (ploti eq 0) THEN vertical_sl, files, type[tct], linestyle = linestyle0,outplot = outplot, title = title, xtitle = '',ytitle = ytitle,colors = colors0,yrange = yrange0,xrange=xrange,maxradii=maxradii[0],thickness=thickness ELSE vertical_sl, files, type[tct], linestyle = linestyle1,outplot = outplot, title = title,colors = colors1,yrange = yrange0, legend_name = legend_name1,xrange=xrange,maxradii=maxradii[0],thickness=thickness
                ENDIF ELSE BEGIN
                    IF (ploti eq 0) THEN vertical_sl, files, type[tct], linestyle = linestyle0,outplot = outplot, title = title, xtitle = xtitle, ytitle = ytitle, colors = colors0,yrange = yrange1,legend_name = legend_name0,xrange=xrange,maxradii=maxradii[1],thickness=thickness ELSE vertical_sl, files, type[tct], linestyle = linestyle1,outplot = outplot, title = title, xtitle = xtitle, colors = colors1, yrange = yrange1,xrange=xrange,maxradii=maxradii[1],thickness=thickness
                ENDELSE
            ENDFOR
        ENDFOR
    ENDFOR
ENDFOR
xyouts,[3000,100100,3000,10100],[11000,11000,6000,6000],[textoidl('10^{12}M')+sunsymbol(),textoidl('10^{12}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol()],charsize = 1.05,/device
device,/close
multiplot,/reset
cd,'/astro/users/christensen/code/Dwarfs'
END
