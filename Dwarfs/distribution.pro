PRO DISTRIBUTION, files, type, validT,validD,valid,tdyn,NBINS = NBINS, MAXRADII = maxradii, OUTPLOT = OUTPLOT, TITLE = TITLE, XTITLE = XTITLE, YTITLE = YTITLE, LEGEND_NAME = LEGEND_NAME, YRANGE = YRANGE, COLORS = COLORS, PRINT = PRINT, LINESTYLE = LINESTYLE,THICKNESS = THICKNESS
;This will plot half mass and ninety-percent mass heights of stellar
;or gas disks
msol_unit = 2.362d5
msol = 1.989d33 ;solar mass in g
kpc_unit = 1
kpc = 3.0857d21 ;kpc in cm
mH = 1.673534d-24
molec_mass = 1.72000 ;average mass of an atom in the universe (0.76 + 0.24*4)
grav = 6.67259d-8 ;in cgs
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
loadct,39
dens_conv = msol*msol_unit/(kpc_unit*kpc)^3/mH*molec_mass ;convers system density into atoms per cm^3

IF NOT KEYWORD_SET(nbins) THEN nbins = 20.
nfiles = N_ELEMENTS(files)

FOR i = 0, nfiles-1 DO BEGIN
    rtipsy, files[i], h, g, d, s
    g.dens = g.dens * dens_conv
    nbins = 140

                                ;density
    minT =  3
    maxT =  7
    attributeT = ALOG10(g.tempg)
    cutoffT = alog10(15000)
    xT = findgen(nbins)/nbins*(maxT-minT)+minT
    yT = (histogram(attributeT,max=maxT,min=minT,nbins=nbins))/FLOAT(N_ELEMENTS(attributeT))
    if ((where(attributeT lt cutoffT))[0] eq -1) THEN validT[i] = 0 else validT[i] = TOTAL(g[where(attributeT lt cutoffT)].mass)/Total(g.mass)
                                ;temperature
    minD = -8
    maxD =  6
    attributeD = alog10(g.dens)
    cutoffD = -1.0
    xD = findgen(nbins)/nbins*(maxD-minD)+minD
    yD = (histogram(attributeD,max=maxD,min=minD,nbins=nbins))/FLOAT(N_ELEMENTS(attributeD))
    if ((where(attributeD gt cutoffD))[0] eq -1) THEN BEGIN
        validD[i] = 0 
        tdyn[i] = 0
    endif ELSE begin
        validD[i] = TOTAL(g[where(attributeD gt cutoffD)].mass)/TOTAL(g.mass)
        tdyn[i] = 1.0/SQRT(4.0*!PI*grav*MEAN(g[where(attributeD gt cutoffD)].dens)*mH/molec_mass)       
    ENDELSE

    if ((where(attributeT lt cutoffT AND attributeD gt cutoffD))[0] eq -1) THEN valid[i] = 0 ELSE valid[i] = TOTAL(g[where(attributeT lt cutoffT AND attributeD gt cutoffD)].mass)/TOTAL(g.mass)
    IF (type eq 0) THEN BEGIN
        x = xD
        y = yD
        min = minD
        max = maxD
        cutoff = cutoffD
     ENDIF ELSE BEGIN
        x = xT
        y = yT
        min = minT
        max = maxT
        cutoff = cutoffT
    ENDELSE
    IF KEYWORD_SET(OUTPLOT) THEN BEGIN
        IF (i eq 0) THEN BEGIN
            plot,x,y,title = title, xtitle = xtitle, ytitle = ytitle ,linestyle = linestyle[i],yrange = yrange,color = colors[i],xrange = [min,max],thick= thickness[i]
            oplot,[cutoff,cutoff],yrange
        ENDIF ELSE  oplot,x,y,linestyle = linestyle[i], color = colors[i],thick= thickness[i]
    ENDIF
    print,files[i],linestyle[i],thickness[i],tdyn[i],validD[i],validT[i],valid[i]
ENDFOR
IF (KEYWORD_SET(legend_name)) THEN legend,legend_name,color = colors,linestyle = linestyle,thick=thickness,/top,/right,charsize = 1.05
END

PRO master_distribution
cd,'/astro/net/scratch1/christensen/DwarfResearch/procedures'
;This master plotter will loop through interesting files

;!P.MULTI = [0, 2, 2, 0, 0] 
!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1
!Y.Charsize = 1
!X.style = 1
times = ['050','100','150','200','250','300']
FOR it=0,N_ELEMENTS(times) - 1 DO BEGIN
!P.MULTI = [0, 2, 2, 0, 0]
base = ['../ResMassTests/','../ResSpecTests/1E3R/']
imf = ['k']
type = ['rho','temp']
type_str = ['Log(Density)','Log(Temperature)']
legend_name0 = [textoidl('10^6'),textoidl('10^5'),textoidl('10^4'),textoidl('10^3')]
legend_name0 = [textoidl('10^5'),textoidl('10^4'),textoidl('10^3'),textoidl('10^2')]
legend_name1 = [textoidl('2.5\times10^{-4}R_{vir}'),textoidl('10^{-3}R_{vir}'),textoidl('10^{-2}R_{vir}'),textoidl('4\times10^{-2}R_{vir}')]
res_str0 = ['1E6R','1E5R','1E4R','1E3R']
res_str0 = ['1E5R','1E4R','1E3R','1E2R']
res_str1 = ['50pc','200pc','2kpc','8kpc']
res0 = [1e6,1e5,1e4,1e3]
res0 = [1e5,1e4,1e3,1e2]
res1 = [2.5e-4,1e-3,1e-2,4e-2]
;xtick_specres = [textoidl('')]

;colors0 = [70,70,70,70,70,70]
;linestyle0 = [0,0,2,2,3,3]
colors0 = [70,70,70,70]
linestyle0 = [0,0,2,2]
colors1 = [70,70,70,70]
linestyle1 = [0,0,2,2]
;thickness0 = [2.5,1,2.5,1,2.5,1]
thickness0 = [2.5,1,2.5,1]
thickness1 = [2.5,1,2.5,1]
yrange0 = [[0,0.06],[0,0.1]]
yrange1 = [[0,0.06],[0,0.3]]
mass_str = ['12','10']
tdyn = fltarr(N_ELEMENTS(res_str0))
tdyn_array = fltarr(N_ELEMENTS(mass_str),N_ELEMENTS(base),N_ELEMENTS(res_str0))

valid_array = fltarr(N_ELEMENTS(mass_str),N_ELEMENTS(base),N_ELEMENTS(type)+1,N_ELEMENTS(res_str0))
valid = fltarr(N_ELEMENTS(res_str0))
validT = fltarr(N_ELEMENTS(res_str0))
validD = fltarr(N_ELEMENTS(res_str0))

;set_plot,'ps'
set_plot,'x'
window,1
FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
    FOR tct = 0, N_ELEMENTS(type) - 1 DO BEGIN
;        device,filename = '../results/distribution_'+type[tct]+'_'+imf[imfct]+'.eps',/color,bits_per_pixel = 8,/times
        FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
            FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN
                IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00'+times[it] ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00'+times[it]
                xtitle = type_str[tct]
                ytitle = textoidl('Distribution')
                multiplot
                IF (mct eq 0) THEN BEGIN 
                     IF (ploti eq 0) THEN distribution, files, tct, validT,validD,valid,tdyn,linestyle = linestyle0,/outplot, title = title, xtitle = '',ytitle = ytitle,colors = colors0,yrange = yrange0[*,tct],thickness=thickness0 ELSE distribution, files, tct, validT,validD,valid,tdyn,linestyle = linestyle1,/outplot, title = title,colors = colors1,yrange = yrange0[*,tct], legend_name = legend_name1,thickness=thickness1
                ENDIF ELSE BEGIN
                    IF (ploti eq 0) THEN distribution, files, tct, validT,validD,valid,tdyn,linestyle = linestyle0,/outplot,title = title, xtitle = xtitle, ytitle = ytitle, colors = colors0,yrange = yrange1[*,tct],legend_name = legend_name0,thickness=thickness0 ELSE distribution, files, tct, validT,validD,valid,tdyn,linestyle = linestyle1,/outplot, title = title, xtitle = xtitle, colors = colors1, yrange = yrange1[*,tct],thickness=thickness1
                ENDELSE
                valid_array[mct,ploti,2,*] = validD
                valid_array[mct,ploti,1,*] = validT
                valid_array[mct,ploti,0,*] = valid
                tdyn_array[mct,ploti,*] = tdyn
            ENDFOR
        ENDFOR
        xyouts,[3200,10200,3200,10200],[11000,11000,6000,6000],[textoidl('10^{12}M')+sunsymbol(),textoidl('10^{12}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol()],charsize = 1.05,/device
;        device,/close
        multiplot,/reset
        if (tct eq 0) THEN window,2 ELSE window,3
    ENDFOR

;    mass_str = ['13','12','11','10','9']
;    ploti = 0
;    tct = 0
;    imfct = 0
;    !P.MULTI = [0] 
;    device,filename = '../results/validity_'+imf[imfct]+'.eps',/color,bits_per_pixel = 8,/times
;    FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
;        IF (mct eq 1 OR mct eq 3) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00100'  else files = base[ploti] + res_str0[1:N_ELEMENTS(res_str0)-1] + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00100' 
;        distribution, files, tct, validT, validD, valid,linestyle = linestyle0,outplot = outplot, title = title, xtitle = '',ytitle = ytitle,colors = colors0,yrange = yrange0[*,tct],thickness=thickness0
;        IF(mct eq 1 OR mct eq 3 ) THEN BEGIN
;            plot,res0,valid,yrange = [0,1],ytitle='Fraction of Gas',xtitle = 'Number of Particles',/xlog ,xrange = [500,2e6]
;            oplot,res0,validT,linestyle = 2
;            oplot,res0,validD,linestyle = 3
;            print,mass_str[mct]
;            print,'Resolution   ',res_str0
;            print,'Total:       ',valid
;            print,'Temperature: ',validT
;            print,'Density:     ',validD
;        ENDIF ELSE BEGIN
;            plot,res0[1:N_ELEMENTS(res0)-1],valid,ytitle='Fraction of Gas',xtitle = 'Number of Particles',yrange = [0,1],/xlog,xrange = [500,2e6]
;            oplot,res0[1:N_ELEMENTS(res0)-1],validT,linestyle = 2
;            oplot,res0[1:N_ELEMENTS(res0)-1],validD,linestyle = 3
;            print,mass_str[mct]
;            print,'Resolution   ',res_str0[1:N_ELEMENTS(res0)-1]
;            print,'Total:       ',valid[0:N_ELEMENTS(res0)-2]
;            print,'Temperature: ',validT[0:N_ELEMENTS(res0)-2]
;            print,'Density:     ',validD[0:N_ELEMENTS(res0)-2]
;        ENDELSE
;        legend,[textoidl('T < T_{max} and'),textoidl('\rho > \rho_{min}'),textoidl('T < T_{max}'),textoidl('\rho > \rho_{min}')],linestyle = [0,0,2,3],/top,/right,charsize = 1.05
;        stop
;   ENDFOR

!P.MULTI = [0, 2, 2, 0, 0] 
mass_str = ['12','10']
;    device,filename = '../results/validity_'+imf[imfct]+'.eps',/color,bits_per_pixel = 8,/times
       FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
            FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN
               multiplot
                IF (ploti eq 0) THEN BEGIN
                    if(mct eq 0) THEN BEGIN
                        plot,res0,valid_array[mct,ploti,0,*],yrange = [0,1],ytitle='Fraction of Gas',/xlog ,xrange = [500,2e6]
                        legend,[textoidl('T < T_{max} and'),textoidl('\rho > \rho_{min}'),textoidl('T < T_{max}'),textoidl('\rho > \rho_{min}')],linestyle = [0,0,2,3],/top,/right,charsize = 1.05,color=[250,60,60,60]
                        ENDIF ELSE plot,res0,valid_array[mct,ploti,0,*],ytitle='Fraction of Gas',xtitle = 'Number of Particles',yrange = [0,1],/xlog,xrange = [10,2e6]
                    oplot,res0,valid_array[mct,ploti,1,*],linestyle = 2
                    oplot,res0,valid_array[mct,ploti,2,*],linestyle = 3

                    print,''
                    print,mass_str[mct]
                    print,'Resolution   ',res_str0
                    print,'Total:       ',valid_array[mct,ploti,0,0],valid_array[mct,ploti,0,1],valid_array[mct,ploti,0,2],valid_array[mct,ploti,0,3]
                    print,'Temperature: ',valid_array[mct,ploti,1,0],valid_array[mct,ploti,1,1],valid_array[mct,ploti,1,2],valid_array[mct,ploti,1,3]
                    print,'Density:     ',valid_array[mct,ploti,2,0],valid_array[mct,ploti,2,1],valid_array[mct,ploti,2,2],valid_array[mct,ploti,2,3]
                    print,'Dynamical Time',tdyn_array[mct,ploti,0],tdyn_array[mct,ploti,1],tdyn_array[mct,ploti,2],tdyn_array[mct,ploti,3]
                ENDIF ELSE BEGIN
                    if (mct eq 0) THEN plot,res1,valid_array[mct,ploti,0,*],yrange=[0,1],/xlog,xrange=[1e-1,1e-4] ELSE plot,res1,valid_array[mct,ploti,0,*],yrange=[0,1],xtitle = textoidl('Softening Length [R_{vir}]'),/xlog,xrange = [1e-1,1e-4],xtickname=[textoidl('10^{-1}'),textoidl('10^{-2}'),textoidl('10^{-3}'),textoidl('10^{-4}')]
                    oplot,res1,valid_array[mct,ploti,1,*],linestyle = 2
                    oplot,res1,valid_array[mct,ploti,2,*],linestyle = 3
                    print,''
                    print,mass_str[mct]
                    print,'Resolution   ',res_str1
                    print,'Total:       ',valid_array[mct,ploti,0,0],valid_array[mct,ploti,0,1],valid_array[mct,ploti,0,2],valid_array[mct,ploti,0,3]
                    print,'Temperature: ',valid_array[mct,ploti,1,0],valid_array[mct,ploti,1,1],valid_array[mct,ploti,1,2],valid_array[mct,ploti,1,3]
                    print,'Density:     ',valid_array[mct,ploti,2,0],valid_array[mct,ploti,2,1],valid_array[mct,ploti,2,2],valid_array[mct,ploti,2,3]
                    print,'Dynamical Time',tdyn_array[mct,ploti,0],tdyn_array[mct,ploti,1],tdyn_array[mct,ploti,2],tdyn_array[mct,ploti,3]
                ENDELSE
            ENDFOR
        ENDFOR
        xyouts,[3200,10200,3200,10200],[10500,10500,5500,5500],[textoidl('10^{12}M')+sunsymbol(),textoidl('10^{12}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol()],charsize = 1.05,/device
;        device,/close
ENDFOR
stop
ENDFOR
cd,'/astro/users/christensen/code/Dwarfs'
END
