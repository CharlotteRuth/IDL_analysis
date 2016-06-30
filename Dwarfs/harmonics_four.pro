;/astro/users/christensen/code/Dwarfs/FFT_SH/s2kit10/call_semi_fly  12M.binrho fftrho.dat  32

PRO HARMONICS_INPUT, files, output, type, NBINS = NBINS, MAXRADII = maxradii, OUTPLOT = OUTPLOT, TITLE = TITLE, XTITLE = XTITLE, YTITLE = YTITLE, LEGEND_NAME = LEGEND_NAME, YRANGE = YRANGE, COLORS = COLORS, PRINT = PRINT, LINESTYLE = LINESTYLE,THICKNESS = THICKNESS
;This will plot half mass and ninety-percent mass heights of stellar
;or gas disks
msol = 2.362e5
kpc_unit = 1
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
loadct,39
type = ['rho','temp']
type_str = ['Log(Density)','Log(Temperature)']

IF NOT KEYWORD_SET(nbins) THEN nbins = 32
nbins = nbins*2.0
nfiles = N_ELEMENTS(files)

binned0 = fltarr(nbins,nbins)
binned1 = fltarr(nbins,nbins)
db_theta = !PI/nbins
db_phi = 2.0*!PI/nbins

FOR i = 0, nfiles-1 DO BEGIN
    rtipsy, files[i], h, g, d, s
 ;   r = SQRT+ g.y*g.y + g.z*g.z)
    ;Theta goes from 0 to pi, angle off of vertical
 ;   theta = ACOS(g.z/r)
    ;phi goes from 0 to 2pi, angle around vertical
 ;   phi = ATAN(g.y/g.x)
    rect_cord = TRANSPOSE([[g.x],[g.y],[g.z]])
    sph_cord = CV_COORD(FROM_RECT = rect_cord,/TO_SPHERE)
    phi = REFORM(sph_cord[0,*]) +!PI
    theta = REFORM(sph_cord[1,*])+!PI/2.0
    openw,1,output[i]+type[0]
    openw,2,output[i]+type[1]

    variable0 = ALOG10(g.dens)
    min0 = -5.0
    variable1 = ALOG10(g.tempg)
    min1 = 2
 
    FOR theta_ct = 0, (nbins)-1 DO BEGIN
        FOR phi_ct = 0, (nbins)-1 DO BEGIN
            ind = where((theta ge theta_ct*db_theta AND theta lt (theta_ct+1.0)*db_theta) AND (phi ge phi_ct*db_phi AND phi lt (phi_ct+1.0)*db_phi))
            IF (ind[0] ne -1) THEN BEGIN
                    binned0[theta_ct,phi_ct] = MEAN(variable0[ind])
                    binned1[theta_ct,phi_ct] = MEAN(variable1[ind])
                ENDIF ELSE BEGIN
                    binned0[theta_ct,phi_ct] = min0
                    binned1[theta_ct,phi_ct] = min1
                ENDELSE
            printf,1,binned0[theta_ct,phi_ct]
            printf,2,binned1[theta_ct,phi_ct]
        ENDFOR
    ENDFOR
    close,1
    close,2
    print,files[i]
    print,output[i]
    print,binned0
    print,binned1
    print,""
ENDFOR
END

PRO HARMONICS, files, harm, NBINS = NBINS, MAXRADII = maxradii, OUTPLOT = OUTPLOT, TITLE = TITLE, XTITLE = XTITLE, YTITLE = YTITLE, LEGEND_NAME = LEGEND_NAME, YRANGE = YRANGE, COLORS = COLORS, PRINT = PRINT, LINESTYLE = LINESTYLE,THICKNESS = THICKNESS
;This will plot half mass and ninety-percent mass heights of stellar
;or gas disks
msol = 2.362e5
kpc_unit = 1
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
loadct,39

IF NOT KEYWORD_SET(nbins) THEN nbins = 32.
nfiles = N_ELEMENTS(files)

FOR i = 0, nfiles-1 DO BEGIN
    rtipsy, files[i], h, g, d, s
    r = SQRT(g.x*g.x + g.y*g.y + g.z*g.z)
    theta = ACOS(g.y/g.x)
    phi = ATAN(g.z/g.r)
    harm[i] = 0
  ENDFOR
END

PRO MASTER_HARMONICS
cd,'/astro/net/scratch1/christensen/DwarfResearch/procedures'
;This master plotter will loop through interesting files

!P.MULTI = [0, 2, 2, 0, 0] 
!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1.25
!Y.Charsize = 1.25
!X.style = 1

base = ['../ResMassTests/','../ResSpecTests/1E5R/']
imf = ['k']
legend_name0 = [textoidl('10^6'),textoidl('10^5'),textoidl('10^4'),textoidl('10^3'),textoidl('100'),textoidl('50')]
legend_name1 = [textoidl('2.5\times10^{-4}R_{vir}'),textoidl('10^{-3}R_{vir}'),textoidl('10^{-2}R_{vir}'),textoidl('4\times10^{-2}R_{vir}')]
res_str0 = ['1E6R','1E5R','1E4R','1E3R']
res_str1 = ['50pc','200pc','2kpc','8kpc']
res0 = [1e6,1e5,1e4,1e3,100,50]
res1 = [2.5e-4,1e-3,1e-2,4e-2]

colors0 = [0,0,0,0,0,0]
linestyle0 = [0,0,2,2,3,3]
colors1 = [0,0,0,0]
linestyle1 = [0,0,2,2]
thickness0 = [2.5,1,2.5,1,2.5,1]
thickness1 = [2.5,1,2.5,1]
yrange0 = [[0,0.06],[0,0.1]]
yrange1 = [[0,0.06],[0,0.3]]
mass_str = ['12','10']
bin=32.0

harm_array = fltarr(N_ELEMENTS(mass_str),N_ELEMENTS(base),N_ELEMENTS(res_str0))
harm = fltarr(N_ELEMENTS(res_str0))

set_plot,'x'
FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
        FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
            FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN
                IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00300' ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
                IF (ploti eq 0) THEN output = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/'+mass_str[mct]+'M'+'.bin' ELSE output = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.bin'
                harmonics_input, files, output, tct, nbins=bin
        ENDFOR
    ENDFOR
ENDFOR
cd,'/astro/users/christensen/code/Dwarfs'
END

PRO HARMONICS_PLOT
cd,'/astro/net/scratch1/christensen/DwarfResearch/procedures'
;This master plotter will loop through interesting files

!P.MULTI = [0, 2, 2, 0, 0] 
!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1.25
!Y.Charsize = 1.25
!X.style = 1

base = ['/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/','/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/']
imf = ['k']
type = ['temp','rho']
type_str = ['Log(Density)','Log(Temperature)']
legend_name0 = [textoidl('10^6'),textoidl('10^5'),textoidl('10^4'),textoidl('10^3'),textoidl('100'),textoidl('50')]
legend_name1 = [textoidl('2.5\times10^{-4}R_{vir}'),textoidl('10^{-3}R_{vir}'),textoidl('10^{-2}R_{vir}'),textoidl('4\times10^{-2}R_{vir}')]
res_str0 = ['1E6R','1E5R','1E4R','1E3R']
res_str1 = ['50pc','200pc','2kpc','8kpc']
res0 = [1e6,1e5,1e4,1e3]
res1 = [2.5e-4,1e-3,1e-2,4e-2]

colors0 = [0,0,0,0,0,0]
linestyle0 = [0,0,2,2,3,3]
colors1 = [0,0,0,0]
linestyle1 = [0,0,2,2]
thickness0 = [2.5,1,2.5,1,2.5,1]
thickness1 = [2.5,1,2.5,1]
yrange = [[[0,20],[0,17]],[[-16,0],[-18,0]]]
mass_str = ['12','10']
bin=32.0
n_comp = 6

harm_array = fltarr(N_ELEMENTS(mass_str),N_ELEMENTS(base),N_ELEMENTS(res0),n_comp)

set_plot,'x'
;    device,filename =
;    '../results/spherical_harmonics_'+imf[imfct]+'.eps',/color,bits_per_pixel = 8,/times
FOR tct = 0, N_ELEMENTS(type) - 1 DO BEGIN
!P.MULTI = [0, 2, 2, 0, 0] 
       FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
            FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN
                IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_k'+'/'+mass_str[mct]+'M'+'.fft'+type[tct] ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_k/'+res_str1+'.fft'+type[tct]          
                FOR rct = 0, N_ELEMENTS(res0) - 1 DO BEGIN
                 ;   openr,1,files[rct]
                    readcol,files[rct],comp
                    harm_array[mct,ploti,rct,*] = comp[0:n_comp-1]
                    print,where(comp eq MAX(comp[1:N_ELEMENTS(comp)-1]))
                  ;  close,1
                ENDFOR
                yr = [-1,1];reform(yrange[*,mct,tct])
                multiplot
                if ploti eq 0 THEN plot, res0,harm_array[mct,ploti,*,0],/xlog,yrange=yr ELSE plot,res1,harm_array[mct,ploti,*,0],/xlog,yrange=yr
;                FOR compct = 1, n_comp-1 DO  if ploti eq 0 THEN oplot, res0,harm_array[mct,ploti,*,compct],linestyle = compct+1 ELSE oplot,res1,harm_array[mct,ploti,*,compct],linestyle = compct+1
                if ploti eq 0 THEN oplot, res0,harm_array[mct,ploti,*,3],linestyle = 2 ELSE oplot,res1,harm_array[mct,ploti,*,3],linestyle = 2
                print,MINMAX(harm_array[mct,ploti,*,*])
                stop
            ENDFOR
        ENDFOR
    ENDFOR
;device,/close
multiplot, /reset
multiplot, /default
!p.multi = 0
cd,'/astro/users/christensen/code/Dwarfs'
END
