;/astro/users/christensen/code/Dwarfs/FFT_SH/s2kit10/call_semi_fly  12M.binrho fftrho.dat  32

PRO HARMONICS_INPUT, files, output, NBINS = NBINS, MAXRADII = maxradii, OUTPLOT = OUTPLOT, TITLE = TITLE, XTITLE = XTITLE, YTITLE = YTITLE, LEGEND_NAME = LEGEND_NAME, YRANGE = YRANGE, COLORS = COLORS, PRINT = PRINT, LINESTYLE = LINESTYLE,THICKNESS = THICKNESS
;This will plot half mass and ninety-percent mass heights of stellar
;or gas disks
msol = 2.362e5
kpc_unit = 1
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
type = ['rho','temp','mass']
loadct,39
IF NOT KEYWORD_SET(nbins) THEN nbins = 16
size = nbins*2.0
nfiles = N_ELEMENTS(files)

binned0 = fltarr(size,size)
binned1 = fltarr(size,size)
binned2 = fltarr(size,size)
db_theta = !PI/size
db_phi = 2.0*!PI/size
close,/all
levels = 75
colors0 = 2.0*(findgen(levels)+1.0)
colors1 = 244.0 - 2.0*findgen(levels)
colors2 = 2.0*(findgen(levels)+1.0)+75.0

db_theta = !PI/size
theta_array = findgen(size + 1)*db_theta
area = COS(theta_array - !PI/2.0)

FOR i = 0, nfiles-1 DO BEGIN
    rtipsy, files[i], h, g, d, s
    disk_ind = WHERE (g.x*g.x + g.y*g.y lt 121 AND ABS(g.z) lt 8)
    if (disk_ind[0] eq -1) THEN disk_ind = fltarr(N_ELEMENTS(g.x) - 1)
    xc = TOTAL(g[disk_ind].x*g[disk_ind].mass)/TOTAL(g[disk_ind].mass)
    yc = TOTAL(g[disk_ind].y*g[disk_ind].mass)/TOTAL(g[disk_ind].mass)
    zc = TOTAL(g[disk_ind].z*g[disk_ind].mass)/TOTAL(g[disk_ind].mass)
;    print,xc,yc,zc
    rect_cord = TRANSPOSE([[g.x - xc],[g.y - yc],[g.z - zc]])
    sph_cord = CV_COORD(FROM_RECT = rect_cord,/TO_SPHERE)
    phi = REFORM(sph_cord[0,*]) +!PI
    theta = REFORM(sph_cord[1,*])+!PI/2.0
    openw,1,output[i]+type[0]
    openw,2,output[i]+type[1]
    openw,3,output[i]+type[2]
    halo_ind = where(kpc_unit*sph_cord[2,*] gt 0)
    variable0 = ALOG10(g[halo_ind].dens)
    min0 = -5.0
    max0 = 2
    variable1 = ALOG10(g[halo_ind].tempg)
    min1 = 2
    max1 = 7.2
    variable2 = g[halo_ind].mass
    min2 = 0
    
    empty = fltarr(2,size*size)
    ct_empty = 0
    FOR theta_ct = 0, (size)-1 DO BEGIN
        FOR phi_ct = 0, (size)-1 DO BEGIN
            ind = where((theta ge theta_ct*db_theta AND theta lt (theta_ct+1.0)*db_theta) AND (phi ge phi_ct*db_phi AND phi lt (phi_ct+1.0)*db_phi))
            IF (ind[0] ne -1) THEN BEGIN
                binned0[phi_ct,theta_ct] = MEAN(variable0[ind])
                binned1[phi_ct,theta_ct] = MEAN(variable1[ind])
                binned2[phi_ct,theta_ct] = 2.0*TOTAL(variable2[ind])/(area[theta_ct] + area[theta_ct+1])
            ENDIF ELSE BEGIN
                binned0[phi_ct,theta_ct] = min0
                binned1[phi_ct,theta_ct] = -1
                binned2[phi_ct,theta_ct] = -1
            ENDELSE
        ENDFOR
    ENDFOR

    ind = where(binned1 eq -1)
    if ind[0] ne -1 then begin
        min0 = min(binned0[where(binned2 ne -1)])
        binned0[ind] = min0
        mean1 = MEAN(binned1[where(binned2 ne -1)])
        binned1[ind] = mean1
        min2 = min(binned2[where(binned2 ne -1)])
        binned2[ind] = min2
    endif    

    FOR theta_ct = 0, (size)-1 Do BEGIN
        FOR phi_ct = 0, (size)-1 DO BEGIN
            printf,1,binned0[phi_ct,theta_ct]
            printf,2,binned1[phi_ct,theta_ct]
            printf,3,binned2[phi_ct,theta_ct]
        ENDFOR
    ENDFOR
    close,1
    close,2
    close,3
    print,files[i]

;    help,binned0
;    help,where(binned0 eq min0)     
    contour,binned0,findgen(size)*db_phi,findgen(size)*db_theta,ystyle=1,xstyle=1,yrange=[0,!PI],xrange=[0,2.0*!PI],ytitle="Theta",xtitle="Phi",/FILL,C_COLOR=colors0,title='Density '+files[i],min_value=min0,max_value=max0,nlevels = levels
    wait,0.5
;    stop
    contour,binned1,findgen(size)*db_phi,findgen(size)*db_theta,ystyle=1,xstyle=1,yrange=[0,!PI],xrange=[0,2.0*!PI],ytitle="Theta",xtitle="Phi",/FILL,C_COLOR=colors1,title='Temperature '+files[i],min_value=min1,max_value=max1,nlevels = levels
    wait,0.5
;    stop
    contour,binned2,findgen(size)*db_phi,findgen(size)*db_theta,ystyle=1,xstyle=1,yrange=[0,!PI],xrange=[0,2.0*!PI],ytitle="Theta",xtitle="Phi",/FILL,C_COLOR=colors2,title='Mass '+files[i],min_value=min2,nlevels = levels
    wait,0.5
;    stop
ENDFOR
END

PRO HARMONICS_INPUT_MASTER
cd,'/astro/net/scratch1/christensen/DwarfResearch/procedures'
;This master program will loop through interesting files and prepare
;them to be a spherical harmonic FFT

!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1.25
!Y.Charsize = 1.25
!X.style = 1

base = ['/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/','/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/']
imf = ['k']
res_str0 = ['1E6R','1E5R','1E4R','1E3R']
res_str1 = ['50pc','200pc','500pc','1kpc','2kpc','8kpc']
iter = [['13M','12M','11M','10M','9M'],['1E5R','1E4R','1E3R','0','0']]
mass = ['13','12','10','11','9']
mass_str = ['12','10']
nbins=16.0
res_str0 = ['1E6R','1E5R','1E4R','1E3R']
res_str1 = ['50pc','200pc','500pc','1kpc','2kpc','8kpc']
res0 = [1e6,1e5,1e4,1e3]
res1 = [2e-2,1e-2,5e-3,2.5e-3,1e-3,2.5e-4]

harm_array = fltarr(N_ELEMENTS(mass_str),N_ELEMENTS(base),N_ELEMENTS(res_str0))
harm = fltarr(N_ELEMENTS(res_str0))

set_plot,'x'
FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
    FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
        FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN
            IF ploti eq 0 THEN nlines = N_ELEMENTS(mass) ELSE nlines = N_ELEMENTS(res0)
            if ploti eq 0 then endloop = N_ELEMENTS(mass) ELSE endloop = 3
            FOR i = 0, endloop - 1 DO BEGIN
                if (ploti eq 0) THEN BEGIN
                    if (i ne 1 AND i ne 3) then res_str0_select = res_str0[1:N_ELEMENTS(res_str0)-1] else res_str0_select = res_str0 
                    files = base[ploti] + res_str0_select + '/' +iter[i,0]+'_'+imf[imfct]+'/o'+iter[i,0]+'_1.00300' 
                    output = base[ploti] + res_str0_select + '/' +iter[i,0]+'_'+imf[imfct]+'/'+iter[i,0]+'.bin' 
                ENDIF ELSE BEGIN
                    files = base[ploti] +iter[i,1]+ '/'+mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
                    output = base[ploti] +iter[i,1]+ '/'+mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.bin'
                ENDELSE
                harmonics_input, files, output, nbins=nbins
            ENDFOR
        ENDFOR
    ENDFOR
ENDFOR
cd,'/astro/users/christensen/code/Dwarfs'
END

PRO HARMONICS_PLOT
;This master plotter will loop through interesting files

;!P.MULTI = [0, 2, 2, 0, 0] 
!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1.25
!Y.Charsize = 1.25
!X.style = 1
kpc_unit = 1
loadct,39

base = ['/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/','/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/']
outplot = [['/astro/net/scratch1/christensen/DwarfResearch/results/SH_massRes_','/astro/net/scratch1/christensen/DwarfResearch/results/SH_massRes_'],['/astro/net/scratch1/christensen/DwarfResearch/results/SH_spec12_','/astro/net/scratch1/christensen/DwarfResearch/results/SH_spec10_']]
imf = ['k']
type = ['rho','temp','mass']
type_str = ['Log(Density)','Log(Temperature)','Mass']
res_str0 = ['1E3R', '1E4R', '1E5R','1E6R']
res_str1 = ['50pc','200pc','500pc','1kpc','2kpc','8kpc']
res0 = [1e3,1e4,1e5,1e6]
res1 = [2.5e-4,1e-3,2.5e-3,5e-3,1e-2,4e-2]
iter = [['13M','12M','11M','10M','9M'],['1E3R','1E4R','1E5R','0','0']]
mass = ['13','12','10','11','9']
mass_str = ['12','10']
symbols = [-2,-4,-5,-6,-7]
linestyle = [0,0,2,2,3]
thickness = [3,1,3,1,3]
xtick = [[textoidl('10^{3}'),textoidl('10^{4}'),textoidl('10^{5}'),textoidl('10^{6}')],[textoidl('10^{-1}'),textoidl('10^{-2}'),textoidl('10^{-3}'),textoidl('10^{-4}')]]
yrange = [[[0,20],[0,17]],[[-16,0],[-18,0]]]
xrange = [[500,2e6],[0.1,0.0001]]
xtitle = ['Number of Particles','Softening Length ' + textoidl('[R_{vir}]')]
ytitle = textoidl('Component')
legend_name0 = [textoidl('10^{13}')+'M'+sunsymbol(),textoidl('10^{12}')+'M'+sunsymbol(),textoidl('10^{11}')+'M'+sunsymbol(),textoidl('10^{10}')+'M'+sunsymbol(),textoidl('10^{9}')+'M'+sunsymbol()]
legend_name1 = [textoidl('10^3'),textoidl('10^4'),textoidl('10^5')]
;legend_name1 = [textoidl('2.5\times10^{-4}R_{vir}'),textoidl('10^{-3}R_{vir}'),textoidl('10^{-2}R_{vir}'),textoidl('4\times10^{-2}R_{vir}')]
mass_str = ['12','10']
nbins=16.0
size = nbins*2
n_comp = 6
comp_plot = 4

db_theta = !PI/size
db_phi = 2.0*!PI/size
phi_array= findgen(size)*db_phi
theta_array = findgen(size + 1)*db_theta
harm_array = fltarr(N_ELEMENTS(mass),N_ELEMENTS(res1),n_comp)
area = COS(theta_array - !PI/2.0)
binned = fltarr(size,size)
binned_i = fltarr(size,size)
levels = 75.0
min = [-5.0,2,5e-10]
max = [1,10,1e-9]
colors = [[2.0*(findgen(levels)+1.0)],[244.0 - 2.0*findgen(levels)],[2.0*(findgen(levels)+1.0)+75.0]]

;comps = [0,1,2,3,4,8,16,241,242,228]
lm_pairs = fltarr(nbins*nbins*2.0,2)
for im = -1.0*(nbins-1), nbins-1 DO BEGIN
    FOR il = abs(im),nbins-1 DO BEGIN
        IF (im ge 0) THEN ind = im*nbins - ((im*(im - 1))/2.0) + (il - im) ELSE ind = ((nbins - 1)*(nbins + 2)/2.0) + 1.0 + (((nbins - 1.0) + im)*(nbins + im)/2.0) + (il - abs(im))
        lm_pairs[ind,*] = [il, im]
;        if ((where(ind eq comps))[0] ne -1) THEN print,ind,'L: ',il,'M: ',im
    ENDFOR
ENDFOR

set_plot,'x'
; set_plot,'ps'
;    device,filename =
;    '../results/spherical_harmonics_'+imf[imfct]+'.eps',/color,bits_per_pixel = 8,/times

FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN                       ;Loops through IMF 
    FOR tct = 0, N_ELEMENTS(type) - 1 DO BEGIN                    ;Loops through temp and density     
        FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN ;Loops twice (the mass resolution is repeated)
            FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN           ;Loops for mass and spatial resolution
                outfile = outplot[mct,ploti]+imf[imfct]+'_'+type[tct]+'.eps'
          ;      device,filename=outfile,/color,bits_per_pixel= 8,/times
                IF ploti eq 0 THEN nlines = N_ELEMENTS(mass) ELSE nlines = N_ELEMENTS(res0)-1
                FOR linesct = 0, nlines - 1 DO BEGIN              ;Loops through the lines for each plot
                    IF (ploti eq 0) THEN BEGIN                    ;Mass Resolution
                        files = base[ploti] + res_str0 + '/' +iter[linesct,0]+'_'+imf[imfct]+'/'+iter[linesct,0]+'.fft'+type[tct] 
                        tipsy = base[ploti] + res_str0 + '/' +iter[linesct,0]+'_'+imf[imfct]+'/o'+iter[linesct,0]+'_1.00300' 
                        IF (linesct eq 1 OR linesct EQ 3) THEN endloop = N_ELEMENTS(res0) ELSE endloop = N_ELEMENTS(res0) - 1
                        xaxes = res0[0:endloop-1]
                    ENDIF ELSE BEGIN                              ;Spatial Resolution
                        files = base[ploti] +iter[linesct,1]+ '/'+mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.fft'+type[tct]
                        tipsy = base[ploti] +iter[linesct,1]+ '/'+mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
                        endloop = N_ELEMENTS(res1)
                        xaxes = res1[0:endloop-1]
                    ENDELSE
                   
                    FOR i = 0, endloop - 1 DO BEGIN     ;Loops through x-axes (masses for mass resolution and particle #'s for spatical resolution)
                        print,tipsy[i]
                        rtipsy, tipsy[i], h, g, d, s
                        disk_ind = WHERE (g.x*g.x + g.y*g.y lt 121 AND ABS(g.z) lt 8)
                        if (disk_ind[0] eq -1) THEN disk_ind = fltarr(N_ELEMENTS(g.x) - 1)
                        xc = TOTAL(g[disk_ind].x*g[disk_ind].mass)/TOTAL(g[disk_ind].mass)
                        yc = TOTAL(g[disk_ind].y*g[disk_ind].mass)/TOTAL(g[disk_ind].mass)
                        zc = TOTAL(g[disk_ind].z*g[disk_ind].mass)/TOTAL(g[disk_ind].mass)
                        rect_cord = TRANSPOSE([[g.x - xc],[g.y - yc],[g.z - zc]])                       
                        sph_cord = CV_COORD(FROM_RECT = rect_cord,/TO_SPHERE)
                        phi = REFORM(sph_cord[0,*]) +!PI
                        theta = REFORM(sph_cord[1,*])+!PI/2.0
                        halo_ind = where(kpc_unit*sph_cord[2,*] gt 0)
                        if (tct eq 0) then variable = ALOG10(g.dens)
                        if (tct eq 1) then variable = ALOG10(g.tempg)
                        if (tct eq 2) then variable = g[halo_ind].mass 
                        FOR theta_ct = 0, (size)-1 DO BEGIN
                            FOR phi_ct = 0, (size)-1 DO BEGIN
                                ind = where((theta ge theta_ct*db_theta AND theta lt (theta_ct+1.0)*db_theta) AND (phi ge phi_ct*db_phi AND phi lt (phi_ct+1.0)*db_phi))
                                IF (ind[0] ne -1) THEN BEGIN
                                    if tct eq 2 then  binned_i[phi_ct,theta_ct] = 2.0*TOTAL(variable[ind])/(area[theta_ct] + area[theta_ct+1]) else binned_i[phi_ct,theta_ct] = MEAN(variable[ind]) 
                                ENDIF ELSE binned_i[phi_ct,theta_ct] = -1
                            ENDFOR
                        ENDFOR
                        ind = where(binned_i eq -1)
                        if ind[0] ne -1 then begin
                            if (tct eq 1) THEN sub = MEAN(binned_i[where(binned_i ne -1)]) else sub = min(binned_i[where(binned_i ne -1)])
                            binned_i[ind] = sub
                        endif 
                        if (tct eq 2) THEN BEGIN
                            max[tct] = MAX(binned_i) - (MAX(binned_i) - Min(binned_i))/10.0 
                            min[tct] = MIN(binned_i) + (MAX(binned_i) - Min(binned_i))/10.0
                        ENDIF
                            
                        print,files[i]
                        window,0
                        contour,binned_i,findgen(size)*db_phi,findgen(size)*db_theta,ystyle=1,xstyle=1,yrange=[0,!PI],xrange=[0,2.0*!PI],ytitle="Theta",xtitle="Phi",/FILL,C_COLOR=colors[*,tct],title=files[i],nlevel = levels;,min_value=min[tct],max_value=max[tct]
                        readcol,files[i],comp,/silent
                        harm_array[linesct,i,*] = comp[0:n_comp-1]/SQRT(TOTAL(comp*comp))
                        sorted = REVERSE((sort(abs(comp))))
                        print,comp[sorted[0:5]]/SQRT(TOTAL(comp*comp)),sorted[0:5]
                        FOR i_theta = 0, size-1 DO BEGIN
                            binned[*,i_theta] = comp[sorted[0]]*SPHER_HARM(theta_array[i_theta], phi_array, lm_pairs[sorted[0],0], lm_pairs[sorted[0],1])
                            for i_comp = 1, 30 do  binned[*,i_theta] = binned[*,i_theta] +comp[sorted[i_comp]]*SPHER_HARM(theta_array[i_theta], phi_array, lm_pairs[sorted[i_comp],0], lm_pairs[sorted[i_comp],1])
                        ENDFOR
                        window,1                   
                        contour,binned,findgen(size)*db_phi,findgen(size)*db_theta,ystyle=1,xstyle=1,yrange=[0,!PI],xrange=[0,2.0*!PI],ytitle="Theta",xtitle="Phi",/FILL,C_COLOR=colors[*,tct],title=files[i],nlevels = levels;,min_value=min[tct],max_value=max[tct]
                   
                    ENDFOR
            ;        IF (linesct EQ 0) THEN plot,xaxes,harm_array[0,0:endloop-1,comp_plot],/xlog,xrange = xrange[*,ploti],yrange=[-0.1,0.28],ytitle = textoidl('Y^0_4')+' Component',xtitle = xtitle[ploti],xtickname=xtick[*,ploti],thick=thickness[linesct],linestyle=linestyle[linesct],psym = symbols[linesct]  ELSE oplot,xaxes,harm_array[linesct,0:endloop-1,comp_plot],thick=thickness[linesct],linestyle=linestyle[linesct],psym = symbols[linesct]
                ENDFOR
          ;      if ploti eq 0 then legend,legend_name0,linestyle = linestyle,thick=thickness,psym=symbols,/top,/left else BEGIN
          ;          IF mct eq 0 THEN BEGIN 
          ;              legend,legend_name1,linestyle = linestyle[0:2],thick=thickness[0:2],psym=symbols[0:2],/top,/right
          ;              xyouts,[3500],[10500],[textoidl('10^{12}M')+sunsymbol()],/device,charsize = 1.5 
          ;          ENDIF ELSE BEGIN
          ;              legend,legend_name1,linestyle = linestyle[0:2],thick=thickness[0:2],psym=symbols[0:2],/bottom,/right
          ;              xyouts,[3500],[10500],[textoidl('10^{10}M')+sunsymbol()],/device,charsize = 1.5
    ;                ENDELSE
    ;            ENDELSE
       ;         device,/close
                ;reform(yrange[*,mct,tct])
;                multiplot
;                if ploti eq 0 THEN plot, res0,harm_array[mct,ploti,*,0],/xlog,yrange=yr ELSE plot,res1,harm_array[mct,ploti,*,0],/xlog,yrange=yr
;                FOR compct = 1, n_comp-1 DO  if ploti eq 0 THEN oplot, res0,harm_array[mct,ploti,*,compct],linestyle = compct+1 ELSE oplot,res1,harm_array[mct,ploti,*,compct],linestyle = compct+1
;                if ploti eq 0 THEN oplot, res0,harm_array[mct,ploti,*,3],linestyle = 2 ELSE oplot,res1,harm_array[mct,ploti,*,3],linestyle = 2
;                print,MINMAX(harm_array[mct,ploti,*,*])
;                stop
            ENDFOR
        ENDFOR
    ENDFOR
ENDFOR
;device,/close
multiplot, /reset
multiplot, /default
!p.multi = 0
cd,'/astro/users/christensen/code/Dwarfs'
END
