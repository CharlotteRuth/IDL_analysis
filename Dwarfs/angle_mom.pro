PRO angle_mom,files,resolutions,xtit = xtit, ytit = ytit,iter = iter, color = color, linestyle = linestyle, psym = psym
msol = 2.362e5
kpc_unit = 1
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
;!X.MARGIN = [10,2]
;!Y.MARGIN = [7,2]

nfiles = N_ELEMENTS(files)-1
yvec = FLTARR(nfiles+1)
FOR i = nfiles,0,-1 DO BEGIN
    PRINT,I
    print,files[i]
    rtipsy, files[i], h, g, d, s
    cofm, h, g, d, s, /pot
    particles = gal_align( h, g, d, s, rbox = 1, lunit = 1)
    p = particles.g
    rvec = fltarr(N_ELEMENTS(particles.g),3)
    vvec = fltarr(N_ELEMENTS(particles.g),3)
    rvec[*,0] = p.x*p.mass
    rvec[*,1] = p.y*p.mass
    rvec[*,2] = p.z*p.mass

    vvec[*,0] = p.vx
    vvec[*,1] = p.vy
    vvec[*,2] = p.vz  
    
    lvec = crossp_multi(rvec, vvec)
    lz_tot = TOTAL(abs(lvec[*,2]))
    l_tot = TOTAL(SQRT(lvec[*,0]^2 + lvec[*,1]^2 +lvec[*,2]^2))
    yvec[i] = lz_tot/l_tot
ENDFOR
;xrange = [resolutions[0],resolutions[N_ELEMENTS(resolutions)-1]]
if (iter eq 0) then plot,resolutions,yvec,ytitle = ytit, xtitle = xtit,/xlog,yrange = [0,1]
oplot,resolutions,yvec,color = color, linestyle = linestyle,psym = psym
print,yvec
END

PRO master_angle_mom
;This master plotter will loop through interesting files
;!P.MULTI = [0, 2, 2, 0, 0] 
!P.thick = 1.5
!X.Charsize = 1.25
!Y.Charsize = 1.25
!P.CHARSIZE = 1.25
;!X.MARGIN = [12,2]
;!Y.MARGIN = [7,2]
loadct,39
set_plot,'x'

base = ['../ResMassTests/','../ResSpecTests/']
outplot = [['angle_mom_massRes_','angle_mom_massRes_'],['angle_mom_spec12_','angle_mom_spec10_']]
mass_str = ['12','10']
imf = ['k']
res_str0 = ['5E1R','1E2R','1E3R','1E4R','1E5R']
res_str1 = ['50pc','200pc','500pc','1kpc','2kpc','8kpc']
resolution0 = [50,100,1000,10000,100000]
resolution1 = [2e-2,1e-2,5e-3,2.5e-3,1e-3,2.5e-4]
iter = [['13M','12M','11M','10M','9M'],['1E5R','1E4R','1E3R','1E2R','5E1R']]
mass_str = ['12','10']
xtitle = ['Number of Particles','Softening Length ' + textoidl('[R_{vir}]')]
ytitle = textoidl('L_z / L_{tot}')
legend_name0 = [textoidl('10^{13}')+'M'+sunsymbol(),textoidl('10^{12}')+'M'+sunsymbol(),textoidl('10^{11}')+'M'+sunsymbol(),textoidl('10^{10}')+'M'+sunsymbol(),textoidl('10^{9}')+'M'+sunsymbol()]
legend_name1 =  ['100K','10K','1K','100','50']

colors0 = [0,0,100,100,0]
colors1 = [0,0,100,100,0]
linestyle0 = [0,2,0,2,1]
linestyle1 = [0,2,0,2,1]
psym0 = [-2,-4,-5,-6,-7]
psym1 = [-2,-4,-5,-6,-7]

FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
    FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
        FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN
            if((ploti eq 0 and mct eq 0) or ploti eq 1) then begin
                set_plot,'ps'
                device,filename ='../results/'+outplot[mct,ploti]+imf[imfct]+'.eps',/color,bits_per_pixel= 8
                print,'../results/'+outplot[ploti,mct]+imf[imfct]+'.eps'
            endif
            FOR i = 0, N_ELEMENTS(iter[*,ploti]) - 1 DO BEGIN
                IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' +iter[i,0]+'_'+imf[imfct]+'/o'+iter[i,0]+'_1.00300' ELSE files = base[ploti] +iter[i,1]+ '/'+mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
                print,files
                IF (mct eq 0) THEN BEGIN 
                    IF (ploti eq 0) THEN angle_mom, files, resolution0, xtit = xtitle[0], ytit = ytitle,iter = i,linestyle = linestyle0[i], psym = psym0[i],color = colors0[i] ELSE angle_mom, files, resolution1, xtit = xtitle[1], iter = i,ytit = ytitle,linestyle = linestyle1[i], psym = psym1[i],color = colors1[i]
                ENDIF ELSE BEGIN
                    IF (ploti eq 1) THEN angle_mom, files, resolution1, xtit = xtitle[1], ytit = ytitle, linestyle = linestyle1[i], iter = i, psym = psym1[i],color = colors1[i]
                ENDELSE
            ENDFOR
            if(ploti eq 0 and mct eq 0) then legend,legend_name0,color = colors0,linestyle = linestyle0,psym = psym0 else begin
                if (ploti eq 1) then legend,legend_name1, color = colors1, linestyle = linestyle1,psym = psym1,/bottom
            endelse
        ENDFOR
        if((ploti eq 0 and mct eq 0) or ploti eq 1) then device,/close
    ENDFOR
ENDFOR
END
