;
;
;
PRO ellipticity,files,res,particle_type = particle_type, symbol = symbol,oplot = oplot, color = color, linestyle = linestyle,xtitle = xtitle,xrange = xrange,cut = cut,thickness=thickness,xtick=xtick
;res = [1e5,1e4,1e3,1e2,5e1]
msol = 2.362e5
kpc_unit = 1
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
loadct,39
thick = [1,2]
print,particle_type
IF NOT KEYWORD_SET(nbins) THEN nbins = 500.
IF NOT KEYWORD_SET(maxr) THEN maxr = 200.
IF NOT KEYWORD_SET(maxz) THEN maxz = 50.
IF NOT KEYWORD_SET(psym) THEN psym = 4.
;IF NOT KEYWORD_SET(particle_type) THEN particle_type = 1.

slr = findgen(nbins)
slz = findgen(nbins)
nfiles = N_ELEMENTS(files)
e = fltarr(nfiles)
r = findgen(nbins + 1)*maxr/nbins
z = findgen(nbins + 1)*maxz/nbins
mr = fltarr(nbins)
mz = fltarr(nbins)

FOR i = 0, nfiles-1 DO BEGIN
    print,files[i],res[i]
    rtipsy, files[i], h, g, d, s
    cofm, h, g, d, s, /pot
    particles = gal_align( h, g, d, s, rbox = 1, lunit = 1)
    h = particles.h
    g = particles.g
    d = particles.d
    s = particles.s
    IF (particle_type eq 0) THEN p = s ELSE p = g
    IF (N_ELEMENTS(p) ge 10) THEN BEGIN
        r = SQRT(p.x^2 + p.y^2)
        r_sort = SORT(r)
        r = r[r_sort]
        pr = p[r_sort]
        mr = pr.mass

        z = ABS(p.z)
        z_sort = SORT(z)
        z = z[z_sort]
        pz = p[z_sort]
        mz = pz.mass

        FOR ip = 1L, N_ELEMENTS(p)-1 DO BEGIN
            mr[ip] = mr[ip] + mr[ip - 1]
            mz[ip] = mz[ip] + mz[ip - 1]        
        ENDFOR
    
        mr = mr/TOTAL(p.mass)
        mz = mz/TOTAL(p.mass)
;    window,1
;    plot,r,mr
;    oplot,[0,500],[0.9,0.9]
;    window,2
;    plot,z,mz
        oplot,[0,100],[0.9,0.9]    
        temp = MIN(ABS(mr - 0.9),r90_ind)
        temp = MIN(ABS(mz - 0.9),z90_ind)
        r90_ind_first = r90_ind - 1
        r90_ind_last = r90_ind + 1
        z90_ind_first = z90_ind - 1
        z90_ind_last = z90_ind + 1
        r90 = SPLINE(mr[r90_ind_first:r90_ind_last],r[r90_ind_first:r90_ind_last],[0.9])
        z90 = SPLINE(mz[z90_ind_first:z90_ind_last],z[z90_ind_first:z90_ind_last],[0.9])
;    stop
        e[i] = SQRT((r90^2 - z90^2)/r90^2)
        print,''
        print,'Ellipticity: ',z90,r90,e[i],particle_type
        print,''
    ENDIF
;    stop
ENDFOR
!X.style = 1
!Y.style = 1
good = where(e ne 0)

IF KEYWORD_SET(cut) THEN BEGIN
    IF NOT KEYWORD_SET(oplot) THEN plot,res[2:N_ELEMENTS(files)-1],e[2:N_ELEMENTS(files)-1],psym = symbol,xtitle = xtitle,ytitle = 'e',/xlog,yrange = [0,1],linestyle = linestyle, xrange = xrange;,ytitle = textoidl('z_{90}/r_{90}')
    IF (particle_type eq 1) THEN oplot,res,e,psym = symbol,color = color,linestyle = linestyle,thick=thickness ELSE oplot,res[2:N_ELEMENTS(files)-1],e[2:N_ELEMENTS(files)-1],psym = symbol,color = color, linestyle = linestyle,thick = thickness
ENDIF ELSE BEGIN
    IF NOT KEYWORD_SET(oplot) THEN plot,res[good],e[good],psym = symbol,xtitle = xtitle,ytitle = 'e',/xlog,yrange = [0.5,1],linestyle = linestyle, xrange = xrange,thick=thickness,xtickname=xtick
    oplot,res[good],e[good],psym = symbol,color = color,linestyle = linestyle,thick=thickness
ENDELSE
;stop
END

PRO ellipticity_master
;This master plotter will loop through interesting files
;!P.MULTI = [0, 2, 2, 0, 0] 
!P.thick = 1.5
!P.CHARSIZE = 1.5
!X.Charsize = 1.35
!Y.Charsize = 1.35
!Y.Style = 1
!X.style = 1

;!X.MARGIN = [12,2]
;!Y.MARGIN = [7,2]
mass = ['9','10','11','12','13']
symbols = [-2,-4,-5,-6,-7]
symbols1 = [-2,-4,-5,-6]
linestyle = [0,0,2,2,3]
linestyle1 = [0,0,2,2]
color = [0,0,0,0,0]
color1 = [0,0,0,0]
thickness0 = [3,1,3,1,3]
thickness1 = [3,1,3,1]
particles = [0];,1]
xtick = [textoidl('10^{-1}'),textoidl('10^{-2}'),textoidl('10^{-3}'),textoidl('10^{-4}')]

base = ['/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/','/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/']
outplot = [['../results/ellip_massRes_','../results/ellip_massRes_'],['../results/ellip_spec12_','../results/ellip_spec10_']]
mass_str = ['12','10']
imf = ['k']
res_str0 = ['5E1R','1E2R','1E3R','1E4R','1E5R','1E6R']
res_str1 = ['8kpc','2kpc','1kpc','500pc','200pc','50pc']
resolution0 = [50,100,1000,10000,100000,1000000]
resolution1 = [2e-2,1e-2,5e-3,2.5e-3,1e-3,2.5e-4]
iter = [['13M','12M','11M','10M','9M'],['1E5R','1E4R','1E3R','1E2R','5E1R']]
mass_str = ['12','10']
xtitle = ['Number of DM Particles','Softening Length ' + textoidl('[R_{vir}]')]
ytitle = textoidl('e')
legend_name0 = [textoidl('10^{13}')+'M'+sunsymbol(),textoidl('10^{12}')+'M'+sunsymbol(),textoidl('10^{11}')+'M'+sunsymbol(),textoidl('10^{10}')+'M'+sunsymbol(),textoidl('10^{9}')+'M'+sunsymbol()]
legend_name1 =  [textoidl('10^5'),textoidl('10^4'),textoidl('10^3'),textoidl('10^2')]
xrange0 = [500,2e6]
xrange1 = [0.1,0.0001]

;colors0 = [0,0,100,100,0]
;colors1 = [0,0,100,100,0]
;linestyle0 = [0,2,0,2,1]
;linestyle1 = [0,2,0,2,1]
;psym0 = [-2,-4,-5,-6,-7]
;psym1 = [-2,-4,-5,-6,-7]

FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
    FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
        FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN
            filename = '~/Scratch1/DwarfResearch/results/'+outplot[mct,ploti]+imf[imfct]+'.eps'
            print,filename
            set_plot,'x'
            set_plot,'ps'
            device,filename = filename,/color,bits_per_pixel= 8,/times
            !P.FONT = 0
            if ploti eq 0 then endloop = N_ELEMENTS(mass) - 1 ELSE endloop = 3
            FOR i = 0, endloop DO BEGIN
                FOR ip = 0, N_ELEMENTS(particles) - 1 DO BEGIN
;                files = ['../ResMassTests/1E5R/'+mass[i]+'M_k/o'+mass[i]+'M_1.00300','../ResMassTests/1E4R/'+mass[i]+'M_k/o'+mass[i]+'M_1.00300','../ResMassTests/1E3R/'+mass[i]+'M_k/o'+mass[i]+'M_1.00300','../ResMassTests/1E2R/'+mass[i]+'M_k/o'+mass[i]+'M_1.00300','../ResMassTests/5E1R/'+mass[i]+'M_k/o'+mass[i]+'M_1.00300']
                    IF (ploti eq 0) THEN BEGIN 
                        files = base[ploti] + res_str0 + '/' +iter[i,0]+'_'+imf[imfct]+'/o'+iter[i,0]+'_1.00300' 
                        IF(i ne 1 and i ne 3) THEN BEGIN
                            resolution0prime = resolution0[0:4]
                            filesprime = files[0:4]
                        ENDIF ELSE BEGIN
                            resolution0prime = resolution0
                            filesprime = files
                        ENDELSE
                        print,filesprime
                        IF (i eq 0 AND particles[ip] eq 0) THEN ellipticity,filesprime,resolution0prime,symbol = symbols[i],particle_type = particles[ip],color = color[i], linestyle = linestyle[i],xtit = xtitle[ploti],xrange = xrange0,thickness=thickness0[i] ELSE ellipticity,filesprime,resolution0prime,symbol = symbols[i],/oplot,particle_type = particles[ip],color = color[i], linestyle = linestyle[i],thickness=thickness0[i]
                    ENDIF ELSE BEGIN
                        files = base[ploti] +iter[i,1]+ '/'+mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
                        print,files
                        IF (i eq 0 AND particles[ip] eq 0) THEN ellipticity,files,resolution1,symbol = symbols[i],particle_type = particles[ip],color = color[i], linestyle = linestyle[i],xtit = xtitle[ploti],xrange = xrange1,thickness=thickness1[i],xtick=xtick ELSE ellipticity,files,resolution1,symbol = symbols[i],/oplot,particle_type = particles[ip],color = color[i], linestyle = linestyle[i],thickness=thickness1[i],xtick=xtick
                    ENDELSE
                ENDFOR
            ENDFOR
             if(ploti eq 0) then legend,legend_name0,color = color,linestyle = linestyle,psym = symbols,/bottom,/right,charsize = 1.5,thick=thickness0,pspacing = 1.5*1.5  else begin
                if (ploti eq 1) then BEGIN 
                    legend,legend_name1, color = color1, linestyle = linestyle1,psym = symbols1,/bottom,charsize = 1.5,thick=thickness1,/right,pspacing = 1.5*1.5
                    IF mct eq 0 THEN xyouts,[3500],[10500],[textoidl('10^{12}M')+sunsymbol()],/device,charsize = 1.5 ELSE xyouts,[3500],[10500],[textoidl('10^{10}M')+sunsymbol()],/device,charsize = 1.5
                ENDIF
            endelse
;            stop
            device,/close
        ENDFOR
    ENDFOR
ENDFOR
END
