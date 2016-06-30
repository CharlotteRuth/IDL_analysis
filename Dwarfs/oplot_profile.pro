PRO master_oplot_profile,sph = sph
;This master plotter will loop through all interesting files
loadct,39
!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1.0
!Y.Charsize = 1.0
!X.STYLE = 1
!Y.STYLE = 1
;!P.FONT = 0
;!X.OMARGIN = 0
;!Y.OMARGIN = 0

base = ['/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/','/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/']
imf = ['k']
legend_name0 = [textoidl('10^6'),textoidl('10^5'),textoidl('10^4'),textoidl('10^3')]
res_str0 = ['1E6R','1E5R','1E4R','1E3R']
colors0 = [0,0,0,0]
linestyle0 = [0,0,2,2]
legend_name1 = [textoidl('2.5\times10^{-4} R_{vir}'),textoidl('10^{-3} R_{vir} '),textoidl('10^{-2} R_{vir} '),textoidl('4\times10^{-2} R_{vir}')]
res_str1 = ['50pc','200pc','2kpc','8kpc']
colors1 = [0,0,0,0]
linestyle1 = [0,0,2,2]
thickness = [3,1,3,1]
mass_str = ['12','10']

type = ['s','g','d']
type_str = ['Stellar','Gas','Dark']
disk_height = 15
eps10= [11,44,440,1760]/1000.0
eps12 = [51.5,206,2060,8240]/1000.0
eps = [[eps12],[eps10]]
;[12][10]
;[s,g,d]                           ;[[1e2,1e10],[1e-3,1e-1],[1e-2,1e0]];[[1e2,9e8],[1e-6,9e-3],[1e-4,7e-2]]
IF(keyword_set(sph)) THEN yrange = [[[1e9,7e10],[1e-3,1e-1],[5e-1,5e3]],[[1e8,7e9],[1e-6,9e-3],[1e-2,1e2]]] $
                     ELSE yrange = [[[1e2,1e12],[1e-1,1e4],[1e0,1e6]],[[1e2,9e11],[1e-1,9e1],[1e-2,9e3]]] 
rmaxes = [100,1000]


FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
    FOR tct = 0, N_ELEMENTS(type) - 1 DO BEGIN
;        set_plot,'x'
        set_plot,'ps'
        if (keyword_set(sph)) then device,filename = '/astro/net/scratch1/christensen/DwarfResearch/results/'+type[tct]+'_profile_'+imf[imfct]+'_sph.eps',/color,bits_per_pixel = 8,/times else device,filename = '/astro/net/scratch1/christensen/DwarfResearch/results/'+type[tct]+'_profile_'+imf[imfct]+'.eps',/color,bits_per_pixel = 8,/times 
        !P.MULTI = [0]         
        !P.MULTI = [0,2,2,0,0] 
        !P.FONT = 0
        FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
            FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN
                IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00300' ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
                xtitle = 'Radius [kpc]'
                print,files

                IF(keyword_set(sph)) THEN BEGIN
                    outplot = '../results/' + mass_str[mct] + 'M_' + type[tct] + 'Spec_profileSph_mres_' + imf[imfct] + '.eps'               
                    print,outplot
                    IF (tct eq 0) THEN ytitle = textoidl('\rho_{star}')+" [M"+sunsymbol()+" kpc!u-3!n]"
                    IF (tct eq 1) THEN ytitle = textoidl('\rho_{gas}')+" [M"+sunsymbol()+" pc!u-3!n]"
                    IF (tct eq 2) THEN ytitle = textoidl('\rho_{dark}')+" [M"+sunsymbol()+" pc!u-3!n]"
                    multiplot
 ;                   oplot_profile, files, type[tct], outplot = outplot, title = title, xtitle = xtitle, ytitle = ytitle, legend_name = legend_name,yrange = yrange[*,tct,mct],colors = colors, /sph,rmax = rmaxes[mct],/PRINT,tickness = thickness
                    IF(mct eq 0) THEN BEGIN
                        IF (ploti eq 0) THEN oplot_profile,files,type[tct],linestyle = linestyle0,outplot = outplot,ytitle = ytitle, legend_name = legend_name0, colors = colors0, yrange = yrange[*,tct,mct],thickness=thickness, /sph $
                                        ELSE oplot_profile,files,type[tct],linestyle = linestyle1,outplot = outplot,                           eps = eps[*,mct], colors = colors1, yrange = yrange[*,tct,mct],thickness = thickness, /sph
                    ENDIF ELSE BEGIN
                        IF (ploti eq 0) THEN oplot_profile,files,type[tct],linestyle = linestyle0,outplot = outplot,ytitle = ytitle, xtitle = xtitle, colors = colors0, yrange = yrange[*,tct,mct], thickness=thickness, /sph $
                                        ELSE oplot_profile,files,type[tct],linestyle = linestyle1,outplot = outplot,                 xtitle = xtitle, colors = colors1, yrange = yrange[*,tct,mct], thickness = thickness, legend_name = legend_name1,eps = eps[*,mct], /sph
;YTICKNAME = [textoidl('10^2'),textoidl('10^4'),textoidl('10^6'),textoidl('10^8'),textoidl('10^{10}')]
                    ENDELSE
                ENDIF ELSE BEGIN
                    outplot = '../results/' + mass_str[mct] + 'M_' + type[tct] + 'profile_mres_' + imf[imfct] + '.eps'
                    print,outplot
                    IF (tct eq 0) THEN ytitle = textoidl('\Sigma_{star}')+" [M"+sunsymbol()+" kpc!u-2!n]"
                    IF (tct eq 1) THEN ytitle = textoidl('\Sigma_{gas}')+" [M"+sunsymbol()+" pc!u-2!n]"
                    IF (tct eq 2) THEN ytitle = textoidl('\Sigma_{dark}')+" [M"+sunsymbol()+" pc!u-2!n]"
                    multiplot
                    IF(mct eq 0) THEN BEGIN
                        IF (ploti eq 0) THEN oplot_profile,files,type[tct],linestyle = linestyle0,outplot = outplot,ytitle = ytitle, legend_name = legend_name0, colors = colors0, yrange = yrange[*,tct,mct],thickness = thickness $
                                        ELSE oplot_profile,files,type[tct],linestyle = linestyle1,outplot = outplot,                           eps = eps[*,mct], colors = colors1, yrange = yrange[*,tct,mct],thickness = thickness
                    ENDIF ELSE BEGIN
                        IF (ploti eq 0) THEN oplot_profile,files,type[tct],linestyle = linestyle0,outplot = outplot,ytitle = ytitle, xtitle = xtitle, colors = colors0, yrange = yrange[*,tct,mct], thickness = thickness $
                                        ELSE oplot_profile,files,type[tct],linestyle = linestyle1,outplot = outplot,                 xtitle = xtitle, colors = colors1, yrange = yrange[*,tct,mct], thickness = thickness, legend_name = legend_name1, eps = eps[*,mct] 
;YTICKNAME = [textoidl('10^2'),textoidl('10^4'),textoidl('10^6'),textoidl('10^8'),textoidl('10^{10}')]
                    ENDELSE
                ENDELSE
            ENDFOR
        ENDFOR
        xyouts,[3000,10100,3000,10100],[11100,11100,6000,6000],[textoidl('10^{12}M')+sunsymbol(),textoidl('10^{12}M')+sunsymbol(),textoidl('10^{10}M'+sunsymbol()),textoidl('10^{10}M'+sunsymbol())],charsize = 1.05,/device
        multiplot,/reset
        device,/close
;        stop
    ENDFOR
ENDFOR
END

PRO oplot_profile, files, comp, outplot = outplot, title = title, xtitle = xtitle, ytitle = ytitle, legend_name = legend_name, yrange = yrange,print = print, colors = colors, sph=sph, rmax=rmax,linestyle = linestyle,eps = eps, thickness=thickness, _EXTRA = _EXTRA
;This will plot a series of profiles over each other
msol = 2.362e5
kpc_unit = 1

nbins = 225
disk_height = 15
IF (keyword_set(sph)) THEN nfiles = N_ELEMENTS(files) -1 ELSE nfiles = 2
nfiles = N_ELEMENTS(files)-1
;print,files
;print,comp,yrange
;print,''
FOR i = nfiles,0,-1 DO BEGIN
    print,files[i]
    rtipsy, files[i], h, g, d, s
    cofm, h, g, d, s, /pot
    particles = gal_align( h, g, d, s, rbox = 1, lunit = 1)
    h = particles.h
    g = particles.g
    d = particles.d
    s = particles.s
    rmax = 3e-1 ;20
    rmin = 3e-2; 0.1
    nbins = 20; 100
    IF (keyword_set(sph)) THEN BEGIN
        IF (comp eq 's') THEN profile = prof(s,'star',h.time,nbins = nbins,rmin = rmin,rmax = rmax,/sph)                
        IF (comp eq 'd') THEN profile = prof(d,'dark',h.time,nbins = nbins,rmin = rmin,rmax = rmax,/sph)         
        IF (comp eq 'g') THEN profile = prof(g,'gas' ,h.time,nbins = nbins,rmin = rmin,rmax = rmax,/sph)      
        profile.rho = profile.rho*msol ; Msol/kpc^3
        IF (comp eq 'd' OR comp eq 'g') THEN profile.rho = profile.rho*1e-9 ; Msol/pc^3
    ENDIF ELSE BEGIN
        s = s[where(abs(s.z*kpc_unit) lt 15)]
        d = d[where(abs(s.z*kpc_unit) lt 15)]
        g = g[where(abs(s.z*kpc_unit) lt 15 AND g.tempg lt 15000)]
        IF (comp eq 's') THEN profile = prof(s,'star',h.time,nbins = nbins,rmin = rmin,rmax = rmax)                
        IF (comp eq 'd') THEN profile = prof(d,'dark',h.time,nbins = nbins,rmin = rmin,rmax = rmax)         
        IF (comp eq 'g') THEN profile = prof(g,'gas' ,h.time,nbins = nbins,rmin = rmin,rmax = rmax)      
        profile.rho = profile.rho*msol ; Msol/kpc^2
        IF (comp eq 'd' OR comp eq 'g') THEN profile.rho = profile.rho*1e-6 ; Msol/pc^2
    ENDELSE
    if NOT (keyword_set(eps)) THEN  begin
        print,comp,MAX(profile.rho),' Prof: '
        print,profile.rbins
        print,profile.rho
        result = linfit(alog10(profile.rbins),alog10(profile.rho))
        print,"fit: ",result
        print,' '
    endif
;    stop
    IF (i eq nfiles) THEN rmax = MAX(profile.rbins)
    IF (i eq nfiles) THEN plot,profile.rbins,profile.rho,/ylog,xtitle = xtitle, ytitle = ytitle, title = title,xrange=[rmin,rmax],yrange = yrange,linestyle = linestyle[i],thick=thickness[i],/xlog;,xrange=[0.1,1]
    oplot,profile.rbins,profile.rho,linestyle = linestyle[i],thick=thickness[i]
    if (keyword_set(eps)) THEN oplot,[0.1,0.1+eps[i]],[10.0^(6.5-i),10.0^(6.5-i)],linestyle = linestyle[i],thick = thickness[i]
ENDFOR
IF (keyword_set(legend_name)) THEN $
  IF (N_ELEMENTS(colors) gt 3) THEN legend,legend_name,linestyle = linestyle,/right,charsize = 1.0, thick=thickness,pspacing = 1.5*1.25 ELSE legend,legend_name,/right,linestyle = linestyle,charsize = 1.0, thick=thickness, pspacing = 1.5*1.25
END
