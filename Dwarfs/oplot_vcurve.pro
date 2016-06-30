PRO master_oplot_vcurve,sph = sph
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
linestyle0 = [0,0,2,1]
legend_name1 = [textoidl('2.5\times10^{-4} R_{vir}'),textoidl('10^{-3} R_{vir} '),textoidl('10^{-2} R_{vir} '),textoidl('4\times10^{-2} R_{vir}')]
res_str1 = ['50pc','200pc','2kpc','8kpc']
colors1 = [0,0,0,0]
linestyle1 = [0,0,2,1]
thickness = [5,1,5,1]
mass_str = ['12','10']

type = ['s','d','g']
type_str = ['Stellar','Gas']
disk_height = 15
eps10= [11,44,440,1760]/1000.0
eps12 = [51.5,206,2060,8240]/1000.0
eps = [[eps12],[eps10]]
yrange = [[0,700],[0,100]] ; spherical
xrange = [[0,11.99],[0,12]] ; spherical
rmaxes = [100,1000]


FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
    FOR tct = 0, 0 DO BEGIN ;N_ELEMENTS(type) - 1 DO BEGIN
        set_plot,'x'
        set_plot,'ps'
        device,filename = '/astro/net/scratch1/christensen/DwarfResearch/results/'+type[tct]+'_vcurve_'+imf[imfct]+'.eps',/color,bits_per_pixel = 8,/times   
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
                    ytitle = 'Velocity'
                    multiplot
                    oplot_vcurve, files, type[tct], outplot = outplot, title = title, xtitle = xtitle, ytitle = ytitle, legend_name = legend_name,colors = colors, /sph,rmax = rmaxes[mct],/PRINT,tickness = thickness,yrange = yrange[*,mct]
                ENDIF ELSE BEGIN
                    outplot = '../results/' + mass_str[mct] + 'M_' + type[tct] + 'profile_mres_' + imf[imfct] + '.eps'
                    print,outplot
                    ytitle = 'Velocity'
                    multiplot
                    IF(mct eq 0) THEN BEGIN
                        IF (ploti eq 0) THEN oplot_vcurve,files,type[tct],linestyle = linestyle0,outplot = outplot,ytitle = ytitle, legend_name = legend_name0, colors = colors0,thickness=thickness,yrange = yrange[*,mct],xrange=xrange[*,ploti] ELSE oplot_vcurve,files,type[tct],linestyle = linestyle1,outplot = outplot, colors = colors1,thickness = thickness,yrange = yrange[*,mct],xrange=xrange[*,ploti];,eps = eps[*,mct]
                    ENDIF ELSE BEGIN
                        IF (ploti eq 0) THEN oplot_vcurve,files,type[tct],linestyle = linestyle0,outplot = outplot,ytitle = ytitle, xtitle = xtitle, colors = colors0,  thickness=thickness,yrange = yrange[*,mct],xrange=xrange[*,ploti] ELSE oplot_vcurve,files,type[tct],linestyle = linestyle1,outplot = outplot,xtitle = xtitle, colors = colors1, legend_name = legend_name1,thickness = thickness,yrange = yrange[*,mct],xrange=xrange[*,ploti];,eps = eps[*,mct]
;YTICKNAME = [textoidl('10^2'),textoidl('10^4'),textoidl('10^6'),textoidl('10^8'),textoidl('10^{10}')]
                    ENDELSE
                ENDELSE
            ENDFOR
        ENDFOR
        xyouts,[3000,10100,3000,10100],[11100,11100,6000,6000],[textoidl('10^{12}M')+sunsymbol(),textoidl('10^{12}M')+sunsymbol(),textoidl('10^{10}M'+sunsymbol()),textoidl('10^{10}M'+sunsymbol())],charsize = 1.05,/device
        multiplot,/reset
        device,/close
    ENDFOR
ENDFOR
END

PRO oplot_vcurve, files, comp, outplot = outplot, title = title, xtitle = xtitle, ytitle = ytitle, legend_name = legend_name, yrange = yrange,print = print, colors = colors, sph=sph, rmax=rmax,linestyle = linestyle,eps = eps, thickness=thickness, _EXTRA = _EXTRA,xrange=xrange
;This will plot a series of profiles over each other
msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 2.362e5
dKpcUnit = 1
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70

nbins = 225
disk_height = 15
IF (keyword_set(sph)) THEN nfiles = N_ELEMENTS(files) -1 ELSE nfiles = 2
nfiles = N_ELEMENTS(files)-1
FOR i = nfiles,0,-1 DO BEGIN
    print,files[i]
    rtipsy, files[i], h, g, d, s
    cofm, h, g, d, s, /pot
    particles = gal_align( h, g, d, s, rbox = 1, lunit = 1)
    h = particles.h
    g = particles.g
    d = particles.d
    s = particles.s

    distancesg = SQRT(g.x^2.0 + g.y^2.0 + g.z^2.0)
    distancess = SQRT(s.x^2.0 + s.y^2.0 + s.z^2.0)
    distancesd = SQRT(d.x^2.0 + d.y^2.0 + d.z^2.0)
    distances = [distancesg,distancess,distancesd]
    massg = g.mass
    masss = s.mass
    massd = d.mass
    mass = [massg,masss,massd]
    maxd = MAX(distances)
    mind = MIN(distances)

    bind = 40.0/dKpcUnit/nbins
    bins = (findgen(nbins)+1)*bind
    tmass = fltarr(nbins)
    density = fltarr(nbins)

    FOR ib = 0, nbins-1 DO BEGIN
        ind = where(distances LE bins[ib])
        if (ind[0] ne -1)then tmass[ib] = TOTAL(mass[ind])else tmass[ib] = 0
    ENDFOR
    velocitycurve = SQRT(tmass*dMsolUnit*msol*grav/bins/dkpcUnit/kpc)
    print,MINMAX(velocitycurve)
    IF (i eq nfiles) THEN plot,bins*dKpcUnit,velocitycurve,xtitle = xtitle, ytitle = ytitle, title = title,linestyle = linestyle[i],thick=thickness[i],yrange = yrange,xrange = xrange
    oplot,bins*dKpcUnit,velocitycurve,linestyle = linestyle[i],thick=thickness[i]
    if (keyword_set(eps)) THEN BEGIN
        dy = (yrange[1]-yrange[0])/20.0
        ystart = (yrange[1]-yrange[0])/2 + 2.5*dy/5.0
        yvals = findgen(N_ELEMENTS(eps))*dy + ystart        
        oplot,[0,eps[i]],[yvals[i],yvals[i]],linestyle = linestyle[i],thick = thickness[i]
        STOP
    ENDIF
ENDFOR
IF (keyword_set(legend_name)) THEN $
  IF (N_ELEMENTS(colors) gt 3) THEN legend,legend_name,linestyle = linestyle,/right,charsize = 1.0, thick=thickness,pspacing = 1.5*1.25 ELSE legend,legend_name,/right,linestyle = linestyle,charsize = 1.0, thick=thickness, pspacing = 1.5*1.25

END
