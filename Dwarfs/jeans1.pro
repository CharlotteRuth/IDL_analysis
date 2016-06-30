PRO gasden_sfr, files, xtitle = xtitle, ytitle = ytitle, title = title, outplot = outplot, legend_name = legend_name, colors = colors,times = times,yrange = yrange
;This program generates a plot like Fig 7 from Saitoh et al 2008 using
;prof.pro
msol = 2.362e5
kpc_unit = 1
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
loadct,39
;colors = [20,70,100]

;IF (keyword_set(print)) THEN BEGIN 
;    set_plot,'ps'
;    print,"plot"
;    device,filename = outplot,/color,bits_per_pixel = 8
;ENDIF ELSE set_plot,'x'
nbins = 225
disk_height = 15
linex = [2.5,7]
liney = [1e-3,1e-1]
symbol = indgen(n_ELEMENTS(times))*2+2
FOR i = N_ELEMENTS(files) - 1,0,-1 DO BEGIN
    print,files[i]
    FOR itime = 0, N_ELEMENTS(times) - 1 DO BEGIN
        filename = files[i] + times[itime]
        rtipsy, filename, h, g, d, s
        cofm, h, g, d, s, /pot
        particles = gal_align( h, g, d, s, rbox = 1, lunit = 1)
        h = particles.h
        g = particles.g
        d = particles.d
        s = particles.s
        s = s[where(abs(s.z*kpc_unit) lt disk_height)]
        d = d[where(abs(s.z*kpc_unit) lt disk_height)]
        g = g[where(abs(s.z*kpc_unit) lt disk_height AND g.tempg lt 15000)]
        sprofile = prof(s,'star',h.time,nbins = 4,rmin = 0,rmax = 10) 
        gprofile = prof(g,'gas' ,h.time,nbins = 4,rmin = 0,rmax = 10) 
        sprofile.rho = sprofile.rho*msol ; Msol/kpc^2  
        sprofile.sfr = sprofile.sfr*msol ; Msol/kpc^2/yr  
        gprofile.rho = gprofile.rho*msol ; Msol/kpc^2
        gprofile.rho = gprofile.rho*1e-6 ; Msol/pc^2
        print,filename
        print,gprofile.rho
        print,sprofile.sfr
        IF (itime eq 0) THEN plot,linex,liney,title = title, xtitle = xtitle, ytitle = ytitle,/ylog,/xlog,yrange = yrange,xrange = [0.01,100] 
        oplot,gprofile.rho,sprofile.sfr,psym = symbol[itime],color = colors[i]
    ENDFOR
;    stop
ENDFOR  
IF (keyword_set(legend_name)) THEN legend,legend_name,psym = symbol,/left,/top
IF (keyword_set(print)) THEN device,/close
END

PRO master_gasden_sfr
;This master plotter will loop through all interesting files
base = '../ResMassTests/'
legend_name = ['2.0 Gyr','2.5 Gyr','3.0 Gyr']
res_str = ['1E5R','1E4R','1E3R']
times = ['200','250','300']
colors = [70,150,240]

;base = '../ResSpecTests/1E5R/'
;legend_name = ['0.5 Gyr','1.0 Gyr','1.5 Gry', '2.0 Gyr','2.5 Gyr','3.0 Gyr']
;res_str = ['50pc','200pc','500pc','1kpc','2kpc','8kpc']
;times = ['050','100','150','200','250','300']
;colors = [20,70,100,150,200,240]

imf = ['k','ms']
mass_str = ['12','10']
yrange = [[1e-4,1e1],[1e-6,1e-1]]

;!P.MULTI = [0, 2, 2, 0, 0] 
!P.CHARSIZE = 1.25
;!X.MARGIN = [12,2]
;!Y.MARGIN = [7,2]
base = ['../ResMassTests/','../ResSpecTests/1E5R/']
res_str0 = ['1E5R','1E4R','1E3R']
res_str1 = ['50pc','200pc','500pc','1kpc','2kpc','8kpc']
colors0 = [70,150,240]
colors1 = [20,70,100,150,200,240]
;yrange0 = [0,8]
;yrange1 = [0,2]
;set_plot,'ps'
;device,filename = '../results/s_Spec_scaleHeight_k.eps',/color,bits_per_pixel = 8
;set_plot,'x'


FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
;    set_plot,'x'
    !P.MULTI = [0]
    multiplot,/reset
;    set_plot,'x'
    set_plot,'ps'
    device,filename = '../results/gasdenVsSfr_'+imf[imfct]+'.eps',/color,bits_per_pixel=8
    !P.MULTI = [0,2,2,0,0]
    FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
        FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN        
 ;       files = base + res_str + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00'
            IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00' ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00'
;            stop
            outplot = '../results/' + mass_str[mct] + 'M_' + 'sfrVsden_' + imf[imfct] + '.eps'
    ;    title = 'SFR vs Gas Density, '+textoidl('10^{'+mass_str[mct]+'}') + ' M' +sunsymbol()+ ' Halo'
            ytitle = textoidl('\Sigma_{SFR}')+" [M"+sunsymbol()+" yr!u-1!n kpc!u-2!n]"
            xtitle = textoidl('\Sigma_{gas}')+" [M"+sunsymbol()+" pc!u-2!n]"            
;        gasden_sfr, files, outplot = outplot, title = title, xtitle = xtitle, ytitle = ytitle, legend_name = legend_name,colors = colors,times = times,yrange=yrange[*,mct]

            multiplot
            IF (mct eq 0) THEN BEGIN 
                IF (ploti eq 0) THEN gasden_sfr, files, outplot = outplot, xtitle = '', ytitle = ytitle, legend_name = legend_name,colors = colors0,times = times,yrange=yrange[*,mct]  ELSE gasden_sfr, files, outplot = outplot, xtitle = '', ytitle = '', colors = colors1,times = times,yrange=yrange[*,mct]  
            ENDIF ELSE BEGIN
                IF (ploti eq 0) THEN gasden_sfr, files, outplot = outplot, xtitle = xtitle, ytitle = ytitle,colors = colors0,times = times,yrange=yrange[*,mct]  ELSE gasden_sfr, files, outplot = outplot, xtitle = xtitle, ytitle = '',colors = colors1,times = times,yrange=yrange[*,mct]     
            ENDELSE
        ENDFOR
    ENDFOR
    device,/close
ENDFOR

END
