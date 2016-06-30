pro checkLW_graph_master
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23

loadct,39
dir = '/astro/net/scratch2/christensen/MolecH/11M/Cloud_1e3/'
file = 'o11M_1.onestar.00001'
file = 'o11M_1.00005'
dMsolUnit       = 2.362e5
dKpcUnit        = 1
;dInitStarMass   = 17.2043010752688

dir = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_zsol/rad.psource/'
file = 'MW_disk.'
step = '00010'
dMsolUnit       = 1.36e17
dKpcUnit        = 1e5
;dInitStarMass   = 2.535817803095e-15   ;3.34632e-14

dir = '/astro/net/scratch2/christensen/MolecH/11M/Cloud_1e5/'
file = 'o11M_00300.'
step = '00001'
dMsolUnit       = 2.362e5
dKpcUnit        = 1d0
key = ['Box','Moment','Point Source','Smoothed Source','Actual']
readcol,dir+file+step+".lw",lw_all,/silent
;lw = lw_all[1:h.ngas]*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc
;lw_star = lw_all[N_ELEMENTS(lw_all) - h.nstar:N_ELEMENTS(lw_all) - 1]*1d30

;readcol,dir+file+step+".lw_calc",lw_calc_all,/silent
;lw_calc = lw_calc_all[1:h.ngas]*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc

;readcol,dir+file+'mom.'+step+".lw",lw_mom_all,/silent
;lw_mom = lw_mom_all[1:h.ngas]*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc

;readcol,dir+file+'psource.'+step+".lw",lw_psource_all,/silent
;lw_psource = lw_psource_all[1:h.ngas]*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc

;readcol,dir+file+'ssource.'+step+".lw",lw_ssource_all,/silent
;lw_ssource = lw_ssource_all[1:h.ngas]*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc

;o11M_00300.00001
cx = -4.692374e-01 
cy = 3.408022e-01 
cz = -6.159824e-01

;h277.cosmo50cmb.1536g1MBWKBH.00001
;cx = 1.742388e-06 
;cy = 7.682247e-07 
;cz = -7.504132e-07

;dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h277.cosmo50cmb.1536g1MBWKBH/'
;file = 'h277.cosmo50cmb.1536g1MBWKBH.'
;step = '00001'
;dMsolUnit       = 1.84793e16
;dKpcUnit        = 50000.

dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/'
tfile = dir + 'h516.cosmo25cmb.2304g14HBWK.00512.halo.1.std'
lwfiles = dir + ['h516.cosmo25cmb.2304g14HBWK.00512.halo.1.lw_calc_sigma','h516.cosmo25cmb.2304g14HBWK.00512.halo.1.lw']
dKpcUnit =  25000.
dMsolUnit = 2.310e15
outplot = '~/h516.cosmo25cmb.paper'
densityunit =  dMsolUnit * gm_per_msol * 5.9753790e+23/dKpcUnit^3/cm_per_kpc^3
rhocut = 0.01/densityunit
checkLW_graph,tfile,lwfiles,dKpcUnit,rhocut = rhocut;,outplot = outplot
end

pro checkLW_graph,tfile,lwfiles,dKpcUnit,outplot = outplot,keys = keys,center = center,color = color,rhocut = rhocut,changa = changa
cmPerKpc = 3.08568025d21
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps' 
    nbins=100.0
    linestyles = [0,2]
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25
;    !p.font=0 
    cb_charsize = 0.75
    l_charsize = 0.75
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    IF NOT KEYWORD_SET(color) THEN BEGIN
        loadct,0
        color = fltarr(N_ELEMENTS(lwfiles) - 1 )
    ENDIF ELSE BEGIN
        loadct,39
        if color[0] eq 1 then color = findgen(N_ELEMENTS(lwfiles))*240/N_ELEMENTS(files)
    ENDELSE
    device,filename=outplot+'_lwcheck.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2
ENDIF ELSE BEGIN
    set_plot,'x'
    nbins=100.0
    linestyles = [0,2]
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    cb_charsize = 1.0
    l_charsize = 1.0
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    IF NOT KEYWORD_SET(color) THEN BEGIN
        color = fltarr(N_ELEMENTS(lwfiles) - 1 ) + 255 
    ENDIF ELSE BEGIN
        loadct,39
        if color[0] eq 1 then color = findgen(N_ELEMENTS(lwfiles))*254/N_ELEMENTS(files)
    ENDELSE
    window,0,xsize = 712,ysize = 392
ENDELSE

rtipsy,tfile,h,g,d,s
IF KEYWORD_SET(center) THEN BEGIN
    cx = center[0]
    cy = center[1]
    cz = center[2]
ENDIF ELSE BEGIN
    cx = 0
    cy = 0
    cz = 0
ENDELSE
s.x = (s.x - cx)*dKpcUnit
s.y = (s.y - cy)*dKpcUnit
s.z = (s.z - cz)*dKpcUnit
g.x = (g.x - cx)*dKpcUnit
g.y = (g.y - cy)*dKpcUnit
g.z = (g.z - cz)*dKpcUnit
d.x = (d.x - cx)*dKpcUnit
d.y = (d.y - cy)*dKpcUnit
d.z = (d.z - cz)*dKpcUnit

rstar = SQRT(s.x*s.x + s.y*s.y + s.z*s.z) ;*dKpcUnit
rgas  = SQRT(g.x*g.x + g.y*g.y + g.z*g.z) ;*dKpcUnit
rdark = SQRT(d.x*d.x + d.y*d.y + d.z*d.z) ;*dKpcUnit
readarr,lwfiles[0],h,lw_star_comp,part = 'star',/ascii
readarr,lwfiles[0],h,lw_comp,part = 'gas',/ascii
stop
IF keyword_set(changa) THEN BEGIN
    lw_comp = 10^lw_comp/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc
    lw_star_comp = 10^lw_star_comp
ELSE BEGIN
    lw_comp = lw_comp*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc
    lw_star_comp = lw_star_comp*1d30
ENDELSE
plot,[0,10],[1, 1],/ylog,xtitle = 'Radius [kpc]',ytitle = 'Simulated LW Flux / Actual LW Flux',xrange = [0,8.0],yrange = [0.0001, 100]

FOR i = 1, N_ELEMENTS(lwfiles) - 1 DO BEGIN    
    readarr,lwfiles[i],h,lw,part = 'gas',/ascii
    lw = lw*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc
    IF KEYWORD_SET(rhocut)  THEN oplot,rgas[where(g.dens gt rhocut)],lw[where(g.dens gt rhocut)]/lw_comp[where(g.dens gt rhocut)],psym = 3,color = color[i - 1] ELSE oplot,rgas,lw/lw_comp,psym = 3,color = color[i - 1]
   
 ;   window,0
;    oplot,rgas,lw/lw_comp,/xlog,/ylog,psym = 3,xrange = [0.1,10],yrange = [1e6,1e12],xtitle = 'Radius [kpc]',ytitle = 'LW Flux [photons s^-1 cm^-2]',title = file
 
;   oplot,rgas,lw,psym = 3,color = 240
 ;   oplot,rgas,lw_calc,psym = 3
 ;   oplot,rgas,lw_mom,psym = 3,color = 120
 ;   oplot,rgas,lw_psource,psym = 3,color = 60
 ;   oplot,rgas,lw_ssource,psym = 3,color = 200
 ;   legend,['Box','Moment','Point Source','Smoothed Source','Actual'],linestyle = [1,1,1,1,1],color = [240,120,60,200,0]
 ;   stop
;    window,1
;    plot,lw,lw_calc,/xlog,/ylog,psym = 3,xtitle = 'Code',ytitle = 'Actual',xrange = [1e6,1e12],yrange = [1e6,1e12],title = file
;    oplot,lw,lw_calc,psym = 3,color = 240
;    oplot,lw_mom,lw_calc,psym = 3,color = 120
;    oplot,lw_psource,lw_calc,psym = 3,color = 60
;    oplot,lw_ssource,lw_calc,psym = 3,color = 200
;    oplot,[1e6,1e12],[1e6,1e12]
;    legend,['Box','Moment','Point Source','Smoothed Source'],linestyle = [1,1,1,1],color = [240,120,60,200]
    
ENDFOR
IF KEYWORD_SET(key) THEN legend,keys,color = color,/right,linestyle = fltarr([N_ELEMENTS(lwfiles - 1)]),/bottom,charsize = l_charsize
IF KEYWORD_SET(outplot) THEN BEGIN
    device,/close
    device,filename=outplot+'_lwcheckhist.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2
ENDIF
stop

FOR i = 1, N_ELEMENTS(lwfiles) - 1 DO BEGIN    
    readarr,lwfiles[i],h,lw,part = 'gas',/ascii
    lw = lw*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc
    IF KEYWORD_SET(rhocut) THEN  yhist = histogram(lw[where(g.dens gt rhocut)]/lw_comp[where(g.dens gt rhocut)],min = -4, max = 4,nbins = 800,locations = xhist)/FLOAT(N_ELEMENTS(lw[where(g.dens gt rhocut)])) ELSE yhist = histogram(lw/lw_comp,min = -4, max = 4,nbins = 800,locations = xhist)/FLOAT(N_ELEMENTS(lw))
    IF i eq 1 then plot, xhist,yhist,xtitle = 'log(Simulated LW Flux / Actual LW Flux)',xrange = [-4,2]
    oplot,xhist,yhist,color = color[i - 1]
stop
ENDFOR
IF KEYWORD_SET(key) THEN legend,keys,color = color,/right,linestyle = fltarr([N_ELEMENTS(lwfiles - 1)]),/bottom,charsize = l_charsize
IF KEYWORD_SET(outplot) THEN device,/close
stop
END
