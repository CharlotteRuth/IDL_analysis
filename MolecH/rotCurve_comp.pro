
pro rotCurve_comp,dt = dt, last=last,first=first
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
sec_per_year = 31556926

;rotCurve,dt=20,last=300,first=20 
;rotCurve,dt = 1,first = 5, last = 5

prefix = "/astro/net/scratch2/christensen/MolecH/12M/"
;prefix = "/astro/net/scratch2/christensen/MolecH/MolecCloud/"

;base = "Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_soft/MW_disk"
;base = "Disk_Iso_1e5/Metal_Cooling_H2_UV_2grad_soft/MW_disk" 
;base = "CoolingCloud/CoolingCloud"
;base = "Metal_H2/12M"
;base2 = "Metal_H2_UV/12M"
;base2 = "Metal_noH2/12M"
;base2 = "CoolingCloud_UV/CoolingCloud"
;base2 = "Disk_Iso_1e5/Metal_Cooling_H2_UV_2grad_soft/MW_disk" 
;base2 = "Disk_Iso_1e5/Metal_Cooling_H2_UV_2grad/MW_disk"
;base2 = "Disk_Iso_1e5/Metal_Cooling/MW_disk"
base = "Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_/MW_disk"
;base2 = "Disk_Iso_1e5/Metal_Cooling_UV_solar_soft/MW_disk"
;base = "Disk_Iso_1e6/MW_disk"
base2 = "Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_soft/MW_disk"
;base2 = "Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_soft_nofit/MW_disk"
;base = "CoolingCloud_UV/CoolingCloud"
;base2 = "CoolingCloud_UV/CoolingCloud"
;base = '1E4R/1E4R'
;base2 = '1E2R/1E2R'

;legend_t = "No UV"
;legend2_t = "UV"
;legend_t = "Small Softening"
;legend2_t = "Original"
;legend_t = "Sim with H2"
;legend2_t =  "Sim with No H2"
legend_t = '1E6 Particles'
legend2_t =  '1E5 Particles'

;legend_t = "Extended Metal Table"
;legend2_t = "Original Metal Table"
;legend_t = '1e4 Particles'
;legend2_t = '1e2 Particles'

outfile = "/astro/net/scratch2/christensen/MolecH/results/UV_solar_soft_H2_noH2_1E5R"
outfile = "/astro/net/scratch2/christensen/MolecH/results/UV_solar_soft_H2_noH2_1E6R"
;outfile = "/astro/net/scratch2/christensen/MolecH/results/UV_solar_soft_exttable_"
;outfile = "/astro/net/scratch2/christensen/MolecH/results/molecCloud"

;msol_per_sysmass = 2.362d5
;kpc_per_syslength = 1d0
 msol_per_sysmass = 1.36e17
 kpc_per_syslength = 1e5
dDelta = 0.003
 ;msol2_per_sysmass = 1.36e17
 ;kpc2_per_syslength = 1e5
 ;msol_per_sysmass = 2.310e15
 ;kpc_per_syslength = 25000.
;sol_per_sysmass = 232510
;pc_per_syslength = 1.0

dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3
vel_convert = 71.0 * (kpc_per_syslength / 1000.0)/2.894405
timeunit = SQRT((cm_per_kpc*kpc_per_syslength)^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year

show_den = 1
show_temp = 1
show_phase = 1
show_frac = 1
show_rot = 1
show_prof = 1
show_sfr = 0
show_column = 0
outprint = 0

;************************************************Specific for Simulation**********************************************
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 2.5
!P.CHARSIZE = 1.25
IF outprint THEN BEGIN
    set_plot,'ps' 
    nbins=50.0
    loadct,0
    colors = [5,80,5,80]
    linestyles = [0,0,2,2]
    thicks = [3,1,3,1]
ENDIF ELSE BEGIN
    set_plot,'x'
    loadct,39
    nbins=50.0
    colors = [240,200,70,140]
    linestyles = [0,0,2,2]
    thicks = [2,2,2,2]
ENDELSE

IF (NOT KEYWORD_SET(last)) THEN last = 10
IF (NOT KEYWORD_SET(dt)) THEN dt = 1
IF (NOT KEYWORD_SET(first)) THEN start = 1 ELSE start = first  ;10/dt + 1
FOR i = start/dt, last/dt DO BEGIN 
    IF outprint THEN loadct,0 ELSE loadct,39
    if (i*dt lt 10) THEN step = '0000'+STRTRIM(i*dt,2) ELSE BEGIN
        if (i*dt lt 100) THEN step = '000'+STRTRIM(i*dt,2) ELSE step = '00'+STRTRIM(dt*i,2)
    ENDELSE
    print,base+"."+step+".HI"
    print,base2+"."+step+".HI"
    rtipsy,prefix+base+"."+step,h,g,d,s
    rtipsy,prefix+base2+"."+step,h_2,g2,d2,s2
    IF (show_column) THEN BEGIN
        readcol,prefix+base+"."+step+'.smoothlength',smoothlengths
        readcol,prefix+base+"."+step+'.column',column
        smoothlengths = smoothlengths[1:h.ngas]*kpc_per_syslength*cm_per_kpc
;*********************************Checking Smoothing Lengths****************************************
;        set_plot,'ps'
;        device,filename=outfile+'smooth_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times
;        y = histogram(alog10(smoothlengths),locations = x,min = 1.5,max = 4,nbins = 50) ;,min = 20,max = 22,nbins = 50)
;        plot,x,y,psym = 10,xtitle = "Log(h) [pc]"
;        oplot,[alog10(MEAN(smoothlengths)),alog10(MEAN(smoothlengths))],[0,500],linestyle = 2
;        device,/close
        IF outprint THEN device,filename=outfile+'smooth_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times ELSE window,8
        smoothl_hist = histogram(alog10(smoothlengths),locations = x,min = 1.5,max = 4,nbins = 50)
        column_hit = histogram(alog10(column),locations = x,min = 1.5,max = 4,nbins = 50)
        mach = smoothlengths/column
        IF outprint THEN device,/close ELSE set_plot,'x'
    ENDIF
;************************************************************************************


    readcol,prefix+base+"."+step+".HI",HI_prime,/silent
    readcol,prefix+base+"."+step+".H2",H2_prime,/silent
    HI_frac = HI_prime[1:h.ngas]
    H2_frac = H2_prime[1:h.ngas]
    total_H = fltarr(N_ELEMENTS(HI_frac))
    total_He = fltarr(N_ELEMENTS(HI_frac))
    total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
    IF(comple[0] ne -1) THEN total_He[comple] = (-0.446*(g[comple].zmetal - 0.1)/0.9 + 0.446)/4.0 ;
    total_H = 1.0 - total_He*4.0 - g.zmetal ;
    rho = alog10(g.dens * dens_convert * total_H) ;H Atoms per CC 
    t = alog10(g.tempg)
    v = g.vx*vel_convert
    HI = HI_prime[1:h.ngas]*g.mass*msol_per_sysmass
    H2 = H2_prime[1:h.ngas]*g.mass*msol_per_sysmass

    readcol,prefix+base2+"."+step+".HI",HI_prime2,/silent
    readcol,prefix+base2+"."+step+".H2",H2_prime2,/silent
    HI2_frac = HI_prime2[1:h_2.ngas]
    H22_frac = H2_prime2[1:h_2.ngas]
    total_H2 = fltarr(N_ELEMENTS(HI2_frac))
    total_He2 = fltarr(N_ELEMENTS(HI2_frac))
    total_He2[where(g2.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g2[where(g2.zmetal le 0.1)].zmetal)/4.0
    IF(comple[0] ne -1) THEN total_He2[comple] = (-0.446*(g2[comple].zmetal - 0.1)/0.9 + 0.446)/4.0 ;
    total_H2 = 1.0 - total_He2*4.0 - g2.zmetal ; 
    rho2 = alog10(g2.dens * dens_convert * total_H2)
    t2 = alog10(g2.tempg)
    v2 = g2.vx*vel_convert
    HI2 = HI_prime2[1:h_2.ngas]*g2.mass*msol_per_sysmass
    H22 = H2_prime2[1:h_2.ngas]*g2.mass*msol_per_sysmass

;********************************* Rotation Curve ***********************************************************
    IF (show_rot) THEN BEGIN
        minv = -300.
        maxv = 300.
        range = maxv - minv
        v_bin = findgen(nbins)*range/(nbins) - range/2.0
        y_H2 = weighted_HISTOGRAM(v,input=H2,binsize=range/nbins,min=minv,max=maxv)
        y_HI = weighted_HISTOGRAM(v,input=HI,binsize=range/nbins,min=minv,max=maxv)
        y_HI2 = weighted_HISTOGRAM(v2,input=HI2,binsize=range/nbins,min=minv,max=maxv)
        y_H22 = weighted_HISTOGRAM(v2,input=H22,binsize=range/nbins,min=minv,max=maxv)
        y_H = weighted_HISTOGRAM(v,input=g.mass*0.765,binsize=range/nbins,min=minv,max=maxv)
        y_H_2 = weighted_HISTOGRAM(v2,input=g2.mass*0.765,binsize=range/nbins,min=minv,max=maxv)
        IF outprint THEN device,filename=outfile+'line_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times ELSE window,2     
        plot,v_bin,y_HI,psym=10,xtitle="Velocity",ytitle="Gas Mass",title='"Line Profiles" -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years'
        oplot,v_bin,y_HI,linestyle=linestyles[0],psym=10,color=colors[0],thick = thicks[0]
        oplot,v_bin,y_HI2,linestyle=linestyles[1],psym=10,color=colors[1],thick = thicks[1]
        oplot,v_bin,y_H2,linestyle=linestyles[2],psym=10,color=colors[2],thick = thicks[2]
        oplot,v_bin,y_H22,linestyle=linestyles[3],psym=10,color=colors[3],thick = thicks[3]
        legend,["HI, "+legend_t,"HI, "+legend2_t,"H2, "+legend_t,"H2, "+legend2_t],linestyle=linestyles,color=colors,thick = thicks
        IF outprint THEN device,/close
    ENDIF

;********************************* Radial Profile ***********************************************************
    IF (show_prof) THEN BEGIN
        minr = 0
        maxr = 50.
        range = maxr - minr
        r = sqrt(g.x*g.x + g.y*g.y)*kpc_per_syslength
        r2 = sqrt(g2.x*g2.x + g2.y*g2.y)*kpc_per_syslength
        r_bin = findgen(nbins)*range/(nbins)
        dr = range/(nbins)
        area = (r_bin + dr)*(r_bin + dr) - r_bin*r_bin
        y_HI = weighted_HISTOGRAM(r,input=HI,binsize=range/nbins,min=minr,max=maxr)
        y_H2 = weighted_HISTOGRAM(r,input=H2,binsize=range/nbins,min=minr,max=maxr)
        y_HI2 = weighted_HISTOGRAM(r2,input=HI2,binsize=range/nbins,min=minr,max=maxr)
        y_H22 = weighted_HISTOGRAM(r2,input=H22,binsize=range/nbins,min=minr,max=maxr)
        y_H = weighted_HISTOGRAM(r,input=g.mass*0.765,binsize=range/nbins,min=minr,max=maxr)
        y_H_2 = weighted_HISTOGRAM(r2,input=g2.mass*0.765,binsize=range/nbins,min=minr,max=maxr)
        IF outprint THEN device,filename=outfile+'profile_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times ELSE window,0
        plot,r_bin,y_HI/area,psym=10,xtitle="Radius",ytitle="Gas Mass Surface Density",title='Gas Surface Density Profile -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years',/ylog,yrange = [1e2,3e7]
        oplot,r_bin,y_HI/area,psym=10,linestyle=linestyles[0],color=colors[0],thick = thicks[0]
        oplot,r_bin,y_HI2/area,psym=10,linestyle=linestyles[1],color=colors[1],thick = thicks[1]
        oplot,r_bin,y_H2/area,psym=10,linestyle=linestyles[2],color=colors[2],thick = thicks[2]
        oplot,r_bin,y_H22/area,psym=10,linestyle=linestyles[3],color=colors[3],thick = thicks[3]
        legend,["HI, "+legend_t,"HI, "+legend2_t,"H2, "+legend_t,"H2, "+legend2_t],linestyle=linestyles,color=colors,thick = thicks,/right
        IF outprint THEN device,/close
        stop
    ENDIF

;********************************* Temperature Distro ***********************************************************
    IF (show_temp) THEN BEGIN
        mint= 1
        maxt = 4.5
        range = maxt - mint
        t_bin = findgen(nbins)*range/(nbins) + mint
        y_H2 = weighted_HISTOGRAM(t,input=H2,binsize=range/nbins,min=mint,max=maxt)
        y_HI = weighted_HISTOGRAM(t,input=HI,binsize=range/nbins,min=mint,max=maxt)
        y_HI2 = weighted_HISTOGRAM(t2,input=HI2,binsize=range/nbins,min=mint,max=maxt)
        y_H22 = weighted_HISTOGRAM(t2,input=H22,binsize=range/nbins,min=mint,max=maxt)
        y_H = weighted_HISTOGRAM(t,input=g.mass*0.765,binsize=range/nbins,min=mint,max=maxt)
        y_H_2 = weighted_HISTOGRAM(t2,input=g2.mass*0.765,binsize=range/nbins,min=mint,max=maxt)
        IF outprint THEN device,filename=outfile+'temp_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times ELSE window,1      
        plot,t_bin,y_HI,psym=10,xtitle="Temperature",ytitle="Gas Mass",title='Temperature Distribution -- '  + strtrim(i*dt*dDelta*timeunit,2) + ' years';,yrange=[0,20000]
        oplot,t_bin,y_HI,psym=10,linestyle=linestyles[0],color=colors[0],thick = thicks[0]
        oplot,t_bin,y_HI2,psym=10,linestyle=linestyles[1],color=colors[1],thick = thicks[1]
        oplot,t_bin,y_H2,psym=10,linestyle=linestyles[2],color=colors[2],thick = thicks[2]
        oplot,t_bin,y_H22,psym=10,linestyle=linestyles[3],color=colors[3],thick = thicks[3]
;        oplot,t_bin,y_H,linestyle=0,psym=10,color=100
;        oplot,t_bin,y_H_2,linestyle=2,psym=10,color=160
        legend,["HI, "+legend_t,"HI, "+legend2_t,"H2, "+legend_t,"H2, "+legend2_t],linestyle=linestyles,color=colors,thick = thicks
        IF outprint THEN device,/close        
    ENDIF
 
;********************************* Density Distro ***********************************************************
    IF (show_den) THEN BEGIN    
        minrho = -5
        maxrho = 4
        range = maxrho - minrho
        rho_bin = findgen(nbins)*range/(nbins) + minrho
        y_H2 = weighted_HISTOGRAM(rho,input=H2,binsize=range/nbins,min=minrho,max=maxrho)
        y_HI = weighted_HISTOGRAM(rho,input=HI,binsize=range/nbins,min=minrho,max=maxrho)
        y_HI2 = weighted_HISTOGRAM(rho2,input=HI2,binsize=range/nbins,min=minrho,max=maxrho)
        y_H22 = weighted_HISTOGRAM(rho2,input=H22,binsize=range/nbins,min=minrho,max=maxrho)
        y_H = weighted_HISTOGRAM(rho,input=g.mass*0.765,binsize=range/nbins,min=minrho,max=maxrho)
        y_H_2 = weighted_HISTOGRAM(rho2,input=g2.mass*0.765,binsize=range/nbins,min=minrho,max=maxrho)
        IF outprint THEN device,filename=outfile+'dens_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times ELSE window,3
        plot,rho_bin,y_HI,psym=10,xtitle="Density",ytitle="Gas Mass",title='Density Distribution -- '  + strtrim(i*dt*dDelta*timeunit,2) + ' years'
        oplot,rho_bin,y_HI,psym=10,linestyle=linestyles[0],color=colors[0],thick = thicks[0]
         oplot,rho_bin,y_HI2,psym=10,linestyle=linestyles[1],color=colors[1],thick = thicks[1]
        oplot,rho_bin,y_H2,psym=10,linestyle=linestyles[2],color=colors[2],thick = thicks[2]
        oplot,rho_bin,y_H22,psym=10,linestyle=linestyles[3],color=colors[3],thick = thicks[3]
        legend,["HI, "+legend_t,"HI, "+legend2_t,"H2, "+legend_t,"H2, "+legend2_t],linestyle=linestyles,color=colors,thick = thicks
        IF outprint THEN device,/close
    ENDIF

;********************************* H2 Fraction ***********************************************************
    IF (show_frac) THEN BEGIN
        IF outprint THEN set_plot,'ps' ELSE window,6
        IF outprint THEN device,filename=outfile+'frac_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times
;        plot,g.dens*dens_convert*total_H,H2_frac/(total_H/2.0),yrange=[0.0,1.1],xrange=[1,1e4],/xlog,psym=3,xstyle=1,ystyle=1,xtitle=textoidl('n_H (cm^{-3})'),ytitle=textoidl('f_{H_2}'),symsize = 0.5,title='Fraction H2 -- '+strtrim(step);/((2.0*H2_frac + HI_frac)/2)
;        oplot,g.dens*dens_convert*total_H,H2_frac/(total_H/2.0),psym=3,color=70,symsize =0.5;/((2.0*H2_frac + HI_frac)/2)
;       oplot,[1e-4,1e4],[1,1]
;stop
        contour_plus,alog10(g.dens*dens_convert*total_H),H2_frac/(total_H/2.0),threshold=3000,nlevels=10,XBINSIZE=0.1,YBINSIZE=0.05,YMIN=0,YMAX=1.0, $
                   xtitle = textoidl("LOG(n_H) (cm^{-3})"),ytitle = textoidl('f_{H_2}'),xrange = [-1,4],yrange = [0, 1.1],title='Fraction H2 -- '+ strtrim(i*dt*dDelta*timeunit,2) + ' years'
       oplot,[-4,4],[1,1]
;        oplot,g2.dens*dens_convert*total_H2 ,H22_frac/(total_H2/2.0),psym=3,color=140,symsize = 0.5
;        legend,[legend_t,legend2_t],linestyle=[1,1],color=[100,140],/right,/bottom
        IF outprint THEN device,/close
    ENDIF


   IF (show_column) THEN BEGIN
        IF outprint THEN set_plot,'ps' ELSE window,7
        IF outprint THEN device,filename=outfile+'frac_columnMass'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times
;        plot,g.dens*dens_convert*total_H,H2_frac/(total_H/2.0),yrange=[0.0,1.1],xrange=[1,1e4],/xlog,psym=3,xstyle=1,ystyle=1,xtitle=textoidl('n_H (cm^{-3})'),ytitle=textoidl('f_{H_2}'),symsize = 0.5,title='Fraction H2 -- '+strtrim(step);/((2.0*H2_frac + HI_frac)/2)
;        oplot,g.dens*dens_convert*total_H,H2_frac/(total_H/2.0),psym=3,color=70,symsize =0.5;/((2.0*H2_frac + HI_frac)/2)
;       oplot,[1e-4,1e4],[1,1]
;stop
;        contour_plus,alog10(g.dens*dens_convert*total_H*smoothlengths),H2_frac/(total_H/2.0),threshold=3000,nlevels=10,XBINSIZE=0.1,YBINSIZE=0.05,YMIN=0,YMAX=1.0, $
;                   xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [0, 1.1],title='Fraction H2 -- '+ strtrim(i*dt*dDelta*timeunit,2) + ' years'
;       oplot,[-4,4],[1,1]

        contour_plus,alog10(g.dens*dens_convert*total_H*column),H2_frac/(total_H/2.0),threshold=3000,nlevels=10,XBINSIZE=0.1,YBINSIZE=0.05,YMIN=0,YMAX=1.0, $
                   xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [0, 1.1],title='Fraction H2 -- '+ strtrim(i*dt*dDelta*timeunit,2) + ' years'
;       oplot,[-4,4],[1,1]

;        oplot,g2.dens*dens_convert*total_H2 ,H22_frac/(total_H2/2.0),psym=3,color=140,symsize = 0.5
;        legend,[legend_t,legend2_t],linestyle=[1,1],color=[100,140],/right,/bottom
        IF outprint THEN device,/close
    ENDIF

;********************************* SFR Diagram ***********************************************************
    IF (show_sfr) THEN BEGIN
  ;      tform=s.tform*timeunit
  ;      smass=s[ind].mass*msol_per_sysmass
  ;      tform2=s2[ind].tform*timeunit
  ;      smass2=s2[ind].mass*msol_per_sysmass
 ;       set_plot,'ps'
        IF outprint THEN device,filename=outfile+'sfr_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times ELSE window,4
        sfr,s,massu=msol_per_sysmass,time=timeunit,OVERPLOT=0,title = 'SFH'
        sfr,s,massu=msol_per_sysmass,time=timeunit,OVERPLOT=1,linestyle=linestyles[2],color=colors[2],thick = thicks[2]
        sfr,s2,massu=msol_per_sysmass,time=timeunit,OVERPLOT=1,gmass=1,linestyle=linestyles[3],color=colors[3],thick = thicks[3]
        legend,[legend_t,legend2_t],linestyle=linestyles[2:3],color=colors[2:3],thick = thicks[2:3]
        IF outprint THEN device,/close

    ENDIF

;********************************* Phase Diagram ***********************************************************
    IF (show_phase) THEN BEGIN
        loadct,39
        IF outprint THEN device,filename=outfile+'phase_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times ELSE window,5
;        plot,g.dens*dens_convert,g.tempg,psym = 3,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",xrange = [10^(-4),10^(4.5)],yrange = [1e1, 1e7],/xlog,/ylog,title='Phase Diagram -- '+strtrim(step)
;        oplot,g.dens*dens_convert,g.tempg,psym = 3,color=240
;        oplot,g2.dens*dens_convert,g2.tempg,psym = 3,color=200
;        legend,[legend_t,legend2_t],linestyle=[0,0],color=[240,200]

        colors_pd = FIX(H2_frac/(total_H/2.0)*206)+50
        plot,g.dens*dens_convert*total_H,g.tempg,psym = 3,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",xrange = [10^(-4),10^(4.5)],yrange = [1e1, 1e7],/xlog,/ylog,title='Phase Diagram '+legend_t+' -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years'
        for ii=0LL,LONG64(N_ELEMENTS(t)) - 1 do oplot,[g[ii].dens*dens_convert*total_H[ii],g[ii].dens*dens_convert*total_H[ii]],[g[ii].tempg,g[ii].tempg],psym = 3,color=colors_pd[ii]
;        contour_plus,alog10(g.dens*dens_convert),alog10(g.tempg),weight = H2_frac/(total_H/2.0),threshold = 0.1,xtitle = textoidl("LOG(n_H) (cm^{-3})"),ytitle = "LOG(T) (K)",xrange = [-4,4.5],yrange = [1, 7],title='Phase Diagram -- '+strtrim(step),levels = levels
        IF outprint THEN device,/close

;        IF outprint THEN device,filename=outfile+'phase2_'+'_'+strtrim(step)+'.eps',/color,bits_per_pixel= 8,/times ELSE window,7
;        colors_pd = FIX(H22_frac/(total_H2/2.0)*206)+50
;        plot,g2.dens*dens_convert*total_H2,g2.tempg,psym = 3,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",xrange = [10^(-4),10^(4.5)],yrange = [1e1, 1e7],/xlog,/ylog,title='Phase Diagram ' +legend2_t+' -- '+ strtrim(i*dt*dDelta*timeunit,2) + ' years'
;        for ii=0LL,LONG64(N_ELEMENTS(g2.dens)) - 1 do oplot,[g2[ii].dens*dens_convert*total_H2[ii],g2[ii].dens*dens_convert*total_H2[ii]],[g2[ii].tempg,g2[ii].tempg],psym = 3,color=colors_pd[ii]
 ;      contour_plus,alog10(g2.dens*dens_convert),alog10(g2.tempg),weight = H22_frac/(total_H2/2.0),threshold = 0.1,xtitle = textoidl("LOG(n_H) (cm^{-3})"),ytitle = "LOG(T) (K)",xrange = [-4,4.5],yrange = [1, 7],title='Phase Diagram -- '+strtrim(step),levels = levels
;        IF outprint THEN device,/close
    ENDIF



    print,'Temperature 1: ',MEAN(g.tempg),MINMAX(g.tempg)
    print,'Temperature 2: ',MEAN(g2.tempg),MINMAX(g2.tempg)
    print,'Density 1: ',MEAN(rho),MINMAX(rho)
    print,'Density 2: ',MEAN(rho2),MINMAX(rho2)
    print,'H2 1: ',MEAN(H2_frac/(total_H/2.0)),MINMAX(H2_frac/(total_H/2.0)),TOTAL(MEAN(H2_frac/(total_H/2.0)))
    print,'H2 2: ',MEAN(H22_frac/(total_H/2.0)),MINMAX(H22_frac/(total_H/2.0)),TOTAL(MEAN(H22_frac/(total_H/2.0)))
    print,' '
    wait,0.5
    stop
endfor
stop
end
