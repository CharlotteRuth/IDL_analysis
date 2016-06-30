;.r /astro/users/christensen/code/MolecH/hist_2d_weighted.pro
;.r /astro/users/christensen/code/MolecH/contour_plus.pro
;.r /astro/users/christensen/code/MolecH/phase_diagram.pro

;phase_diagram,'h516.cosmo25cmb.3072g3bwK.00512.1.std',25000.0,2.310e15

;phase_diagram,'/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Standard/MW_disk.00001','/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Standard/phase_dia.ps'
;phase_diagram,'/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling/MW_disk.00001','/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling/phase_dia.ps'
;phase_diagram,'/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_H2/MW_disk.00001','/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_H2/phase_dia.ps'
;phase_diagram,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.768g1HBWK/steps/h516.cosmo25cmb.768g1HBWK.00408.dir/h516.cosmo25cmb.768g1HBWK.00408',25000.,2.310e15,outfile= '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.768g1HBWK/steps/h516.cosmo25cmb.768g1HBWK.00408.dir/phase_dia.ps'
;phase_diagram,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.768g1MBWK_jill/steps/h516.cosmo25cmb.768g1MBWK.00408.dir/h516.cosmo25cmb.768g1MBWK.00408',25000.,2.310e15,outfile= '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.768g1MBWK_jill/steps/h516.cosmo25cmb.768g1MBWK.00408.dir/phase_dia.ps'

;-0.776722      7.22328

PRO phase_diagram_h516
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/"
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/"
base = "h516.cosmo25cmb.1536g3HBWK"
msol_per_sysmass       = 2.310e15
kpc_per_syslength        = 25000.
last = 300; 168;512
dt = 12
start = 0; 168;512 
steps = [15,30,45,60,75,90,105,120,135,150,165,180,195,210,216,228,240,252,264,276,288,300,312,324,336,348,360,372]
;steps = [216,228,240,252,264,276,288,300]
;steps = [288,300]



nsteps = last/dt
nsteps = N_ELEMENTS(steps) - 1
cd,prefix
FOR i = start/dt, nsteps DO BEGIN 
    step = i*dt
    step = steps[i]
    if (step lt 10) THEN step = '0000'+STRTRIM(step,2) ELSE BEGIN
        if (step lt 100) THEN step = '000'+STRTRIM(step,2) ELSE step = '00'+STRTRIM(step,2)
    ENDELSE
    filename = base+"/steps/"+base+"."+step+".dir/"+base + "." + step
    print,filename
    phase_diagram,filename,base+'/'+base + '.starlog',kpc_per_syslength,msol_per_sysmass,/molecH,/outfile,title = "Step "+step,maxsfr = 0.2
ENDFOR

END

PRO phase_diagram_h603
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/"
base = "h603.cosmo50cmb.1536g3HBWK"
msol_per_sysmass       = 1.84793e16
kpc_per_syslength        = 50000.
last = 300; 168;512
dt = 12
start = 0; 168;512 
steps = [37,46,48,60,61,72,84,96,108,120,128,132,144,156,168,180,192,204,216,227,228,240,252,264,271,276,288,300]
steps = [312,324,336,348,360,372,384,396]
;steps = [216,228,240,252,264,276,288,300]
;steps = [288,300]

nsteps = last/dt
nsteps = N_ELEMENTS(steps) - 1
cd,prefix
FOR i = start/dt, nsteps DO BEGIN 
    step = i*dt
    step = steps[i]
    if (step lt 10) THEN step = '0000'+STRTRIM(step,2) ELSE BEGIN
        if (step lt 100) THEN step = '000'+STRTRIM(step,2) ELSE step = '00'+STRTRIM(step,2)
    ENDELSE
    filename = base+"/steps/"+base+"."+step+".dir/"+base + "." + step
    print,filename
    phase_diagram,filename,base+'/'+base + '.starlog',kpc_per_syslength,msol_per_sysmass,/molecH,/outfile,title = "Step "+step,maxsfr = 6
ENDFOR

END

pro phase_diagram,infile,slfile,dunit,massunit,outfile = outfile,molecH = molecH,title = title,maxsfr = maxsfr,redshift = redshift,smooth = smooth
;This will plot the phase diagram for the gas in a given simulation
cm_per_kpc = 3.0857d21
gm_per_solarmass = 1.989d33
amu_per_gm = 6.022142d23
molec_weight = (0.76*1 + 0.24*4.0)
dens_convert =  massunit*gm_per_solarmass*amu_per_gm/molec_weight/cm_per_kpc^3/dunit^3
k = 1.38d-16
mp = 1.67d-24

close,/all
IF NOT KEYWORD_SET (maxsfr) THEN maxsfr = 1
IF KEYWORD_SET(outfile) THEN BEGIN
    set_plot,'ps' 
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.MARGIN = [4,2]
    !Y.THICK=4
    !p.charsize=0.70
 ;   !x.charsize=2.25
 ;   !y.charsize=2.25
;    !p.font=0 
    !X.STYLE = 1
    !Y.STYLE = 1
;    loadct,39
ENDIF ELSE BEGIN
    set_plot,'x'
    !P.CHARTHICK=1.5 ;1.25
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5
    !X.STYLE = 1
    !Y.STYLE = 1  
    loadct,39
ENDELSE
!p.multi = [0,2,2]
;multiplot

IF (KEYWORD_SET(redshift)) THEN BEGIN ;redshift  = 0
    ind = STRSPLIT(infile,'/')
    base = STRMID(infile,0,ind[N_ELEMENTS(ind) - 1])
    spawn,'ls '+base+'*.AHF_halos',result
    segments = strsplit(result,'.',/extract)
    redshift = segments[N_ELEMENTS(segments) - 3] + '.'  +segments[N_ELEMENTS(segments) - 2]
    redshift = STRMID(redshift,1,STRLEN(redshift) - 1)
ENDIF ELSE redshift = 0
a = 1.0/(1.0 + redshift)

rtipsy,infile,h,g,d,s
a = h.time
;g = g[where(SQRT(g.x*g.x + g.y*g.y) LT 8.62309e-05/2 AND 2e-06 GT ABS(g.z))]
g.x = g.x * dunit*a
g.y = g.y * dunit*a
g.z = g.z * dunit*a
g.mass = g.mass*massunit
IF h.nstar ne 0 THEN BEGIN
    s.x = s.x * dunit*a
    s.y = s.y * dunit*a
    s.z = s.z * dunit*a
    s.mass = s.mass*massunit
    s.tform = s.tform*SQRT((dunit*3.086d21)^3/(6.67d-8*massunit*1.99d33))/(3600.*24.*365.24)
ENDIF ELSE s = {mass: 0.,x: 0.,y : 0., z:0.,vx:0.,vy:0.,vz:0.,metals:0.,tform:0.,eps: 0.,phi: 0.}
maxtime = 8.5e9

total_H = fltarr(h.ngas)
total_He = fltarr(h.ngas)
total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
IF(comple[0] ne -1) THEN total_He[comple] = (-0.446*(g[comple].zmetal - 0.1)/0.9 + 0.446)/4.0 ;
total_H = 1.0 - total_He*4.0 - g.zmetal
IF KEYWORD_SET(molecH) THEN BEGIN
    readarr,infile+".HI",h,HI_prime,/ascii
    readarr,infile+".H2",h,H2_prime,/ascii
    HI_frac = HI_prime[0:h.ngas - 1]
    H2_frac = H2_prime[0:h.ngas - 1]   
ENDIF ELSE H2_frac = g.dens*0
if Keyword_set(smooth) then BEGIN
    readarr,smooth,h,dens,/ascii,part = 'gas'
    density = dens*dens_convert*total_H/a^3
ENDIF  ELSE density = g.dens*dens_convert*total_H/a^3

;plot,g.dens,g.tempg,psym = 1,xtitle = 'Density/Critical Density',ytitle = 'Temperature (K)'
;plot,g.tempg,g.tempg*density*k/mp,psym = 1,xtitle = 'Temperature',ytitle = 'Pressure',title = 'Phase Diagram'


;x = density
;y = g.tempg
;xrange = MAX(x) - MIN(X)
;yrange = MAX(y) - MIN(y)
;contour_plus,x,y,xbinsize = xrange/30.,ybinsize = yrange/30.,threshold = 5,nlevels = 5,xtitle = 'Density/Critical Density',ytitle = 'Temperature (K)',title = 'Phase Diagram',psym = 1
;stop
;x = ALOG10(x)
;y = ALOG10(y)
;xrange = MAX(x) - MIN(X)
;yrange = MAX(y) - MIN(y)
IF NOT KEYWORD_SET(outfile) THEN BEGIN
    plot,density,g.tempg,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",psym = 3,yrange = [1,1e6],xrange = [1e-6,1e3],/ylog,/xlog,title = title,ystyle = 1
;contour_plus,x,y,xbinsize = xrange/40.,ybinsize = yrange/40.,threshold = 2.5e5,nlevels = 50,xtitle = 'LOG(Density/Critical Density)',ytitle = 'LOG(Temperature K)',title = 'Phase Diagram',psym = 3,yrange = [1,6],xrange = [-6,3],c_color = ((INDGEN(50)+1)*200/50)+55 
    IF KEYWORD_SET(molecH) THEN BEGIN
        logfrac = alog10(H2_frac/(total_H/2.0))
        minlogfrac = -6         ;MIN(alog10(H2_frac/(total_H/2.0)))
        maxlogfrac = -0.0518840 ;MAX(alog10(H2_frac/(total_H/2.0))) 
        colors_pd = FIX((logfrac - minlogfrac)/(maxlogfrac - minlogfrac)*254)
        colors_pd[where(logfrac lt minlogfrac)] = 255
        for ii=0LL,LONG64(N_ELEMENTS(y)) - 1 do oplot,[density[ii],density[ii]],[g[ii].tempg,g[ii].tempg],psym = 3,color=colors_pd[ii]
;        colorbar,maxrange = maxlogfrac,minrange = minlogfrac,format = '(F10.1)'
    ENDIF
;oplot,[2,2],[1,4.5]
;oplot,[0,0],[1,4.5]
;oplot,[-6,3],[4.0,4.0]

;    multiplot
    IF KEYWORD_SET(molecH) THEN data = rstarlog(slfile,/MOLECULARH) ELSE data = rstarlog(slfile)
    data.timeform = data.timeform*SQRT((dunit*3.086d21)^3/(6.67d-8*massunit*1.99d33))/(3600.*24.*365.24)
    data = data[where(data.timeform le MAX(s.tform))]
    rhoform = data.rhoform*dens_convert*total_H
    plot,density,g.tempg,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",psym = 3,yrange = [1,1e6],xrange = [1e-6,1e3],/ylog,/xlog
    oplot,rhoform,data.tempform,psym = 3
    mintime =0             
    colors_pd = FIX((data.timeform - mintime)/(maxtime - mintime)*254)
    IF (where(data.timeform lt mintime))[0] ne -1 THEN colors_pd[where(data.timeform lt mintime)] = 1
    for ii=0LL,LONG64(N_ELEMENTS(rhoform)) - 1 do oplot,[rhoform[ii],rhoform[ii]],[data[ii].tempform,data[ii].tempform],psym = 3,color=colors_pd[ii]

;    multiplot
    plot,density,H2_frac/(total_H/2.0),xtitle = textoidl("n_H (cm^{-3})"),ytitle = textoidl('f_{H_2}'),xrange = [1e0,1e3],yrange = [0, 1.1],psym=3,/xlog
    oplot,[-4,4],[1,1]

;    multiplot
    sfr,s,massunit = 1,timeunit = 1,binsize = 1e8,yrange = [0,maxsfr],xstyle = 1,xrange = [0,maxtime/1e9]
stop
multiplot, /default
!p.multi = 0

    ntemp = histogram(alog10(data.tempform),nbins = 500,locations = x,min = alog10(20),max = alog10(2e4))
    ntempcum = ntemp
    for i = 1, N_ELEMENTS(x)-1 DO BEGIN
        ntempcum[i] = 1.0*TOTAL(ntemp[0:i])
    ENDFOR
    plot,x,ntempcum,xtitle = 'Log Temperature',ytitle = 'Fraction of Stars at or below Temperature'
stop

;set_plot,'ps'
;device,filename = infile+'.SFtemp.eps',/color,bits_per_pixel= 8,/times
;    plot,x,ntempcum,xtitle = 'Log Temperature',ytitle = 'Fraction of Stars at or below Temperature'
;device,/close
;stop

;set_plot,'ps'
;device,filename = infile+'.SFphase.eps',/color,bits_per_pixel= 8,/times
    x = alog10(density)
    y = alog10(g.tempg)
    xrange = 10 + 1
    yrange = 6 + 1
;    contour_plus,x,y,xbinsize = xrange/100.,ybinsize = yrange/100.,threshold = 400000,nlevels = 4,levels = [400000,600000,800000,1000000],xtitle = 'Log Density',ytitle = 'Log Temperature [K]',title = 'Phase Diagram',psym = 3,yrange = [1,7],xrange = [-6,4]
;    oplot,alog10(rhoform),alog10(data.tempform),psym = 3,color = 60
    plot,density,g.tempg,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",psym = 3,yrange = [1,1e7],xrange = [1e-6,1e4],/ylog,/xlog
    oplot,rhoform,data.tempforexit
m,psym = 3,color = 60
;device,/close
    stop

ENDIF ELSE BEGIN
    device,filename = infile+'.pd.eps',/color,bits_per_pixel= 8,/times,ysize=6,xsize=8,/inch
 ;   print,'Phase D',infile+'.pd.eps'
    plot,density,g.tempg,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",psym = 3,yrange = [1,1e6],xrange = [1e-6,1e3],/ylog,/xlog,title = title
 ;   print,'plotted'
;contour_plus,x,y,xbinsize = xrange/40.,ybinsize = yrange/40.,threshold = 2.5e5,nlevels = 50,xtitle = 'LOG(Density/Critical Density)',ytitle = 'LOG(Temperature K)',title = 'Phase Diagram',psym = 3,yrange = [1,6],xrange = [-6,3],c_color = ((INDGEN(50)+1)*200/50)+55 
    IF KEYWORD_SET(molecH) THEN BEGIN
        logfrac = alog10(H2_frac/(total_H/2.0))
        minlogfrac = -6         ;MIN(alog10(H2_frac/(total_H/2.0)))
        maxlogfrac = -0.0518840 ;MAX(alog10(H2_frac/(total_H/2.0))) 
        colors_pd = FIX((logfrac - minlogfrac)/(maxlogfrac - minlogfrac)*254)
        colors_pd[where(logfrac lt minlogfrac)] = 1
        for ii=0LL,LONG64(N_ELEMENTS(y)) - 1 do oplot,[density[ii],density[ii]],[g[ii].tempg,g[ii].tempg],psym = 3,color=colors_pd[ii]
;        colorbar,maxrange = maxlogfrac,minrange = minlogfrac,format = '(F10.1)'
    ENDIF

;    multiplot
 ;   print,'Star log Plot'
    plot,density,g.tempg,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",psym = 3,yrange = [1,1e6],xrange = [1e-6,1e3],/ylog,/xlog
    data = rstarlog(slfile)
    data.timeform = data.timeform*SQRT((dunit*3.086d21)^3/(6.67d-8*massunit*1.99d33))/(3600.*24.*365.24)
    IF (where(data.timeform le MAX(s.tform)))[0] ne -1 THEN BEGIN
        data = data[where(data.timeform le MAX(s.tform))]
        rhoform = data.rhoform*dens_convert*total_H
        oplot,rhoform,data.tempform,psym = 3
        mintime =0             
        colors_pd = FIX((data.timeform - mintime)/(maxtime - mintime)*254)
        IF (where(data.timeform lt mintime))[0] ne -1 THEN colors_pd[where(data.timeform lt mintime)] = 1
        for ii=0LL,LONG64(N_ELEMENTS(rhoform)) - 1 do oplot,[rhoform[ii],rhoform[ii]],[data[ii].tempform,data[ii].tempform],psym = 3,color=colors_pd[ii]
    ENDIF

;    multiplot
 ;   print,'abundance'
    plot,density,H2_frac/(total_H/2.0),xtitle = textoidl("n_H (cm^{-3})"),ytitle = textoidl('f_{H_2}'),xrange = [1e0,1e3],yrange = [0, 1.1],psym=3,/xlog
    oplot,[-4,4],[1,1]

;    multiplot
;    print,'SFH'
    sfr,s,massunit = 1,timeunit = 1,binsize = 1e8,yrange = [0,maxsfr],xstyle = 1,xrange = [0,maxtime/1e9]
    device,/close
    
    ntemp = histogram(alog10(data.tempform),nbins = 100,x =locations,min = alog10(20),max = alog10(2e4))
    ntempcum = ntemp
    for i = 1, N_ELEMENTS(histogram)-1 DO ntempcum[i] = TOTAL(ntemp[0:i])
    plot,locations,ntempcum,xtitle = 'Temperature',ytitle = 'Fraction of Stars at or below Temperature'
    stop
ENDELSE
print,'Plots done'
;contour_plus,x,y,xbinsize = xrange/30.,ybinsize = yrange/30.,threshold = 15,nlevels = 5,xtitle = 'LOG(Density/Critical Density)',ytitle = 'LOG(Temperature K)',title = 'Phase Diagram',psym = 3,/loglevel
;stop

;multiplot, /reset
;multiplot, /default
!p.multi = 0
stop
END
