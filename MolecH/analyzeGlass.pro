PRO analyzeGlass,dt = dt, last=last,first=first

prefix = "/astro/net/scratch1/christensen/MolecH/ShockTube/st_glass/"
;base = "glass_highd/glass_highd"
base = "glass_highd_nometal/glass_highd"
outfile = "/astro/net/scratch2/christensen/MolecH/results/"



msol_per_sysmass       = 99278160.
kpc_per_syslength      = 1.0
dDelta                 = 0.0021128488
outprint = 0

!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 2.5
!P.CHARSIZE = 1.25
IF outprint THEN BEGIN
    set_plot,'ps' 
;    nbins=50.0
    loadct,0
;    colors = [5,10]
;    linestyles = [0,2]
;    thicks = [3,2]
;    !P.CHARTHICK=4
;    !X.THICK=4
;    !Y.THICK=4
;    !p.charsize=1.8
;    !p.font=0 
ENDIF ELSE BEGIN
    set_plot,'x'
    loadct,39
;    nbins=50.0
;    colors = [240,240]
;    linestyles = [0,2]
;    thicks = [2,2]
ENDELSE

cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3
vel_convert = 71.0 * (kpc_per_syslength / 1000.0)/2.894405
sec_per_year = 31556926
timeunit = SQRT((cm_per_kpc*kpc_per_syslength)^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year
IF (NOT KEYWORD_SET(last)) THEN last = 300
IF (NOT KEYWORD_SET(dt)) THEN dt = 1
IF (NOT KEYWORD_SET(first)) THEN start = 1 ELSE start = first  ;10/dt + 1

time     = fltarr((last - start)/dt + 1)
meantemp = fltarr((last - start)/dt + 1)
meanH2_frac = fltarr((last - start)/dt + 1)
iarray   = 0

FOR i = start/dt, last/dt DO BEGIN 
    IF outprint THEN loadct,0 ELSE loadct,39
    if (i*dt lt 10) THEN step = '0000'+STRTRIM(i*dt,2) ELSE BEGIN
        if (i*dt lt 100) THEN step = '000'+STRTRIM(i*dt,2) ELSE step = '00'+STRTRIM(dt*i,2)
    ENDELSE
    
    readcol,prefix+base+"."+step+".HI",HI_prime,/silent
    readcol,prefix+base+"."+step+".H2",H2_prime,/silent
    rtipsy,prefix+base+"."+step,h,g,d,s
    
    HI_frac = HI_prime[1:h.ngas]
    H2_frac = H2_prime[1:h.ngas]
    readcol,prefix+base+"."+step+'.mach_sheer',mach
    column = 1.0/mach[1:h.ngas]*kpc_per_syslength*cm_per_kpc
    readcol,prefix+base+"."+step+'.smoothlength',smoothlengths_L
    smoothlengths = smoothlengths_L[1:h.ngas]*kpc_per_syslength*cm_per_kpc    
    total_H = fltarr(N_ELEMENTS(HI_frac))
    total_He = fltarr(N_ELEMENTS(HI_frac))
    total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
    IF(comple[0] ne -1) THEN total_He[comple] = (-0.446*(g[comple].zmetal - 0.1)/0.9 + 0.446)/4.0 ;
    total_H = 1.0 - total_He*4.0 - g.zmetal ;
    rho = alog10(g.dens * dens_convert * total_H) ;H Atoms per CC 
    t = alog10(g.tempg)
    v = g.vx*vel_convert
    HI = HI_frac*g.mass*msol_per_sysmass
    H2 = H2_frac*g.mass*msol_per_sysmass

    time[iarray]     = i*dt*dDelta*timeunit
    meantemp[iarray] = MEAN(g.tempg)
    meanH2_frac[iarray] = MEAN(H2_frac/(total_H/2.0))
    iarray = iarray +1
    window,0
    plot,g.dens*dens_convert*total_H,g.tempg,psym = 3,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",xrange = [10,100],yrange = [1e1, 1e3],/xlog,/ylog,xstyle = 1
 ;   stop
ENDFOR
window,1
plot,time,meantemp,xtitle = 'Time [yrs]',ytitle = 'Mean Temperature [K]'
window,2
plot,time,meanH2_frac,xtitle = 'Time [yrs]',ytitle = 'Fraction of H that is H2'
stop
END
