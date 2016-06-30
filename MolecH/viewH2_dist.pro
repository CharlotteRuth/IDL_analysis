PRO VIEWH2_DIST,dt = dt, last=last,first=first
;rotCurve,dt=20,last=300,first=20 
;rotCurve,dt = 1,first = 5, last = 5

;prefix = "/astro/net/scratch2/christensen/MolecH/12M/"
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/"
cosmo = 1
halo = 1

;base = "Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_HIshield/MW_disk"
;base = "Disk_Iso_1e6/MW_disk"
;base = "Disk_Iso_1e6g2/MW_disk"
;base = "h516.cosmo25cmb.768g5bwK/h516.cosmo25cmb.768g5bwK"
;base = "h603.cosmo50cmb.768g5bwK/h603.cosmo50cmb.768g5bwK"
base = "h516.cosmo25cmb.1536g2HBWK/h516.cosmo25cmb.1536g2HBWK"

;legend_t = '1E5 Particles'
;legend_t = '1E6 Particles'
;legend_t = '1E6 Particles, g2'
;legend_t = 'h516'
legend_t = 'h516g2'
;legend_t = 'h603'

;outfile = "/astro/net/scratch2/christensen/MolecH/results/UV_solar_soft_1E6Rg2"
;outfile = "/astro/net/scratch2/christensen/MolecH/results/UV_solar_HIshield_1E5R"
;outfile = "/astro/net/scratch2/christensen/MolecH/results/h516"
outfile = "/astro/net/scratch2/christensen/MolecH/results/h516g2"
;outfile = "/astro/net/scratch2/christensen/MolecH/results/h603"

; msol_per_sysmass = 1.36e17
; kpc_per_syslength = 1e5
 msol_per_sysmass       = 2.310e15
 kpc_per_syslength        = 25000.

outprint = 1 ;true if you want to make a hard copy

!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 2.5
!P.CHARSIZE = 1.25
IF outprint THEN BEGIN
    set_plot,'ps' 
    nbins=50.0
    loadct,3
    colors = [5,80]
    linestyles = [0,2]
    thicks = [3,1]
    !P.THICK=4
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.8
    !p.font=0 
ENDIF ELSE BEGIN
    set_plot,'x'
    loadct,22
    nbins=50.0
    colors = [240,200]
    linestyles = [0,2]
    thicks = [2,2]
ENDELSE
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3
vel_convert = 71.0 * (kpc_per_syslength / 1000.0)/2.894405
dDelta = 0.003
sec_per_year = 31556926
timeunit = SQRT((cm_per_kpc*kpc_per_syslength)^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year
IF (NOT KEYWORD_SET(last)) THEN last = 300
IF (NOT KEYWORD_SET(dt)) THEN dt = 1
IF (NOT KEYWORD_SET(first)) THEN start = 1 ELSE start = first  ;10/dt + 1

FOR it = start/dt, last/dt DO BEGIN 
    IF outprint THEN loadct,0 ELSE loadct,39
    if (it*dt lt 10) THEN step = '0000'+STRTRIM(it*dt,2) ELSE BEGIN
        if (it*dt lt 100) THEN step = '000'+STRTRIM(it*dt,2) ELSE step = '00'+STRTRIM(dt*it,2)
    ENDELSE

    if cosmo then begin
        filenametip = prefix+base+'.'+step+'.'+STRTRIM(halo,2)+'.std'
        filename = prefix+base+"."+step
        grpfile = filename+'.amiga.grp'
        gtpfile = filename+'.amiga.gtp'
        statfile = filename+'.amiga.stat'
        grp = read_lon_array(grpfile)

        print,filename+".HI"

        readcol,filename+".HI",HI_prime,/silent
        readcol,filename+".H2",H2_prime,/silent
        rtipsy,filename,h,g,d,s

        ind = where(grp eq halo,comp=indcomp)
        inds = ind[where(ind ge h.ngas+h.ndark)]-h.ngas-h.ndark
        indg = ind[where(ind lt h.ngas)]
        indd = ind[where(ind ge h.ngas and ind lt h.ngas+h.ndark)]-h.ngas

        s = s[inds]
        g = g[indg]
        d = d[indd] 
        
       ; rtipsy,filenametip,h,g,d,s
        HI_frac = HI_prime[1:h.ngas]
        HI_frac = HI_frac[indg]
        H2_frac = H2_prime[1:h.ngas]
        H2_frac = H2_frac[indg]
rtipsy, gtpfile, h1,g1,d1,s1

cx= s1[0].x
cy= s1[0].y
cz= s1[0].z

g.x = g.x - cx
g.y = g.y - cy
g.z = g.z - cz
d.x = d.x - cx
d.y = d.y - cy
d.z = d.z - cz
s.x =  s.x - cx
s.y =  s.y - cy
s.z =  s.z - cz
dist_units=kpc_per_syslength*h.time ; multiplies the distand units, which you fed it, by the expansion 
limit=5*h.time/dist_units            ;use stars within this limit (physical kpc/dist_units)
align,g,d,s,limit

    ENDIF else begin
        print,base+"."+step+".HI"
        
        readcol,prefix+base+"."+step+".HI",HI_prime,/silent
        readcol,prefix+base+"."+step+".H2",H2_prime,/silent
        rtipsy,prefix+base+"."+step,h,g,d,s
        
        HI_frac = HI_prime[1:h.ngas]
        H2_frac = H2_prime[1:h.ngas]
        
    endelse
;*************************************************************

    g.x = g.x*kpc_per_syslength
    g.y = g.y*kpc_per_syslength
    g.z = g.z*kpc_per_syslength
    total_H = fltarr(N_ELEMENTS(HI_frac))
    total_He = fltarr(N_ELEMENTS(HI_frac))
    total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
    IF(comple[0] ne -1) THEN total_He[comple] = (-0.446*(g[comple].zmetal - 0.1)/0.9 + 0.446)/4.0 ;
    total_H = 1.0 - total_He*4.0 - g.zmetal ;
    rho = alog10(g.dens * dens_convert * total_H) ;H Atoms per CC 
    HI = HI_frac*g.mass*msol_per_sysmass
    H2 = H2_frac*g.mass*msol_per_sysmass

    H2_frac_norm = H2_frac/(total_H/2.0)
    
    H2_frac_norm_hist = histogram(alog10(H2_frac_norm),locations = x,nbins = 100)
;    window,0
;    plot,x,H2_frac_norm_hist,psym = 10

    nlevels = 21
    xbinsize = 0.05
    ybinsize = 0.05
    xmax = 12                  ;200
    xmin = -12                 ;-200
    ymax = 12                  ;400
    ymin = -12                 ;-400
    zmax = 12                  ;400
    zmin = -12                 ;-400

    threshold = -5
    maxdist = 0
    levels = (findgen(nlevels)*(maxdist-threshold)/(nlevels - 1.0) + threshold) ;*scale[ct]
;    threshold = MIN(levels)
    print,levels
    ind = where(alog10(H2_frac_norm) gt threshold)
    weights = alog10(H2_frac_norm[ind])
    x = g[ind].x
    y = g[ind].y
    z = g[ind].z

 ;   window,1
    colors = FIX((weights-threshold)/(maxdist - threshold)*205)+50
;    window,0
;    plot,rho,t,psym = 3,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",xrange = [10^(-6),10^(4.5)],yrange = [1e1, 1e7],/xlog,/ylog,title='Phase Diagram -- '+strtrim(step)
;    for i=0LL,LONG64(N_ELEMENTS(t)) - 1 do oplot,[rho[i],rho[i]],[t[i],t[i]],psym = 3,color=colors[i]
  
    device,/color,bits_per_pixel=8,filename=outfile+'logH2dist_point'+step+'.eps',/times,ysize=8,xsize=8,/inch
    loadct,39
    plot,x,y,psym = 3,xrange=[xmin,xmax],yrange = [ymin,ymax],symsize = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]'
    for i=0LL,LONG64(N_ELEMENTS(x)) - 1 do oplot,[x[i],x[i]],[y[i],y[i]],psym = 3,color=colors[i],symsize = 1
    device,/close
;stop
    device,/color,bits_per_pixel=8,filename=outfile+'logH2dist'+step+'.eps',/times,ysize=8,xsize=8,/inch
    loadct,3
    contour_plus,x,y,weight = weights,xmin = xmin, xmax = xmax, ymin = ymin, ymax  = ymax, xbinsize = xbinsize, ybinsize = ybinsize, threshold = threshold,levels = levels,xtitle = 'X [kpc]',ytitle = 'Y [kpc]'
    device,/close
;    window,2
    nlevels = 21
    xbinsize = 0.05
    ybinsize = 0.05
    xmax = 12                  ;200
    xmin = -12                 ;-200
    ymax = 12                  ;400
    ymin = -12                 ;-400
    zmax = 12                  ;400
    zmin = -12                 ;-400

    threshold = 0.05
    maxdist = 1
    levels = (findgen(nlevels)*(maxdist-threshold)/(nlevels - 1.0) + threshold) ;*scale[ct]
;    threshold = MIN(levels)
    print,levels
    ind = where(H2_frac_norm gt threshold)
    weights = H2_frac_norm[ind]
    colors = FIX(H2_frac_norm[ind]*205)+50
    x = g[ind].x
    y = g[ind].y
    z = g[ind].z

    device,/color,bits_per_pixel=8,filename=outfile+'H2dist_point'+step+'.eps',/times,ysize=8,xsize=8,/inch
    loadct,39
    plot,x,y,psym = 3,xrange=[xmin,xmax],yrange = [ymin,ymax],symsize = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]'
    for i=0LL,LONG64(N_ELEMENTS(x)) - 1 do oplot,[x[i],x[i]],[y[i],y[i]],psym = 3,color=colors[i],symsize = 1
    device,/close
;stop
    device,/color,bits_per_pixel=8,filename=outfile+'H2dist'+step+'.eps',/times,ysize=8,xsize=8,/inch
    loadct,3
    contour_plus,x,y,weight = weights,xmin = xmin, xmax = xmax, ymin = ymin, ymax  = ymax, xbinsize = xbinsize, ybinsize = ybinsize, threshold = threshold,levels = levels,xrange=[xmin,xmax],yrange = [ymin,ymax],xtitle = 'X [kpc]',ytitle = 'Y [kpc]'
    device,/close
  ;  stop
ENDFOR
;device,/close
END
