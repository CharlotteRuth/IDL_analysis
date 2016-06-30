pro rotCurve_12M
prefix = "/astro/net/scratch2/christensen/MolecH/12M/"

;base = "Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_HIshield/MW_disk"
;base = "Disk_Iso_1e6/MW_disk"
base = "Disk_Iso_1e6g2/MW_disk"

;legend_t = '1E6 Particles'
;legend_t = '1E6 Particles, g2'

;outfile ="/astro/net/scratch2/christensen/MolecH/results/UV_solar_soft_1E6Rg2"
;outfile ="/astro/net/scratch2/christensen/MolecH/results/UV_solar_soft_1E6R"
;outfile = "/astro/net/scratch2/christensen/MolecH/results/UV_solar_soft_1E5R"
;outfile = "/astro/net/scratch2/christensen/MolecH/results/UV_solar_soft_1E6R_11M"
;outfile = "/astro/net/scratch2/christensen/MolecH/results/UV_solar_HIshield_1E5R"

;msol_per_sysmass = 1.36e17
;kpc_per_syslength = 1e5

IF (NOT KEYWORD_SET(last)) THEN last = 300
IF (NOT KEYWORD_SET(dt)) THEN dt = 1
IF (NOT KEYWORD_SET(first)) THEN start = 1 ELSE start = first  ;10/dt + 1

FOR i = start/dt, last/dt DO BEGIN 
    if (i*dt lt 10) THEN step = '0000'+STRTRIM(i*dt,2) ELSE BEGIN
        if (i*dt lt 100) THEN step = '000'+STRTRIM(i*dt,2) ELSE step = '00'+STRTRIM(dt*i,2)
    ENDELSE
    filenametip = prefix+base+'.'+step+'.'+STRTRIM(cosmo,2)+'.std'
    filename = prefix+base+"."+step
ENDFOR

end

pro rotCurve_h516
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/"

;base = "h516.cosmo25cmb.768g5bwK/h516.cosmo25cmb.768g5bwK"
base = "h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.00512.dir/h516.cosmo25cmb.1536g2HBWK"
base = "h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir/h516.cosmo25cmb.1536g1MBWK"
base = ["h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.00512.dir/h516.cosmo25cmb.1536g2HBWK","h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir/h516.cosmo25cmb.1536g1MBWK"]
base = ["h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00168.dir/h516.cosmo25cmb.1536g1MBWK", $
        "h516.cosmo25cmb.1536g2HBWK_nocool/steps/h516.cosmo25cmb.1536g2HBWK.00168.dir/h516.cosmo25cmb.1536g2HBWK",$
        "h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.00168.dir/h516.cosmo25cmb.1536g2HBWK"]
base = ["h516.cosmo25cmb.768g1HBWK/steps/h516.cosmo25cmb.768g1HBWK.00408.dir/h516.cosmo25cmb.768g1HBWK", $
        "h516.cosmo25cmb.768g1MBWK_jill/steps/h516.cosmo25cmb.768g1MBWK.00408.dir/h516.cosmo25cmb.768g1MBWK"]
base = ["h516.cosmo25cmb.1536g1MBWK/h516/h516.cosmo25cmb.1536g1MBWK",$
        "h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.00512.dir/h516.cosmo25cmb.1536g2HBWK", $
        "h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00513.dir/h516.cosmo25cmb.1536g3HBWK"]
base = ["h516.cosmo25cmb.3072g3HBWK/steps/h516.cosmo25cmb.3072g3HBWK.00396.00012.dir/h516.cosmo25cmb.3072g3HBWK.00396.00012.halo.1"]

keys = ["h516","h516 H2 no cool","h516 H2 cooling"]
keys = ['h516 H2', 'h516 noH2, M']
keys = ['h516 noH2','h516 H2','h516 H2+SF']
keys = ['With H2']

outfile = "/astro/net/scratch2/christensen/MolecH/results/h516g2cool"
outfile = "/astro/net/scratch2/christensen/MolecH/results/h516g4noH2"
outfile = "/astro/net/scratch2/christensen/MolecH/results/h516_noH2_H2cool_H2"
outfile = "/astro/net/scratch2/christensen/MolecH/results/h516g768"
outfile = "/astro/net/scratch2/christensen/MolecH/results/h516_noH2_H2_H2SF"
outfile = "/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g3HBWK/steps/h516.cosmo25cmb.3072g3HBWK.00396.00012.dir/h516_H2"

cosmo = 1
cosmo= 6

msol_per_sysmass       = 2.310e15
kpc_per_syslength        = 25000.

last = 408; 168;512
dt = 1
start = 408; 168;512 

FOR i = start/dt, last/dt DO BEGIN 
    if (i*dt lt 10) THEN step = '0000'+STRTRIM(i*dt,2) ELSE BEGIN
        if (i*dt lt 100) THEN step = '000'+STRTRIM(i*dt,2) ELSE step = '00'+STRTRIM(dt*i,2)
    ENDELSE
    filename = prefix+base+"."+step
ENDFOR

filename = prefix+base ;+"."+['00512','00512','00513']+'.halo.1'
rotcurve,[filename],msol_per_sysmass,kpc_per_syslength,keys =keys,/color;,outfile = outfile,/color
end

pro rotCurve_h603
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g2HBWK/"
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/"
base = ["steps/h603.cosmo50cmb.1536g1BWK.00512.dir/h603.cosmo50cmb.1536g1BWK.00512.halo.1", $
        "h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512.halo.1", $
        "h603.cosmo50cmb.2304g5bwK/h603.cosmo50cmb.2304g5bwK.00512.halo.1"]
;base = ["steps/h603.cosmo50cmb.1536g1BWK.00512.dir/h603.cosmo50cmb.1536g1BWK.00512.halo.1"]
base = ["h603.cosmo50cmb.1536g3HBWK/steps/h603.cosmo50cmb.1536g3HBWK.00120.dir/h603.cosmo50cmb.1536g3BWK.00120.halo.1"]
keys = ["H2","no H2, low threshold","no H2, high threshold"]
keys = ["H2"]

outfile = "/astro/net/scratch2/christensen/MolecH/results/h603g452"
outfile = "/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g3HBWK/steps/h603.cosmo50cmb.1536g3HBWK.00120.dir/h603.1536g3HBWK.00120"
msol_per_sysmass       = 1.84793e16
kpc_per_syslength      = 50000.

filename = prefix + base
rotcurve,[filename],msol_per_sysmass,kpc_per_syslength,keys = keys,/color;,outfile = outfile,/color
end

pro leroy_twins
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/"
base = ["h603.cosmo50cmb.2304g3HBWK_00504/steps_column/h603.cosmo50cmb.2304g3HBWK.00504.00020.dir/h603.cosmo50cmb.2304g3HBWK.00504.00020.halo.1", $
"h516.cosmo25cmb.2304g3HBWK_00504/steps/h516.cosmo25cmb.2304g3HBWK.00504.00030.dir/h516.cosmo25cmb.2304g3HBWK.00504.00030.halo.1"]
base = ["h603.cosmo50cmb.2304g3HBWK_00504/steps_column/h603.cosmo50cmb.2304g3HBWK.00504.00020.dir/h603.cosmo50cmb.2304g3HBWK.00504.00020.halo.1"]
base = ["h603.cosmo50cmb.2304g3HBWK_00504/steps_ssource/h603.cosmo50cmb.2304g3HBWK.00504.00018.dir/h603.cosmo50cmb.2304g3HBWK.00504.00018.halo.1"]

keys = ["h603","h516"]
keys = "h603"
outfile = "/astro/net/scratch2/christensen/MolecH/results/twins"
outfile = "/astro/net/scratch2/christensen/MolecH/results/h603ssource"
msol_per_sysmass = [1.84793e16,2.310e15 ]
msol_per_sysmass = 1.84793e16
kpc_per_syslength = [50000.,25000]
kpc_per_syslength = 50000.

filename = prefix + base
leroy,[filename],msol_per_sysmass,kpc_per_syslength,keys = keys,/color;,outfile = outfile
end

pro leroy_res
prefix = "/astro/net/scratch2/christensen/MolecH/11M/"
;prefix = "/astro/net/nbody1/christensen/MolecH/MWHR/"
base = [  $
         "Disk_Iso_1e4/largeStar/MW_disk",$
         "Disk_Iso_1e5/largeStar/MW_disk" $
       ]
;base = [ $
;       "gas_merger0.1_split_s16n10/steps/gas_merger0.1_split_s16n10.00060.dir/gas_merger0.1_split_s16n10", $
;       "gas_merger0.1_single/steps/gas_merger0.1_single.00060.dir/gas_merger0.1_single" $
;]

keys = ['1e4 GP','1e5 GP' ];k
;keys = ['Split','Single' ]

outfile = "/astro/net/scratch2/christensen/MolecH/results/11Mdisk_res_"
; outfile = "/astro/net/scratch2/christensen/MolecH/results/gas_merger_res"

msol_per_sysmass         = 1.36e17
kpc_per_syslength        = 1e5

;msol_per_sysmass         = 2.3262e5
;kpc_per_syslength        = 1.0

last = 3;512
dt = 1
start = 3;512 

;last = 60;512
;dt = 1
;start = 60;512

FOR i = start/dt, last/dt DO BEGIN 
    if (i*dt lt 10) THEN step = '0000'+STRTRIM(i*dt,2) ELSE BEGIN
        if (i*dt lt 100) THEN step = '000'+STRTRIM(i*dt,2) ELSE step = '00'+STRTRIM(dt*i,2)
    ENDELSE
    filename = prefix+base+"."+step
    leroy,[filename],msol_per_sysmass,kpc_per_syslength,keys = keys,outfile = outfile + step
ENDFOR
END



pro leroy,files,dMsolUnit,dKpcUnit, cosmo = cosmo, outfile = outfile,keys = keys,color = color
;rotCurve,dt=20,last=300,first=20 
;rotCurve,dt = 1,first = 5, last = 5

;*********************************** Units ****************************************
hubble = 73.0
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
dDelta = 0.003
sec_per_year = 31556926
molec_weight = (0.76*1 + 0.24*4.0)

;************************************ Plot Format ************************************
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outfile) THEN BEGIN
    set_plot,'ps' 
    nbins=100.0
    linestyles = [0,2]
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5
 ;   !x.charsize=2.25
 ;   !y.charsize=2.25
;    !p.font=0 
ENDIF ELSE BEGIN
 ;   !P.CHARSIZE = 1.25
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    set_plot,'x'
    nbins=100.0
    linestyles = [0,2]
ENDELSE
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    colors = (findgen(N_ELEMENTS(files)) + 1)*240/N_ELEMENTS(files)
    thicks = fltarr(N_ELEMENTS(files)) + 2
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(N_ELEMENTS(files)) + 1)*10.0 + 5.0
    thicks = (findgen(N_ELEMENTS(files)) + 1)*6/N_ELEMENTS(files) - 1
ENDELSE

IF NOT (KEYWORD_SET(legend_t)) THEN legend_t = ''
names = files

;*********************************** Loading Sims ***********************************
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    IF (N_ELEMENTS(dKpcUnit) GT 1) THEN kpc_per_syslength = dKpcUnit[ifile] ELSE kpc_per_syslength = dKpcUnit
    IF (N_ELEMENTS(dMsolUnit) GT 1) THEN msol_per_sysmass = dMsolUnit[ifile] ELSE msol_per_sysmass = dMsolUnit
    dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3
    vel_convert = hubble * (kpc_per_syslength / 1000.0)/2.894405
    timeunit = SQRT((cm_per_kpc*kpc_per_syslength)^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year

    filename = files[ifile]
    names[ifile] = (STRSPLIT(filename,'/',/extract))[N_ELEMENTS(STRSPLIT(filename,'/')) - 1]
    IF KEYWORD_SET(cosmo) THEN BEGIN
        filenametip = filename + '.halo.'+STRTRIM(cosmo,2)+'.std'
        grpfile = filename+'.amiga.grp'
        grp = read_lon_array(grpfile)
        readcol,filename+".HI",HI_prime,/silent
        readcol,filename+".H2",H2_prime,/silent
        IF (FILE_TEST(filenametip)) THEN BEGIN  ;If the halo has already been selected (saves lots of time)
            rtipsy,filename,h,g,d,s,/justhead
            ind = where(grp eq cosmo,comp=indcomp)
            indg = ind[where(ind lt h.ngas)]
            rtipsy,filenametip,hhalo,g,d,s 
        ENDIF ELSE BEGIN  ;If the halo hasn't been selected, then select and rotate it
            rtipsy,filename,h,g,d,s
            ind = where(grp eq cosmo,comp=indcomp)
            inds = ind[where(ind ge h.ngas+h.ndark)]-h.ngas-h.ndark
            indg = ind[where(ind lt h.ngas)]
            indd = ind[where(ind ge h.ngas and ind lt h.ngas+h.ndark)]-h.ngas
            s = s[inds]
            g = g[indg]
            d = d[indd]
            rtipsy, filename+'.amiga.gtp', h1,g1,d1,s1       
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
            limit=5*h.time/dist_units ;use stars within this limit (physical kpc/dist_units)                    
            align,g,d,s,limit
        ENDELSE
        HI_frac = (HI_prime[1:h.ngas])[indg]
        H2_frac = (H2_prime[1:h.ngas])[indg]      
        readcol,filename+'.shear',mach
        column = (1.0/mach[1:h.ngas]*kpc_per_syslength*cm_per_kpc)[indg]
        readcol,filename+'.smoothlength',smoothlengths_L
        smoothlengths = (smoothlengths_L[1:h.ngas]*kpc_per_syslength*cm_per_kpc)[indg]
        length = MAX([[column],[smoothlengths]],DIMENSION = 2)
    ENDIF ELSE BEGIN
        rtipsy,filename,h,g,d,s
        print,"file read"
        readcol,filename+".HI",HI_prime,/silent
        readcol,filename+".H2",H2_prime,/silent
        readcol,filename+".lw",lw_prime,/silent
        HI_frac = HI_prime[1:h.ngas]
        H2_frac = H2_prime[1:h.ngas]
        lw = lw_prime[1:h.ngas]
        readcol,filename+'.shear',mach,/silent
        column = 1.0/mach[1:h.ngas]*kpc_per_syslength*cm_per_kpc
        readcol,filename+'.smoothlength',smoothlengths_L,/silent
        smoothlengths = smoothlengths_L[1:h.ngas]*kpc_per_syslength*cm_per_kpc
        length = MIN([[column],[smoothlengths]],DIMENSION = 2)        
        print,"arrays read"
    ENDELSE

    total_H = fltarr(N_ELEMENTS(HI_frac))
    total_He = fltarr(N_ELEMENTS(HI_frac))
    total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
    IF(comple[0] ne -1) THEN total_He[comple] = (-0.446*(g[comple].zmetal - 0.1)/0.9 + 0.446)/4.0 ;
    total_H = 1.0 - total_He*4.0 - g.zmetal 
    HI = HI_frac*g.mass*msol_per_sysmass
    H2 = H2_frac*g.mass*msol_per_sysmass
    marray = g.mass*amu_per_gm*gm_per_msol
    hmarray = (HI+2.0*H2)*amu_per_gm*gm_per_msol
    h2marray = (2.0*H2)*amu_per_gm*gm_per_msol
    rho = alog10(g.dens * dens_convert * total_H) ;H Atoms per CC 
    t = alog10(g.tempg)
    v = g.vx*vel_convert
    r = sqrt(g.x*g.x + g.y*g.y)*kpc_per_syslength


;********************************* H2 Fraction Column ***********************************************************
    readcol,'/astro/users/christensen/code/MolecH/Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
    readcol,'/astro/users/christensen/code/MolecH/Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'
    IF (KEYWORD_SET(outfile)) THEN device,filename=outfile+STRTRIM(ifile,2)+'frac_columnMass2.eps',/color,bits_per_pixel= 8,/times ELSE window,3
   plot,alog10(g.dens*dens_convert*total_H*length),H2_frac/(total_H/2.0), xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),title = names[ifile],yrange = [1e-6,1.0],/ylog,psym=3,xstyle=1,ystyle=1,symsize =0.5,xrange = [19,23] 
    
    IF 1 THEN BEGIN
        ind = where(H2_frac/total_H ge 5e-7)
        dens_ind = alog10(g[ind].dens*dens_convert*total_H[ind]*length[ind])
        frac_ind = H2_frac[ind]/(total_H[ind]/2.0) 
        
        colorarr = alog10(lw[ind])
        min_colorarr = 14 ;MIN(colorarr)
        max_colorarr = MAX(colorarr)
        colorarr[where(colorarr lt min_colorarr)] = min_colorarr
                                ;   colorarr[where(colorarr gt max_colorarr)] = max_colorarr
        colors_pd = FIX((colorarr - min_colorarr)/(max_colorarr - min_colorarr)*254)
        for ii=0LL,LONG64(N_ELEMENTS(g[ind].dens)) - 1 do $
          oplot,[dens_ind[ii],dens_ind[ii]], [frac_ind[ii],frac_ind[ii]],psym = 3,color=colors_pd[ii]
;        position = [0.88, 0.1, 0.95, 0.9]
;        colorbar,maxrange = max_colorarr,minrange = min_colorarr,/vertical,position = position,format = '(F10.1)'
;        colorarr_hist = histogram(colorarr,locations = x,nbins = 100,min = min_colorarr,max = max_colorarr)
;        plot,colorarr_hist,x,psym = 10,position = position,/NOERASE,XRANGE = [0,MAX(colorarr_hist)],YRANGE = [min_colorarr, max_colorarr],xstyle = 1,ystyle = 1,xticks = 1,yticks = 1,XTICKFORMAT='(A1)'
    ENDIF
    IF (KEYWORD_SET(outfile)) THEN BEGIN
        device,/close
    ENDIF
        
    minr = 0
    maxr = 10.
    range = maxr - minr
    r_bin = findgen(nbins)*range/(nbins)
    y_H2 = weighted_HISTOGRAM(r,input=H2,binsize=range/nbins,min=minr,max=maxr)
    y_HI = weighted_HISTOGRAM(r,input=HI,binsize=range/nbins,min=minr,max=maxr)
    dr = range/(nbins)
    area = ((r_bin + dr)*(r_bin + dr) - r_bin*r_bin)*1e6
    IF iFile eq 0 THEN prof_y = fltarr(2,N_ELEMENTS(files),N_ELEMENTS(y_H2))
    prof_y[0,iFile,*] = y_H2/area
    prof_y[1,iFile,*] = y_HI/area        

    range = 10
    xmin  = -1.0*range/kpc_per_syslength
    xmax  =  1.0*range/kpc_per_syslength 
    ymin  = -1.0*range/kpc_per_syslength 
    ymax  =  1.0*range/kpc_per_syslength 
    delta = 0.4/kpc_per_syslength 
    nx = (xmax - xmin)/delta
    ny = (ymax - ymin)/delta    
    grid_sd = fltarr(nx,ny)
    grid_frac = grid_sd
    grid_z = grid_sd
    grid_rad = grid_sd
    xarray = findgen(nx + 1)*delta + xmin
    yarray = findgen(ny + 1)*delta + ymin
    FOR ix = 0, nx - 1 DO BEGIN
        FOR iy = 0, ny - 1 DO BEGIN
            ind = where(g.x gt xarray[ix] AND g.x lt xarray[ix + 1] AND g.y gt yarray[iy] AND g.y lt yarray[iy + 1])
            if ind[0] NE -1 THEN BEGIN
                grid_frac[ix,iy] = TOTAL(h2marray[ind])/TOTAL(hmarray[ind])
                grid_sd[ix,iy] = TOTAL(hmarray[ind])/cm_per_kpc/cm_per_kpc/kpc_per_syslength/kpc_per_syslength/delta/delta 
                grid_z[ix,iy] = TOTAL(g[ind].zmetal*g[ind].mass)/TOTAL(g[ind].mass)
                grid_rad[ix,iy] = SQRT((xarray[ix]+xarray[ix + 1])^2/4.0 + (yarray[iy]+yarray[iy + 1])^2/4.0)
          ;      oplot,[grid_sd[ix,iy],grid_sd[ix,iy]],[grid_frac[ix,iy],grid_frac[ix,iy]],psym = 3,color = FIX(grid_z[ix,iy]/zsol*240)
            ENDIF ELSE BEGIN
                grid_frac[ix,iy] = 0
                grid_sd[ix,iy] = 0
                grid_z[ix,iy] = 0
            ENDELSE
        ENDFOR
    ENDFOR    
    oplot,alog10(grid_sd),grid_frac,psym = 1
        position = [0.88, 0.1, 0.95, 0.9]
        colorbar,maxrange = max_colorarr,minrange = min_colorarr,/vertical,position = position,format = '(F10.1)'
        colorarr_hist = histogram(colorarr,locations = x,nbins = 100,min = min_colorarr,max = max_colorarr)
        plot,colorarr_hist,x,psym = 10,position = position,/NOERASE,XRANGE = [0,MAX(colorarr_hist)],YRANGE = [min_colorarr, max_colorarr],xstyle = 1,ystyle = 1,xticks = 1,yticks = 1,XTICKFORMAT='(A1)' 
        grid_rad = grid_rad *kpc_per_syslength 

        plot,grid_rad,grid_frac,psym = 3,/ylog
    print,'Temperature (Mean, Min, Max): ',MEAN(g.tempg),MINMAX(g.tempg)
    print,'Density (Mean, Min, Max): ',MEAN(rho),MINMAX(rho)
    print,'H2 (Mean, Min, Max): ',MEAN(H2_frac/(total_H/2.0)),MINMAX(H2_frac/(total_H/2.0))
    print,'Total Mass in H2 + HI: ',TOTAL((HI_frac + 2.0*H2_frac)/(total_H)*g.mass*msol_per_sysmass)   
    print,'Total Mass in H2: ',TOTAL(H2_frac/(total_H/2.0)*g.mass*msol_per_sysmass)
    print,'Percent Mass in Cold ISM: ',TOTAL(H2_frac/(total_H/2.0)*g.mass*msol_per_sysmass)/TOTAL(g.mass*msol_per_sysmass)*100.0
    print,' '
    IF NOT (KEYWORD_SET(outfile)) THEN stop
ENDFOR
IF NOT (KEYWORD_SET(keys)) THEN keys = names

;********************************* Radial Profile ***********************************************************
IF (KEYWORD_SET(outfile)) THEN BEGIN
    device,filename=outfile+'profile2_.eps',/color,bits_per_pixel= 8,/times 
    IF KEYWORD_SET(COLOR) THEN BEGIN
        plot,r_bin,prof_y[1,0,*],psym=10,xtitle="Radius [kpc]",yrange = [1,MAX(prof_y[1,*,*])],ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' pc^2')+"]" ,/ylog,thick = thicks[0]
;           oplot,r_bin,prof_y[1,0,*],psym=10,color=colors[0],thick = thicks[0]
    ENDIF ELSE plot,r_bin,prof_y[1,0,*],psym=10,color=colors[0],xtitle="Radius [kpc]",yrange = [1,MAX(prof_y[1,*,*])],ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' pc^2')+"]" ,/ylog 
ENDIF ELSE BEGIN 
    window,5
    plot,r_bin,prof_y[1,0,*],psym=10,xtitle="Radius [kpc]",yrange = [1,1000],ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' pc^2')+"]" ,/ylog,thick = thicks[0] ;,title='Gas Surface Density Profile -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years',/ylog,yrange = [1e2,3e7]
ENDELSE
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    oplot,r_bin,prof_y[1,iFile,*],psym=10,linestyle=linestyles[0],color=colors[iFile],thick = thicks[iFile]
    oplot,r_bin,prof_y[0,iFile,*],psym=10,linestyle=linestyles[1],color=colors[iFile],thick = thicks[iFile]
ENDFOR
legend,["HI "+legend_t,"H2 "+legend_t],linestyle=linestyles,/right,charsize=1.2 ;,color=colors,thick = thicks
IF (N_ELEMENTS(files) gt 1) THEN legend,keys,color = colors,thick = thicks,/left ,/bottom ,charsize = 0.75,linestyle = fltarr([N_ELEMENTS(files)])
IF (KEYWORD_SET(outfile)) THEN device,/close

stop
END
