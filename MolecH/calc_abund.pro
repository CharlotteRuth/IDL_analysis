pro calc_abund,dt = dt, last=last, first = first

IF (NOT KEYWORD_SET(last)) THEN last = 3000
IF (NOT KEYWORD_SET(dt)) THEN dt = 1
IF (NOT KEYWORD_SET(first))THEN first = 1;10/dt + 1

;base = "/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/"
;base = "/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_H2_UV_solar/MW_disk" ;calc_abund,dt = 5,first = 10, last = 30
;base = "/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/Metal_Cooling_H2_UV_solar_soft_noSF/MW_disk"
;base = "/astro/net/scratch2/christensen/MolecH/12M/Disk_Iso_1e5/MW_disk"  ;abund,dt = 1,first = 10, last = 10
;base = "/astro/net/scratch1/christensen/MolecH/12M/CoolingCloud/CoolingCloud"
;base = "/astro/net/scratch1/christensen/MolecH/12M/H2_cooling/12M"

base = "/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g2HBWK/h516.cosmo25cmb.1536g2HBWK"
!Y.STYLE = 1
!X.STYLE = 1
loadct,10

cm_per_kpc = 3.0857d21
gm_per_msol = 1.989d33
amu_per_gm = 6.022142d23
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)

;msol_per_sysmass = 1.36e17
;kpc_per_syslength = 1e5
msol_per_sysmass       = 2.310e15
kpc_per_syslength        = 25000.

lengthunit =  kpc_per_syslength*cm_per_kpc;system length unit in cm (=1kpc)
massunit   = msol_per_sysmass*gm_per_msol;system mass unit in gm
dens_convert =  msol_per_sysmass * gm_per_msol/kpc_per_syslength^3/cm_per_kpc^3 ;massunit*amu_per_gm/molec_weight/lengthunit^3 ;Converts to grams/cm^3

FOR is = first/dt, last/dt DO BEGIN 
    if (is*dt lt 10) THEN step = '0000'+STRTRIM(is*dt,2) ELSE BEGIN
        if (is*dt lt 100) THEN step = '000'+STRTRIM(is*dt,2) ELSE step = '00'+STRTRIM(dt*is,2)
    ENDELSE
    readcol,base+"."+step+".HI",HI_prime,/silent
    readcol,base+"."+step+".H2",H2_prime,/silent
    readcol,base+"."+step+".mach_sheer",mach_prime,/silent    
    rtipsy,base+"."+step,h,g,d,s
;    readcol,base + "MW_disk_solar_soft.HI",HI_prime,/silent
;    readcol,base + "MW_disk_solar_soft.H2",H2_prime,/silent
;    rtipsy,base + 'MW_disk_solar_soft.std',h,g,d,s
    HI = HI_prime[1:h.ngas]
    H2 = H2_prime[1:h.ngas]
    mach = mach_prime[1:h.ngas]
    t = g.tempg
    rho = g.dens*dens_convert
    zmetal = g.zmetal
    smooth = cm_per_kpc*kpc_per_syslength/mach

    total_H = fltarr(N_ELEMENTS(HI))
    total_He = fltarr(N_ELEMENTS(HI))
    total_He[where(zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*zmetal[where(zmetal le 0.1)])/4.0
    IF(comple[0] ne -1) THEN total_He[comple] = (-0.446*(zmetal[comple] - 0.1)/0.9 + 0.446)/4.0 ;
    total_H = 1.0 - total_He*4.0 - zmetal ; 

    print,'Computed Range'
    print,'HI: ',MINMAX(HI/total_H)
    print,'H2: ',MINMAX(H2*2/total_H)
    print,'LOG(HI): ',aLOG10(MINMAX(HI/total_H))
    print,'LOG(H2): ',ALOG10(MINMAX(H2*2/total_H))

    all = fltarr(9,N_ELEMENTS(t))
    for i=0LL,LONG64(N_ELEMENTS(t)) - 1 do  BEGIN
        all[*,i] = LTE_abund(t[i],g[i].dens*dens_convert,zmetal[i],smooth[i])
    ENDFOR

    Y_HI = REFORM(all[0,*])
    Y_H2 = REFORM(all[1,*])    
    print,''
    print,'LTE Range'
    print,'HI: ',MINMAX(Y_HI/total_H)
    print,'H2: ',MINMAX(Y_H2*2/total_H)
    print,'LOG(HI): ',aLOG10(MINMAX(Y_HI/total_H))
    print,'LOG(H2): ',ALOG10(MINMAX(Y_H2*2/total_H))
    HI_prime[1:h.ngas] = Y_HI
    H2_prime[1:h.ngas] = Y_H2
    openw,1,base+"."+step+".HI_LTE"
    openw,2,base+"."+step+".H2_LTE"
;    openw,1,base + "MW_disk_solar_soft.HI_LTE"
;    openw,2,base + "MW_disk_solar_soft.H2_LTE"
    for i=0LL,LONG64(N_ELEMENTS(HI_prime)) - 1 do  BEGIN
        printf,1,HI_prime[i]
        printf,2,H2_prime[i]
    ENDFOR
    close,1
    close,2

    rho = g.dens*dens_convert*total_H* 5.9753790e+23
    rhoCol = g.dens*dens_convert*total_H*smooth* 5.9753790e+23
    rho_HI = g.dens*dens_convert*Y_HI*smooth* 5.9753790e+23
    rho_H2 = g.dens*dens_convert*Y_H2*smooth* 5.9753790e+23
    N_H = 10^(findgen(100)*5./100. + 18)
    N_H2 = N_H/2.0
    x = N_H2/5d14
    omega_H2 = 0.2
    sigma_d = 2d-21

    colors = FIX(Y_H2/(total_H/2.0)*205)+50
;    window,0
;    plot,rho,t,psym = 3,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",xrange = [10^(-6),10^(4.5)],yrange = [1e1, 1e7],/xlog,/ylog,title='Phase Diagram -- '+strtrim(step)
;    for i=0LL,LONG64(N_ELEMENTS(t)) - 1 do oplot,[rho[i],rho[i]],[t[i],t[i]],psym = 3,color=colors[i]
    window,2
    plot,rho,Y_H2/(total_H/2.0),yrange=[0.0,1],xrange=[1e-1,1e4],/xlog,psym=3,xstyle=1,ystyle=1,xtitle=textoidl('n_H (cm^{-3})'),ytitle=textoidl('f_{H_2}'),symsize = 0.5,title='Fraction H2 -- '+strtrim(step) ;/((2.0*H2_frac + HI_frac)/2)
    oplot,N_H/1d19,exp(-1.0*sigma_d*(N_H + 2.0*N_H))*(1 - omega_H2)/(1 + 2.0*x)^2 + omega_H2/SQRT(1 + 2.0*x)*exp(-0.00085*SQRT(1 + 2.0*x))

    window,3
    plot,rho,REFORM(all[7,*]),psym = 3,/xlog,yrange=[1,0],xrange = [1e-1,1e4],title = 'Dust Sheilding'    
    oplot,N_H/1d19,exp(-1.0*sigma_d*(N_H))

    window,1
    plot,rho,REFORM(all[6,*]),psym = 3,/xlog,yrange=[1,0],xrange = [1e-1,1e4],title = 'Self Sheilding'
    oplot,(N_H2*2 + N_H)/1d19,(1 - omega_H2)/(1 + x)^2 + omega_H2/SQRT(1 + x)*exp(-0.00085*SQRT(1 + x))

    window,4
    plot,rhoCol,REFORM(all[7,*]),psym = 3,/xlog,yrange=[1,0],title = 'Dust Sheilding'
    oplot,N_H,exp(-1.0*sigma_d*(N_H))

    window,5
    plot,rhoCol,REFORM(all[6,*]),psym = 3,/xlog,yrange=[1,0],title = 'Self Sheilding'    
    oplot,(N_H2*2 + N_H),(1 - omega_H2)/(1 + x)^2 + omega_H2/SQRT(1 + x)*exp(-0.00085*SQRT(1 + x))

    window,6
    plot,rhoCol,Y_H2/(total_H/2.0),yrange=[1e-6,1],/xlog,psym=3,xstyle=1,ystyle=1,xtitle=textoidl('n_H (cm^{-3})'),ytitle=textoidl('f_{H_2}'),symsize = 0.5,title='Fraction H2 -- '+strtrim(step),/ylog,xrange=[1e18,1e22]
    stop
endfor
end
