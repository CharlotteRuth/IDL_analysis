
PRO replicate_H2frac_master
set_plot,'x'
set_plot,'ps'
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
!P.CHARTHICK=4
!X.THICK=4
!Y.THICK=4

dir = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_repl/'
device,filename=dir+'H2frac_metallicty.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 25,xoffset =  2,yoffset =  2
files = dir + ['z0.1_lw3e7_cp10/','z0.3_lw3e7_cp10/','z1_lw3e7_cp10/'] + 'MW_disk.00010'
;replicate_H2frac,files
device,/close
;stop

device,filename=dir+'H2frac_lw.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 25,xoffset =  2,yoffset =  2
files = dir + ['z1_lw3e6_cp10/','z1_lw3e7_cp10/','z1_lw3e8_cp10/'] + 'MW_disk.00010'
;replicate_H2frac,files
device,/close
;stop

device,filename=dir+'H2frac_clump.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 25,xoffset =  2,yoffset =  2
files = dir + ['z1_lw3e7_cp1/','z1_lw3e7_cp10/','z1_lw3e7_cp100/'] + 'MW_disk.00010'
;replicate_H2frac,files
device,/close
;stop

dir = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_repl/'
files = dir + ['z1_lw3e7_cp10/'] + 'MW_disk.00010'
device,filename=dir+'H2frac_metallicty_dat.eps',/color,bits_per_pixel= 8,/times,xsize = 9*N_ELEMENTS(files),ysize = 25,xoffset =  2,yoffset =  2
replicate_H2frac,files

readcol,'/astro/users/christensen/code/MolecH/Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
readcol,'/astro/users/christensen/code/MolecH/Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'
oplot,nhw,10^logfH2w,psym = 1,color = 240
oplot,nhg,10^logfH2g,psym = 4,color = 240
;legend,['Simulated Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [3,4,1],color=[0,240,240],/right,/bottom,charsize = 1.2
device,/close

END



PRO replicate_H2frac,files
hubble = 73.0
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
dDelta = 0.003
sec_per_year = 31556926
molec_weight = (0.76*1 + 0.24*4.0)
dMsolUnit = 1.36e17
dKpcUnit = 1e5
dens_convert =  dMsolUnit * gm_per_msol * 5.9753790e+23/dKpcUnit^3/cm_per_kpc^3
loadct,39
!p.multi = [0,N_ELEMENTS(files),3]
FOR i = 0,N_ELEMENTS(files)-1 DO BEGIN
    rtipsy,files[i],h,g,d,s   
    readarr,files[i]+".HI",h,HI_prime,/ascii  
    readarr,files[i]+".H2",h,H2_prime,/ascii 
    readarr,files[i]+'.shear',h,mach,/ascii
    column = (1.0/mach*dKpcUnit*cm_per_kpc)
    readarr,files[i]+".lw",h,lw_prime,/ascii
    lw = lw_prime[0:h.ngas - 1]*1d30/cm_per_kpc^2/dKpcUnit^2
    total_H = fltarr(N_ELEMENTS(HI_prime))
    total_He = fltarr(N_ELEMENTS(HI_prime))
    total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
    HI_frac = HI_prime[0:h.ngas - 1]
    H2_frac = H2_prime[0:h.ngas - 1] 

    total_H = 1.0 - total_He*4.0 - g.zmetal 
    HI = HI_frac*g.mass*dMsolUnit
    H2 = H2_frac*g.mass*dMsolUnit
    rho = alog10(g.dens * dens_convert * total_H)
    plot,g.dens*dens_convert*total_H,g.tempg,xtitle = textoidl("LOG(n_H) (cm^{-3})"),ytitle = textoidl('Temperature [K]'),xrange = [1e-2,1e3],yrange = [1e1, 1e6],psym=3,xstyle=1,ystyle=1,symsize = 0.5,/ylog,/xlog,thick = 3,charsize = 2.0
ENDFOR

FOR i = 0,N_ELEMENTS(files) -1 DO BEGIN
    rtipsy,files[i],h,g,d,s   
    readarr,files[i]+".HI",h,HI_prime,/ascii  
    readarr,files[i]+".H2",h,H2_prime,/ascii 
    readarr,files[i]+'.shear',h,mach,/ascii
    column = (1.0/mach*dKpcUnit*cm_per_kpc)
    readarr,files[i]+".lw",h,lw_prime,/ascii
    lw = lw_prime[0:h.ngas - 1]*1d30/cm_per_kpc^2/dKpcUnit^2
    total_H = fltarr(N_ELEMENTS(HI_prime))
    total_He = fltarr(N_ELEMENTS(HI_prime))
    total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
    HI_frac = HI_prime[0:h.ngas - 1]
    H2_frac = H2_prime[0:h.ngas - 1] 

    total_H = 1.0 - total_He*4.0 - g.zmetal 
    HI = HI_frac*g.mass*dMsolUnit
    H2 = H2_frac*g.mass*dMsolUnit
    rho = alog10(g.dens * dens_convert * total_H)  

    plot,alog10(g.dens*dens_convert*total_H),alog10(H2_frac/(total_H/2.0)),xtitle = textoidl("LOG(n_H) (cm^{-3})"),ytitle = textoidl('f_{H_2}'),xrange = [-2,3],yrange = [-6, 0],psym=3,xstyle=1,ystyle=1,symsize = 0.5,thick = 3,charsize = 2.0

    ind = where(H2_frac/total_H ge 5e-7)
    if 0 THEN BEGIN
        dens_ind = alog10(g[ind].dens*dens_convert*total_H[ind])
        frac_ind = alog10(H2_frac[ind]/(total_H[ind]/2.0))    
        colorarr = alog10(lw[ind])
        min_colorarr = MIN(colorarr)
        max_colorarr = MAX(colorarr)
        min_colorarr = 3        ;25.0;-4.5
        max_colorarr = 9        ;30.0;2
        colors_pd = FIX((colorarr - min_colorarr)/(max_colorarr - min_colorarr)*254)
        if (where(colorarr lt min_colorarr))[0] ne -1 THEN colors_pd[where(colorarr lt min_colorarr)] = 255
        if (where(colorarr gt max_colorarr))[0] ne -1 THEN colors_pd[where(colorarr gt max_colorarr)] = 255
        for ii=0LL,LONG64(N_ELEMENTS(g[ind].dens)) - 1 do $
          oplot,[dens_ind[ii],dens_ind[ii]], [frac_ind[ii],frac_ind[ii]],psym = 3,color=colors_pd[ii]
        if KEYWORD_SET(outfile) THEN position = [0.88, 0.14, 0.95, 0.92] ELSE position = [0.88, 0.1, 0.95, 0.9]     
        colorbar,maxrange = max_colorarr,minrange = min_colorarr,/vertical,position = position,format = '(F10.1)'
        colorarr_hist = histogram(colorarr,locations = x,nbins = 100,min = min_colorarr,max = max_colorarr)
        plot,colorarr_hist,x,psym = 10,position = position,/NOERASE,XRANGE = [0,MAX(colorarr_hist)],YRANGE = [min_colorarr, max_colorarr],xstyle = 9,ystyle = 9,xticks = 1,yticks = 1,XTICKFORMAT='(A1)',YTICKFORMAT='(A1)'
    ENDIF
;     plot,alog10(g.dens*dens_convert*total_H*smoothlengths),H2_frac/(total_H/2.0), xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),title = names[ifile],psym=3,xstyle=1,ystyle=1,symsize =0.5,xrange = [19,23],yrange = [1e-6,1.0],/ylog
ENDFOR
FOR i = 0,N_ELEMENTS(files)-1 DO BEGIN
    rtipsy,files[i],h,g,d,s   
    readarr,files[i]+".HI",h,HI_prime,/ascii  
    readarr,files[i]+".H2",h,H2_prime,/ascii 
    readarr,files[i]+'.shear',h,mach,/ascii
    column = (1.0/mach*dKpcUnit*cm_per_kpc)
    readarr,files[i]+".smoothlength",h,smoothlengths_L,/ascii
    smoothlengths = smoothlengths_L[0:h.ngas - 1]*dKpcUnit*cm_per_kpc
    readarr,files[i]+".lw",h,lw_prime,/ascii
    lw = lw_prime[0:h.ngas - 1]*1d30/cm_per_kpc^2/dKpcUnit^2    
    column = smoothlengths

    total_H = fltarr(N_ELEMENTS(HI_prime))
    total_He = fltarr(N_ELEMENTS(HI_prime))
    total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
    HI_frac = HI_prime[0:h.ngas - 1]
    H2_frac = H2_prime[0:h.ngas - 1] 

    total_H = 1.0 - total_He*4.0 - g.zmetal 
    HI = HI_frac*g.mass*dMsolUnit
    H2 = H2_frac*g.mass*dMsolUnit
    rho = alog10(g.dens * dens_convert * total_H)  

;     plot,alog10(g.dens*dens_convert*total_H),alog10(H2_frac/(total_H/2.0)),xtitle = textoidl("LOG(n_H) (cm^{-3})"),ytitle = textoidl('f_{H_2}'),xrange = [-2,3],yrange = [-6, 1],psym=3,xstyle=1,ystyle=1,symsize = 0.5
    plot,alog10(g.dens*dens_convert*total_H*column),H2_frac/(total_H/2.0), xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),psym=3,xstyle=1,ystyle=1,symsize =0.5,xrange = [19,23],yrange = [1e-6,1.0],/ylog,thick = 3,charsize = 2.0
    
    ind = where(H2_frac/total_H ge 5e-7)
    if 0 THEN BEGIN
        dens_ind = alog10(g[ind].dens*dens_convert*total_H[ind]*column[ind])
        frac_ind = H2_frac[ind]/(total_H[ind]/2.0)
        colorarr = alog10(lw[ind])
        min_colorarr = MIN(colorarr)
        max_colorarr = MAX(colorarr)
        min_colorarr = 3        ;25.0;-4.5
        max_colorarr = 9        ;30.0;2
        colors_pd = FIX((colorarr - min_colorarr)/(max_colorarr - min_colorarr)*254)
        if (where(colorarr lt min_colorarr))[0] ne -1 THEN colors_pd[where(colorarr lt min_colorarr)] = 255
        if (where(colorarr gt max_colorarr))[0] ne -1 THEN colors_pd[where(colorarr gt max_colorarr)] = 255
        for ii=0LL,LONG64(N_ELEMENTS(g[ind].dens)) - 1 do $
          oplot,[dens_ind[ii],dens_ind[ii]], [frac_ind[ii],frac_ind[ii]],psym = 3,color=colors_pd[ii]
        if KEYWORD_SET(outfile) THEN position = [0.88, 0.14, 0.95, 0.92] ELSE position = [0.88, 0.1, 0.95, 0.9]     
        colorbar,maxrange = max_colorarr,minrange = min_colorarr,/vertical,position = position,format = '(F10.1)'
        colorarr_hist = histogram(colorarr,locations = x,nbins = 100,min = min_colorarr,max = max_colorarr)
        plot,colorarr_hist,x,psym = 10,position = position,/NOERASE,XRANGE = [0,MAX(colorarr_hist)],YRANGE = [min_colorarr, max_colorarr],xstyle = 9,ystyle = 9,xticks = 1,yticks = 1,XTICKFORMAT='(A1)',YTICKFORMAT='(A1)'
        ENDIF
ENDFOR
END
