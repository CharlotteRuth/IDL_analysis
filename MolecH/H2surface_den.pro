FUNCTION read_dencube_fits, infile, header, kpcunit=kpcunit, noscale = noscale
if not keyword_set(infile) then begin
    print,"Syntax:"
    print,""
    print,"read_cube_fits(infile, [header, mrdh])"
    print,""
    return,0.
endif

array = mrdfits(infile, 0, mrdh, /silent)
finite = where(FINITE(array),complement = nonfinite)
IF (nonfinite[0] ne -1) then array[nonfinite] = 0
;array = double(10^(double(array)*1.0))
header = {cubedenheader, $
          naxis: 0, $
          naxis1: 0, $
          naxis2: 0, $
          crval1: 0.0, $
          cdelt1: 0.0, $
          crpix1: 0.0, $
          crval2: 0.0, $
          cdelt2: 0.0, $
          crpix2: 0.0, $
          bscale: 0.0, $
          bzero: 0.0, $
          blank: 0.0, $
          kpcunit: 0.0, $
          munit: 0.0 $
         }
header.naxis = sxpar(mrdh, 'NAXIS')
header.naxis1 = sxpar(mrdh, 'NAXIS1')
header.naxis2 = sxpar(mrdh, 'NAXIS2')

header.crval1 = sxpar(mrdh, 'CRVAL1')
header.cdelt1 = sxpar(mrdh, 'CDELT1')
header.crpix1 = sxpar(mrdh, 'CRPIX1')

header.crval2 = sxpar(mrdh, 'CRVAL2')
header.cdelt2 = sxpar(mrdh, 'CDELT2')
header.crpix2 = sxpar(mrdh, 'CRPIX2')

header.bscale = sxpar(mrdh, 'BSCALE')
header.bzero = sxpar(mrdh, 'BZERO')
header.blank = sxpar(mrdh, 'BLANK')

IF NOT KEYWORD_SET(noscale) THEN BEGIN
; scale the array
    array = array * header.bscale + header.bzero
    array = 10.0^(float(array))

; set blanks to 0.
    vals = where(array EQ 10^(float(header.blank * header.bscale + header.bzero)))
    IF (vals[0] NE -1) THEN array[vals] = 0.0
ENDIF

IF keyword_set(kpcunit) THEN BEGIN
    header.crval1 = header.crval1 * kpcunit
    header.crval2 = header.crval2 * kpcunit
    header.cdelt1 = header.cdelt1 * kpcunit
    header.cdelt2 = header.cdelt2 * kpcunit
ENDIF

RETURN,array
END

;prefix = '/home/christensen/Storage2/UW/'
;prefix = '/astro/net/scratch2/christensen/'

;dir = prefix + 'MolecH/11M/Disk_Iso_1e5/largeStar/'
;file = 'MW_disk.00003' 

;dir = '/astro/net/nbody1/christensen/MolecH/MWHR/12M_hr'
;file = '12M_hr.00800'
;dKpcUnit        = 1
;dMsolUnit       = 2.362e5

;---------------- h516 --------------------
;dKpcUnit        = 25000.
;dMsolUnit       = 2.310e15

;dir = prefix + 'MolecH/Cosmo/h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.00512.dir/'
;file = 'h516.cosmo25cmb.1536g2HBWK.00512.halo.1'

;dir = prefix + 'MolecH/Cosmo/h516.cosmo25cmb.3072g3HBWK/steps/h516.cosmo25cmb.3072g3HBWK.00396.00012.dir'
;file = 'h516.cosmo25cmb.3072g3HBWK.00396.00012.halo.1'

;dir = prefix + 'MolecH/Cosmo/h516.cosmo25cmb.3072g3HBWK_00504/steps/h516.cosmo25cmb.3072g3HBWK.00504.00002.dir'
;file = 'h516.cosmo25cmb.3072g3HBWK.00504.00002.halo.1'

;dir = prefix + 'MolecH/Cosmo/h516.cosmo25cmb.2304g3HBWK_00504/steps/h516.cosmo25cmb.2304g3HBWK.00504.00002.dir'
;file = 'h516.cosmo25cmb.2304g3HBWK.00504.00002.halo.1'

;dir = prefix + 'MolecH/Cosmo/h516.cosmo25cmb.2304g3HBWK_00504/steps/h516.cosmo25cmb.2304g3HBWK.00504.00030.dir'
;file = 'h516.cosmo25cmb.2304g3HBWK.00504.00030.halo.1'

;dir = prefix + 'MolecH/Cosmo/h516.cosmo25cmb.3072g3HBWK/steps/h516.cosmo25cmb.3072g3HBWK.00396.00012.dir'
;file = 'h516.cosmo25cmb.3072g3HBWK.00396.00012.halo.1'

;dir = prefix + 'MolecH/Cosmo/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir'
;file = 'h516.cosmo25cmb.1536g3HBWK.00512.halo.1'

;dir = prefix + 'MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir'
;file = 'h516.cosmo25cmb.3072g14HBWK.00492.halo.1'

;----------------- h799 -----------------

;----------------- h603 -----------------
;dMsolUnit       = 1.84793e16
;dKpcUnit        = 50000.

;file = 'h603.cosmo50cmb.2304g3HBWK.00504.00002.halo.1'
;dir = prefix + 'MolecH/Cosmo/h603.cosmo50cmb.2304g3HBWK_00504/steps/h603.cosmo50cmb.2304g3HBWK.00504.00002.dir'

;----------------- h986 -----------------

;H2surface_den,dir,file,tipsyfile=file,dKpcUnit = dKpcUnit, dMsolUnit= dMsolUnit,/outplot
;H2surface_den,dir,file,dKpcUnit = dKpcUnit, dMsolUnit= dMsolUnit,/outplot
PRO H2surface_den,dir,file,tipsyfile=tipsyfile,dKpcUnit = dKpcUnit, dMsolUnit = dMsolUnit,outplot = outplot
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790d+23
molec_weight = (0.76*1 + 0.24*4.0)
loadct,39
set_plot,'x'
f_H = 0.764

IF NOT KEYWORD_SET(dKpcUnit) THEN dKpcUnit = 50000.
IF NOT KEYWORD_SET(dMsolUnit) THEN dMsolUnit = 1.36e17

cd,dir
cubeH2 = read_dencube_fits(file+'.H2.arr.90.fits',headerH2)*2.0
;headerH2.CDELT1 = headerH2.CDELT1*h.time
;headerH2.CDELT2 = headerH2.CDELT2*h.time
pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*cm_per_kpc*cm_per_kpc ;cm^2
;pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*1e6  ;pc^2
cubeH2_surface_den = transpose(cubeH2)/pix_areaH2*gm_per_msol*amu_per_gm;*f_H ;amu/cm^2
;cubeH2_surface_den = transpose(cubeH2)/pix_areaH2;*f_H ;Msun/pc^2

;cubeHI = read_dencube_fits(file+'.HI.arr.45.fits',headerHI)
;pix_areaHI = headerHI.CDELT1*headerHI.CDELT2*cm_per_kpc*cm_per_kpc

;cubeHI = read_dencube_fits(file+'.90.smoothed.mom0.fits',headerHI,/noscale)
cubeHI = read_dencube_fits(file+'.90.cube.mom0.fits'    ,headerHI,/noscale)
cubeHI_surface_den = cubeHI;*f_H ;amu/cm^2
;cubeHI_surface_den = cubeHI/gm_per_msol/amu_per_gm*cm_per_kpc/1000*cm_per_kpc/1000;*f_H ;Msun/pc^2

;cubemetal = read_dencube_fits(file+'.zmetal.fits',headermetal)
;cubemetal = cubemetal - 0.000680577
;cubemetal[where(cubemetal lt 0)] = 0.0250000
;cubemetal = cubemetal - 1.23200e-05
;cubemetal[where(cubemetal lt 0)] = 0
loadct,0
range = 24.0
xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1
yaxes = (findgen(headerHI.NAXIS2) - headerHI.NAXIS1/2.0)*headerHI.CDELT2
window,0,xsize = 600,ysize = 600
;cubeH2_surface_den[where(cubeH2_surface_den lt 3e12)] = 3e12
;cubeHI_surface_den[where(cubeHI_surface_den lt 3e19)] = 3e19
contour,alog10(cubeH2_surface_den),xaxes,yaxes,/fill,nlevels = 60,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = 22.0,min_value = 12.0,title = textoidl('H_2 Surface Density')
colorbar,range = [12,22],/vertical,divisions = 3
window,1,xsize = 600,ysize = 600
contour,alog10(cubeHI_surface_den),xaxes,yaxes,/fill,nlevels = 60,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = 22.0,min_value = 19.0,title = textoidl('HI Surface Density')
colorbar,range = [19,22],/vertical,divisions = 3
stop

window,0,xsize = 600,ysize = 600
contour,alog10(cubeH2_surface_den/gm_per_msol/amu_per_gm*cm_per_kpc/1000*cm_per_kpc/1000),xaxes,yaxes,/fill,nlevels = 60,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = 2,min_value = -2,title = textoidl('H2 Surface Density')
colorbar,range = [0,2],/vertical,divisions = 3

window,1,xsize = 600,ysize = 600
contour,alog10(cubeHI_surface_den/gm_per_msol/amu_per_gm*cm_per_kpc/1000*cm_per_kpc/1000),xaxes,yaxes,/fill,nlevels = 60,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = 2,min_value = -2,title = textoidl('HI Surface Density')
colorbar,range = [0,2],/vertical,divisions = 3

IF KEYWORD_SET(outplot) THEN BEGIN
    squareplot,filename = 'H2surface.ps'
    contour,alog10(cubeH2_surface_den),xaxes,yaxes,/fill,nlevels = 60,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = 22.0,min_value = 12.0,title = textoidl('H_2 Surface Density')
    colorbar,range = [12,22],/vertical,position = [0.86,0.14,0.93,0.90],/right,divisions = 5
    squareplot,/close
    squareplot,filename = 'HIsurface.ps'
    contour,alog10(cubeHI_surface_den/1.98892d33*3.08568021d18*3.08568021d18/6.022142d23),xaxes,yaxes,/fill,nlevels = 60,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = 2,min_value = 0,title = textoidl('HI Surface Density')
    colorbar,range = [0,2],/vertical,position = [0.86,0.14,0.93,0.90],/right,divisions = 3
    squareplot,/close
ENDIF
stop

loadct,39
distance = fltarr(headerHI.NAXIS1,headerHI.NAXIS2)
FOR ix = 0, headerHI.NAXIS1 - 1 DO $
  FOR iy = 0, headerHI.NAXIS2 - 1 DO distance[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy])
ind = where(distance lt 4)
fraction = 2.0*cubeH2_surface_den/(cubeHI_surface_den + 2.0*cubeH2_surface_den)
window,2
;plot,cubeHI_surface_den[ind]+2.0*cubeH2_surface_den[ind],fraction[ind],psym = 3,/ylog,xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.0],xrange = [1e19,1e24],xstyle=1,ystyle=1,symsize =0.5,/xlog
plot, (cubeHI_surface_den     +2.0*cubeH2_surface_den     )*7.9493849e-21,fraction,     psym = 3,xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.0],xrange = [1e19,1e24]*7.9493849e-21,xstyle=1,ystyle=1,symsize =0.5,/xlog,/ylog
oplot,(cubeHI_surface_den[ind]+2.0*cubeH2_surface_den[ind])*7.9493849e-21,fraction[ind],psym = 3,color = 240 
oplot,[1e17,1e25],[1,1]
stop
;oplot,nhw,10^logfH2w,psym = 1,color = 240
;oplot,nhg,10^logfH2g,psym = 4,color = 240
;legend,['Simulated Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [3,4,1],color=[0,240,240],/right,/bottom,charsize = 1.2

IF 0 THEN BEGIN
    zmetal = 0.025
    readcol,'/astro/users/christensen/code/MolecH/Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
    readcol,'/astro/users/christensen/code/MolecH/Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'
    oplot,10^nhw,10^logfH2w,psym = 1,color = 240
    oplot,10^nhg,10^logfH2g,psym = 4,color = 240
    legend,['Simulated Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [3,4,1],color=[0,240,240],/right,/bottom,charsize = 1.2
ENDIF

IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps'
    device,filename = dir+file+'frac_column.ps',/color,bits_per_pixel= 8,/times 
    fraction = cubeH2_surface_den/(cubeHI_surface_den + cubeH2_surface_den)
    plot,cubeHI_surface_den+cubeH2_surface_den,fraction,psym = 3,xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),xrange = [1e19,24],xstyle=1,ystyle=1,symsize =0.5,yrange = [1e-6,1.0],/ylog,/xlog
    oplot,[1e17,1e25],[1,1]
    
    IF 0 THEN BEGIN
        zmetal = 0.025
        readcol,'/astro/users/christensen/code/MolecH/Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
        readcol,'/astro/users/christensen/code/MolecH/Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'
        oplot,10^nhw,10^logfH2w,psym = 1,color = 240
        oplot,10^nhg,10^logfH2g,psym = 4,color = 240
        legend,['Simulated Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [3,4,1],color=[0,240,240],/right,/bottom,charsize = 1.2
    ENDIF
;sigma_d = 2d-21                 ;4d-21
;zsol = 0.0177                   ;0.025
;    columnd = findgen(100)*15.0/100. + 10
;    d_shield = 1.0 - exp(-1.0*sigma_d*zmetal/zsol*(10.0^columnd))
;    oplot,columnd,d_shield,color = 120
    device,/close
    set_plot,'x'
ENDIF
print,'Total Hydrogen Mass: ',TOTAL(cubeH2 + cubeHI),' Solar Masses'

IF KEYWORD_SET(tipsyfile) THEN BEGIN 
    rtipsy,tipsyfile,h,g,d,s
    readcol,tipsyfile+'.H2',h2,/silent
    readcol,tipsyfile+'.HI',hI,/silent 
    readcol,tipsyfile+'.shear',mach,/silent
  ;  readcol,tipsyfile+'.smoothlength',smoothlengths_L,/silent
  ;  sprad = sqrt(g.x*g.x + g.y*g.y + g.z*g.z)
  ;  ind = sprad lt 0.000044
    dens_convert =  dMsolUnit*gm_per_msol*amu_per_gm/dKpcUnit^3/cm_per_kpc^3

    HI = HI[1:h.ngas]
    H2 = H2[1:h.ngas] 
    print,'Total Hydrogen Mass: ',TOTAL((HI + 2.0*H2)*g.mass*dMsolUnit)*f_H,' Solar Masses'
    column = 1.0/mach[1:h.ngas]*dKpcUnit*cm_per_kpc  
   ; smoothlengths = smoothlengths_L[1:h.ngas]*dKpcUnit*cm_per_kpc
    rho = g.dens*dens_convert*(HI + 2.0*H2)

    oplot,g.dens*dens_convert*(HI + 2.0*H2)*column,2.0*H2/(HI + 2.0*H2),psym = 3,color = 50
 
    omega_H2 = 0.035            ;0.2
    sigma_d = 2d-21             ;4d-21
    zsol = 0.0177               ;0.025
    x = (g.dens*dens_convert*(HI + 2.0*H2)*column)/5d14
    zmetal = TOTAL(g.zmetal*g.mass*(HI + 2.0*H2))/TOTAL(g.mass*(HI + 2.0*H2))
    print,'Mean Metallicity: ',zmetal/zsol,' Zsol'
    h_shield = 1 - (1 - omega_H2)/(1 + x)^2 + omega_H2/SQRT(1 + x)*exp(-0.00085*SQRT(1 + x))
    d_shield = 1.0 - exp(-1.0*sigma_d*zmetal/zsol*(g.dens*dens_convert*(HI + 2.0*H2)*column))
;    oplot,(g.dens*dens_convert*(HI + 2.0*H2)*column)/surf_den_convert,d_shield*h_shield,psym = 3,color = 150
;    oplot,g.dens*dens_convert*(HI + 2.0*H2)*column,d_shield*h_shield,psym = 3,color = 150
    
    IF (0) THEN BEGIN
        set_plot,'ps'
        device,filename = dir+file+'frac_column.ps',/color,bits_per_pixel= 8,/times 
        plot,alog10(g.dens*dens_convert*(HI + 2.0*H2)*column),2.0*H2/(HI + 2.0*H2),psym = 3,xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),xrange = [17,24],xstyle=1,ystyle=1,symsize =0.5,yrange = [1e-6,1.0],/ylog
        device,/close
        set_plot,'x'
    ENDIF

    window,3
    surf_den_convert = 1.258075e20
    plot,(g.dens*dens_convert*(HI + 2.0*H2)*column)/surf_den_convert,2.0*H2/(HI + 2.0*H2),psym = 3,xtitle = textoidl("N_{HI} + 2N_{H_2} (M_{\odot} pc^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [0,1.0],xrange = [1,1e4],xstyle=1,ystyle=1,symsize =0.5,/xlog
    krumholz,MEAN(g.zmetal),x,y
    oplot,x,y,color = 240
    krumholz,0.0177*10.0,x,y
    oplot,x,y,color = 240,linestyle = 3 
    krumholz,0.0177/10.0,x,y
    oplot,x,y,color = 240,linestyle = 2
    krumholz,0.0177/100.0,x,y
    oplot,x,y,color = 240,linestyle = 1
    oplot,(g.dens*dens_convert*(HI + 2.0*H2)*column)/surf_den_convert,d_shield*h_shield,psym = 3,color = 150
    d_shield = 1.0 - exp(-1.0*sigma_d*(g.dens*dens_convert*(HI + 2.0*H2)*column))
    oplot,(g.dens*dens_convert*(HI + 2.0*H2)*column)/surf_den_convert,d_shield*h_shield,psym = 3,color = 180
    stop

    window,4
  ;  smoothl_hist = histogram(alog10(smoothlengths),locations = x,min = 14,max = 24,nbins = 50)
    dens_ind = where(rho gt 0)
    column_hit = histogram(alog10(column),locations = x,min = 14,max = 24,nbins = 50)
    columnden_hit = histogram(alog10(g.dens*dens_convert*(HI + 2.0*H2)*column),locations = x,min = 14,max = 24,nbins = 50)
    IF (dens_ind[0] ne -1) THEN column_hit_dens = histogram(alog10(column[dens_ind]),locations = x,min = 14,max = 24,nbins = 50)
 ;   plot,x,smoothl_hist,psym = 10,xtitle = "Log(h) [cm]"
    oplot,x,columnden_hit,linestyle = 1,psym = 10
    oplot,x,column_hit,linestyle = 2,psym = 10
    IF (dens_ind[0] ne -1) THEN oplot,x,column_hit_dens,linestyle = 3,psym = 10
    legend,["Smoothing Length","Column Length","Column Density","High Density Column Length"],linestyle=[0,2,1,3],charsize = 1.2

ENDIF
stop

END
