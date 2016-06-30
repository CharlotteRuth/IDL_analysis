PRO H2frac_particle_master
dir = ['/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.2304g3HBWK_00504/steps/h603.cosmo50cmb.2304g3HBWK.00504.00002.dir','/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g3HBWK_00504/steps/h516.cosmo25cmb.2304g3HBWK.00504.00002.dir']
file =['h603.cosmo50cmb.2304g3HBWK.00504.00002.halo.1','h516.cosmo25cmb.2304g3HBWK.00504.00002.halo.1']

dir = ['/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.2304g3HBWK_00504/steps/h603.cosmo50cmb.2304g3HBWK.00504.00002.dir','/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g3HBWK_00504/steps/h516.cosmo25cmb.2304g3HBWK.00504.00030.dir']
file =['h603.cosmo50cmb.2304g3HBWK.00504.00002.halo.1','h516.cosmo25cmb.2304g3HBWK.00504.00030.halo.1']

dir =['/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.2304g3HBWK_00504/steps/h603.cosmo50cmb.2304g3HBWK.00504.00002.dir','/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g3HBWK/steps/h516.cosmo25cmb.3072g3HBWK.00396.00012.dir']
file =['h603.cosmo50cmb.2304g3HBWK.00504.00002.halo.1','h516.cosmo25cmb.3072g3HBWK.00396.00012.halo.1']

dir =['/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.2304g3HBWK_00504/steps_ssource.J/h603.cosmo50cmb.2304g3HBWK.00504.00002.dir']
file = ['h603.cosmo50cmb.2304g3HBWK.00504.00002.halo.1']

dKpcUnit = ['50000.','25000.']
dMsolUnit = ['1.84793e16','2.310e15']

dir=['/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/', $
     '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_lowCstar/h516.cosmo25cmb.1536g3HBWK_lowCstar.00456.dir/']
file = ['h516.cosmo25cmb.1536g3HBWK.00512.halo.1','h516.cosmo25cmb.1536g3HBWK_lowCstar.00456.halo.1']
dKpcUnit = ['25000.','25000.']
dMsolUnit = ['2.310e15','2.310e15']
key = ['H2','H2, low Cstar']

H2frac_particle,dir,file,dKpcUnit,dMsolUnit,key

END


PRO H2frac_particle,dir,file,dKpcUnit,dMsolUnit,key,outplot = outplot
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790d+23
molec_weight = (0.76*1 + 0.24*4.0)
loadct,39
f_H = 0.764
colors = (findgen(N_ELEMENTS(file)) + 1)*240/N_ELEMENTS(file)
if (KEYWORD_SET(outplot)) then begin
    !Y.STYLE = 1
    !X.STYLE = 1
    !P.THICK = 2.0
    !P.CHARTHICK=2.0
    !P.charsize=2.0
    set_plot,'ps' 
    device,filename = '~/h516_H2_noH2_SK.eps',/color,bits_per_pixel= 8,/times,xsize = 8,ysize = 8,/inch,xoffset =  2,yoffset =  2 
endif else begin
    set_plot,'x'
    window,0,xsize = 600,ysize = 600
endelse

xsigmalow = findgen(300)/100 - 1.
xsigma = 10.0^(findgen(400)/100 + 1.) ;xsigmalow
ysigma=2.5e-4*xsigma^1.4
ysigma1=2.5e-4*xsigma^(1.4 + 0.15)
ysigma2=2.5e-4*xsigma^(1.4 - 0.15)
ysigmalow = xsigmalow*2.4 - 5.0
range = 4*(fltarr(N_ELEMENTS(file))+1)
distance2 = 4*(fltarr(N_ELEMENTS(file))+1)
delta = 0.750/dKpcUnit

IF 0 THEN BEGIN
    FOR i = 0, N_ELEMENTS(file) - 1 DO BEGIN
        cubeH2 = read_dencube_fits(dir[i]+file[i]+'.H2.lr.fits',headerH2)
        pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*3.08568021d21*3.08568021d21
        cubeH2_surface_den = cubeH2/pix_areaH2*gm_per_msol*amu_per_gm*f_H
    
        cubeHI = read_dencube_fits(dir[i]+file[i]+'.HI.lr.fits',headerHI)
        pix_areaHI = headerHI.CDELT1*headerHI.CDELT2*3.08568021d21*3.08568021d21
        cubeHI_surface_den = cubeHI/pix_areaHI*gm_per_msol*amu_per_gm*f_H

        xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1
        yaxes = (findgen(headerHI.NAXIS2) - headerH2.NAXIS1/2.0)*headerH2.CDELT2
    
        distance = fltarr(headerHI.NAXIS1,headerHI.NAXIS2)
        FOR ix = 0, headerHI.NAXIS1 - 1 DO $
          FOR iy = 0, headerHI.NAXIS2 - 1 DO distance[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy])
        ind = where(distance lt distance2[i])
        fraction = 2.0*cubeH2_surface_den/(cubeHI_surface_den + 2.0*cubeH2_surface_den)
        IF i eq 0 THEN plot,cubeHI_surface_den[ind]+2.0*cubeH2_surface_den[ind],fraction[ind],psym = 3,/ylog,xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.0],xrange = [1e19,1e24],xstyle=1,ystyle=1,symsize =0.5,/xlog
        oplot,cubeHI_surface_den[ind]+2.0*cubeH2_surface_den[ind],fraction[ind],psym = 3 ;,color = colors[i]
    ENDFOR
ENDIF
   
FOR i = 0, N_ELEMENTS(file) - 1 DO BEGIN
    rtipsy,dir[i]+file[i],h,g,d,s
    dens_convert =  dMsolUnit[i]*gm_per_msol*amu_per_gm/dKpcUnit[i]^3/cm_per_kpc^3
    readarr,dir[i]+file[i]+'.HI',h,hI,part = 'gas',/ascii
    IF (FILE_TEST(dir[i]+file[i]+".H2")) THEN BEGIN
        readarr,dir[i]+file[i]+'.H2',h,h2,part = 'gas',/ascii
    ENDIF ELSE h2 = fltarr(N_ELEMENTS(hI))
    marray = g.mass*amu_per_gm*gm_per_msol*dMsolUnit[i]
    hmarray = (HI+2.0*H2)*g.mass*amu_per_gm*gm_per_msol*dMsolUnit[i]
    h2marray = (2.0*H2*g.mass)*amu_per_gm*gm_per_msol*dMsolUnit[i]
;     column = 1.0/mach*dKpcUnit[i]*cm_per_kpc
    rho = alog10(g.dens * dens_convert * (HI + 2.0*H2))
    ind = where(g.tempg lt 1e4)
    zmetal = TOTAL(g[ind].zmetal*g[ind].mass)/TOTAL(g[ind].mass)
    print,'Total Hydrogen Mass: ',TOTAL((2.0*H2)*g.mass*dMsolUnit)*f_H,' Solar Masses'
    print,'Total Hydrogen Mass: ',TOTAL((HI + 2.0*H2)*g.mass*dMsolUnit)*f_H,' Solar Masses'
    print,'Mean Metallicity: ',zmetal/0.0177 ,' Zsol'
    zsol = 0.0177
    
    IF i eq 0 THEN BEGIN
        plot,[1e19,1e24],[0,1],psym = 3,/ylog,xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.0],xrange = [1e19,1e23],xstyle=1,ystyle=1,symsize =0.5,/xlog
    ENDIF
    xmin  = -1.0*range[i]/dKpcUnit[i]
    xmax  =  1.0*range[i]/dKpcUnit[i] 
    ymin  = -1.0*range[i]/dKpcUnit[i] 
    ymax  =  1.0*range[i]/dKpcUnit[i]  
    nx = (xmax - xmin)/delta[i]
    ny = (ymax - ymin)/delta[i]    
    grid_sd = fltarr(nx,ny)
    grid_frac = grid_sd
    grid_z = grid_sd
    grid_sfr = grid_sd
    xarray = findgen(nx + 1)*delta[i] + xmin
    yarray = findgen(ny + 1)*delta[i] + ymin
    deltat = 100.e6
    timeunit = 3.872d3*SQRT((dKpcUnit[i]*3.0856805d21)^3/(dMsolUnit[i]*1.98892d33))/31556926.0
    tform=s.tform*(timeunit)[0]
    tcurrent = 1e10 ;max(tform)
    tclip = tcurrent - deltat
    FOR ix = 0, nx - 1 DO BEGIN
        FOR iy = 0, ny - 1 DO BEGIN
            ind = where(g.x gt xarray[ix] AND g.x lt xarray[ix + 1] AND g.y gt yarray[iy] AND g.y lt yarray[iy + 1])
            inds = where(s.x gt xarray[ix] AND s.x lt xarray[ix + 1] AND s.y gt yarray[iy] AND s.y lt yarray[iy + 1] AND tform gt tclip)           
            if ind[0] NE -1 THEN BEGIN
                grid_frac[ix,iy] = TOTAL(h2marray[ind])/TOTAL(hmarray[ind])
                grid_sd[ix,iy] = TOTAL(hmarray[ind])/cm_per_kpc/cm_per_kpc/dKpcUnit[i]/dKpcUnit[i]/delta[i]/delta[i]  
                grid_z[ix,iy] = TOTAL(g[ind].zmetal*g[ind].mass)/TOTAL(g[ind].mass)
                oplot,[grid_sd[ix,iy],grid_sd[ix,iy]],[grid_frac[ix,iy],grid_frac[ix,iy]],psym = 2,color = colors[i]
            ENDIF ELSE BEGIN
                grid_frac[ix,iy] = 0
                grid_sd[ix,iy] = 0
                grid_z[ix,iy] = 0
            ENDELSE
            if inds[0] NE -1 THEN BEGIN
                grid_sfr[ix,iy] = TOTAL(s[inds].mass*dMsolUnit[i])/dKpcUnit[i]/dKpcUnit[i]/delta[i]/delta[i]/deltat
            ENDIF ELSE BEGIN
                grid_sfr[ix,iy] = 0
            ENDELSE
          ENDFOR
    ENDFOR    
    stop
    IF 0 THEN BEGIN
        openr,1,dir[i] + '/H2clouds.txt'
        nlines = 5
        clouds = fltarr(3,nlines)
        readf,1,clouds
        FOR ic = 0, nlines - 1 DO BEGIN
            distance = SQRT((g.x - clouds[0,ic])*(g.x - clouds[0,ic]) + (g.y - clouds[1,ic])*(g.y - clouds[1,ic]))
            ind = where(distance lt clouds[2,ic]/2.0)
            frac = TOTAL(h2marray[ind])/TOTAL(hmarray[ind])
            sden = TOTAL(hmarray[ind])/cm_per_kpc/cm_per_kpc/dKpcUnit[i]/dKpcUnit[i]/(clouds[2,ic]/2.0)/(clouds[2,ic]/2.0)
            oplot,[sden,sden],[frac,frac],psym = 4
        ENDFOR
    ENDIF

;    plot,alog10(grid_sd*cm_per_kpc*cm_per_kpc/6.022142d23/gm_per_msol/1000.0/1000.0),alog10(grid_sfr),psym = 2,xrange = [-1,2.5],yrange = [-6,0],xtitle = 'Log Gas Surface Den',ytitle= 'Log SFR Surface Den'
;    plot,alog10(grid_frac*grid_sd*cm_per_kpc*cm_per_kpc/6.022142d23/gm_per_msol/1000.0/1000.0),alog10(grid_sfr),psym = 2,xrange = [-1,2.5],yrange = [-6,0],xtitle = 'Log H2 Surface Den',ytitle= 'Log SFR Surface Den'
ENDFOR
legend,key,psym = 2,color = colors
if (KEYWORD_SET(outplot)) then device,/close
END
