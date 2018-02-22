;This plots the resolved schmidt law, a la Bigel

PRO master_resSchmidtLaw
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/"
dir = prefix + ['h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/',$
        'h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir/',$
        'h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir/']

file = ["h516.cosmo25cmb.1536g3HBWK.00512.halo.1","h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.halo.1","h516.cosmo25cmb.1536g6MbwK.00512.halo.1"]
msol_per_sysmass       =  [2.310e15,2.310e15,2.310e15]
kpc_per_syslength        = [25000.,25000.,25000.]
keys = ['Standard','No H2 SF','Fabio'] 
;resSchmidtLaw,dir,file,kpc_per_syslength,msol_per_sysmass,keys,outplot='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/h516_512_h2sf_SK_nonobs.eps'

cd,prefix
base0 = "h516.cosmo25cmb.1536g3HBWK"
base1 = "h516.cosmo25cmb.1536g6MbwK"
base2 = "h516.cosmo25cmb.1536g6HBWK.jeans.prev"
base3 = "h516.cosmo25cmb.1536g3HBWK" 
steps = [12,24,36,48,60,72,84,96,108,120,132,180,192,204,216,228,240,252,264,276,288,300,312,324,336,348,360,372,384,396,408,420,432,444,456,464]
step0 = ['00195','00240','00324','00384','00408','00456','00480','00512']
step1 = ['00192','00240','00328','00384','00406','00455','00480','00512'];192
step2 = ['00144','00192','00276','00336','00360','00408','00432','00464'];144
step3 = ['00192','00240','00324','00384','00408','00456','00480','00512']

nsteps = N_ELEMENTS(step0) - 1

cd,prefix
start = 0
dt = 1
keys = 'fracLim'
keys = ['no H2','H2','H2 SF','MC SF']
msol_per_sysmass       =  (fltarr(N_ELEMENTS(keys))+ 1)*2.310e15
kpc_per_syslength        =  (fltarr(N_ELEMENTS(keys))+ 1)*25000.


FOR i = start/dt, nsteps DO BEGIN 
;    step = i*dt
;    step = steps[i]
;    if (step lt 10) THEN step = '0000'+STRTRIM(step,2) ELSE BEGIN
;        if (step lt 100) THEN step = '000'+STRTRIM(step,2) ELSE step = '00'+STRTRIM(step,2)
;    ENDELSE
;    filename0 = base+"/Jeans_oldLW/steps/"+base+"."+step+".dir/"+base + "." + step+'.halo.1';+STRTRIM(halo,2)

    dir0 = prefix + base0 + "/steps/"+base0+"."+step0[i]+".dir/"
    file0 = base0 + "." + step0[i]+'.halo.1'
    dir1 = prefix + base1 + "/steps/"+base1+"."+step1[i]+".dir/"
    file1 = base1 + "." + step1[i]+'.halo.1'
    dir2 = prefix + "h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/"+base2+"."+step2[i]+".dir/"
    file2 = base2 + "." + step2[i]+'.halo.1' ;+STRTRIM(halo,2)
    dir3 = prefix + base3 + "/steps_noH2SF/"+base3+"_noH2SF."+step3[i]+".dir/"
    file3 = base0 + "_noH2SF." + step3[i]+'.halo.1'
    dir = [dir1,dir3,dir0,dir2]
    file = [file1,file3,file0,file2]
    colors = [30,90,120,240]
    
;    compSF,filename,slfile,dMsolUnit,dKpcUnit,key =key,color = [30,90,120,240],outplot = 'stellarprof.'+step3[i]+'.eps';,redshift = redshift[i];,outfile = outfile,/color
;    stop
    resSchmidtLaw,dir,file,kpc_per_syslength,msol_per_sysmass,keys,colors = colors;,outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516_SK_nonobs'+step0[i]+'.eps'
    stop
ENDFOR

END

;resSchmidtLaw,'/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK.00512.dir/','h603.cosmo50cmb.3072gs1MbwK.00512.halo.1',50000,1.84793e16,outplot= '/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK.00512.dir/SK.eps'
PRO resSchmidtLaw,dir,file,dKpcUnit,dMsolUnit,key,outplot = outplot,colors = colors
set_plot,'ps'
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790d+23
molec_weight = (0.76*1 + 0.24*4.0)
loadct,39
f_H = 0.764
IF NOT KEYWORD_SET(colors) THEN colors  = (findgen(N_ELEMENTS(file)) + 1)*240/N_ELEMENTS(file)
;colors = [240,160,80] ;270 - (alog10(s.tform) - min)/(MAX(alog10(s.tform)) - min)*254

!p.multi = 0
if (KEYWORD_SET(outplot)) then begin
    !P.CHARTHICK=4.0
    !X.THICK=4
    !Y.THICK=4
    !P.charsize=1.5
    !x.charsize=1.0;0.75
    !y.charsize=1.0;0.75
    set_plot,'ps' 
    device,filename = outplot,/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2 
endif else begin
    set_plot,'x'
    window,0,xsize = 600,ysize = 600
endelse

xsigmalow = findgen(300)/100 - 1.
xsigma = 10.0^(findgen(400)/100 + 1.)
ysigma=2.5e-4*xsigma^1.4
ysigma1=2.5e-4*xsigma^(1.4 + 0.15)
ysigma2=2.5e-4*xsigma^(1.4 - 0.15)
ysigmalow = xsigmalow*2.4 - 5.0
range = 11.25*(fltarr(N_ELEMENTS(file))+1)
distance2 = 4*(fltarr(N_ELEMENTS(file))+1)
delta = 0.750;/dKpcUnit

IF 0 THEN BEGIN
    FOR i = 0, N_ELEMENTS(file) - 1 DO BEGIN
        cd,dir[i]
        rtipsy,file[i],h,g,d,s
        g.x = g.x * dKpcUnit[i]*h.time
        g.y = g.y * dKpcUnit[i]*h.time
        g.z = g.z * dKpcUnit[i]*h.time
        g.mass = g.mass*dMsolUnit[i]
        s.x = s.x * dKpcUnit[i]*h.time
        s.y = s.y * dKpcUnit[i]*h.time
        s.z = s.z * dKpcUnit[i]*h.time
        s.mass = s.mass*dMsolUnit[i]
        s.tform = s.tform*SQRT((dKpcUnit[i]*3.086d21)^3/(6.67d-8*dMsolUnit[i]*1.99d33))/(3600.*24.*365.24)
                                ;   prof_s = prof(s,'star',MAX(s.tform),nbins = 50,rmax = 25)
        cuml_star = findgen(50)
        radii = (findgen(50) + 1)/50.0*25.0
        FOR ibin = 0, 49 DO BEGIN
            cuml_star[ibin] = TOTAL(s[where(SQRT(s.x*s.x + s.y*s.y + s.z*s.z) lt radii[ibin])].mass)
        ENDFOR
        IF (i EQ 0) THEN plot,radii,cuml_star,/ylog,xtitle = 'Radius [kpc]',ytitle = 'Mass of Interior Stars',xrange = [0,25]
        oplot,radii,cuml_star,color = colors[i]
        oplot,[0,25],[cuml_star[49]*0.75,cuml_star[49]*0.75],color = colors[i],linestyle = 2
    ENDFOR
ENDIF

FOR i = 0, N_ELEMENTS(file) - 1 DO BEGIN
    cd,dir[i]
    IF 0 THEN BEGIN
        cubeH2 = read_dencube_fits(file[i]+'.H2.lr.fits',headerH2)
        pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*3.08568021d21*3.08568021d21
        cubeH2_surface_den = cubeH2/pix_areaH2*gm_per_msol*amu_per_gm*f_H
    
        cubeHI = read_dencube_fits(file[i]+'.HI.lr.fits',headerHI)
        pix_areaHI = headerHI.CDELT1*headerHI.CDELT2*3.08568021d21*3.08568021d21
        cubeHI_surface_den = cubeHI/pix_areaHI*gm_per_msol*amu_per_gm*f_H

        xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1
        yaxes = (findgen(headerHI.NAXIS2) - headerH2.NAXIS1/2.0)*headerH2.CDELT2
    
        distance = fltarr(headerHI.NAXIS1,headerHI.NAXIS2)
        FOR ix = 0, headerHI.NAXIS1 - 1 DO $
          FOR iy = 0, headerHI.NAXIS2 - 1 DO distance[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy])
        ind = where(distance lt distance2[i])
        fraction = 2.0*cubeH2_surface_den/(cubeHI_surface_den + 2.0*cubeH2_surface_den)
    ENDIF
    rtipsy,file[i],h,g,d,s
    g.x = g.x *h.time*dKpcUnit[i]
    g.y = g.y *h.time*dKpcUnit[i]
    g.z = g.z *h.time*dKpcUnit[i]
    s.x = s.x *h.time*dKpcUnit[i]
    s.y = s.y *h.time*dKpcUnit[i]
    s.z = s.z *h.time*dKpcUnit[i]
;    dens_convert =  dMsolUnit[i]*gm_per_msol*amu_per_gm/dKpcUnit[i]^3/cm_per_kpc^3/h.time/h.time/h.time
    readarr,file[i]+'.HI',h,hI,part = 'gas',/ascii
    IF (FILE_TEST(file[i]+".H2")) THEN BEGIN
        readarr,file[i]+'.H2',h,h2,part = 'gas',/ascii
    ENDIF ELSE h2 = fltarr(N_ELEMENTS(hI))
    marray = g.mass*amu_per_gm*gm_per_msol*dMsolUnit[i]
    hmarray = (HI+2.0*H2)*g.mass*dMsolUnit[i] ;*amu_per_gm*gm_per_msol
    h2marray = (2.0*H2*g.mass)*dMsolUnit[i] ;*amu_per_gm*gm_per_msol
 ;   rho = alog10(g.dens * dens_convert * (HI + 2.0*H2))
    ind = where(g.tempg lt 1e4)
    zmetal = TOTAL(g[ind].zmetal*g[ind].mass)/TOTAL(g[ind].mass)
    print,'Total H2 Mass: ',TOTAL((2.0*H2)*g.mass*dMsolUnit)*f_H,' Solar Masses'
    print,'Total HI + H2 Mass: ',TOTAL((HI + 2.0*H2)*g.mass*dMsolUnit)*f_H,' Solar Masses'
    print,'Mean Metallicity: ',zmetal/0.0177 ,' Zsol'
;    stop
    xmin  = -1.0*range[i]
    xmax  =  1.0*range[i]
    ymin  = -1.0*range[i] 
    ymax  =  1.0*range[i] 
    nx = (xmax - xmin)/delta
    ny = (ymax - ymin)/delta    
    grid_sd = fltarr(nx,ny)
    grid_frac = grid_sd
    grid_z = grid_sd
    grid_sfr = grid_sd
    xarray = findgen(nx + 1)*delta + xmin
    yarray = findgen(ny + 1)*delta + ymin
    deltat = 100.e6 ;Appearently, in Bigiel, this is based on the FUV flux, so OB stars with have a lifetime of 10^8 years
    deltat = 50.e6
    
    timeunit = 3.872d3*SQRT((dKpcUnit[i]*3.0856805d21)^3/(dMsolUnit[i]*1.98892d33))/31556926.0
    tform=s.tform*(timeunit)[0]
    tcurrent = max(tform) ;1e10
    tclip = tcurrent - deltat

;    window,0
;    inds = where(tform gt tclip)
;    plot,s[inds].x*dkpcUnit[0],s[inds].y*dKpcUnit[0],psym = 3,xrange = [-3,3],yrange = [-3,3]

;    window,1
    if (i eq 0) then begin
            plot,[-5,-5],[-5,-5],xrange = [-7,5],yrange = [-5,3],xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]",xstyle = 1,ystyle = 1
            oplot,[1,1],[-5,3],linestyle = 1
            oplot,xsigmalow,ysigmalow
            oplot,alog10(xsigma),alog10(ysigma)
    ENDIF


;    stop
    FOR ix = 0, nx - 1 DO BEGIN
        FOR iy = 0, ny - 1 DO BEGIN
            ind = where(g.x gt xarray[ix] AND g.x lt xarray[ix + 1] AND g.y gt yarray[iy] AND g.y lt yarray[iy + 1])
            inds = where(s.x gt xarray[ix] AND s.x lt xarray[ix + 1] AND s.y gt yarray[iy] AND s.y lt yarray[iy + 1] AND tform gt tclip)           
            if ind[0] NE -1 THEN BEGIN
                grid_frac[ix,iy] = TOTAL(h2marray[ind])/TOTAL(hmarray[ind])
                grid_sd[ix,iy] = TOTAL(hmarray[ind])/delta/delta/1000.0/1000.0;/cm_per_kpc/cm_per_kpc/dKpcUnit[i]/dKpcUnit[i]/delta[i]/delta[i]  
                grid_z[ix,iy] = TOTAL(g[ind].zmetal*g[ind].mass)/TOTAL(g[ind].mass)
            ENDIF ELSE BEGIN
                grid_frac[ix,iy] = 0
                grid_sd[ix,iy] = 0
                grid_z[ix,iy] = 0
            ENDELSE
            if inds[0] NE -1 THEN BEGIN
                grid_sfr[ix,iy] = TOTAL(s[inds].mass*dMsolUnit[i])/delta/delta/deltat
            ENDIF ELSE BEGIN
                grid_sfr[ix,iy] = 0
            ENDELSE
            oplot,[alog10(grid_sd[ix,iy]),alog10(grid_sd[ix,iy])],[alog10(grid_sfr[ix,iy]),alog10(grid_sfr[ix,iy])],psym = 2,color = colors[i]
;            oplot,[alog10(grid_sd[ix,iy]*cm_per_kpc*cm_per_kpc/6.022142d23/gm_per_msol/1000.0/1000.0),alog10(grid_sd[ix,iy]*cm_per_kpc*cm_per_kpc/6.022142d23/gm_per_msol/1000.0/1000.0)], $
;                  [alog10(grid_sfr[ix,iy]),alog10(grid_sfr[ix,iy])],psym = 2,color = colors[i]
 ;           print,alog10(grid_sd[ix,iy]),alog10(grid_sfr[ix,iy])
 ;           stop
        ENDFOR
    ENDFOR    
ENDFOR
legend,key,psym = 2,color = colors
if (KEYWORD_SET(outplot)) then device,/close
stop
END
