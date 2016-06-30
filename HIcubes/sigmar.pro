;Comparing to Tamburro et al 2009

pro sigmar_master 
;filename = 'h603.cosmo50cmb.2304g5bwK.00512'
;cd,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG'

IF KEYWORD_SET(vetinari) THEN BEGIN
   prefix = '~/Data/MolecH/Cosmo/' 
   outfile = '~/Figures/h516.cosmo25cmb.paper_'
ENDIF ELSE BEGIN
   prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
   outfile = '~/plots/h516.cosmo25cmb.paper_'   
ENDELSE

cd,prefix+'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir'
filename = 'h516.cosmo25cmb.2304g14HBWK.00512'
pfile = 'h516.cosmo25cmb.2304g14HBWK.param'
;outplot = prefix+'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/h516.cosmo25cmb.2304g14HBWK.00512.halo.1_' + strtrim(res,2)
;schmidtlaw_res_obs,filename,pfile,res = res,/useH2,outplot = outplot

cd,prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir'
filename = 'h516.cosmo25cmb.3072g1MBWK.00492'
pfile = 'h516.cosmo25cmb.3072g1MBWK.param'
;outplot = prefix+'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1_' + strtrim(res,2)

cd,prefix+'/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir'
filename = 'h516.cosmo25cmb.3072g14HBWK.00512'
pfile = '../../h516.cosmo25cmb.3072g14HBWK.param'
;outplot = prefix+'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1_' + strtrim(res,2)

sigmar,filename,pfile,/outplot ;= outplot
end

pro sigmar,filename,pfile,keys = keys, colors = colors, thicks = thicks, linestyle= linestyle,label = label,ctables = ctables,outplot = outplot,camera = camera, intq = intq,center = center,halo_str = halo_str,sfr = sfr
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
speed_o_light = 299792458 

n = N_ELEMENTS(filename)
formatplot,outplot = outplot
cb_charsize = 0.75
IF KEYWORD_SET(outplot) THEN BEGIN
    fg = 0
ENDIF ELSE BEGIN
    fg = 255
ENDELSE
!p.multi = 0
!X.MARGIN = [12,12]


angle = !pi/2 ;!pi/4
nbins = 60
IF NOT KEYWORD_SET(camera) THEN extno = 14 ELSE extno = camera
IF NOT KEYWORD_SET(intq) THEN extno_int = 32 ELSE extno_int = intq
IF NOT KEYWORD_SET(center) THEN center = [0,0] ELSE center = center
IF NOT KEYWORD_SET(halo_str) THEN halo_str = 1

formatplot,outplot = outplot
IF KEYWORD_SET(outplot) THEN     device,filename=outplot,/color,bits_per_pixel= 8,/times,ysize = 12,xsize = 18,xoffset =  2,yoffset =  2 ELSE window,0
IF KEYWORD_SET(outplot) THEN     fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(outplot) THEN     bgcolor = 255 ELSE bgcolor = 0
IF KEYWORD_SET(colors) THEN BEGIN
    loadct,39
    if colors[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = colors
    if NOT keyword_set(ctables) then ctables = [39,39,39]
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT KEYWORD_SET(symbols) THEN symbols = (fltarr(n)+4)
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    if NOT keyword_set(ctables) then ctables = [0,0,0]
    IF NOT KEYWORD_SET(thicks) THEN thicks = (findgen(n) + 1)*6/n - 1
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = REVERSE(findgen(n)*2)   
    IF NOT KEYWORD_SET(symbols) THEN symbols = (findgen(n)+2)*2
ENDELSE

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
lb_linestyle = intarr(n)

FOR ifile = 0, n - 1 DO BEGIN
    units = tipsyunits(pfile[ifile])
    massunit = units.massunit
    lengthunit = units.lengthunit
    r25 = opticalRadii(filename = filename[ifile],verbose = 0,extno = extno[ifile], center = center[*,ifile])
    maxr = r25*4
    radii_bins = findgen(nbins + 1)*maxr/nbins

;*************** Tipsy Data ***********************
    rtipsy,filename[ifile] + '.halo.1.std',h,g,d,s
    g.x = g.x*units.lengthunit*h.time
    g.y = g.y*units.lengthunit*h.time
    g.z = g.z*units.lengthunit*h.time
    g.mass = g.mass*massunit
    readarr,filename[ifile] + '.halo.1.HI',h,HI,/ascii,part = 'gas'
    
;*************** Rotation ****************
;angle0 = ATAN(s.z/s.y)
;angle0[where(s.y lt 0)] = angle0[where(s.y lt 0)]+!PI
;s1 = s
;Ryz = SQRT(s.y*s.y + s.z*s.z)
;s1.y = Ryz*COS(angle0 + angle)
;s1.z = Ryz*SIN(angle0 + angle)
;s1.x = s1.x
;angle = 1.0*!PI/4.0
;angle0 = ATAN(g.z/g.y)
;angle0[where(g.y lt 0)] = angle0[where(g.y lt 0)]+!PI
;g1 = g
;Ryz = SQRT(g.y*g.y + g.z*g.z)
;g1.y = Ryz*COS(angle0 + angle)
;g1.z = Ryz*SIN(angle0 + angle)
;g1.x = g1.x

;***************** Gas Cubes **********************
    filebase = filename + '.halo.1.cube'
    filebase = filename + '.halo.1.smoothed'
    filebase = filename + '.halo.1.thermal.smoothed'
    filebase = filename[ifile] + '.halo.1.90.smoothed'
    print,filebase
    cubeHI_surface_den = read_dencube_fits(filebase + '.mom0.fits',headerHI,/noscale)
    cubeHI_sigma = read_dencube_fits(filebase + '.mom2.fits',headerHI,/noscale) ;/amu_per_gm/gm_per_msol*cm_per_kpc*cm_per_kpc/1d6
    xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1
    yaxes = (findgen(headerHI.NAXIS2) - headerHI.NAXIS1/2.0)*headerHI.CDELT2
    radius = fltarr(N_ELEMENTS(xaxes),N_ELEMENTS(yaxes))
    FOR ix = 0,N_ELEMENTS(xaxes) - 1 DO $
      FOR iy = 0,N_ELEMENTS(yaxes) - 1 DO $
      radius[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy]/SIN(angle)/SIN(angle))
    
;window, 0, xsize=700,ysize=570
    loadct, 39
    nl2 = 6                     ; HACK ADDED BY AMS
    cvec = dindgen(nl2) / nl2 * 250 + 20
    nlevels = 10.0
    levs2 = 30.0*(findgen(nlevels))/(nlevels - 1) 
    titl = 'HI Velocity Dispersion (km/s)'
;contour, cubeHI_sigma, xaxes, yaxes, levels=levs2, xtitle='kpc', ytitle='kpc', $
;  title=titl, position=[0.1,0.1,0.75,0.90],/cell_fill, xstyle=1,ystyle=1,$
;  c_colors=cvec,xrange = [-maxr,maxr],yrange = [-maxr,maxr]
;contour,radius, xaxes, yaxes,levels = radii_bins*r25,/overplot,nlevels = N_ELEMENTS(radii_bins*r25),c_linestyle = 0,c_thick = 0.5
;colorbar, /right, /vertical, position=[0.1,0.80,0.90,0.85], $
;  format='(%"%7.3g")', range=[levs2[0], levs2[n_elements(levs2)-1]],$
;  divisions=6, bottom=cvec[0],ncolors=cvec[n_elements(cvec)-1] - cvec[0]


;**************** Sunrise **********************
    IF KEYWORD_SET(sfr) THEN BEGIN
        sfilename = filename + '.1/broadband.fits'
        fuv_num = 0
        mic24_num = 15
        fuv_lambda = 1.5180502483701E-07
        mic24_lambda = 2.38717045592424E-05
        extno = 15
        sr_surface_den = mrdfits(sfilename,extno,sr_head)
        
;Sunrise units are in W/m/n^2/sr
;Bigiel Units are in MJy*str^(-1)
        fuv_surface_den = sr_surface_den[*,*,fuv_num]*fuv_lambda*fuv_lambda/speed_o_light/1e-20 
        mic24_surface_den = sr_surface_den[*,*,mic24_num]*mic24_lambda*mic24_lambda/speed_o_light/1e-20
        SFR_surface_den = (3.2e-3*mic24_surface_den + 8.1e-2*fuv_surface_den)/1e6 ;EQ 1, Bigiel et al, 2008; then converted to per pc^2 rather than kpc^2
        fuv_surface_den = ROTATE(fuv_surface_den,5)
        mic24_surface_den = ROTATE(mic24_surface_den,5)
        SFR_surface_den_log = ALOG10(ROTATE(SFR_surface_den,5))
        sr_range = 24           ;headerHI.NAXIS1*headerHI.CDELT1
        sr_NAXIS = 480.0
        xaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0)*sr_range/sr_NAXIS
        yaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0)*sr_range/sr_NAXIS
        deltaS_SR = sr_range/sr_NAXIS*sr_range/sr_NAXIS
        radius_sr = fltarr(N_ELEMENTS(xaxes_sr),N_ELEMENTS(yaxes_sr))
        FOR ix = 0,N_ELEMENTS(xaxes_sr) - 1 DO $
          FOR iy = 0,N_ELEMENTS(yaxes_sr) - 1 DO $
          radius_sr[ix,iy] = SQRT(xaxes_sr[ix]*xaxes_sr[ix] + yaxes_sr[iy]*yaxes_sr[iy]/SIN(angle)/SIN(angle))
        
;window, 1, xsize=700,ysize=570
        loadct, 39
        nl2 = 6                 ; HACK ADDED BY AMS
        cvec = dindgen(nl2) / nl2 * 250 + 20
        nlevels = 10.0
        levs2 = 5.0*(findgen(nlevels))/(nlevels - 1) - 10.0 
        titl = 'SFR from FUV and 24 micron M_sol/yr/pc/pc'
;contour, SFR_surface_den_log, xaxes_sr, yaxes_sr, xtitle='kpc', ytitle='kpc', $
;  title=titl, position=[0.1,0.1,0.75,0.90],/cell_fill, xstyle=1,ystyle=1,$
;  c_colors=cvec,xrange = [-maxr,maxr],yrange = [-maxr,maxr],levels=levs2
;contour,radius, xaxes, yaxes,levels = radii_bins*r25,/overplot,nlevels = N_ELEMENTS(radii_bins*r25),c_linestyle = 0,c_thick = 0.5
;colorbar, /right, /vertical, position=[0.1,0.80,0.90,0.85], $
;  format='(%"%7.3g")', range=[levs2[0], levs2[n_elements(levs2)-1]],$
;  divisions=6, bottom=cvec[0],ncolors=cvec[n_elements(cvec)-1] - cvec[0]

;**************** Gas Dispersion & SFR
        sfr = fltarr(nbins)
    ENDIF
    sigmar = fltarr(nbins)
    sigmar_norm = sigmar
    radii_x = findgen(nbins)*maxr/nbins + maxr/nbins/2.0

    for i = 0 , nbins - 1 DO BEGIN 
        ind = where(radius ge radii_bins[i] and radius lt radii_bins[i + 1])
        IF radii_bins[0] ne -1 THEN BEGIN
        sigmar[i] = MEAN(cubeHI_sigma[ind])
        sigmar_norm[i] = TOTAL(cubeHI_sigma[ind]*cubeHI_surface_den[ind])/TOTAL(cubeHI_surface_den[ind])
        ELSE BEGIN
            sigmar[i] = 0
            sigmar_norm[i] = 0
        ENDELSE
        IF KEYWORD_SET(sfr) THEN BEGIN
            indsr = where(radius_sr ge radii_bins[i] and radius_sr lt radii_bins[i + 1])
            sfr[i] = MEAN(SFR_surface_den_log[indsr])
        ENDIF
;    stop
    ENDFOR
    loadct,ctables[ifile]
    if ifile eq 0 THEN plot,radii_x/r25,sigmar,xrange = [0,4],yrange = [0,25],psym = -5,xtitle = textoidl('r/r_{25}'),ytitle = textoidl(' \sigma_{HI} [km s^{-1}]'),/nodata
    oplot,radii_x/r25,sigmar,thick = thicks[iFile],linestyle = linestyle[iFile],color = colors[iFile]

;oplot,radii_x/r25,sigmar_norm,psym = -5
    IF KEYWORD_SET(SFR) THEN BEGIN
        if ifile eq 0 THEN plot,radii_x/r25,sigmar,xrange = [0,4],yrange = [0,25],psym = -5,xtitle = 'r/r_25',ytitle = 'sigma_HI [km s^-1]',YSTYLE=8,/nodata
        oplot,radii_x/r25,sigmar,thick = thicks[iFile],linestyle = linestyle[iFile],color = colors[iFile]
        AXIS,YAXIS = 1,YRANGE = [-10.5,-6.8],ytitle = 'log SFR [M_sun yr^-1 pc^-2]',/save,ystyle = 1
        oplot,radii_x/r25,sfr,psym = -6,color = 50
        legend,['Sigma','SFR'],psym = [-5,-6],color = [fg,50],/right
;stop
    ENDIF
    IF KEYWORD_SET(keys) THEN BEGIN
        l_keys = lb_keys
        l_keys[iFile] = keys[iFile]
        l_thicks = lb_thicks
        l_thicks[iFile] = thicks[iFile]
        l_color = lb_color
        l_color[iFile] = colors[iFile]
        l_linestyle = lb_linestyle
        l_linestyle[iFile] = l_linestyle[iFile]
        legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/right,box = 0
    ENDIF
ENDFOR
    IF keyword_set(outplot) THEN device,/close

END
