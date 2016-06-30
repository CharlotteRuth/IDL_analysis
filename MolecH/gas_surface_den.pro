pro gas_surface_den,filename,outplot = outplot,useH2 = useH2,verbose = verbose,range = range,color = color, halo_str = halo_str,nocolorbar = nocolorbar,boxside = boxside
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
speed_o_light = 299792458 ;m per sec
molec_weight = (0.76*1 + 0.24*4.0)

IF NOT KEYWORD_SET(boxside) THEN boxside = 1; 2: bottom right, 3: top right, 4: top left
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps' 
    nbins=100.0
    linestyles = [0,2]
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25 
    !X.MARGIN = [12,5]
    !Y.MARGIN = [6,4]
    cb_charsize = 0.75
    black = 0
    white = 255
ENDIF ELSE BEGIN
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    set_plot,'x'
    nbins=100.0
    linestyles = [0,2]
    !X.MARGIN = [12,5]
    !Y.MARGIN = [6,4]
    cb_charsize = 0.75
    white = 255
    black = 0
ENDELSE
!p.multi = 0
IF NOT KEYWORD_SET(halo_str) THEN halo_str = '1'
blank = [' ',' ',' ',' ',' ']

loadct,0
rtipsy,filename,h,g,d,s,/justhead
;----------------- Read in gas cubes ----------------------
;cubeHI_surface_den = read_dencube_fits(filename +'.halo.1.90.cube.mom0.fits',headerHI,/noscale)
cubeHI_surface_den = read_dencube_fits(filename + '.halo.' + halo_str + '.90.smoothed.mom0.fits',headerHI,/noscale);*2.0
;/amu_per_gm/gm_per_msol*cm_per_kpc*cm_per_kpc/1d6
headerHI.CDELT1 = headerHI.CDELT1*h.time
headerHI.CDELT2 = headerHI.CDELT2*h.time
IF NOT (where(~finite(cubeHI_surface_den)))[0] EQ -1 THEN cubeHI_surface_den[where(~finite(cubeHI_surface_den))] = 0
IF keyword_set( useH2) THEN BEGIN
    cubeH2 = read_dencube_fits(filename+'.halo.' + halo_str + '.H2.arr.90.fits',headerH2)*2.0;*gm_per_msol*amu_per_gm
    headerH2.CDELT1 = headerH2.CDELT1*h.time
    headerH2.CDELT2 = headerH2.CDELT2*h.time
    pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*1e6;*cm_per_kpc*cm_per_kpc
    cubeH2_surface_den = TRANSPOSE(cubeH2)/pix_areaH2
    cubeH2_surface_den_plot = cubeH2_surface_den;/gm_per_msol/amu_per_gm*cm_per_kpc*cm_per_kpc
ENDIF

xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1 + headerHI.CDELT1/2
yaxes = (findgen(headerHI.NAXIS2) - headerHI.NAXIS1/2.0)*headerHI.CDELT2 + headerHI.CDELT1/2
dxcell = headerHI.cdelt1;/h.time
areapc = 1e6*dxcell*dxcell
areakpc = dxcell*dxcell
minrange = headerHI.CRVAL1
IF NOT KEYWORD_SET(range) THEN range = -1.0*minrange
radius = fltarr(headerHI.NAXIS1,headerHI.NAXIS2)
angle = !PI/4
FOR ix = 0,N_ELEMENTS(xaxes) - 1 DO $
  FOR iy = 0,N_ELEMENTS(yaxes) - 1 DO $
  radius[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy]/SIN(angle)/SIN(angle))

;------------------- Plot in amu/cc --------------
IF (KEYWORD_SET(outplot)) THEN device,filename=outplot+'_HI_H2sd2.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE BEGIN
    window,1,xsize = 800,ysize = 800
ENDELSE
position =  [0.85, 0.15, 0.9, 0.9]
mingas = 1e-2
maxgas = 1e2
nlevels = 255
userlevels = indgen(nlevels)*(alog10(maxgas) - alog10(mingas))/nlevels + alog10(mingas)
cubeHI_surface_den_plot = cubeHI_surface_den/gm_per_msol/amu_per_gm*cm_per_kpc/1000*cm_per_kpc/1000
ind = where(cubeHI_surface_den_plot lt mingas)
IF ind[0] NE -1 THEN cubeHI_surface_den_plot[ind] = mingas*1.01
contour,alog10(cubeHI_surface_den_plot),xaxes,yaxes,$
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,max_value = alog10(maxgas),min_value = alog10(mingas),levels = userlevels,xtickname = blank, ytickname = blank ;,title = textoidl('HI Surface Density [M')+sunsymbol()+textoidl(' pc^{-2}]')
;xtitle = 'X [kpc]',ytitle = 'Y [kpc]'
IF (boxside EQ 1 OR boxside EQ 2) $
  THEN axis,xaxis = 0,xtitle = 'X [kpc]' $
  ELSE axis,xaxis = 1,xtitle = 'X [kpc]'
IF (boxside EQ 1 OR boxside EQ 3) $
  THEN axis,yaxis = 0,ytitle = 'Y [kpc]' $
  ELSE axis,yaxis = 1,ytitle = 'Y [kpc]'
IF keyword_set(useH2) THEN BEGIN
    IF KEYWORD_SET(color) THEN loadct,39
    contour,alog10(cubeH2_surface_den_plot),xaxes,yaxes,/overplot,levels = [-2,0,2],color = 240;,c_linestyle = [1,4,0]
    loadct,0
ENDIF
IF NOT KEYWORD_SET(nocolorbar) THEN  colorbar,range = [alog10(mingas),alog10(maxgas)],/vertical,divisions = 4,position = position,color = white
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop

IF 0 THEN BEGIN
    IF (KEYWORD_SET(outplot)) THEN device,filename=outplot+'_H2sd.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,1,xsize = 800,ysize = 800
    loadct,0
    mingas = 1e-2
    maxgas = 1e2
    userlevels = indgen(nlevels)*(alog10(maxgas) - alog10(mingas))/nlevels + alog10(mingas)
    cubeH2_surface_den_plot = cubeH2_surface_den/gm_per_msol/amu_per_gm*cm_per_kpc/1000*cm_per_kpc/1000
    ind = where(cubeH2_surface_den_plot lt mingas)
    IF ind[0] NE -1 THEN cubeH2_surface_den_plot[ind] = mingas*1.01
;;cubeH2_surface_den_plot[where(NOT FINITE(cubeH2_surface_den_plot))] = 1e-2
    contour,alog10(cubeH2_surface_den_plot),xaxes,yaxes,$
      /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,max_value = alog10(maxgas),min_value = alog10(mingas),levels = userlevels ;,title = textoidl('H_2 Surface Density [M')+sunsymbol()+textoidl(' pc^{-2}]')
;,xtitle = 'X [kpc]',ytitle = 'Y [kpc]'
    
    IF NOT KEYWORD_SET(nocolorbar) THEN colorbar,range = [alog10(mingas),alog10(maxgas)],/vertical,divisions = 4,position = position,color = white
    IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
endif
if not keyword_set(color) then loadct,0
end
