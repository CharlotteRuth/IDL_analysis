
function enclosedL_R,filename = filename,outplot = outplot,verbose = verbose,extno = extno, sr_range = sr_range,F_num = F_num,zero_point = zero_point,center = center,halo_str = halo_str, level = level
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
speed_o_light = 299792458 ;m per sec
molec_weight = (0.76*1 + 0.24*4.0)
ster_to_arcs2 = 4.25451703e10 
Fsun = 3.196e-10 ; erg/s/cm/cm Flux from Sun @ 10pc
Msun = 4.76 ; absolute mag of sun

IF NOT keyword_set(center) THEN center = [0,0] ELSE center = center
IF NOT keyword_set(halo_str) THEN halo_str = '1'
IF NOT keyword_set(level) THEN level = 0.9 ;90%

IF keyword_set(verbose) THEN BEGIN
    loadct,0
    nbins=100.0
    linestyles = [0,2]
    IF keyword_set(outplot) THEN formatplot,/outplot ELSE formatplot
ENDIF

;****************** Sunrise Cube **************************\
;IF keyword_SET(filename) THEN sfilename = filename+'.' + halo_str + '/broadband.fits' ELSE sfilename = 'broadband.fits'
IF keyword_SET(filename) THEN spawn,'ls ' +  filename+'.' + halo_str + '/*broadband.fits',sfilenames $ 
ELSE spawn,'ls *broadband.fits',sfilenames
sfilename = sfilenames[0]
print,sfilename

IF NOT keyword_set(extno) THEN extno = 14 ;Faceon
IF NOT keyword_set(F_num) THEN F_num = 4 ;SDSS r band
sr_surface_den = mrdfits(sfilename,extno,sr_head)
dxy = sr_head[where(strcmp('CD1_1',sr_head,5))]
dxy = double((strsplit(dxy,' ',/extract))[2])
sr_NAXIS = (size(sr_surface_Den))[1] ;500 ;480.0
IF NOT keyword_set(sr_range) THEN sr_range = dxy*sr_NAXIS ;headerHI.NAXIS1*headerHI.CDELT1

xaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0 - center[0])*sr_range/sr_NAXIS + sr_range/sr_NAXIS/2.0
yaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0 - center[1])*sr_range/sr_NAXIS + sr_range/sr_NAXIS/2.0
print,minmax(xaxes_sr)
print,minmax(yaxes_sr)
radius = fltarr(n_elements(xaxes_sr),n_elements(yaxes_sr))
FOR ix = 0,n_elements(xaxes_sr) - 1 DO $
  FOR iy = 0,n_elements(yaxes_sr) - 1 DO $
  radius[ix,iy] = SQRT(xaxes_sr[ix]*xaxes_sr[ix] + yaxes_sr[iy]*yaxes_sr[iy])

surface_den = sr_surface_den[*,*,F_num]/ster_to_arcs2*1e-7 ;Converting W/m/m^2/sr to ergs/s/A/cm^2/arcsec^2

;Set up circles of increasing radii 
dr = FLOAT(sr_range)/sr_NAXIS*4.0/2.0 ;0.25  ;difference in radii
radius_arr = (findgen(sr_range/dr/2.0) + 1)*dr  ;array of radii
surface_den_enclosed = fltarr(sr_range/dr/2.0) ;Set up an array to hold the surface brightnesses
FOR ir = 0, n_elements(radius_arr) - 1 DO BEGIN
  ind = where(radius LT radius_arr[ir])  ;select smaller than the given distance
  surface_den_enclosed[ir] = TOTAL(surface_den(ind))
ENDFOR 
surface_den_enclosed = surface_den_enclosed/surface_den_enclosed[n_elements(radius_arr) - 1]


IF keyword_set(verbose) THEN BEGIN
    window,1 ;plot average surface brightness against radii
    plot,radius_arr,surface_den_enclosed,xtitle = 'Radius [kpc]',ytitle = 'Normalized Cumulative Flux Distribution',yrange = [0,1],xrange = [0,20],title = filename;,yrange = [30,MIN(m_B_surface_den_annuli)]
    oplot,[0,20],[level,level],linestyle = 2
ENDIF

;iords = where(surface_den_enclosed gt mag + 1)
temp = min(abs(surface_den_enclosed - level),iord) ;Find the annuli in which the surface brightness is closest to 26 and limit everything to there
IF (n_elements(surface_den_enclosed) GT iord + 3 AND iord - 3 GE 0) THEN BEGIN
    surface_den_enclosed_spline = surface_den_enclosed[iord - 3: iord + 3]
    radius_arr_spline = radius_arr[iord - 3: iord + 3]
ENDIF ELSE BEGIN
    IF (n_elements(surface_den_enclosed) LT iord + 3) THEN BEGIN
        surface_den_enclosed_spline = surface_den_enclosed[iord - 3: n_elements(surface_den_enclosed) - 1]
        radius_arr_spline = radius_arr[iord - 3: n_elements(surface_den_enclosed) - 1]
    ENDIF ELSE BEGIN
        surface_den_enclosed_spline = surface_den_enclosed[0: iord + 3]
        radius_arr_spline = radius_arr[0: iord + 3]
    ENDELSE
ENDELSE
order  = sort(surface_den_enclosed_spline) ;sort by brightness
rlimit = spline(surface_den_enclosed_spline[order],radius_arr_spline[order],level)
temp   = radius_arr_spline[order]
IF rlimit GT temp[n_elements(temp) - 1] THEN rlimit = radius_arr[iord]

IF keyword_set(verbose) THEN BEGIN
    oplot,[rlimit,rlimit],[0,1],linestyle = 2

    window,0,xsize = 800,ysize = 800
;    print,minmax(B_surface_Den[where(B_surface_den ne 0)])
;    print,minmax(m_B_surface_den[where(m_B_surface_den ne 0)])
    range = minmax(alog10(surface_den[where(surface_den NE 0)]))
    max = ceil(range[1])
    min = floor(range[0])
    contour,alog10(surface_den),xaxes_sr,yaxes_sr,min = range[0],max = range[1],xrange = [-7,7],yrange = [-7,7],/fill,nlevels = 245
    contour,radius,xaxes_sr,yaxes_sr,/overplot,levels = [1,2,3,4,5,6,7,8,9,10,11,12],color = 100
    contour,radius,xaxes_sr,yaxes_sr,/overplot,levels = [0,rlimit]    
ENDIF

RETURN,rlimit
END

