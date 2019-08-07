
function opticalRadii,filename = filename,outplot = outplot,verbose = verbose,extno = extno, sr_range = sr_range,B_num = B_num,zero_point = zero_point,center = center,halo_str = halo_str
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
speed_o_light = 299792458 ;m per sec
molec_weight = (0.76*1 + 0.24*4.0)
arcs2_per_ster = 4.25451703e10 
Fsun = 3.196e-10 ; erg/s/cm/cm Flux from Sun @ 10pc
Msun = 4.76 ; absolute mag of sun

IF NOT keyword_set(center) THEN center = [0,0] ELSE center = center
IF NOT keyword_set(halo_str) THEN halo_str = '1'

IF keyword_set(verbose) THEN BEGIN
    loadct,0
    !p.multi = 0
    !Y.STYLE = 1
    !X.STYLE = 1
    !P.THICK = 3.5
    nbins=100.0
    linestyles = [0,2]
    formatplot,outplot = outplot
ENDIF

;****************** Sunrise Cube **************************\
IF keyword_SET(filename) THEN sfilename = filename+'.' + halo_str + '/broadband.fits' ELSE sfilename = 'broadband.fits'
;mfilename = 'mcrx.fits'

temp = mrdfits(sfilename,2,header)
;fovstr = strsplit(header[13],' ',/extract) ;FOV in kpc
fovstr = strsplit(header[14],' ',/extract) ;FOV in rad
;fov = float(fovstr[3]) 
fov = float(fovstr[2])          
;diststr = strsplit(header[12],' ',/extract)
diststr = strsplit(header[17],' ',/extract)
distance = float(diststr[3])    ;Camera distance in kpc
fov = fov*distance

IF NOT keyword_set(extno) THEN extno = 14 ;Faceon
IF NOT keyword_set(B_num) THEN B_num = 21
sr_surface_den = mrdfits(sfilename,extno,sr_head)
dxy = sr_head[where(strcmp('CD1_1',sr_head,5))]
dxy = double((strsplit(dxy,' ',/extract))[2])
sr_naxis = (size(sr_surface_den))[1] ;500 ;480.0
IF NOT keyword_set(sr_range) THEN sr_range = dxy*sr_naxis ;headerHI.naxis1*headerHI.CDELT1

filter = 14 ;13
filters= mrdfits(sfilename,filter)
binsize = 0.05 ;kpc
angular_width = tan(binsize/distance) ;angular width of bins in radians
angular_area = angular_width*angular_width*arcs2_per_ster ;angular area in sterradians
xaxes_sr = (findgen(sr_naxis) - sr_naxis/2.0 - center[0])*sr_range/sr_naxis + sr_range/sr_naxis/2.0
yaxes_sr = (findgen(sr_naxis) - sr_naxis/2.0 - center[1])*sr_range/sr_naxis + sr_range/sr_naxis/2.0
print,minmax(xaxes_sr)
print,minmax(yaxes_sr)
radius = fltarr(n_elements(xaxes_sr),n_elements(yaxes_sr))
FOR ix = 0,n_elements(xaxes_sr) - 1 DO $
  FOR iy = 0,n_elements(yaxes_sr) - 1 DO $
  radius[ix,iy] = sqrt(xaxes_sr[ix]*xaxes_sr[ix] + yaxes_sr[iy]*yaxes_sr[iy])

;ergs_watt = 1e7 ;Watts/erg
;m_A = 1e-10 ;m/Angstrom
;m_cm = 1e-2 ;m/cm
;IF NOT keyword_set(zero_point) THEN zero_point = zeropoint_flux(filters[B_num].filter)
;B_surface_den = sr_surface_den[*,*,B_num]*ergs_watt*m_A*(m_cm)^2/arcs2_per_ster/zero_point ;Converting W/m/m^2/sr to ergs/s/A/cm^2/arcsec^2
;m_B_surface_den = -2.5*alog10(B_surface_den/zero_point) ;Convert flux into a magnitude
Lunit = filters[B_num].ewidth_lambda/filters[B_num].ewidth_nu ;for conversion between W/m^2/m and W/m^2/Hz
Fo=3.631e-23 ;zero point of AB magnitude system for every filter, in W/m^2/Hz
zero_point = Fo/Lunit; W/m/m^2
B_surface_den = sr_surface_den[*,*,B_num]/arcs2_per_ster ;  W/m/m^2/arcsec^2
m_B_surface_den = -2.5*alog10(B_surface_den/zero_point) ;Convert flux into a magnitude
m_B_surface_den[where(~FINITE(m_B_surface_den) )] = 0

;Set up circles of increasing radii 
dr = float(sr_range)/sr_naxis*4.0/2.0 ;0.25  ;difference in radii
radius_arr = findgen(sr_range/dr/2.0 + 1)*dr  ;array of radii
B_surface_den_annuli = fltarr(sr_range/dr/2.0) ;Set up an array to hold the surface brightnesses
area_annuli = fltarr(sr_range/dr/2.0) ;Array to hold the areas
FOR ir = 0, n_elements(radius_arr) - 2 DO BEGIN
  ind = where(radius ge radius_arr[ir] AND radius lt radius_arr[ir + 1])  ;select pixels at the appropriate distance
  area_annuli[ir] = n_elements(ind)  ;area of annuli
  B_surface_den_annuli[ir] = TOTAL(B_surface_den(ind))/area_annuli[ir] ;divide by area
ENDFOR 
m_B_surface_den_annuli = -2.5*alog10(B_surface_den_annuli/zero_point) ;convert to magnitude
IF (where(~FINITE(m_B_surface_den_annuli)))[0] ne -1 THEN m_B_surface_den_annuli[where(~FINITE(m_B_surface_den_annuli))] = 0

mag = 25 ;6
IF keyword_set(verbose) THEN BEGIN
    window,1 ;plot average surface brightness against radii
    plot,radius_arr,m_B_surface_den_annuli,xtitle = 'Radius [kpc]',ytitle = 'Surface Brightness',yrange = [30,22],xrange = [0,20],title = filename;,yrange = [30,MIN(m_B_surface_den_annuli)]
    oplot,[0,20],[mag,mag],linestyle = 2
ENDIF

iords = where(m_B_surface_den_annuli gt mag + 1)
;temp = MIN(ABS(m_B_surface_den_annuli - 26),iord) ;Find the annuli in which the surface brightness is closest to 26 and limit everything to there
m_B_surface_den_annuli = m_B_surface_den_annuli[0:iords[0]]
radius_arr = radius_arr[0:iords[0]]

temp = MIN(ABS(m_B_surface_den_annuli - mag),iord) ;Find the annuli in which the surface brightness is closest to 25
;select region surrounding that radii
IF (n_elements(m_B_surface_den_annuli) gt iord + 3 AND iord - 3 ge 0) THEN BEGIN
    m_B_surface_den_annuli_spline = m_B_surface_den_annuli[iord - 3:iord + 3] 
    radius_arr_spline = radius_arr[iord - 3:iord + 3] ;select region of radius array
ENDIF ELSE BEGIN
    if (n_elements(m_B_surface_den_annuli) le iord + 3) THEN BEGIN
        m_B_surface_den_annuli_spline = m_B_surface_den_annuli[iord - 3:n_elements(m_B_surface_den_annuli)-1] 
        radius_arr_spline = radius_arr[iord - 3:n_elements(m_B_surface_den_annuli)-1]
    ENDIF ELSE BEGIN
        m_B_surface_den_annuli_spline = m_B_surface_den_annuli[0:iord + 3] 
        radius_arr_spline = radius_arr[0:iord + 3] 
    ENDELSE
ENDELSE
order = sort(m_B_surface_den_annuli_spline) ;sort by brightness
r25 = SPLINE(m_B_surface_den_annuli_spline[order],radius_arr_spline[order],mag) ;find r25
temp = radius_arr_spline[order]
if r25 gt temp[n_elements(temp) - 1] then r25 = radius_arr[iord]

IF keyword_set(verbose) THEN BEGIN
    oplot,[r25,r25],[30,20],linestyle = 2

    window,0,xsize = 800,ysize = 800
;    print,minmax(B_surface_Den[where(B_surface_den ne 0)])
;    print,minmax(m_B_surface_den[where(m_B_surface_den ne 0)])
    range = minmax(m_B_surface_den[where(m_B_surface_den ne 0)])
    max = CEIL(range[1])
    min = FLOOR(range[0])
    contour,m_B_surface_den,xaxes_sr,yaxes_sr,min = min,max = max,xrange = [-7,7],yrange = [-7,7],/fill,nlevels = 245
    contour,radius,xaxes_sr,yaxes_sr,/overplot,levels = [1,2,3,4,5,6,7,8,9,10,11,12],color = 100
    contour,radius,xaxes_sr,yaxes_sr,/overplot,levels = [0,r25]    
;    stop
ENDIF

;aux_data = mrdfits(mfilename,27,header)
;lb_slice_den = aux_data[*,*,6]*1.05026322e-36  ;6,7  1 (watt / kpc) / kpc = 1.05026322 × 10-36 ((ergs / s) / cm) / cm
;lb_slice_den = (Msun - 2.5*alog10(lb_slice_den/angular_area/Fsun) + 5*alog10(distance) - 5)
;lb_slice_den[where( ~FINITE(lb_slice_den) )] = 0
;print,minmax(lb_slice_den[where(lb_slice_den ne 0)])
;range = minmax(lb_slice_den[where(lb_slice_den ne 0)])
;min = CEIL(range[1])
;max = FLOOR(range[0])
;lb_slice_den = -1.0*lb_slice_den
;contour,lb_slice_den,xaxes_sr,yaxes_sr,min = -1.0*min,max = -1.0*max,xrange = [-10,10],yrange = [-10,10],/fill,nlevels = 245
;contour,lb_slice_den,xaxes_sr,yaxes_sr,levels = [25],/overplot

;window,1
;plot,

;fuv_surface_den = sr_surface_den[*,*,fuv_num]*fuv_lambda*fuv_lambda/speed_o_light/1e-20 
;mic24_surface_den = sr_surface_den[*,*,mic24_num]*mic24_lambda*mic24_lambda/speed_o_light/1e-20
;SFR_surface_den = 3.2e-3*mic24_surface_den + 8.1e-2*fuv_surface_den ;EQ 1, Bigiel et al, 2008

;stop
;IF keyword_SET(filename) THEN cd,'..'
return,r25
END

pro opticalRadiiMaster 
;filename = 'h603.cosmo50cmb.2304g5bwK.00512'
;cd,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG'

res = 4.0

cd,'/astro/net/nbody1/abrooks/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.00512/'
filename = 'h516.cosmo25cmb.3072g1MBWK.00512'
outplot = '~/opticaltest.eps'
;pfile = '../h516.cosmo25cmb.3072g1MBWK.param'
r25 = opticalRadii(filename,/verbose);,outplot = outplot
stop

cd,'/astro/net/nbody1/abrooks/h799.cosmo25cmb.3072g1MBWK/h799.cosmo25cmb.3072g1MBWK.00512/'
filename = 'h799.cosmo25cmb.3072g1MBWK.00512'
r25 = opticalRadii(filename,/verbose,sr_range = 50);,outplot = outplot
stop

cd,'/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.00512/'
filename = 'h603.cosmo50cmb.3072gs1MbwK.00512'
r25 = opticalRadii(filename,/verbose,sr_range = 50);,outplot = outplot
stop

cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464.dir'
filename = 'h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464'
outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464.dir/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464_' + strtrim(res,2)
r25 = opticalRadii(filename,/verbose);,outplot = outplot

cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir'
filename = 'h516.cosmo25cmb.1536g3HBWK.00512'
pfile = '../../h516.cosmo25cmb.1536g3HBWK.param'
outplot ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512_' + strtrim(res,2)
r25 = opticalRadii(filename,/verbose);,outplot = outplot

cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir'
filename = 'h516.cosmo25cmb.1536g6MbwK.00512'
pfile = '../../h516.cosmo25cmb.1536g6MbwK.param'
outplot ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir/h516.cosmo25cmb.1536g6MbwK.00512_' + strtrim(res,2)
r25 = opticalRadii(filename,/verbose);,outplot = outplot

cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir'
filename = 'h516.cosmo25cmb.1536g3HBWK_noH2SF.00512'
pfile = '../../h516.cosmo25cmb.1536g3HBWK.param'
outplot ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512_' + strtrim(res,2)
r25 = opticalRadii(filename,/verbose);,outplot = outplot

end
