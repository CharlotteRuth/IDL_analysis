
function eighty_flux,filename,extno = extno,band = band,filtername = filtername,verbose = verbose,center = center,halo_str = halo_str
arcs2_per_ster = 4.25451703e10 

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
sfilename = filename+'.' + halo_str + '/broadband.fits'
temp = mrdfits(sfilename,2,header)
fovstr = strsplit(header[13],' ',/extract)
fov = float(fovstr[3]) ;FOV in kpc
diststr = strsplit(header[12],' ',/extract)
distance = float(diststr[3]) ;Camera distance in kpc

IF NOT keyword_set(band) THEN band = 13
IF NOT keyword_set(extno) THEN extno = 14 ;Faceon
IF NOT keyword_set(filtername) THEN filtername = "i_SDSS.res"
filters = mrdfits(sfilename,band,filt_head)
filtind = 5 ;where(strcmp(filters.filter,filtername))
sr_surface_den = mrdfits(sfilename,extno,sr_head)
dxy = sr_head[where(strcmp('CD1_1',sr_head,5))]
dxy = double((strsplit(dxy,' ',/extract))[2])
sr_naxis = (size(sr_surface_den))[1] ;500 ;480.0
IF NOT keyword_set(sr_range) THEN sr_range = dxy*sr_naxis ;headerHI.naxis1*headerHI.CDELT1

;filter = 13
;filters= mrdfits(sfilename,filter)
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
;IF NOT keyword_set(zero_point) THEN zero_point = zeropoint_flux(filters[filtind].filter)
;B_surface_den = sr_surface_den[*,*,filtind]*ergs_watt*m_A*(m_cm)^2/arcs2_per_ster/zero_point ;Converting W/m/m^2/sr to ergs/s/A/cm^2/arcsec^2
;m_B_surface_den = -2.5*alog10(B_surface_den/zero_point) ;Convert flux into a magnitude
Lunit = filters[filtind].ewidth_lambda/filters[filtind].ewidth_nu ;for conversion between W/m^2/m and W/m^2/Hz
Fo=3.631e-23 ;zero point of AB magnitude system for every filter, in W/m^2/Hz
zero_point = Fo/Lunit; W/m/m^2
surface_den = sr_surface_den[*,*,filtind]/arcs2_per_ster ;  W/m/m^2/arcsec^2

;histogramp,radius,weight = surface_den,/cum,/normalize
hist = weighted_histogram(radius,weight = surface_den,/cum,locations = xhist,nbins = 100)
hist = hist/max(hist)
r = spline(hist,xhist,[0.8])

IF keyword_set(verbose) THEN BEGIN
    window,0,xsize = 400,ysize = 400
;    print,minmax(B_surface_Den[where(B_surface_den ne 0)])
;    print,minmax(m_B_surface_den[where(m_B_surface_den ne 0)])
    range = minmax(surface_den[where(surface_den ne 0)])
    max = ceil(range[1])
    min = floor(range[0])
    contour,alog10(surface_den),xaxes_sr,yaxes_sr,min = -15,xrange = [-10,10],yrange = [-10,10],/fill,nlevels = 245
;    contour,radius,xaxes_sr,yaxes_sr,/overplot,levels = [1,2,3,4,5,6,7,8,9,10,11,12],color = 100
    contour,radius,xaxes_sr,yaxes_sr,/overplot,levels = [0,r]

    window,1,xsize = 600,ysize = 400
    plot,xhist,hist
    oplot,[r,r],[0,1],linestyle = 2
    oplot,[0,100],[0.8,0.8],linestyle = 2
ENDIF
RETURN,r25
END
