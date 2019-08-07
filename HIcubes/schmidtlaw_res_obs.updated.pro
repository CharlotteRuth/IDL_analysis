
pro schmidtlaw_res_obs_master 
;filename = 'h603.cosmo50cmb.2304g5bwK.00512'
;cd,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG'
spawn,'hostname',hostname
IF hostname EQ 'vetinari' THEN BEGIN
   prefix = '~/Data/MolecH/Cosmo/' 
   outfile = '~/Figures/h516.cosmo25cmb.paper_'
ENDIF ELSE BEGIN
   IF hostname EQ 'quirm.math.grinnell.edu' THEN BEGIN
      prefix = '/home/christenc/Storage/Cosmo/'
      outfile = '/home/christenc/Figures/h516.cosmo25cmb.3072g1HBWKS.000068'
   ENDIF ELSE BEGIN
      prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
      outfile = '~/plots/'
   ENDELSE
ENDELSE

res = 0.75
;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir'
filename = 'h516.cosmo25cmb.1536g3HBWK.00512'
pfile = '../../h516.cosmo25cmb.1536g3HBWK.param'
outplot ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512_' + strtrim(res,2)
;schmidtlaw_res_obs,filename,pfile,res = res,outplot = outplot,/useH2 
;stop

;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir'
filename = 'h516.cosmo25cmb.1536g6MbwK.00512'
pfile = '../../h516.cosmo25cmb.1536g6MbwK.param'
outplot ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir/h516.cosmo25cmb.1536g6MbwK.00512_' + strtrim(res,2)
;res = 1.0
;schmidtlaw_res_obs,filename,pfile,res = res,outplot = outplot

;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir'
filename = 'h516.cosmo25cmb.1536g3HBWK_noH2SF.00512'
pfile = '../../h516.cosmo25cmb.1536g3HBWK.param'
outplot ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512_' + strtrim(res,2)
;res = 1.0
;schmidtlaw_res_obs,filename,pfile,res = res,/useH2,outplot = outplot

;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464.dir'
filename = 'h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464'
pfile = '../../h516.cosmo25cmb.1536g6HBWK.jeans.prev.param'
outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464.dir/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464_' + strtrim(res,2)
;res = 1.0
;schmidtlaw_res_obs,filename,pfile,res = res,/useH2,outplot = outplot

;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir'
filename = 'h516.cosmo25cmb.1536g1MBWK.00512'
pfile = '../../h516.cosmo25cmb.1536g1MBWK.param'
outplot = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir/h516.cosmo25cmb.1536g1MBWK.00512.halo.1_' + strtrim(res,2)
;schmidtlaw_res_obs,filename,pfile,res = res,outplot = outplot

;cd,prefix+'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir'
filename = 'h516.cosmo25cmb.1536g14HBWK.00512'
pfile = '../../h516.cosmo25cmb.1536g14HBWK.param'
outplot = prefix+'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK.00512.halo.1_' + strtrim(res,2)
;schmidtlaw_res_obs,filename,pfile,res = res,/useH2,outplot = outplot

;cd,prefix+'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir'
filename = 'h516.cosmo25cmb.2304g14HBWK.00512'
pfile = 'h516.cosmo25cmb.2304g14HBWK.param'
outplot = prefix+'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/h516.cosmo25cmb.2304g14HBWK.00512.halo.1_' + strtrim(res,2)
;schmidtlaw_res_obs,filename,pfile,res = res,/useH2,outplot = outplot

;cd,prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir'
filename = 'h516.cosmo25cmb.3072g1MBWK.00492'
pfile = 'h516.cosmo25cmb.3072g1MBWK.param'
outplot = prefix+'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1_' + strtrim(res,2)
;schmidtlaw_res_obs,filename,pfile,res = res,useH2 = 0,/halpha,angle = 45,extno = 15, rotateAngle = 7;,/verbose;outplot = outplot

;cd,prefix+'/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir'
filename = 'h516.cosmo25cmb.3072g14HBWK.00492'
pfile = '../../h516.cosmo25cmb.3072g14HBWK.param'
outplot = prefix+'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir/h516.cosmo25cmb.3072g14HBWK.00492.halo.1_' + strtrim(res,2)
;schmidtlaw_res_obs,filename,pfile,res = res,/useH2,/halpha,angle = 45,extno = 15, rotateAngle = 7;,outplot = outplot

filename = 'h516.cosmo25cmb.3072g1HBWKS.000068'
pfile = 'h516.cosmo25cmb.3072g1HBWKS.param'

outplot = outfile + '/h516.cosmo25cmb.3072g1HBWKS.000068.halo.1_' + strtrim(res,2)
;schmidtlaw_res_obs,filename,pfile,res = res,useH2 = 1,/halpha,angle = 45,extno = 16, rotateAngle = 7,/verbose

cd,prefix + '/h516.cosmo25cmb/h516.cosmo25cmb.2304g1HBWKS'
filename = 'h516.cosmo25cmb.2304g1HBWKS.000512'
pfile = 'h516.cosmo25cmb.2304g1HBWKS.param'
outplot = outfile + '/h516.cosmo25cmb.2304g1HBWKS.000512.halo.1_' + strtrim(res,2)

cd,prefix + '/h516.cosmo25cmb/h516.cosmo25cmb.2304g1HBWK'
filename = 'h516.cosmo25cmb.2304g1HBWK.000512'
pfile = 'h516.cosmo25cmb.2304g1HBWK.param'
outplot = outfile + '/h516.cosmo25cmb.2304g1HBWK.000512.halo.1_' + strtrim(res,2)

cd,prefix + '/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWKS'
filename = 'h516.cosmo25cmb.3072g1HBWKS.000512'
pfile = 'h516.cosmo25cmb.3072g1HBWKS.param'
outplot = outfile + '/h516.cosmo25cmb.3072g1HBWKS.000512.halo.1_' + strtrim(res,2)
schmidtlaw_res_obs,filename,pfile,res = res,useH2 = 1,angle = 45,extno = 16, rotateAngle = 0,/verbose;,outplot = outplot
stop

cd,prefix + '/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWK'
filename = 'h516.cosmo25cmb.3072g1HBWK.000512'
pfile = 'h516.cosmo25cmb.3072g1HBWK.param'
outplot = outfile + '/h516.cosmo25cmb.3072g1HBWK.000512.halo.1_' + strtrim(res,2)

schmidtlaw_res_obs,filename,pfile,res = res,useH2 = 1,angle = 45,extno = 16, rotateAngle = 0,/verbose ;,outplot = outplot              ;,/Halpha

end

pro schmidtlaw_res_obs,filename,pfile,res = res,outplot = outplot,overplot = overplot,useH2 = useH2,verbose = verbose, angle = angle, extno = extno,center = center,rotateAngle = rotateAngle,range = range,halo_str = halo_str, Halpha = Halpha,formatthick = formatthick
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
speed_o_light = 299792458 ;m per sec
molec_weight = (0.76*1 + 0.24*4.0)
units = tipsyunits(pfile)
massunit = units.massunit
lengthunit = units.lengthunit
IF keyword_set(Halpha) THEN use24mm = 1 ELSE use24mm = 1 ;If not based on H-alpha, choice as to whether to use Bigiel (FUV + 24mm) or Madau and Dickinson 2014 (FUV)

IF NOT keyword_set(angle) THEN BEGIN 
    angle_str = '90'
    angle = !PI/2
ENDIF ELSE BEGIN 
    angle_str = STRTRIM(FIX(angle),2)
    angle = !PI*angle/180.0
ENDELSE
IF NOT keyword_set(extno) THEN extno = 14 ELSE extno = extno
IF NOT keyword_set(center) THEN center = [0,0] ELSE center = center
IF NOT keyword_set(rotateAngle) THEN rotateAngle = 0 ;3
IF NOT keyword_set(range) THEN range = 12.0
IF NOT keyword_set(halo_str) THEN halo_str = '1'
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF keyword_set(outplot) THEN BEGIN
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
    !Y.MARGIN = [6,4]*sin(angle)
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
verbose = 1
;*************** Tipsy Data ***********************
if keyword_set(verbose) THEN BEGIN
    rtipsy,filename + '.halo.' + halo_str + '.std',h,g,d,s
    s.x = s.x*units.lengthunit
    s.y = s.y*units.lengthunit
    s.z = s.z*units.lengthunit
    g.x = g.x*units.lengthunit
    g.y = g.y*units.lengthunit
    g.z = g.z*units.lengthunit
    g.mass = g.mass*massunit
    readarr,filename + '.halo.' + halo_str + '.HI',h,HI,/ascii,part = 'gas'
    IF keyword_set( useH2) THEN BEGIN
        readarr,filename + '.halo.' + halo_str + '.H2',h,H2,/ascii,part = 'gas'
    ENDIF
    z = (h.time)-1.
    tcurrent=13.7e9-wmap3_lookback(-1.0*z) ;tcurrent=13.73d9
    deltat = 100.e6
    s.tform = s.tform*units.timeunit
    tcurrent = max(s.tform)     ;1e10
    tclip = tcurrent - deltat

;*************** Rotation ****************
   angle0 = ATAN(s.z/s.x)
   angle0[where(s.x lt 0)] = angle0[where(s.x lt 0)]+!PI
;    angle0 = ATAN(s.z/s.y)
;    angle0[where(s.y lt 0)] = angle0[where(s.y lt 0)]+!PI
    s1 = s
    s1.x = s1.x*h.time
    s1.y = s1.y*h.time
    s1.z = s1.z*h.time    
    if angle ne !PI/2 THEN BEGIN 
       Ryz = SQRT(s.x*s.x + s.z*s.z)
       s1.x = Ryz*COS(angle0 + angle)
;        Ryz = SQRT(s.y*s.y + s.z*s.z)
;        s1.y = Ryz*COS(angle0 + angle)
        s1.z = Ryz*SIN(angle0 + angle)
    ENDIF 
    temp = s1.x
    s1.x = s1.y
    s1.y = temp
    
   angle0 = ATAN(g.z/g.x)
   angle0[where(g.x lt 0)] = angle0[where(g.x lt 0)]+!PI
;    angle0 = ATAN(g.z/g.y)
;    angle0[where(g.y lt 0)] = angle0[where(g.y lt 0)]+!PI    
    g1 = g
    g1.x = g1.x*h.time
    g1.y = g1.y*h.time
    g1.z = g1.z*h.time 
    if angle ne !PI/2 THEN BEGIN 
       Ryz = SQRT(g.x*g.x + g.z*g.z)
       g1.x = Ryz*COS(angle0 + angle)
;        Ryz = SQRT(g.y*g.y + g.z*g.z)
;        g1.y = Ryz*COS(angle0 + angle)
        g1.z = Ryz*SIN(angle0 + angle)
    ENDIF 

    temp = g1.x
    g1.x = g1.y
    g1.y = temp

;        s1 = s
;        s1.y = s1.y*h.time 
;        s1.z = s1.z*h.time                         
;        s1.x = s1.x*h.time
;        g1 = g
;        g1.y =     g1.y*h.time                       
;        g1.z = g1.z*h.time                        
;        g1.x = g1.x*h.time
;        angle = !PI/2.0
ENDIF ELSE rtipsy,filename + '.halo.' + halo_str + '.std',h,g,d,s,/justhead

;***************** Gas Cubes **********************
cubeHI_surface_den = read_dencube_fits(filename + '.halo.' + halo_str + '.' + angle_str + '.smoothed.mom0.fits',headerHI,/noscale);amu/cm^2
;cubeHI_surface_den = read_dencube_fits(filename + '.halo.' + halo_str + '.90.smoothed.mom0.fits',headerHI,/noscale);/amu_per_gm/gm_per_msol*cm_per_kpc*cm_per_kpc/1d6
;cubeHI_surface_den = read_dencube_fits(filename + '.halo.' + halo_str + '.90.cube.mom0.fits',headerHI,/noscale)/2.0;/amu_per_gm/gm_per_msol*cm_per_kpc*cm_per_kpc/1d6
headerHI.CDELT1 = headerHI.CDELT1;*h.time
headerHI.CDELT2 = headerHI.CDELT2;*h.time
IF NOT (where(~finite(cubeHI_surface_den)))[0] EQ -1 THEN cubeHI_surface_den[where(~finite(cubeHI_surface_den))] = 0
IF keyword_set( useH2) THEN BEGIN
    cubeH2 = read_dencube_fits(filename+'.halo.' + halo_str + '.H2.arr.' + angle_str + '.fits',headerH2)*2.0*gm_per_msol*amu_per_gm
    headerH2.CDELT1 = headerH2.CDELT1*h.time
    headerH2.CDELT2 = headerH2.CDELT2*h.time
    pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*cm_per_kpc*cm_per_kpc
;cubeHI_surface_den = TRANSPOSE(cubeHI)/pix_areaH2*gm_per_msol*amu_per_gm
    cubeH2_surface_den = TRANSPOSE(cubeH2)/pix_areaH2;/h.time/h.time ;AMU/cm^2
ENDIF
;cubeGas_den = (read_cube_fits(filename+'.halo.' + halo_str + '.cube.all.fits',head,kpcunit = lengthunit,munit = massunit))/pix_areaH2;*gm_per_msol*amu_per_gm
cubeHI_surface_den_prof = cubeHI_surface_den

;****************** Sunrise Cube **************************
sfilename = filename + '.' + halo_str + '/broadband.fits'
;extno = 13
fuv_num = 0
mic24_num = 15
fuv_lambda = 1.5180502483701E-07 ;m
mic24_lambda = 2.38717045592424E-05 ;m
;mic24_dlambda = 5.23145335670255E-06;m
Halpha_lambda = 6.56d-7 ;m
;fuv_nu = speed_o_light/fuv_lambda ; Hz
;mic24_nu = speed_o_light/mic24_lambda
;extno = 15 ;45 degree angle
;extno = 14 ; Face on looking down z axes with x up and y to right
sr_surface_den = mrdfits(sfilename,extno,sr_head)
dxy = sr_head[where(strcmp('CD1_1',sr_head,5))]
dxy = double((strsplit(dxy,' ',/extract))[2])
sr_NAXIS = (size(sr_surface_Den))[1] ;500 ;480.0
IF NOT keyword_set(sr_range) THEN sr_range = dxy*sr_NAXIS ;headerHI.NAXIS1*headerHI.CDELT1
;sr_fov = 24
;Sunrise units are in W/m/m^2/sr

;Bigiel Units are in MJy*str^(-1) (Jansky = 1e-26 W/Hz/m^2,MegaJansky = 1e6 Jy, MJy = 1e-20 W/Hz/m^2)
mic24_surface_den_unrotated = sr_surface_den[*,*,mic24_num]*mic24_lambda*mic24_lambda/speed_o_light/1e-20
fuv_surface_den_unrotated   = sr_surface_den[*,*,fuv_num]*fuv_lambda*fuv_lambda/speed_o_light/1e-20

;Madau and Dickinson 2014 units are ergs s^-1 Hz^-1
fuv_surface_den_unrotatedW  = sr_surface_den[*,*,fuv_num]*9.5214226d38*1d7*fuv_lambda*fuv_lambda/speed_o_light ;in ergs/s/kpc^2/Hz for

IF keyword_set(Halpha) THEN BEGIN
                                ;Bolatto units are in ergs s^-1 kpc^-2
   IF NOT file_test(filename + '.' + halo_str + '/Halpha.fits') THEN makeLineImage,filename + '.' + halo_str + '/mcrx.fits',extno,6.54e-7,6.583e-7,filename + '.' + halo_str + '/Halpha.fits'
   Halpha_surface_den_unrotatedW = mrdfits(filename + '.' + halo_str + '/Halpha.fits')*9.5214226d38*1d7 ;*Halpha_lambda
;   stop
    mic24_surface_den_unrotatedW  = sr_surface_den[*,*,mic24_num]*9.5214226d38*1d7*mic24_lambda

    SFR_surface_den_unrotated = 5.3d-42*(Halpha_surface_den_unrotatedW + 0.031*mic24_surface_den_unrotatedW) ;Bolatto
    SFR_surface_den_unrotated2 = 3.2e-3*mic24_surface_den_unrotated + 8.1e-2*fuv_surface_den_unrotated ;Bigiel
 ENDIF ELSE BEGIN
    IF (use24mm EQ 1) THEN BEGIN
       SFR_surface_den_unrotated = 3.2e-3*mic24_surface_den_unrotated + 8.1e-2*fuv_surface_den_unrotated ;EQ 1, Bigiel et al, 2008
       SFR_surface_den_unrotated2 = 1.3d-28*FUV_surface_den_unrotatedW ;From Madau and Dickinson 2014 where FUVrebin is the luminosity in units of ergs s^-1 Hz^-1. Eq 10 for solar metallicity
    ENDIF ELSE BEGIN
       SFR_surface_den_unrotated2 = 3.2e-3*mic24_surface_den_unrotated + 8.1e-2*fuv_surface_den_unrotated ;EQ 1, Bigiel et al, 2008
       SFR_surface_den_unrotated = 1.3d-28*FUV_surface_den_unrotatedW ;From Madau and Dickinson 2014 where FUVrebin is the luminosity in units of ergs s^-1 Hz^-1. Eq 10 for solar metallicity       
    ENDELSE
ENDELSE
    
;histogramp,alog10(SFR_surface_den_unrotated),nbins = 100,min = -6,max = -0.5
;histogramp,alog10(SFR_surface_den_unrotated2),nbins = 100, min = -6,max = -0.5, color = 100, /overplot
;histogramp,alog10(5.3d-42*0.031*mic24_surface_den_unrotatedW),nbins = 100,min = -6,max = -0.5
;histogramp,alog10(3.2e-3*mic24_surface_den_unrotated),nbins = 100, min = -6,max = 0.5, color = 100, /overplot
;histogramp,alog10(5.3d-42*Halpha_surface_den_unrotatedW),nbins = 100,min = -6,max = -0.5
;histogramp,alog10(8.1e-2*fuv_surface_den_unrotated),nbins = 100,min = -6,max = 0.5, color = 100, /overplot

;plot,5.3d-42*0.031*mic24_surface_den_unrotatedW,3.2e-3*mic24_surface_den_unrotated,psym=3,/xlog,/ylog,xrange = [1e-14,1e-1],yrange = [1e-14,1e-1]
;plot,5.3d-42*     Halpha_surface_den_unrotatedW,8.1e-2*  fuv_surface_den_unrotated,psym=3,/xlog,/ylog,xrange = [1e-14,1e-1],yrange = [1e-14,1e-1]
;plot,SFR_surface_den_unrotated,SFR_surface_den_unrotated2,psym = 3, /xlog,/ylog,xrange = [1e-8,1e-1],yrange = [1e-8,1e-1]
;oplot,[1e-14,1e-1],[1e-14,1e-1]

;fuv_surface_den = ROTATE(fuv_surface_den_unrotated,7) ;used with 45 degree angle
;mic24_surface_den = ROTATE(mic24_surface_den_unrotated,7);used with 45 degree angle
;SFR_surface_den = ROTATE(SFR_surface_den_unrotated,7);used with 45 degree angle
;fuv_surface_den2 = ROTATE(fuv_surface_den_unrotated,1)
;mic24_surface_den2 = ROTATE(mic24_surface_den_unrotated,1)
;SFR_surface_den2 = ROTATE(SFR_surface_den_unrotated,1)
    IF (use24mm EQ 1) THEN $
       FUV_surface_den = ROTATE(FUV_surface_den_unrotated,rotateAngle) ELSE $ ;For Bigiel
          FUV_surface_den = ROTATE(FUV_surface_den_unrotatedW,rotateAngle) ;For Madau and Dickinson 2014
mic24_surface_den = ROTATE(mic24_surface_den_unrotated,rotateAngle)
SFR_surface_den = ROTATE(SFR_surface_den_unrotated,rotateAngle)
SFR_surface_den2 = ROTATE(SFR_surface_den_unrotated2,rotateAngle)
;center = [center[1],-1.0*center[0]] 

;**************** Rebinning and Spatial Axes

;Axes for the HI/H2 surface density
xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1 + headerHI.CDELT1/2
yaxes = (findgen(headerHI.NAXIS2) - headerHI.NAXIS1/2.0)*headerHI.CDELT2 + headerHI.CDELT1/2
dxcell = headerHI.cdelt1
areapc = 1e6*dxcell*dxcell
areakpc = dxcell*dxcell
minrange = headerHI.CRVAL1
radius = fltarr(headerHI.NAXIS1,headerHI.NAXIS2)
FOR ix = 0,N_ELEMENTS(xaxes) - 1 DO $
  FOR iy = 0,N_ELEMENTS(yaxes) - 1 DO $
  radius[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy]/SIN(angle)/SIN(angle))

;Axes for the sunrise image
;sr_range = 24 ;headerHI.NAXIS1*headerHI.CDELT1
;sr_NAXIS = 480.0
dcell_sr = sr_range/sr_NAXIS
;xaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0)*dcell_sr + dcell_sr/2 
;yaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0)*dcell_sr + dcell_sr/2
xaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0 - center[0])*sr_range/sr_NAXIS + sr_range/sr_NAXIS/2.0
yaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0 - center[1])*sr_range/sr_NAXIS + sr_range/sr_NAXIS/2.0

;spatial_axes = findgen(headerHI.NAXIS1)*headerHI.CDELT1 + headerHI.CRVAL1
minrange = sr_range/(-2.0) 
radius_sr = fltarr(sr_NAXIS,sr_NAXIS)
FOR ix = 0,N_ELEMENTS(xaxes_sr) - 1 DO $
  FOR iy = 0,N_ELEMENTS(yaxes_sr) - 1 DO $
  radius_sr[ix,iy] = SQRT(xaxes_sr[ix]*xaxes_sr[ix] + yaxes_sr[iy]*yaxes_sr[iy]/SIN(angle)/SIN(angle))

;Axes rebinned
IF NOT keyword_set(res) THEN res = 0.75 ;kpc
res_st = strtrim(res,2)
ncell = round(sr_range/res)
res = 1.0*sr_range/ncell
xaxes_lr = (findgen(ncell)*res - sr_range/2.0 + res/2.0)
yaxes_lr = (findgen(ncell)*res - sr_range/2.0 + res/2.0)
HIrebin = fltarr(ncell,ncell)
IF keyword_set( useH2) THEN H2rebin = fltarr(ncell,ncell)
FUVrebin = fltarr(ncell,ncell)
mic24rebin = fltarr(ncell,ncell)
SFRrebin = fltarr(ncell,ncell)
SFRrebin2 = fltarr(ncell,ncell)
SFRrebin_par = fltarr(ncell,ncell)
HI_par = fltarr(ncell,ncell)
H2_par = fltarr(ncell,ncell)
gas_sd = fltarr(ncell,ncell)
radius_rebin = fltarr(ncell,ncell)
FOR ix = 0, ncell - 1 DO BEGIN
    FOR iy = 0,ncell - 1 DO BEGIN 
        radius_rebin[ix,iy] = sqrt((ix*res + 0.5*res + minrange)^2 + (iy*res + 0.5*res + minrange)^2/SIN(angle)/SIN(angle))

        indx = where((xaxes GT (ix*res + minrange)) AND (xaxes LE ((ix+1)*res + minrange)))
        indy = where((yaxes GT (iy*res + minrange)) AND (yaxes LE ((iy+1)*res + minrange)))
        IF ((indx[0] EQ -1) OR (indy[0] EQ -1)) THEN HIrebin[ix,iy] = 0 ELSE HIrebin[ix,iy] = MEAN((cubeHI_surface_den[indx,*])[*,indy])
        IF keyword_set( useH2) THEN BEGIN
            IF ((indx[0] EQ -1) OR (indy[0] EQ -1)) THEN H2rebin[ix,iy] = 0 ELSE H2rebin[ix,iy] = MEAN((cubeH2_surface_den[indx,*])[*,indy])
        ENDIF
;        HIrebin[ix,iy] = MEAN(cubeHI_surface_den[ix*dcell : ix*dcell + dcell - 1, iy*dcell :iy*dcell + dcell - 1])
;        IF keyword_set( useH2) THEN H2rebin[ix,iy] = MEAN(cubeH2_surface_den[ix*dcell : ix*dcell + dcell - 1, iy*dcell :iy*dcell + dcell - 1])

        ind_srx = where((xaxes_sr GT (ix*res + minrange)) AND (xaxes_sr LE ((ix+1)*res + minrange)))
        ind_sry = where((yaxes_sr GT (iy*res + minrange)) AND (yaxes_sr LE ((iy+1)*res + minrange)))
        IF(ind_srx[0] NE -1 AND ind_sry[0] NE -1 ) THEN BEGIN
            FUVrebin[ix,iy] =   MEAN((  fuv_surface_den[ind_srx,*])[*,ind_sry])
            mic24rebin[ix,iy] = MEAN((mic24_surface_den[ind_srx,*])[*,ind_sry])
            SFRrebin[ix,iy] =   MEAN((  SFR_surface_den[ind_srx,*])[*,ind_sry])
            SFRrebin2[ix,iy] =   MEAN((  SFR_surface_den2[ind_srx,*])[*,ind_sry])
        ENDIF ELSE BEGIN
            FUVrebin[ix,iy] =   0
            mic24rebin[ix,iy] = 0
            SFRrebin[ix,iy] =   0
            SFRrebin2[ix,iy] =   0
        ENDELSE
;        FUVrebin[ix,iy] =   MEAN(  fuv_surface_den[ix*dcell_sr : ix*dcell_sr + dcell_sr - 1, iy*dcell_sr :iy*dcell_sr + dcell_sr - 1])
;        mic24rebin[ix,iy] = MEAN(mic24_surface_den[ix*dcell_sr : ix*dcell_sr + dcell_sr - 1, iy*dcell_sr :iy*dcell_sr + dcell_sr - 1])
;        SFRrebin[ix,iy] =   MEAN(  SFR_surface_den[ix*dcell_sr : ix*dcell_sr + dcell_sr - 1, iy*dcell_sr :iy*dcell_sr + dcell_sr - 1])

        IF keyword_set(verbose) THEN BEGIN
            inds = where((s1.x gt (res*ix + minrange)) AND (s1.x le (res*(ix+1) + minrange)) AND (s1.y gt (res*iy + minrange)) AND (s1.y lt (res*(iy+1) + minrange)) AND s1.tform gt tclip)
            IF inds[0] NE -1 THEN SFRrebin_par[ix,iy] = TOTAL(s1[inds].mass*massunit)/res/res/deltat $
            ELSE SFRrebin_par[ix,iy] = 0
            indg = where((g1.x gt (res*ix + minrange)) AND (g1.x le (res*(ix+1) + minrange)) AND (g1.y gt (res*iy + minrange)) AND (g1.y lt (res*(iy+1) + minrange)))
            IF indg[0] NE -1 THEN HI_par[ix,iy] = TOTAL(g1[indg].mass*HI[indg])/res/res/1000.0/1000.0 ELSE HI_par[ix,iy] = 0
            IF keyword_set( useH2) THEN BEGIN
                IF indg[0] NE -1 THEN H2_par[ix,iy] = TOTAL(g1[indg].mass*2.0*H2[indg])/res/res/1000.0/1000.0 ELSE H2_par[ix,iy] = 0
                IF indg[0] NE -1 THEN gas_sd[ix,iy] = TOTAL(g1[indg].mass*(HI[indg] + 2.0*H2[indg]))/res/res/1000.0/1000.0 ELSE gas_sd[ix,iy] = 0
            ENDIF ELSE BEGIN
                H2_par[ix,iy] = 0
                IF indg[0] NE -1 THEN gas_sd[ix,iy] = TOTAL(g1[indg].mass*(HI[indg]))/res/res/1000.0/1000.0 ELSE gas_sd[ix,iy] = 0
            ENDELSE
        ENDIF ELSE BEGIN
            SFRrebin_par[ix,iy] = 0
            HI_par[ix,iy] = 0
            H2_par[ix,iy] = 0
            gas_sd[ix,iy] = 0
        ENDELSE
     ENDFOR
 ENDFOR
IF (use24mm EQ 1) THEN $
   SFRrebin_ave = 3.2e-3*mic24rebin + 8.1e-2*FUVrebin ELSE $;Bigiel
      SFRrebin_ave = 1.3d-28*FUVrebin                       ;Madau and Dickinson 2014

;*************** Plotting **************
;IF keyword_set( useH2) THEN cubeH2_surface_den[where(cubeH2_surface_den lt 3e12)] = 3e12
;cubeHI_surface_den[where(cubeHI_surface_den lt 3e19)] = 3e19

IF (keyword_set(verbose)) THEN BEGIN
;    !p.multi = [0,2,2]
    IF (keyword_set(outplot)) THEN device,filename=outplot+'_part.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 800,ysize = 800
    loadct,39
    multiplot,[2,2],gap=0.06,/square,/doxaxis,/doyaxis
    contour,alog10(SFRrebin),xaxes_lr,yaxes_lr,xrange = [-8,8],yrange = [-8,8],nlevels = 60,/fill,xstyle = 1,ystyle = 1,min_value = -8,max_value = 0,xtitle = 'X [kpc]',ytitle = 'Y [kpc]';,title = textoidl('Log \Sigma')
    contour,alog10(SFRrebin_par),xaxes_lr,yaxes_lr,nlevels = 10,/overplot,min_value = -3,max_value = 0

    multiplot,/doxaxis,/doyaxis    
    min = 8                     ;alog10(0.32)
    color = 270 - (alog10(s.tform - 10.0^min))/(alog10(MAX(s.tform - 10^min)))*254
    color = 270 - 254*(alog10(tcurrent + 1e6 - s.tform) - 6)/(alog10(tcurrent) - 6)
    IF (where(color lt 0))[0] ne -1 THEN color[where(color lt 0)] = 0
    plot,s1.x,s1.y,psym = 3,xrange = [-8,8],yrange = [-8,8],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',title = 'Stellar Age'
    FOR i =0L, N_ELEMENTS(s) - 1 DO $
      oplot,[s1[i].x,s1[i].x],[s1[i].y,s1[i].y],psym = 3,color = color[i]
;contour,alog10(SFRrebin_par),xaxes_lr,yaxes_lr,nlevels = 20,/overplot,min_value = -8,max_value = -2
;contour,alog10(SFR_surface_den),xaxes_sr,yaxes_sr,nlevels = 5,/overplot

    multiplot,/doxaxis,/doyaxis    
    contour,alog10(gas_sd),xaxes_lr,yaxes_lr,xrange = [-8,8],yrange = [-8,8],nlevels = 60,/fill,xstyle = 1,ystyle = 1,min_value = -2,max_value = 2, xtitle = 'X [kpc]',ytitle = 'Y [kpc]';title=textoidl('Log \Sigma'),
;contour,alog10(HIrebin + H2rebin),xaxes_lr,yaxes_lr,xrange = [-8,8],yrange = [-8,8],nlevels = 60,/fill,xstyle = 1,ystyle = 1,min_value = 12
    IF keyword_set( useH2) THEN contour,alog10((HIrebin + H2rebin)/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18),xaxes_lr,yaxes_lr,nlevels = 10,/overplot,min_value = -2,max_value = 2 $
    ELSE contour,alog10((HIrebin)/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18),xaxes_lr,yaxes_lr,nlevels = 5,/overplot,min_value = -2,max_value = 2

    multiplot,/doxaxis,/doyaxis    
    min = 2
    IF keyword_set( useH2) THEN color = (alog10(HI*g.mass + H2*g.mass) - min)/(MAX(alog10(HI*g.mass + H2*g.mass)) - min)*254 ELSE color = (alog10(HI*g.mass) - min)/(MAX(alog10(HI*g.mass)) - min)*254
;color = (alog10(HI*g.mass) - min)/(MAX(alog10(HI*g.mass)) - min)*254
    color[where(color lt 0)] = 0
    plot,g1[where(color gt 0)].x,g1[where(color gt 0)].y,psym = 3,xrange = [-8,8],yrange = [-8,8],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',title = 'Hydrogen Gas Mass'
    FOR i =0L, N_ELEMENTS(g) - 1 DO $
      oplot,[g1[i].x,g1[i].x],[g1[i].y,g1[i].y],psym = 3,color = color[i]
    IF keyword_set( useH2) THEN contour,alog10(cubeHI_surface_den + cubeH2_surface_den),xaxes,yaxes,nlevels = 6,/overplot ELSE  contour,alog10(cubeHI_surface_den),xaxes,yaxes,nlevels = 6,/overplot
    IF keyword_set(outplot) THEN device,/close ELSE stop
;    !p.multi = 0
    multiplot,/reset
endif

loadct,0
nlevels = 60
;!p.multi = 0
;position = [0.42, 0.6, 0.45, 0.9]
position =  [0.85, 0.15, 0.9, 0.9]
IF keyword_set( useH2) THEN BEGIN
    IF (keyword_set(outplot)) THEN device,filename=outplot+'_H2sd.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,1,xsize = 800,ysize = 800
    loadct,0
    mingas = 1e-2
    maxgas = 1e2
    userlevels = indgen(nlevels)*(alog10(maxgas) - alog10(mingas))/nlevels + alog10(mingas)
    cubeH2_surface_den_plot = cubeH2_surface_den/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18
    ind = where(cubeH2_surface_den_plot lt mingas)
    IF ind[0] NE -1 THEN cubeH2_surface_den_plot[ind] = mingas*1.01
;;cubeH2_surface_den_plot[where(NOT FINITE(cubeH2_surface_den_plot))] = 1e-2
    contour,alog10(cubeH2_surface_den_plot),xaxes,yaxes,$
      /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),levels = userlevels ;,title = textoidl('H_2 Surface Density [M')+sunsymbol()+textoidl(' pc^{-2}]')
    colorbar,range = [alog10(mingas),alog10(maxgas)],/vertical,divisions = 4,position = position,color = white
    IF keyword_set(outplot) THEN device,/close ELSE stop
endif

IF (keyword_set(outplot)) THEN device,filename=outplot+'_HIsd.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE BEGIN
    window,1,xsize = 800,ysize = 800
ENDELSE
mingas = 1e-2
maxgas = 1e2
userlevels = indgen(nlevels)*(alog10(maxgas) - alog10(mingas))/nlevels + alog10(mingas)
cubeHI_surface_den_plot = cubeHI_surface_den/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18
ind = where(cubeHI_surface_den_plot lt mingas)
IF ind[0] NE -1 THEN cubeHI_surface_den_plot[ind] = mingas*1.01
contour,alog10(cubeHI_surface_den_plot),xaxes,yaxes,$
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),levels = userlevels ;,title = textoidl('HI Surface Density [M')+sunsymbol()+textoidl(' pc^{-2}]')
colorbar,range = [alog10(mingas),alog10(maxgas)],/vertical,divisions = 4,position = position,color = white
IF keyword_set(outplot) THEN device,/close

IF (keyword_set(outplot)) THEN device,filename=outplot+'_SFsd.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE BEGIN
    stop
    window,1,xsize = 800,ysize = 800
ENDELSE
minsf = 1e-7
maxsf = 1e1
userlevels = indgen(nlevels)*(alog10(maxsf) - alog10(minsf))/nlevels + alog10(minsf)
SFR_surface_den_plot = SFR_surface_den
ind = where(SFR_surface_den lt minsf)
IF ind[0] NE -1 THEN SFR_surface_den_plot[ind] = minsf*1.01
contour,alog10(SFR_surface_den_plot),xaxes_sr,yaxes_sr,$
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxsf),min_value = alog10(minsf),levels = userlevels ;,title=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]"
colorbar,range = [alog10(minsf),alog10(maxsf)],/vertical,divisions = 8,position = position,color = white;,position = position
IF keyword_set(outplot) THEN device,/close ELSE stop

r25 = opticalRadii(filename = filename,verbose = verbose,halo_str = halo_str,extno = 15)
;stop
r25 = r25*2 ;7 ;opticalRadii(filename = filename,verbose = verbose,halo_str = halo_str,extno = 15)

mingas = 1e18
maxgas = 1e23
cubeHI_surface_den_cut = cubeHI_surface_den
ind = where(cubeHI_surface_den_cut le mingas)
IF ind[0] NE -1 THEN cubeHI_surface_den_cut[ind] = 1.01*mingas
HIrebin_cut = HIrebin
ind = where(HIrebin_cut le mingas)
IF ind[0] NE -1 THEN HIrebin_cut[ind] = 1.01*mingas
userlevels = indgen(nlevels)*(alog10(maxgas) - alog10(mingas))/nlevels + alog10(mingas)
IF keyword_set( useH2) THEN BEGIN
    cubeH2_surface_den_cut = cubeH2_surface_den
    ind = where(cubeH2_surface_den_cut le mingas) 
    IF ind[0] NE -1 THEN cubeH2_surface_den_cut[ind] = 1.01*mingas
    H2rebin_cut = H2rebin
    ind = where(H2rebin_cut le mingas)
    IF ind[0] NE -1 THEN H2rebin_cut[ind] = 1.01*mingas
    position = [0.415, 0.59, 0.44, 0.90]
    !X.MARGIN = [4,1]
    !Y.MARGIN = [2,1]
    IF (keyword_set(outplot)) THEN device,filename=outplot+'_gas.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,1,xsize = 800,ysize = 800
    multiplot,[2,2],gap=0.06,/square,/doxaxis,/doyaxis

    contour,alog10(cubeHI_surface_den_cut + cubeH2_surface_den_cut),xaxes,yaxes,title = textoidl('Log \Sigma_{HI + H_2} [amu cm^{-2}]'), $
      /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels 
    contour,radius,xaxes,yaxes,/overplot,levels = [0,r25]
    colorbar,range = [alog10(mingas),alog10(maxgas)],/vertical,divisions = 4,position = position,charsize=cb_charsize,color = white
    multiplot,/doxaxis,/doyaxis
    contour,alog10(HIrebin_cut + H2rebin_cut),xaxes_lr,yaxes_lr,title = textoidl('Log \Sigma_{HI + H_2} [amu cm^{-2}]'),$                
      /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels 
    contour,radius_rebin,xaxes_lr,yaxes_lr,/overplot,levels = [0,r25]
;colorbar,range = [19,23],/vertical,divisions = 3
    multiplot,/doxaxis,/doyaxis
    contour,alog10(cubeH2_surface_den_cut),xaxes,yaxes,title = textoidl('Log \Sigma_{H_2} [amu cm^{-2}]'), $                        
      /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels 
;colorbar,range = [12,23],/vertical,divisions = 3
    multiplot,/doxaxis,/doyaxis
    contour,alog10(cubeHI_surface_den_cut),xaxes,yaxes,title = textoidl('Log \Sigma_{HI} [amu cm^{-2}]'),$                         
      /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels 
;colorbar,range = [19,23],/vertical,divisions = 3
    IF keyword_set(outplot) THEN device,/close
    multiplot,/reset
ENDIF ELSE BEGIN
;    !p.multi = [0,2,1]
;if keyword_set(outplot) THEN position = [0.88, 0.14, 0.95, 0.92] ELSE
    position = [0.415, 0.17, 0.44, 0.83]
    !X.MARGIN = [4,1]
    !Y.MARGIN = [2,1]
    IF (keyword_set(outplot)) THEN device,filename=outplot+'_gas.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 9,xoffset =  2,yoffset =  2  ELSE window,1,xsize = 800,ysize = 400
    multiplot,[2,1],gap=0.06,/square,/doxaxis,/doyaxis
    contour,alog10(cubeHI_surface_den_cut),xaxes,yaxes,$
      /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels ;,title = textoidl('Log \Sigma_{H} [amu cm^{-2}]')
    contour,radius,xaxes,yaxes,/overplot,levels = [0,r25]    
    colorbar,range = [alog10(mingas),alog10(maxgas)],/vertical,divisions = 4,position = position,charsize=cb_charsize,color = white
    multiplot,/doxaxis,/doyaxis
    contour,alog10(HIrebin_cut),xaxes_lr,yaxes_lr,$
      /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxgas),min_value = alog10(mingas),xcharsize = 1.0,ycharsize = 1.0, levels = userlevels;title = textoidl('Log \Sigma_{H}[amu cm^{-2}]')
    contour,radius_rebin,xaxes_lr,yaxes_lr,/overplot,levels = [0,r25] 
    IF keyword_set(outplot) THEN device,/close
    multiplot,/reset
ENDELSE

minsf = 1e-7
maxsf = 1e1
userlevels = indgen(nlevels)*(alog10(maxsf) - alog10(minsf))/nlevels + alog10(minsf)
SFR_surface_den_plot = SFR_surface_den
ind = where(SFR_surface_den lt minsf)
IF ind[0] NE -1 THEN SFR_surface_den_plot[ind] = minsf*1.01
SFRrebin_plot = SFRrebin
ind = where(SFRrebin_plot lt minsf)
IF ind[0] NE -1 THEN SFRrebin_plot[ind] = minsf*1.01
fuv_surface_den_plot = fuv_surface_den
ind = where(fuv_surface_den lt minsf/8.1e-2)
IF ind[0] NE -1 THEN fuv_surface_den_plot[ind] = minsf*1.01/8.1e-2
mic24_surface_den_plot = mic24_surface_den
ind = where(mic24_surface_den lt minsf/3.2e-3)
IF ind[0] NE -1 THEN mic24_surface_den_plot[ind] = minsf*1.01/3.2e-3

position = [0.415, 0.59, 0.44, 0.90]
!X.MARGIN = [4,1]
!Y.MARGIN = [2,1]
IF (keyword_set(outplot)) THEN device,filename=outplot+'_star.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 18,xoffset =  2,yoffset =  2  ELSE window,3,xsize = 800,ysize = 800
multiplot,[2,2],gap=0.06,/square,/doxaxis,/doyaxis
contour,alog10(SFR_surface_den_plot),xaxes_sr,yaxes_sr,$
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxsf),min_value = alog10(minsf),xcharsize = 1.0,ycharsize = 1.0,levels = userlevels;,title=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]",
contour,radius_sr,xaxes_sr,yaxes_sr,/overplot,levels = [0,r25]
colorbar,range = [alog10(minsf),alog10(maxsf)],/vertical,divisions = 8,position = position,charsize=cb_charsize,color = white
multiplot,/doxaxis,/doyaxis
contour,alog10(SFRrebin_plot),xaxes_lr,yaxes_lr,$
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxsf),min_value = alog10(minsf),xcharsize = 1.0,ycharsize = 1.0,levels = userlevels;,title=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]",
contour,radius_rebin,xaxes_lr,yaxes_lr,/overplot,levels = [0,r25]
multiplot,/doxaxis,/doyaxis
contour,alog10(fuv_surface_den_plot*8.1e-2),xaxes_sr,yaxes_sr,title = 'FUV',$
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxsf),min_value = alog10(minsf),xcharsize = 1.0,ycharsize = 1.0,levels = userlevels
multiplot,/doxaxis,/doyaxis
contour,alog10(mic24_surface_den_plot*3.2e-3),xaxes_sr,yaxes_sr,title = '24 microns',$
  /fill,nlevels = nlevels,yrange = [-1.0*range,range],xrange = [-1.0*range,range],xstyle = 1,ystyle = 1,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',max_value = alog10(maxsf),min_value = alog10(minsf),xcharsize = 1.0,ycharsize = 1.0,levels = userlevels
IF keyword_set(outplot) THEN device,/close
multiplot,/reset

!X.MARGIN = [12,5]
!Y.MARGIN = [6,4]
!p.multi = 0
loadct,39
IF keyword_set(outplot) THEN  device,filename = outplot + 'resSchmidt.eps' ELSE window,2,xsize = 400,ysize = 400

xsigmalow = findgen(300)/100 - 1.
xsigma = 10.0^(findgen(400)/100 + 1.) ;xsigmalow
ysigma=2.5e-4*xsigma^1.4
ysigma1=2.5e-4*xsigma^(1.4 + 0.15)
ysigma2=2.5e-4*xsigma^(1.4 - 0.15)
ysigmalow = xsigmalow*2.4 - 5.0

HIrebin = HIrebin/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18
IF keyword_set( useH2) THEN  H2rebin = H2rebin/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18
inside = where(radius_rebin lt r25)
IF (inside[0]) ne -1 THEN BEGIN
   IF keyword_set( useH2) THEN plot,alog10(sin(angle)*(HIrebin[inside] + H2rebin[inside])*1.4),alog10(sin(angle)*SFRrebin_par[inside]),psym = 2, xstyle=1, ystyle=1,xrange = [-1,5], yrange = [-5,3] $ ;,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]",xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"
   ELSE plot,alog10(sin(angle)*(HIrebin[inside])),                  alog10(sin(angle)*SFRrebin[inside]),psym = 2, xstyle=1, ystyle=1,xrange = [-1,5], yrange = [-5,3];, xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]" ,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]"

;IF keyword_set( useH2) THEN plot,alog10(sin(angle)*(gas_sd[inside])),alog10(sin(angle)*SFRrebin_par[inside]),psym = 2,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-1,5], yrange = [-5,3], xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]" ELSE $
;                            plot,alog10(sin(angle)*(gas_sd[inside])),alog10(sin(angle)*SFRrebin_par[inside]),psym = 2,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-1,5], yrange = [-5,3], xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]" 

oplot,alog10(xsigma),alog10(ysigma)
oplot,xsigmalow,ysigmalow 
oplot,[1,1],[-5,3],linestyle = 1

finite = where(FINITE(FUVrebin),complement = nfinite)
if ((nfinite[0]) ne -1) then FUVrebin[nfinite] = 0
finite = where(FINITE(mic24rebin),complement = nfinite)
if ((nfinite[0]) ne -1) then mic24rebin[nfinite] = 0
finite = where(FINITE(SFRrebin_par),complement = nfinite)
if ((nfinite[0]) ne -1) then SFRrebin_par[nfinite] = 0
finite = where(FINITE(SFRrebin),complement = nfinite)
if ((nfinite[0]) ne -1) then SFRrebin[nfinite] = 0
IF keyword_set( useH2) THEN BEGIN
    finite = where(FINITE(H2rebin),complement = nfinite)
    if ((nfinite[0]) ne -1) then H2rebin[nfinite] = 0
ENDIF
finite = where(FINITE(HIrebin),complement = nfinite)
if ((nfinite[0]) ne -1) then HIrebin[nfinite] = 0
finite = where(FINITE(HI_par),complement = nfinite)
if ((nfinite[0]) ne -1) then HI_par[nfinite] = 0
finite = where(FINITE(H2_par),complement = nfinite)
if ((nfinite[0]) ne -1) then H2_par[nfinite] = 0
finite = where(FINITE(gas_sd),complement = nfinite)
if ((nfinite[0]) ne -1) then gas_sd[nfinite] = 0

;IF keyword_set( useH2) THEN writecol,'schmidtlaw_res_obs'+strtrim(res,2)+'.dat', $
;reform(FUVrebin,N_ELEMENTS(FUVrebin))*sin(angle), reform(mic24rebin,N_ELEMENTS(FUVrebin))*sin(angle), reform(SFRrebin_par,N_ELEMENTS(FUVrebin))*sin(angle),reform( SFRrebin,N_ELEMENTS(FUVrebin))*sin(angle), $
;reform(HIrebin,N_ELEMENTS(FUVrebin))*sin(angle), reform(H2rebin,N_ELEMENTS(FUVrebin))*sin(angle), reform(gas_sd,N_ELEMENTS(FUVrebin))*sin(angle), reform(HIrebin,N_ELEMENTS(FUVrebin))*sin(angle)+reform(H2rebin,N_ELEMENTS(FUVrebin))*sin(angle),FORMAT = '(8F)' ELSE writecol,'schmidtlaw_res_obs'+strtrim(res,2)+'.dat', $
;reform(FUVrebin,N_ELEMENTS(FUVrebin))*sin(angle), reform(mic24rebin,N_ELEMENTS(FUVrebin))*sin(angle), reform(SFRrebin_par,N_ELEMENTS(FUVrebin))*sin(angle), reform(SFRrebin,N_ELEMENTS(FUVrebin))*sin(angle), $
;reform(HIrebin,N_ELEMENTS(FUVrebin))*sin(angle), fltarr(N_ELEMENTS(HIrebin)), reform(gas_sd,N_ELEMENTS(FUVrebin))*sin(angle), reform(HIrebin,N_ELEMENTS(FUVrebin))*sin(angle),FORMAT = '(8F)'

;FUV        24micron        SFR (particle)        SFR (obs)        HI
;(particle)        HI (obs)        H2 (particle)        H2 (obs)        Gas
;(particle)        Gas (obs)
;stop
IF keyword_set(Halpha) THEN outname = 'schmidtlaw_res_obs_all_Ha'+res_st+'_2.dat' ELSE outname = 'schmidtlaw_res_obs_all'+res_st+'_2.dat'
IF keyword_set( useH2) THEN writecol,outname, $
reform(radius_rebin[inside],N_ELEMENTS(inside)), $
reform(FUVrebin[inside],N_ELEMENTS(inside))*sin(angle),reform(mic24rebin[inside],N_ELEMENTS(inside))*sin(angle), reform(SFRrebin_par[inside],N_ELEMENTS(inside))*sin(angle),reform(SFRrebin[inside],N_ELEMENTS(inside))*sin(angle), $
reform(  HI_par[inside],N_ELEMENTS(inside))*sin(angle),reform(   HIrebin[inside],N_ELEMENTS(inside))*sin(angle), reform(      H2_par[inside],N_ELEMENTS(inside))*sin(angle),reform( H2rebin[inside],N_ELEMENTS(inside))*sin(angle), $
                                     reform(  gas_sd[inside],N_ELEMENTS(inside))*sin(angle),reform(   HIrebin[inside],N_ELEMENTS(inside))*sin(angle)+reform(H2rebin[inside],N_ELEMENTS(inside))*sin(angle),FORMAT = '(11E)' ELSE writecol,outname, $
reform(radius_rebin[inside],N_ELEMENTS(inside)), $   
reform(FUVrebin[inside],N_ELEMENTS(inside))*sin(angle),reform(mic24rebin[inside],N_ELEMENTS(inside))*sin(angle), reform(SFRrebin_par[inside],N_ELEMENTS(inside))*sin(angle),reform(SFRrebin[inside],N_ELEMENTS(inside))*sin(angle), $
reform(  HI_par[inside],N_ELEMENTS(inside))*sin(angle),reform(   HIrebin[inside],N_ELEMENTS(inside))*sin(angle), reform(      H2_par[inside],N_ELEMENTS(inside))*sin(angle),                 fltarr(N_ELEMENTS(inside))           , $
reform(  gas_sd[inside],N_ELEMENTS(inside))*sin(angle),reform(   HIrebin[inside],N_ELEMENTS(inside))*sin(angle),FORMAT = '(11F)'
ENDIF

IF keyword_set(outplot) THEN  device,filename = outplot + 'gasProf.eps' ELSE window,4,xsize = 400,ysize = 400
minr = 0
maxr = range
;range = maxr - minr
nbins = 40.0
r_bin = findgen(nbins + 1)*range/(nbins)
dr = range/nbins
area = ((r_bin + dr)*(r_bin + dr) - r_bin*r_bin)*1e6
prof_H2 = fltarr(nbins)
prof_HI = fltarr(nbins)
xaxes_grid = REBIN(xaxes + 0.0234375,N_ELEMENTS(xaxes),N_ELEMENTS(xaxes))
yaxes_grid = TRANSPOSE(REBIN(xaxes + 0.0234375,N_ELEMENTS(xaxes),N_ELEMENTS(xaxes)))
distance = SQRT((xaxes_grid)*(xaxes_grid) + (yaxes_grid)* (yaxes_grid))
FOR ir = 0, nbins - 1 DO BEGIN
    ind = where(distance le r_bin(ir + 1) AND distance gt r_bin(ir))
    area = !PI*(r_bin(ir + 1)*r_bin(ir + 1) - r_bin(ir)*r_bin(ir))*1e6
    IF keyword_set(useH2) THEN prof_H2[ir] = MEAN(cubeH2_surface_den[ind]/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18) ;TOTAL(cubeH2_surface_den[ind]/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18)/area
    prof_HI[ir] = MEAN(cubeHI_surface_den_prof[ind]/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18) ;TOTAL(cubeHI_surface_den[ind]/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18)/area
ENDFOR
plot,r_bin[0:nbins - 1] + dr/2.0,alog10(prof_HI/2.0*sin(angle)),yrange = [-4,2],xrange = [0,5],xtitle = 'Radius [kpc]';,ytitle = textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"
oplot,[r25,r25],[-5,6],linestyle = 1
IF keyword_set( useH2) THEN BEGIN
    oplot,r_bin[0:nbins - 1] + dr/2.0,alog10(prof_H2*sin(angle)),linestyle = 2
    oplot,[0,10],[1,1],linestyle = 1
    legend,['HI',textoidl('H_2')],linestyle = [0,2],/right
ENDIF
IF keyword_set(outplot) THEN device,/close ELSE stop

;window,4
;plot,HIrebin + H2rebin,gas_sd/(HIrebin + H2rebin),psym = 2,/ylog,/xlog,yrange = [1e-2,1e2],xrange = [1e-7,1e2],xstyle = 1,ystyle = 1
;IF keyword_set(outplot) THEN device,/close
END

;outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.4.00512'
;outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.4.0.25.00512'
;outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.4.1.00512'
;outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.4.2.00512'
;outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.4.4.00512'


;file=['/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.2304g1HBWK/schmidtlaw_res_obs_all_Ha0.750000.dat','/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.2304g1HBWKS/schmidtlaw_res_obs_all_Ha0.750000.dat']
;file =['/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWK/schmidtlaw_res_obs_all0.750000.dat','/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWKS/schmidtlaw_res_obs_all0.750000.dat'] 
;file =['/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWK/schmidtlaw_res_obs_all0.750000.dat','/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWKS/schmidtlaw_res_obs_all_Ha0.750000.dat'] 
;file = ['/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWK/schmidtlaw_res_obs_all0.750000_2.dat','/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWKS/schmidtlaw_res_obs_all0.750000_2.dat']
;outplot = '~/KS_shield_h516_3072'
;key = ['H2','S']
;color = [60,254]
;schmidtlaw_res_obs_master_out,file,color= color,key = key,/dwarf;,outplot= outplot

PRO schmidtlaw_res_obs_master_out,file,outplot = outplot,vetinari = vetinari,color= color, thick = thick,symbols = symbols,key = key,symsize = symsize,contour = contour,multiframe = multiframe,true = true,scaleHe = scaleHe, lowz = lowz,formatthick = formatthick,ctables = ctables,label = label,dwarf = dwarf
n = N_ELEMENTS(file)

;res = '4.00773'
;res = '0.257807'
;res = '1.99215'
;res = '1.00779'
;res = '0.749984'
res = '0.750000'

IF keyword_set(outplot) THEN BEGIN
    formatplot,/outplot,thick = formatthick
    set_plot,'ps' 
    outfile = outplot
    fgcolor = 0
    bgcolor = 255
    IF keyword_set(multiframe) THEN BEGIN
        xsize = 10*n
        ysize = 12
        mxTitSize = 1.5
        mxTitOffset = 2
    ENDIF ELSE BEGIN
        xsize = 18
        ysize = 18
    ENDELSE
ENDIF ELSE BEGIN
    formatplot
    set_plot,'x'
    fgcolor = 255
    bgcolor = 0
    IF keyword_set(multiframe) THEN BEGIN
        xsize = 400*n
        ysize = 475
        mxTitSize = 1.5
        mxTitOffset = 1
    ENDIF ELSE BEGIN
        xsize = 400
        ysize = 400
    ENDELSE
ENDELSE
!p.multi = 0

;IF keyword_set(vetinari) THEN BEGIN
;   prefix = '~/Data/MolecH/Cosmo/' 
;ENDIF ELSE BEGIN
;   prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
;ENDELSE
spawn,'hostname',hostname
IF hostname EQ 'ozma' OR hostname Eq 'quirm.math.grinnell.edu' THEN BEGIN 
   prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' 
   obsprefix = '~/Code/'
ENDIF ELSE BEGIN
   prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
   obsprefix = '~/code/'
ENDELSE

IF keyword_set(color) THEN BEGIN
    loadct,39
    IF color[0] EQ 1 THEN  colors = (findgen(n) + 1)*254/n ELSE colors = color
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n) + 39
    IF NOT keyword_set(thick) THEN $
      IF keyword_set(formatthick) THEN thick = fltarr(n) + 6 ELSE thick = fltarr(n) + 2
    IF NOT keyword_set(symbols) THEN symbols = fltarr(n) + 4
    obscolor = fgcolor
ENDIF ELSE BEGIN
    loadct,0
    colors = fltarr(n) + fgcolor
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n)
    IF NOT keyword_set(thick) THEN $
      IF keyword_set(formatthick) THEN thick = fltarr(n) + 6 ELSE  thick = (fltarr(n) + 1)
    IF NOT keyword_set(symbols) THEN symbols = fltarr(n) + 4 ;[4,15];symbols = (findgen(n)+2)*2
    obscolor = 150
ENDELSE

IF NOT keyword_set(symsize) THEN symsize = fltarr(n) + 2
nlevels = 154 
IF keyword_set(outplot) THEN color_hist = reverse(findgen(nlevels)) ELSE color_hist = findgen(nlevels) + (254 - nlevels)
obssym = 16
obssym2 = 15
obssymsize = 0.8;1.5
obssymsize2 = 1
c_line = [2,5,0]
c_line = [0,0,0]
thresh_points = 5000

;------------------------- Read in the obesrvational data ------------------------
fileobs = obsprefix + 'Datafiles/HIcubes/bigiel10.dat'
readcol,fileobs,type,name1,name2,logHI,e_logHI,logH2,e_logH2,logSFR,e_logSFR,format = '(A,A,A,F,F,F,F,F,F)'
name = name1 + name2
uniqname = name[uniq(name)]
weights = fltarr(n_elements(name)) + 1
;FOR i = 0, n_elements(uniqname) - 1 DO weights[where(name EQ uniqname[i])] = 1/float(n_elements(where(name EQ uniqname[i])))
;stop
fileobs_outr25 = obsprefix + 'Datafiles/HIcubes/bigiel10_outr25.dat'
readcol,fileobs_outr25,type_outr25,name1_outr25,name2_outr25,logHI_outr25,e_logHI_outr25,SFR_outr25,e_SFR_outr25,format = '(A,A,A,F,F,F,F)',/silent
readcol,fileobs_outr25,type_outr25_b,name1_outr25_b,name2_outr25_b,logHI_outr25_b,e_logHI_outr25_b,SFR_outr25_b,format = '(A,A,A,F,F,F)',/silent
name_outr25 = name1_outr25 + name2_outr25
name_outr25_b = name1_outr25_b + name2_outr25_b
uniqname_outr25_b = name_outr25_b[uniq(name_outr25_b)]
weights_outr25 = fltarr(n_elements(name_outr25)) + 1
;FOR i = 0, n_elements(uniqname_outr25) - 1 DO weights_outr25[where(name_outr25 EQ uniqname_outr25[i])] = 1/float(n_elements(where(name_outr25_b EQ uniqname_outr25_b[i])))

logHI_outr25 = logHI_outr25/1.36
e_logHI_outr25 = e_logHI_outr25/1.36
logSFR_outr25 = alog10(SFR_outr25*1e-5)
e_logSFR_outr25 = alog10(e_SFR_outr25*1e-5)
temp = where(type EQ 'Spirals',complement = inddwarf)
IF keyword_set(dwarf) THEN BEGIN
   logHI = logHI[inddwarf]
   logH2 = logH2[inddwarf]
   logSFR = logSFR[inddwarf]
   name1 = name1[inddwarf]
   name2 = name2[inddwarf]
ENDIF
logH2[where(logH2 eq 0)] = -6
logHI[where(logHI eq 0)] = -6
ind = where(logSFR NE 0); Select for gas with some star formation
logSFR = logSFR[ind]
logHI = alog10(10^logHI[ind]/1.36)
logH2 = alog10(10^logH2[ind]/1.36)
ind = where(name2 NE '925')
logSFR_lowz = logSFR[ind]
logHI_lowz = logHI[ind]
logH2_lowz = logH2[ind]

logH =   alog10((10^(logHI) + 10^(logH2)))
loggas = alog10((10^(logHI) + 10^(logH2))*1.36)
logH_lowz   = logH[ind]
loggas_lowz = loggas[ind]

;filebollato = '~/code/Datafiles/HIcubes/1kpcdata_Christensen.dat'
;readcol,filebollato,loggas_bolatto, logSFR_bolatto

min1 = -0.5
max1 =  2.5
min2 = -4.0
min2 = -4.7 ;Include outer disk
max2 = -0.5;0
nbins = 100
;bin1 = (max1 - min1)/nbins                      
;bin2 = (max2 - min2)/nbins                     
;nxbins=FLOOR((max1-min1) / bin1)
;nybins=FLOOR((max2-min2) / bin2)
;xbins=(FINDGEN(nxbins)*bin1)+min1+bin1/2.0
;ybins=(FINDGEN(nybins)*bin2)+min2+bin2/2.0
bin1 = 0.03                 
bin2 = 0.04              
nxbins=FLOOR((max1-min1) / bin1) + 1
nybins=FLOOR((max2-min2) / bin2) + 1
xbins=(FINDGEN(nxbins)*bin1) + min1
ybins=(FINDGEN(nybins)*bin2) + min2
bin1_outr25 = 0.2                 
bin2_outr25 = 0.2              
nxbins_outr25=FLOOR((max1-min1) / bin1_outr25) + 1
nybins_outr25=FLOOR((max2-min2) / bin2_outr25) + 1
xbins_outr25=(FINDGEN(nxbins_outr25)*bin1_outr25) + min1
ybins_outr25=(FINDGEN(nybins_outr25)*bin2_outr25) + min2

IF keyword_set(lowz) THEN BEGIN
    logH_points = logH_lowz
    logHI_points = logHI_lowz
    logH2_points = logH2_lowz
    loggas_points = loggas_lowz
    logSFR_points = logSFR_lowz
ENDIF ELSE BEGIN
    logH_points = logH
    logHI_points = logHI
    logH2_points = logH2
    loggas_points = loggas
    logSFR_points = logSFR
ENDELSE
histH  = HIST_2D_W(logH_points, logSFR_points, weights, min1=min1,max1=max1,min2=min2,max2=max2,bin1 = bin1, bin2 = bin2)
histHI = HIST_2D_W(logHI_points, logSFR_points, weights, min1=min1,max1=max1,min2=min2,max2=max2,bin1 = bin1, bin2 = bin2)
histHI_outr25 = HIST_2D_W(logHI_outr25, logSFR_outr25, weights_outr25, min1=min1,max1=max1,min2=min2,max2=max2,bin1 = bin1_outr25, bin2 = bin2_outr25)
;histHI_outr25 = HIST_2D_W(logHI_outr25[where(type_outr25 eq "Dwarfs")], logSFR_outr25[where(type_outr25 eq "Dwarfs")], weights_outr25[where(type_outr25 eq "Dwarfs")], min1=min1,max1=max1,min2=min2,max2=max2,bin1 = bin1_outr25, bin2 = bin2_outr25)
histH2 = HIST_2D_W(logH2_points, logSFR_points, weights, min1=min1,max1=max1,min2=min2,max2=max2,bin1 = bin1, bin2 = bin2)
histgas= HIST_2D_W(loggas_points,logSFR_points, weights, min1=min1,max1=max1,min2=min2,max2=max2,bin1 = bin1, bin2 = bin2)
IF NOT keyword_set(scaleHe) THEN BEGIN
   hist = histH
   hist_outr25 = histHI_outr25
    logx_points = logH_points
ENDIF ELSE BEGIN 
    hist = histgas 
    logx_points = loggas_points
ENDELSE
IF n_elements(logH_points) LE 500 THEN points = 1 ELSE points = 0

nbins_sim = 50
bin1_sim = (max1 - min1)/nbins_sim                      
bin2_sim = (max2 - min2)/nbins_sim                     
nxbins_sim=FLOOR((max1-min1) / bin1_sim)
nybins_sim=FLOOR((max2-min2) / bin2_sim)
xbins_sim=(FINDGEN(nxbins_sim)*bin1_sim)+min1+bin1_sim/2.0
ybins_sim=(FINDGEN(nybins_sim)*bin2_sim)+min2+bin2_sim/2.0

xsigma = 10.0^(findgen(600)/100 - 1.) ;xsigmalow
ysigma=2.5e-4*xsigma^1.4

;---------------------------- TOTAL GAS -----------------------------------------
IF keyword_set(outplot) THEN  device,filename = outfile + 'resSchmidt_all_MC.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,0,xsize = xsize,ysize = ysize

IF keyword_set(scaleHe) THEN $
  xtitle1 = textoidl('Log \Sigma')+"!lgas!n [M" + sunsymbol() + " pc!u-2!n]" ELSE $
  xtitle1 = textoidl('Log \Sigma')+"!lH!n [M"   + sunsymbol() + " pc!u-2!n]"

loadct,0
IF keyword_set(multiframe) THEN BEGIN
    multiplot,[n,1],/square,mxtitle = xtitle1,mxTitSize = mxTitSize, mxTitOffset = mxTitOffset,mtitle = label
    IF points THEN BEGIN
        plot, logx_points,logSFR_points,xrange = [min1,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1, ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]"
        oplot,logx_points,logSFR_points,color = obscolor, psym = symcat(obssym), symsize = obssymsize 
;       oplot,alog10(loggas_bolatto),alog10(logSFR_bolatto),color = obscolor, psym = symcat(obssym2), symsize = obssymsize2
    ENDIF ELSE $
      contour,hist,xbins,ybins,         xrange = [min1,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1, ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", nlevels = nlevels, c_colors = color_hist, /fill
ENDIF ELSE BEGIN
    IF points THEN BEGIN
        plot, logx_points,logSFR_points,xrange = [min1,max1],yrange = [min2,max2], xstyle = 1,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xtitle=xtitle1,/nodata, title = label
        oplot,logx_points,logSFR_points,color = obscolor, psym = symcat(obssym), symsize = obssymsize
;       oplot,alog10(loggas_bolatto),alog10(logSFR_bolatto),color = obscolor, psym = symcat(obssym2), symsize = obssymsize2
     ENDIF ELSE BEGIN
        contour,hist,xbins,ybins,xrange = [min1,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1, ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]",  nlevels = nlevels, c_colors = color_hist, /fill, title = label, xtitle=xtitle1,min_value = min(hist[where(hist NE 0)]) - 1,max_value = max(hist)
;        contour_levels = contourlevels(hist,[0.25,0.50,0.75,0.90])
;        contour,hist,xbins,ybins,xrange = [min1,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1, ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", /fill, title = label, xtitle=xtitle1, levels = max(hist)*[0.25,0.50,0.75,0.90]
        contour_levels = contourlevels(hist_outr25,[0.25,0.50,0.75,0.90])
        contour,hist_outr25,xbins_outr25,ybins_outr25,levels = max(hist_outr25)*[0.1,0.25,0.50,0.75],/overplot, thick = 4
     ENDELSE
ENDELSE

FOR i = 0, n -1 DO BEGIN
    loadct,ctables[i]
    readcol,file[i],radius,fuv,mic24,SFR_par,SFR_obs,HI_par,HI,H2_par,H2,gas_sd,total_gas
    IF keyword_set(true)    THEN total_gas = gas_sd
    IF keyword_set(scaleHe) THEN total_gas = total_gas*1.36
    ind = where(alog10(SFR_obs) ge -3.8)
    ind_inr25 = where(radius lt max(radius/2) + float(res)/2,complement = ind_outr25)
    IF n_elements(total_gas) GE thresh_points THEN BEGIN
        hist_sim   = HIST_2D(alog10(total_gas), alog10(SFR_obs),min1=min1 + bin1_sim/2,max1=max1 - bin1_sim/2.0,min2=min2,max2=max2 - bin2_sim,bin1 = bin1_sim, bin2 = bin2_sim)
        contour,hist_sim,xbins_sim[0:(size(hist_sim))[1] - 1],ybins_sim[0:(size(hist_sim))[2] - 1],xrange = [min1,max1],yrange = [min2,max2],color = colors[i],/overplot,levels = [1,2,3,4],thick = thick[i],c_linestyle = c_line;levels = [1,3,6]
     ENDIF ELSE  oplot,alog10(total_gas[ind_inr25]),alog10(SFR_obs[ind_inr25]),psym = symcat(symbols[i]),color = colors[i],thick = 3,symsize = symsize[i] ;160,thick = thick[i]
    oplot,alog10(total_gas[ind_outr25]),alog10(SFR_obs[ind_outr25]),psym = symcat(16),color = colors[i],symsize = 0.75 ;symsize[i],thick = thick[i]
    stop
    IF keyword_set(multiframe) THEN multiplot
ENDFOR
IF keyword_set(multiframe) THEN multiplot,/reset
;IF keyword_set(key) THEN legend,[key,'Bigiel et al., 2010','Bolatto et al., 2011'],psym = [symbols,obssym,obssym2],color = [colors,obscolor,obscolor],thick = [thick,1,1]
IF  keyword_set(key) THEN legend,[key,'Bigiel et al. 2010'],                        psym = [symbols,obssym],        color = [colors,obscolor],         thick = [thick,1]
IF keyword_set(outplot) THEN  device,/close ELSE stop

;----------------------------------------------------------- HI ---------------------------------------------------------------------------------------------
IF keyword_set(outplot) THEN  device,filename = outfile + 'resSchmidt_HI.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,1,xsize = xsize,ysize = ysize

loadct,0
IF keyword_set(multiframe) THEN BEGIN
    multiplot,[n,1],/square,mxtitle = textoidl('Log \Sigma')+"!lH!n [M"+sunsymbol()+" pc!u-2!n]",mxTitSize = mxTitSize, mxTitOffset = mxTitOffset, mtitle = label
    contour,histHI,xbins,ybins,xrange = [min1,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", nlevels = nlevels,c_colors = color_hist
ENDIF ELSE contour,histHI,xbins,ybins,xrange = [min1,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xtitle=textoidl('Log \Sigma_{HI} [M')+sunsymbol()+" pc!u-2!n]", nlevels = nlevels,c_colors = color_hist, title = label

IF keyword_set(color) THEN loadct,39
FOR i = 0, n - 1 DO BEGIN 
    readcol,file[i],fuv,mic24,SFR_par,SFR_obs,HI_par,HI,H2_par,H2,gas_sd,total_gas
    IF keyword_set(contour) THEN BEGIN
        hist_sim   = HIST_2D(alog10(HI), alog10(SFR_obs),min1=min1 + bin1_sim/2,max1=max1 - bin1_sim/2.0,min2=min2,max2=max2 - bin2_sim,bin1 = bin1_sim, bin2 = bin2_sim)
        contour,hist_sim,xbins_sim[0:(size(hist_sim))[1] - 1],ybins_sim[0:(size(hist_sim))[2] - 1],xrange = [min1,max1],yrange = [min2,max2],color = colors[i],/overplot,levels = [1,3,6],thick = thick[i],c_linestyle = c_line
    ENDIF ELSE  oplot,alog10(HI),alog10(SFR_obs),psym = symbols[i],color = colors[i],thick = thick[i],symsize = symsize[i]
    IF keyword_set(multiframe) THEN BEGIN
        multiplot
        if i lt n -1 THEN BEGIN
            loadct,0
            contour,histHI,xbins[0:98],ybins[0:99],xrange = [min1 + 0.01,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill, nlevels = nlevels,c_colors = color_hist        
            loadct,39
        ENDIF
    ENDIF
ENDFOR
IF keyword_set(multiframe) THEN multiplot,/reset
IF keyword_set(key) THEN legend,key,psym = symbols,color = colors,thick = thick
IF keyword_set(outplot) THEN  device,/close ELSE stop

;------------------------------------------------------------ H2 ------------------------------------------------------------------------
IF keyword_set(outplot) THEN  device,filename = outfile + 'resSchmidt_H2.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2 ELSE window,2,xsize = xsize,ysize = ysize
loadct,0
IF keyword_set(multiframe) THEN BEGIN
    multiplot,[n,1],/square,mxtitle = textoidl('Log \Sigma_{H_2} [M')+sunsymbol()+" pc!u-2!n]",mxTitSize = mxTitSize, mxTitOffset = mxTitOffset,mtitle = label
    contour,histH2,xbins,ybins,xrange = [min1,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", nlevels = nlevels,c_colors = color_hist
ENDIF ELSE contour,histH2,xbins,ybins,xrange = [-4,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xtitle=textoidl('Log \Sigma_{H_2} [M')+sunsymbol()+" pc!u-2!n]", nlevels = nlevels,c_colors = color_hist,title = label;,xrange = [min1,max1];[-3.5,max1]

IF keyword_set(color) THEN loadct,39
FOR i = 0, n - 1 DO BEGIN 
    readcol,file[i],fuv,mic24,SFR_par,SFR_obs,HI_par,HI,H2_par,H2,gas_sd,total_gas
    IF keyword_set(contour) THEN BEGIN
        hist_sim   = HIST_2D(alog10(H2), alog10(SFR_obs),min1=min1 + bin1_sim/2,max1=max1 - bin1_sim/2.0,min2=min2,max2=max2 - bin2_sim,bin1 = bin1_sim, bin2 = bin2_sim)
        contour,hist_sim,xbins_sim[0:(size(hist_sim))[1] - 1],ybins_sim[0:(size(hist_sim))[2] - 1],xrange = [min1,max1],yrange = [min2,max2],color = colors[i],/overplot,levels = [1,3,6],thick = thick[i],c_linestyle = c_line
    ENDIF ELSE oplot,alog10(H2),alog10(SFR_obs),psym = symcat(symbols[i]),color = colors[i],thick = thick[i],symsize = symsize[i]
    IF keyword_set(multiframe) THEN BEGIN
        multiplot
        if i lt n -1 THEN BEGIN
            loadct,0
            contour,histH2,xbins,ybins,xrange = [min1 + 0.01,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill, nlevels = nlevels,c_colors = color_hist        
            loadct,39
        ENDIF
    ENDIF
ENDFOR
IF keyword_set(key) THEN legend,key,psym = symbols,color = colors,thick = thick
IF keyword_set(outplot) THEN device,/close ELSE stop
END

; device,filename = outfile + 'resSchmidt_all.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 7,xoffset =  2,yoffset =  2
; multiplot,[3,1],/square
; contour,hist,xbins[0:98],ybins[0:99],xrange = [min1,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xtitle=textoidl('Log \Sigma')+"!lH!n [M"+sunsymbol()+" pc!u-2!n]", nlevels = nlevels,c_colors = color_hist
; multiplot
; contour,histHI,xbins[0:98],ybins[0:99],xrange = [-0.49,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill, xtitle=textoidl('Log \Sigma_{HI} [M')+sunsymbol()+" pc!u-2!n]", nlevels = nlevels,c_colors = color_hist
; multiplot
; contour,histH2,xbins[0:98],ybins[0:99],xrange = [-0.49,max1],yrange= [min2,max2], xstyle = 1, ystyle = 1,/fill, xtitle=textoidl('Log \Sigma_{H_2} [M')+sunsymbol()+" pc!u-2!n]",nlevels = nlevels,c_colors = color_hist
; multiplot,/reset
; device,/close
