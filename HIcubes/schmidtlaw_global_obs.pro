;filename = 'h603.cosmo50cmb.2304g5bwK.00512'
;cd,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG'
;outplot ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/globalSK.eps'
;outplot = '~/h516.cosmo25cmb.paper.globalSK.eps'

pro schmidtlaw_global_obs_master,color=color,outplot = outplot,halo_str = halo_str
!p.multi = 0
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    fgcolor  = 0
ENDIF ELSE BEGIN
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    fgcolor = 255
ENDELSE

dir=['/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir', $
     '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir', $
     '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir', $
     '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464.dir']

pfile = ['../../h516.cosmo25cmb.1536g6MbwK.param', $
         '../../h516.cosmo25cmb.1536g3HBWK.param', $
         '../../h516.cosmo25cmb.1536g3HBWK.param', $
         '../../h516.cosmo25cmb.1536g6HBWK.jeans.prev.param']

filename = ['h516.cosmo25cmb.1536g6MbwK.00512', $
            'h516.cosmo25cmb.1536g3HBWK_noH2SF.00512', $
            'h516.cosmo25cmb.1536g3HBWK.00512', $
            'h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464']
useH2 = [0,1,1,1]
key = ['no H2','H2','H2 SF','MC SF']

;-----------------------------------------------------

prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
base = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir',$
        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir']
dir = prefix+base
pfile = ['../../h516.cosmo25cmb.3072g1MBWK.param',$
         '../../h516.cosmo25cmb.2304g14HBWK.param']
filename = ['h516.cosmo25cmb.3072g1MBWK.00492',$
            'h516.cosmo25cmb.2304g14HBWK.00512']
tfile = ['h516.cosmo25cmb.3072g1MBWK.00492.halo.1.std',$
         'h516.cosmo25cmb.2304g14HBWK.00512.halo.1.std']
useH2 = [0,1]
key = ['no '+textoidl('H_2'), textoidl('H_2')]

;-----------------------------------------------------------

prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
base = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir',$
        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir']
dir = prefix+base
pfile = ['../../h516.cosmo25cmb.3072g1MBWK.param',$
         '../../h516.cosmo25cmb.3072g14HBWK.param']
filename = ['h516.cosmo25cmb.3072g1MBWK.00492',$
            'h516.cosmo25cmb.3072g14HBWK.00512']
tfile = ['h516.cosmo25cmb.3072g1MBWK.00492.halo.1.std',$
         'h516.cosmo25cmb.3072g14HBWK.00512.halo.1.std']
useH2 = [0,0]
key = ['no '+textoidl('H_2'), textoidl('H_2')]
key = ['DnoH2','DH2']

;-----------------------------------------------------

data=fltarr(2,N_ELEMENTS(dir))
FOR i = 0, N_ELEMENTS(dir) - 1 DO BEGIN
    data[*,i] = schmidtlaw_global_obs(dir[i],filename[i],pfile[i],useH2 = useH2[i],Halpha=1,tipsyfile = tfile[i])
ENDFOR      

IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  colors = (findgen(N_ELEMENTS(filename)) + 1)*240/N_ELEMENTS(filename) else colors = color
    symbols = fltarr(N_ELEMENTS(filename)) + 4
    obscolor = fgcolor
ENDIF ELSE BEGIN
    loadct,0    
;    colors = (findgen(N_ELEMENTS(filename)) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(files)) + 5
    colors = fltarr(N_ELEMENTS(filename)) + fgcolor
;    symbols = (findgen(N_ELEMENTS(filename))+1)*2
    symbols = (findgen(N_ELEMENTS(filename))+2)*2
    obscolor = 150
ENDELSE
if (KEYWORD_SET(outplot)) then begin
    set_plot,'ps'
    device,filename = outplot,/color,bits_per_pixel= 8,/times,ysize=18,xsize=18,xoffset =  2,yoffset =  2
endif else begin
    set_plot,'x'
    window,1
endelse
xsigma = 10.0^(findgen(600)/100 - 1.) ;xsigmalow
ysigma=2.5e-4*xsigma^1.4

readcol,'~/code/HIcubes/ks98.dat',name,D,logHI,logH2,logH,logSFR,tdyn,co,HI,halpha;,format='(A9D)'
readcol,'~/code/HIcubes/uavpl.dat',logHdarf,logSFRdwarf
plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,2.5], yrange = [-4,-0.5], xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"   
FOR i = 0, N_ELEMENTS(dir) - 1 DO oplot,[alog10(data[0,i]),alog10(data[0,i])],[alog10(data[1,i]),alog10(data[1,i])],psym = symbols[i],color = colors[i] 
oplot,logH,logSFR,psym = 2,color = obscolor
oplot,logHdarf,logSFRdwarf,psym = 4,color = obscolor
oplot,[-0.2,0.2],[-1,-1]
oplot,[0,0],[-0.8,-1.2]
legend,[key,'Kennicutt 98'],psym = [symbols,2],color = [colors,obscolor],/bottom,/right
if (KEYWORD_SET(outplot)) then device,/close else stop

END

function schmidtlaw_global_obs,dir,filename,pfile,useH2 = useH2,tipsyfile = tipsyfile,Halpha = Halpha,verbose = verbose,camera = camera, intq = intq,center = center,halo_str = halo_str, angle = angle,intrinsic = intrinsic
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
speed_o_light = 299792458 ;m per sec
molec_weight = (0.76*1 + 0.24*4.0)
;angle = !pi/4 ;Angle on
;angle = !PI/2 ;Face on
IF NOT KEYWORD_SET(angle) THEN BEGIN 
    angle_str = '90'
    angle = !PI/2
ENDIF ELSE BEGIN 
    angle_str = STRTRIM(FIX(angle),2)
    angle = !PI*angle/180.0
ENDELSE

;extno = 15  ;Angle on
;extno = 14 ; Face on
IF NOT KEYWORD_SET(camera) THEN extno = 14 ELSE extno = camera
IF NOT KEYWORD_SET(intq) THEN extno_int = 32 ELSE extno_int = intq
IF NOT KEYWORD_SET(center) THEN center = [0,0] ELSE center = center
IF NOT KEYWORD_SET(halo_str) THEN halo_str = '1'

cd,dir
units = tipsyunits(pfile)
massunit = units.massunit
lengthunit = units.lengthunit
timeunit = units.timeunit
rtipsy,filename,h,g,d,s,/justhead
r25 = opticalRadii(filename = filename, extno = extno, center = center,halo_str = halo_str,verbose = verbose)
cd,dir

stop
cubeHI_surface_den = read_dencube_fits(filename + '.halo.' + halo_str + '.' + angle_str + '.cube.mom0.fits',headerHI,/noscale)/amu_per_gm/gm_per_msol*cm_per_kpc*cm_per_kpc/1d6 ;converting to Msol/pc
;cubeHI_surface_den = cubeHI_surface_den ;/2.0
headerHI.CDELT1 = headerHI.CDELT1;*h.time
headerHI.CDELT2 = headerHI.CDELT2;*h.time
IF NOT (where(~finite(cubeHI_surface_den)))[0] EQ -1 THEN cubeHI_surface_den[where(~finite(cubeHI_surface_den))] = 0
IF keyword_set(useH2) THEN BEGIN
    cubeH2 = read_dencube_fits(filename+'.halo.' + halo_str + '.H2.arr.' + angle_str + '.fits',headerH2)*2.0
    headerH2.CDELT1 = headerH2.CDELT1*h.time
    headerH2.CDELT2 = headerH2.CDELT2*h.time
    pix_areaH2 = headerH2.CDELT1*headerH2.CDELT2*1e6
;cubeH2_surface_den = TRANSPOSE(cubeH2)/amu_per_gm/gm_per_msol*cm_per_kpc*cm_per_kpc/1d6*2.0*2.0 ;/pix_areaH2*gm_per_msol*amu_per_gm*2.0
    cubeH2_surface_den = TRANSPOSE(cubeH2)/pix_areaH2;*gm_per_msol*amu_per_gm*2.0
ENDIF

xaxes = (findgen(headerHI.NAXIS1) - headerHI.NAXIS1/2.0)*headerHI.CDELT1
yaxes = (findgen(headerHI.NAXIS2) - headerHI.NAXIS1/2.0)*headerHI.CDELT2
radius = fltarr(N_ELEMENTS(xaxes),N_ELEMENTS(yaxes))
deltaS = headerHI.CDELT1*headerHI.CDELT2*1d6 ;In pc^2
FOR ix = 0,N_ELEMENTS(xaxes) - 1 DO $
  FOR iy = 0,N_ELEMENTS(yaxes) - 1 DO $
  radius[ix,iy] = SQRT(xaxes[ix]*xaxes[ix] + yaxes[iy]*yaxes[iy]/SIN(angle)/SIN(angle))
IF KEYWORD_SET(verbose) THEN BEGIN
    window,4,xsize = 800,ysize = 800
    contour,alog10(cubeHI_surface_den),xaxes,yaxes,/fill,nlevels = 512,xrange = [-7,7],yrange = [-7,7]
    contour,radius,xaxes,yaxes,/overplot,levels = [1,2,3,4,5,6,7,8,9,10,11,12],color = 100
    contour,radius,xaxes,yaxes,/overplot,levels = [0,r25]
ENDIF

sfilename = filename + '.' + halo_str + '/broadband.fits'
fuv_num = 0
mic24_num = 15
fuv_lambda = 1.5180502483701E-07
mic24_lambda = 2.38717045592424E-05
filtermags = mrdfits(sfilename,13)
print,filtermags[5].filter,filtermags[5].ab_mag0
print,'g-r',filtermags[3].ab_mag0 - filtermags[4].ab_mag0
sr_surface_den = mrdfits(sfilename,extno,sr_head)
dxy = sr_head[where(strcmp('CD1_1',sr_head,5))]
dxy = double((strsplit(dxy,' ',/extract))[2])
sr_NAXIS = (size(sr_surface_Den))[1] ;500 ;480.0
sr_range = dxy*sr_NAXIS       ;headerHI.NAXIS1*headerHI.CDELT1
;sr_range = 24 ;headerHI.NAXIS1*headerHI.CDELT1
;sr_NAXIS = 480.0
xaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0 - center[0])*sr_range/sr_NAXIS + sr_range/sr_NAXIS/2.0
yaxes_sr = (findgen(sr_NAXIS) - sr_NAXIS/2.0 - center[1])*sr_range/sr_NAXIS + sr_range/sr_NAXIS/2.0
deltaS_SR = sr_range/sr_NAXIS*sr_range/sr_NAXIS
radius_sr = fltarr(N_ELEMENTS(xaxes_sr),N_ELEMENTS(yaxes_sr))
FOR ix = 0,N_ELEMENTS(xaxes_sr) - 1 DO $
  FOR iy = 0,N_ELEMENTS(yaxes_sr) - 1 DO $
  radius_sr[ix,iy] = SQRT(xaxes_sr[ix]*xaxes_sr[ix] + yaxes_sr[iy]*yaxes_sr[iy]/SIN(!PI/4)/SIN(!PI/4))
print,minmax(xaxes_sr)
print,minmax(yaxes_sr)

;Sunrise units are in W/m/n^2/sr
;Bigiel Units are in MJy*str^(-1)
fuv_surface_den = sr_surface_den[*,*,fuv_num]*fuv_lambda*fuv_lambda/speed_o_light/1e-20 
mic24_surface_den = sr_surface_den[*,*,mic24_num]*mic24_lambda*mic24_lambda/speed_o_light/1e-20
SFR_surface_den = 3.2e-3*mic24_surface_den + 8.1e-2*fuv_surface_den ;EQ 1, Bigiel et al, 2008
fuv_surface_den = ROTATE(fuv_surface_den,5)
mic24_surface_den = ROTATE(mic24_surface_den,5)
SFR_surface_den = ROTATE(SFR_surface_den,5)

IF KEYWORD_SET(Halpha) THEN BEGIN
    sfilename = filename + '.' + halo_str + '/mcrx.fits'
    lambda1 = 6.54e-7
    lambda2 = 6.58e-7
    lc1 = 6.4e-7 ;Continuum is defined at these wavelengths in Kennicutt 89, page 1097
    lc2 = 6.8e-7
    lambdacube = mrdfits(sfilename, extno_int)
    lambda = lambdacube.LAMBDA
    indlambda = where(lambda ge lambda1 AND lambda le lambda2)
    temp = MIN(ABS(lc1 - lambda),c1)
    temp = MIN(ABS(lc1 - lambda),c2)
    LLambda = lambdacube.L_lambda_scatter1
    cont = (LLambda(c1) + LLambda(c2))/2.0
;    plot,lambda,llambda,xrange = [6e-7,7e-7]
;    oplot,[6e-7,7e-7],[cont,cont],linestyle = 1
    camera1scat = mrdfits(sfilename, 21)
ENDIF

IF KEYWORD_SET(verbose) OR keyword_set(intrinsic) OR NOT keyword_set(Halpha) THEN BEGIN
    dtime = 100.0e6
    rtipsy,filename+'.halo.' + halo_str + '.std',h,g,d,s
    readarr,filename+'.halo.' + halo_str + '.HI',h,HI,/ascii,part = gas
    IF keyword_set(useH2) THEN  readarr,filename+'.halo.' + halo_str + '.H2',h,H2,/ascii,part = gas
    tcurrent = max(s.tform*timeunit)
    massform = MAX(s.mass)*massunit
    radius_g = SQRT(g.x*g.x + g.y*g.y)*lengthunit*h.time
    radius_g = SQRT(g.x*g.x + g.y*g.y)*lengthunit*h.time
    radius_s = SQRT(s.x*s.x + s.y*s.y)*lengthunit*h.time
    sinterior = where(radius_s le r25 AND s.tform*timeunit GT tcurrent - dtime)
    ginterior = where(radius_g le r25)
    SFR_surface_den_global_true = N_ELEMENTS(sinterior)*massform/(!PI*r25^2*dtime)
    HI_surface_den_global_true = TOTAL(g[ginterior].mass*(HI[ginterior]))*massunit/(!PI*r25^2*1e6)
    IF keyword_set(useH2) THEN H2_surface_den_global_true = TOTAL(g[ginterior].mass*2.0*H2[ginterior])*massunit/(!PI*r25^2*1e6) ELSE H2_surface_den_global_true = 0
    IF keyword_set(useH2) THEN Gas_surface_den_global_true = TOTAL(g[ginterior].mass*(HI[ginterior] + 2.0*H2[ginterior])*massunit)/(!PI*r25^2*1e6) ELSE $
      Gas_surface_den_global_true = TOTAL(g[ginterior].mass*(HI[ginterior])*massunit)/(!PI*r25^2*1e6)
ENDIF ELSE BEGIN
    SFR_surface_den_global_true = 0
    HI_surface_den_global_true = 0
    IF keyword_set(useH2) THEN H2_surface_den_global_true = 0
    Gas_surface_den_global_true = 0
ENDELSE


IF keyword_set(useH2) THEN $
  gas_surface_den_global = TOTAL([cubeHI_surface_den[where(radius le r25)],cubeH2_surface_den[where(radius le r25)]])*deltaS/!PI/r25/r25/1e6 ELSE $
  gas_surface_den_global = TOTAL(cubeHI_surface_den[where(radius le r25)])*deltaS/!PI/r25/r25/1e6
;h.time*h.time?
;Divide deltaS by h.time to compensate for expansion that already counted

; TOTAL([cubeHI_surface_den[where(radius le r25)],cubeH2_surface_den[where(radius le r25)]])/N_ELEMENTS(where(radius le r25))/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18 ELSE $
;  gas_surface_den_global = TOTAL(cubeHI_surface_den[where(radius le r25)])/N_ELEMENTS(where(radius le r25))/gm_per_msol/amu_per_gm*3.08568021d18*3.08568021d18



SFR_surface_den_global_24 = TOTAL(SFR_surface_den[where(radius_sr le r25)])*deltaS_SR/!PI/r25/r25
IF keyword_set(Halpha) THEN SFR_surface_den_global = 7.9e-42*Halpha(filename + '.' + halo_str + '/mcrx.fits',verbose = verbose,extno_int = extno_int)/!PI/r25/r25 ELSE SFR_surface_den_global = SFR_surface_den_global_true
print,'R_25: ',r25
print,' '
print,'True SFR: ',SFR_surface_den_global_true
print,'FIR/24mic: ',SFR_surface_den_global_24*SIN(angle)
print,'Halpha: ',SFR_surface_den_global,SFR_surface_den_global*!PI*r25*r25
print,' '
print,'True HI: ',HI_surface_den_global_true
print,'Obs HI: ',TOTAL(cubeHI_surface_den[where(radius le r25)])*deltaS/!PI/r25/r25/1e6*SIN(angle)
print,' '
IF keyword_set(useH2) THEN BEGIN
    print,'True H2: ',H2_surface_den_global_true
    print,'Obs H2: ',TOTAL(cubeH2_surface_den[where(radius le r25)])*deltaS/!PI/r25/r25/1e6*SIN(angle)
    print,' '
ENDIF
print,'True Gas: ',gas_surface_den_global_true
print,'Obs Gas: ',gas_surface_den_global*SIN(angle)
;print,alog10([gas_surface_den_global_true,SFR_surface_den_global])
IF KEYWORD_SET(verbose) THEN stop
if SFR_surface_den_global lt 0 THEN SFR_surface_den_global = SFR_surface_den_global_true
return,[gas_surface_den_global*SIN(angle),SFR_surface_den_global]
end
