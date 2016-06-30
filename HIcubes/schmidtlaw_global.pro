;filename = 'h603.cosmo50cmb.2304g5bwK.00512'
;cd,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG'
;outplot ='/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/globalSK.eps'
;outplot = '~/h516.cosmo25cmb.paper.globalSK.eps'

;outplot = '/astro/store/student-scratch1/christensen/Aquilla/gas_cc_c5/aq-c-5-sph_sk.g'  
;outplot = '/astro/store/student-scratch1/christensen/Aquilla/gas_met_c5/aq-c-5-sph_sk.g'  

pro schmidtlaw_global_master,color=color,outplot = outplot,verbose = verbose
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
    set_plot,'ps'
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
    set_plot,'x'
ENDELSE

;-----------------------------------------------------------

prefix = '/astro/store/student-scratch1/christensen/Aquilla/gas_cc_c5/'
steps = ['00046','00060','00084','00127','00166','00225','00327','00512']
first = 2

;prefix = '/astro/store/student-scratch1/christensen/Aquilla/gas_met_c5/'
;steps = ['00084','00127','00225','00512']
;;r25 = [1.57,0.89,0.84,4.30]
;first = 0

base = 'aq-c-5-sph' ;+steps
dir = prefix+base+'.'+steps+ '.dir'
pfile = strarr(N_ELEMENTS(steps)) + '/astro/store/student-scratch1/christensen/Aquilla/aq-c-5-sph.param'
filename = base+'.'+steps
tfile = base+'.'+steps + '.halo.1.std'
useH2 = intarr(N_ELEMENTS(steps))
key = steps

;-----------------------------------------------------

data=fltarr(4,N_ELEMENTS(dir))
FOR i = 0,N_ELEMENTS(dir) - 1 DO BEGIN
    data[*,i] = schmidtlaw_global(dir[i],filename[i],pfile[i],useH2 = useH2[i],tipsyfile = tfile[i],verbose = verbose)
ENDFOR      

IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  colors = (findgen(N_ELEMENTS(filename)) + 1)*240/N_ELEMENTS(filename) else colors = color
    symbols = fltarr(N_ELEMENTS(filename)) + 4
    obscolor = fgcolor
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(N_ELEMENTS(filename)) + fgcolor
    symbols =  fltarr(N_ELEMENTS(filename)) + 4 ;(findgen(N_ELEMENTS(filename))+2)*2
    obscolor = 150
    symsizes = (findgen(N_ELEMENTS(filename)) + 1)*0.25 + 0.5
ENDELSE
if (KEYWORD_SET(outplot)) then begin
    set_plot,'ps'
    device,filename = outplot + '.eps',/color,bits_per_pixel= 8,/times,ysize=18,xsize=18
endif else begin
    set_plot,'x'
    window,1
endelse
xsigma = 10.0^(findgen(600)/100 - 1.) ;xsigmalow
ysigma=2.5e-4*xsigma^1.4

readcol,'~/code/HIcubes/ks98.dat',name,D,logHI,logH2,logH,logSFR,tdyn,co,HI,halpha;,format='(A9D)'
;plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,2.5], yrange = [-4,0.5], xtitle=textoidl('Log \Sigma')+"!lHI!n [M"+sunsymbol()+" pc!u-2!n]"  
plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,2.5], yrange = [-4,0.5], xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]" 
;FOR i = 0, N_ELEMENTS(dir) - 1 DO oplot,[alog10(data[0,i]),alog10(data[0,i])],[alog10(data[2,i]),alog10(data[2,i])],psym = symbols[i],color = colors[i],symsize = symsizes[i]
FOR i = 0, N_ELEMENTS(dir) - 1 DO oplot,[alog10(data[1,i]),alog10(data[1,i])],[alog10(data[2,i]),alog10(data[2,i])],psym = symcat(symbols[i]),color = colors[i],symsize = symsizes[i]
;oplot,logH,logSFR,psym = 2,color = obscolor
;oplot,1.4*logH,logSFR,psym = 2,color = obscolor
;oplot,[-0.2,0.2],[-1,-1]
;oplot,[0,0],[-0.8,-1.2]
;legend,[key,'Kennicutt 98'],psym = [symbols,2],color = [colors,obscolor],symsize = [symsizes,1],/bottom,/right
legend,[key],psym = symcat([symbols]),color = [colors],symsize = [symsizes],/bottom,/right
if (KEYWORD_SET(outplot)) then device,/close else stop

if (KEYWORD_SET(outplot)) then device,filename = outplot + '2.eps',/color,bits_per_pixel= 8,/times,ysize=12,xsize=18
;plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,5], yrange = [-3,4], xtitle=textoidl('Log \Sigma')+"!lHI!n [M"+sunsymbol()+" pc!u-2!n]"   
plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,5], yrange = [-3,4], xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]"   
;FOR i = first, N_ELEMENTS(dir) - 1 DO oplot,[alog10(data[0,i]),alog10(data[0,i])],[alog10(data[2,i]),alog10(data[2,i])],psym = symbols[i],color = colors[i],symsize = symsizes[i]
FOR i = first, N_ELEMENTS(dir) - 1 DO oplot,[alog10(data[1,i]),alog10(data[1,i])],[alog10(data[2,i]),alog10(data[2,i])],psym = symcat(symbols[i]),color = colors[i],symsize = symsizes[i]
;oplot,logH,logSFR,psym = 2,color = obscolor
;oplot,logH,logSFR,psym = 2,color = obscolor
;oplot,[-0.2,0.2],[-1,-1]
;oplot,[0,0],[-0.8,-1.2]
;legend,[key,'Kennicutt 98'],psym = [symbols,2],color = [colors,obscolor],symsize = [symsizes,1],/bottom,/right
legend,[key],psym = symcat([symbols]),color = [colors],symsize = [symsizes],/bottom,/right
if (KEYWORD_SET(outplot)) then device,/close else stop

if (KEYWORD_SET(outplot)) then begin
    openw,1,outplot+'.dat'
    printf,1,'Name   ','Radius [kpc]   ','Sigma_HI [M_sun/pc^2]   ','Sigma_gas [M_sun/pc^2]   ','Sigma_SFR [M_sun/kpc^2/yr]   ',format = '(5A)'
    FOR i = 0, N_ELEMENTS(dir) - 1 DO printf,1,base+steps[i],data[3,i],data[0,i],data[1,i],data[2,i]
    close,1
ENDIF

END

function schmidtlaw_global,dir,filename,pfile,useH2 = useH2,tipsyfile = tipsyfile,verbose = verbose,camera = camera, intq = intq,center = center,r25 = r25
scale_scalelength = 2

cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
speed_o_light = 299792458 ;m per sec
molec_weight = (0.76*1 + 0.24*4.0)
angle = !pi/4 ;Angle on
angle = !PI/2 ;Face on

cd,dir
units = tipsyunits(pfile,silent = 1)
massunit = units.massunit
lengthunit = units.lengthunit
timeunit = units.timeunit

dtime = 100.0e6
rtipsy,tipsyfile,h,g,d,s
tcurrent = max(s.tform*timeunit)
print,'Scale factor: ',h.time

g.mass = g.mass*massunit
g.x = g.x*lengthunit*h.time
g.y = g.y*lengthunit*h.time
g.z = g.z*lengthunit*h.time 
radius_g = SQRT(g.x*g.x + g.y*g.y)

s.mass = s.mass*massunit
s.x = s.x*lengthunit*h.time
s.y = s.y*lengthunit*h.time
s.z = s.z*lengthunit*h.time  
radius_s = SQRT(s.x*s.x + s.y*s.y)

;Find the appropriate radii
IF 0 THEN BEGIN
    r25 = scalelength(s,rmin = 1,rmax = 10,verbose = verbose)*scale_scalelength
ENDIF ELSE BEGIN
    ind_gall = where(radius_g lt 30 AND ABS(g.z) lt 3 AND g.tempg le 1e5)
    gall = g[ind_gall]
    radius_gall = radius_g[ind_gall]
    radius_gall = radius_gall[SORT(radius_gall)]
    IF NOT KEYWORD_SET(r25) THEN BEGIN
        r25 = radius_gall[N_ELEMENTS(radius_gall)/2]
        IF KEYWORD_SET(verbose) THEN BEGIN
            histogramp,radius_gall,/normalize
            histogramp,radius_g,linestyle = 2,/overplot,/normalize
            oplot,[r25,r25],[0,1],linestyle = 1
        ENDIF
    ENDIF
ENDELSE

;Selecting the right type of gas
IF 0 THEN BEGIN
    readarr,filename+'.halo.1.HI',h,HI,/ascii,part = gas
    IF keyword_set(useH2) THEN  readarr,filename+'.halo.1.H2',h,H2,/ascii,part = gas
    ginterior = where(radius_g le r25)
    IF keyword_set(useH2) THEN $
      H_surface_den_global_true = TOTAL(g[ginterior].mass*(HI[ginterior] + 2.0*H2[ginterior]))/(!PI*r25^2*1e6) ELSE $
      H_surface_den_global_true = TOTAL(g[ginterior].mass*(HI[ginterior]))/(!PI*r25^2*1e6)
    gas_surface_den_global_true = H_surface_den_global_true*1.4
ENDIF ELSE BEGIN
    ginterior = where(radius_g le r25 AND g.tempg le 1e5)
    gas_surface_den_global_true = TOTAL(g[ginterior].mass)/(!PI*r25*r25*1e6)
    H_surface_den_global_true = gas_surface_den_global_true/1.4
ENDELSE

sinterior = where(radius_s le r25 AND s.tform*timeunit GT tcurrent - dtime)
SFR_surface_den_global_true = N_ELEMENTS(sinterior)*units.istarmass/(!PI*r25*r25*dtime)
stop
print,'R_25: ',r25
print,N_ELEMENTS(where(g.tempg le 1e5)),N_ELEMENTS(ind_gall),N_ELEMENTS(ginterior)
print,'True HI: ', TOTAL(g[ginterior].mass),gas_surface_den_global_true,alog10(gas_surface_den_global_true)
print,'True SFR: ',N_ELEMENTS(sinterior)*units.istarmass,SFR_surface_den_global_true,alog10(SFR_surface_den_global_true)
print,' '
IF KEYWORD_SET(verbose) THEN stop
return,[H_surface_den_global_true,gas_surface_den_global_true,SFR_surface_den_global_true,r25]
end
