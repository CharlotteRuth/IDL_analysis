pro LW_luminosity,wtipsyarr = wtipsyarr,outplot = outplot,readarr = readarr
loadct,39

c = 3d10 ;cm s^-1
sigmaSB = 5.670400d-5  ;erg cm^-2 s^-1 K^-4
nu = 2.99d15  ; s^-1
dnu = 5.88e14 ; s^-1
h = 6.626068d-27 ;cm^2 g s^-1
k = 1.3806503d-16 ;cm^2 g s^-2 K^-1 
L_sol = 3.839d33 ;  erg s^-1
cmperkpc = 3.08568025d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
zsolar = 0.0130215 ;0.0177

dir = '/astro/net/scratch1/abrooks/FABIO/h277.cosmo50cmb.1536g2bwdK/h277.cosmo50cmb.1536g2bwdK.00512/'
file = 'h277.cosmo50cmb.1536g2bwdK.00512'
pfile = '~/REPOSITORY/e12Gals/h277.cosmo50cmb.1536g/h277.cosmo50cmb.1536g.param'
;2.64e-3

dir = '~/REPOSITORY/e12Gals/h285.cosmo50cmb.1536g2bwdK/h285.cosmo50cmb.1536g2bwdK.00512/'
file = 'h285.cosmo50cmb.1536g2bwdK.00512'
pfile = '~/REPOSITORY/e12Gals/h285.cosmo50cmb.1536g2bwdK/h285.cosmo50cmb.1536g2bwdK.param'

dir = '~/REPOSITORY/e12Gals/h258.cosmo50cmb.1536g2bwK/h258.cosmo50cmb.1536g2bwK.00512/'
file = 'h258.cosmo50cmb.1536g2bwK.00512'
pfile = '~/REPOSITORY/e12Gals/h258.cosmo50cmb.1536g2bwK/h258.cosmo50cmb.1536g2bwK.param'

dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/'
file = 'h516.cosmo25cmb.1536g3HBWK.00512.halo'
pfile = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK.param'

dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h277.cosmo50cmb.1536g1MBWKBH/'
file = 'h277.cosmo50cmb.1536g1MBWKBH.mompsource.00001'
pfile = '/astro/net/scratch2/christensen/MolecH/Cosmo/h277.cosmo50cmb.1536g1MBWKBH/h277.cosmo50cmb.1536g1MBWKBH.momssource.param'
;file = 'h277.cosmo50cmb.1536g1MBWKBH.psource2.00001'
;pfile = '/astro/net/scratch2/christensen/MolecH/Cosmo/h277.cosmo50cmb.1536g1MBWKBH/h277.cosmo50cmb.1536g1MBWKBH.psource2.param'

dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g8HBWK.00504/'
file = 'h516.cosmo25cmb.1536g3HBWK.00504.halo.1'
file = 'h516.cosmo25cmb.1536g8HBWK.combgmomcmass3.00504.halo.1'
pfile = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g8HBWK.00504/h516.cosmo25cmb.1536g8HBWK.00504.param'
dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g7HBWK_00480/steps.JeansSF/h516.cosmo25cmb.1536g8HBWK.JeansSF.00408.00032.dir/'
file = 'h516.cosmo25cmb.1536g8HBWK.JeansSF.00408.00032.halo.1'
pfile = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g7HBWK_00480/h516.cosmo25cmb.1536g8HBWK.param'

dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g8HBWK/steps/h516.cosmo25cmb.1536g9HBWK.00408.dir/'
file = 'h516.cosmo25cmb.1536g9HBWK.00408.halo.1'
pfile = '../../h516.cosmo25cmb.1536g9HBWK.param'

dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/'
file = 'h516.cosmo25cmb.2304g14HBWK.00512.halo.1'
pfile = '../../h516.cosmo25cmb.2304g14HBWK.param'

IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps'
    loadct,39
    !P.MULTI = 0
    !P.CHARSIZE = 2             ;3
    !P.THICK=1.5                ;4
    !P.CHARTHICK=1.5            ;4
    !P.font=0
    !X.THICK=1.5                ;4
    !Y.THICK=1.5                ;4
    device,filename = file+'_UVrad.eps',/color,bits_per_pixel= 8,/times,xsize = 6,ysize = 8,/inch
ENDIF ELSE BEGIN
    set_plot,'x'
    loadct,39
    !P.MULTI = 0
    !P.CHARSIZE = 2             ;3
    !P.THICK=1.5                ;4
    !P.CHARTHICK=1.5            ;4
    !P.font=0
    !X.THICK=1.5                ;4
    !Y.THICK=1.5                ;4
    window,0,xsize = 800,ysize = 1200
ENDELSE

units = tipsyunits(pfile,/silent)
dunit = units.lengthunit ;50000.000 ;kpc per sys
munit =  units.massunit  ;1.8479300d16 ;msol per sys
tunit =  units.timeunit ;3.8781500d10; yrs per sys
rtipsy,dir + file,h,g,d,s
mind = 2.50000e-07/dunit
IF ((where(s.x eq 0))[0] ne -1) THEN s[where(s.x eq 0)].x = mind
IF ((where(s.y eq 0))[0] ne -1) THEN s[where(s.y eq 0)].y = mind
IF ((where(s.z eq 0))[0] ne -1) THEN s[where(s.z eq 0)].z = mind
IF keyword_set(readarr) THEN BEGIN
    readarr,dir+file+".lw",h,lw_prime,/ascii
    lw = lw_prime[0:h.ngas - 1]*1d30/cmperkpc^2/dunit^2
ENDIF
IF (units.istarmass eq 0) THEN imass = MAX(s.mass)*munit ELSE imass = units.istarmass ;64021.535; solarmasses
gas_disk = where(SQRT(g.x*g.x + g.y*g.y) < 0.000612841/2.0 AND ABS(g.y) < 3.27241e-05/2.0)

rstar = SQRT(s.x*s.x + s.y*s.y)*dunit
rgas = SQRT(g.x*g.x + g.y*g.y)*dunit
age = (MAX(s.tform) - s.tform)*tunit
LW_lum = LW_lum(age)*imass/2.0
LW_lum_write = dblarr(h.n + 1)
LW_lum_write[0] = h.n
LW_lum_write[h.ngas + h.ndark + 1:h.ngas + h.ndark + h.nstar] = LW_lum
IF KEYWORD_SET(wtipsyarr) THEN BEGIN
    openw,1,dir + file + '.1.LW_lum'
    printf,1,TRANSPOSE(LW_lum_write)
    close,1
ENDIF

!P.MULTI = [0, 1, 3]
X_H = 0.7
z_dust = 0.4
g.dens = g.dens*munit * gm_per_msol * 5.9753790e+23/dunit^3/cmperkpc^3
cross_sec = 2.1d-21
minr = 0
maxr = 8 ;12.5
range = maxr - minr
nbins = 100
r_bin = findgen(nbins)*range/(nbins)
;mass = weighted_HISTOGRAM(rgas,input=g.mass*munit,binsize=range/nbins,min=minr,max=maxr)
dr = range/(nbins)
area = ((r_bin + dr)*(r_bin + dr) - r_bin*r_bin)*1e6
;prof = mass/area
;plot,r_bin,prof
dens = weighted_HISTOGRAM(rgas,input=g.dens,binsize=range/nbins,min=minr,max=maxr)
metal = weighted_HISTOGRAM(rgas,input=g.zmetal,binsize=range/nbins,min=minr,max=maxr)
;dens_metal = weighted_HISTOGRAM(rgas,input=g.dens*g.zmetal*X_H*z_dust*cross_sec,binsize=range/nbins,min=minr,max=maxr)
dens_metal = weighted_HISTOGRAM(rgas,input=g.dens*X_H*g.zmetal/zsolar*cross_sec,binsize=range/nbins,min=minr,max=maxr)
number = weighted_HISTOGRAM(rgas,input=1.0,binsize=range/nbins,min=minr,max=maxr)
plot,r_bin,dens/number,/ylog,ytitle = 'density [cm^-3]',xtitle = 'radius [kpc]'
plot,r_bin,metal/number/zsolar,/ylog,ytitle = 'metallicity [z/zsolar]',xtitle = 'radius [kpc]'
plot,r_bin,dens_metal/number,/ylog,ytitle = 'alpha_LW',xtitle = 'radius [kpc]'
;stop

maxrad = maxr/dunit;0.00016
nel_rad = 10
radii = findgen(nel_rad)/(nel_rad - 1)*maxrad
Lw_flux_rad = radii
LW_flux_rad_alpha = radii
dr = maxrad*dunit/(nel_rad - 1)
colors = findgen(nel_rad)/(nel_rad - 1)*(240 - 40) + 20
alpha = 0.4
dens_metal = weighted_HISTOGRAM(rgas,input=g.dens*g.zmetal*X_H*z_dust*cross_sec,binsize=range/nel_rad,min=0,max=maxr,location = location)
number = weighted_HISTOGRAM(rgas,input=1.0,binsize=range/nel_rad,min=0,max=maxr)
alpha = dens_metal/number

IF KEYWORD_SET(outplot) then begin
    set_plot,'ps'
    device,filename = file+'_UVrad.eps',/color,bits_per_pixel= 8,/times,xsize = 6,ysize = 8,/inch
endif else window,1,xsize = 800,ysize = 1200;window,0
!P.MULTI = [0, 1, 3]
FOR i = 0, nel_rad - 1 DO BEGIN
    distance = dunit*SQRT((s.x - radii[i])*(s.x - radii[i]) + s.y*s.y + s.z*s.z)
    LW_flux = LW_lum/distance/distance/cmperkpc/cmperkpc/4.0/!PI
    y = weighted_histogram(distance,input = LW_flux,locations = locations,min = 0,max = 2.0*maxrad*dunit)
    IF (i eq 0) THEN plot,locations,y,/ylog,xtitle = 'Distance [kpc]',ytitle = 'F[d]',xrange = [0,2.0*maxrad*dunit],xstyle = 1,yrange = [1e-5,1e15]
    oplot,locations,y,color = colors[i]
    flux_weighted_r = TOTAL(distance[where(finite(Lw_flux))]*LW_flux[where(finite(Lw_flux))])/TOTAL(LW_flux[where(finite(Lw_flux))])
    oplot,[flux_weighted_r,flux_weighted_r],[1e-10,1e20],linestyle = 1, color = colors[i]
    print,TOTAL(LW_flux[where(finite(Lw_flux))])/3e10*0.6,' photons cm^-3 at r = ',STRTRIM(radii[i]*dunit,2),' kpc'
    lw_flux_rad[i] = TOTAL(LW_flux[where(finite(Lw_flux))]);/3e10*0.6
ENDFOR
linestyle = fltarr(nel_rad)
names = 'R = ' + STRTRIM(radii*dunit,2) + ' kpc'
legend,names,linestyle = linestyle,color = colors,CHARSIZE = 1,/right
print,''

IF KEYWORD_SET(outplot) then begin
    device,/close
    device,filename = file+'_UVrad_sigma.eps',/color,bits_per_pixel= 8,/times
endif ;else window,1
FOR i = 0, nel_rad - 1 DO BEGIN
;    alpha = MEAN(g[where(ABS(rgas - radii[i]*dunit) lt dr/2.0)].dens*g[where(ABS(rgas - radii[i]*dunit) lt dr/2.0)].zmetal*X_H*z_dust*cross_sec)*cmperkpc
    alpha = MEAN(g[where(ABS(rgas - radii[i]*dunit) lt dr/2.0)].dens*g[where(ABS(rgas - radii[i]*dunit) lt dr/2.0)].zmetal/zsolar*X_H*cross_sec)*cmperkpc
;    print,alpha
;    stop
  ;  alpha = 100.0*alpha
    distance = dunit*SQRT((s.x - radii[i])*(s.x - radii[i]) + s.y*s.y + s.z*s.z)
    LW_flux = LW_lum/distance/distance/cmperkpc/cmperkpc*exp(-1.0*alpha*distance)/4.0/!PI
    y = weighted_histogram(distance,input = LW_flux,locations = locations,min = 0,max = 2.0*maxrad*dunit)
    IF (i eq 0) THEN plot,locations,y,/ylog, xtitle = 'Distance [kpc]',ytitle = 'F[r]',xrange = [0,2.0*maxrad*dunit],xstyle = 1,yrange = [1e-5,1e15]
    oplot,locations,y,color = colors[i]
    flux_weighted_r = TOTAL(distance[where(finite(Lw_flux))]*LW_flux[where(finite(Lw_flux))])/TOTAL(LW_flux[where(finite(Lw_flux))])
    oplot,[flux_weighted_r,flux_weighted_r],[1e-10,1e20],linestyle = 1, color = colors[i]
    print,TOTAL(LW_flux[where(finite(Lw_flux))])/3e10*0.6,' photons cm^-3 at r = ',STRTRIM(radii[i]*dunit,2),' kpc'
    lw_flux_rad_alpha[i] = TOTAL(LW_flux[where(finite(Lw_flux))]);/3e10*0.6
ENDFOR
legend,names,linestyle = linestyle,color = colors,CHARSIZE = 1,/right

inrad = 3.0  ;6.18499e-05
outrad = 10.0  ;0.000131241
height = 0.06
box_age = where( rstar gt inrad AND rstar lt outrad AND ABS(s.z) le height AND age lt 1e8)
box = where( rstar gt inrad AND rstar lt outrad AND ABS(s.z) le height AND age lt 2e9)
volume = !PI*(outrad^2 - inrad^2) * height*dunit^3*(cmperkpc)^3
leg = volume^(1/3.0)/2.0

plot,radii*dunit,LW_flux_rad_alpha,xtitle = 'radius [kpc]',ytitle = 'LW photons cm^-3',/ylog,yrange = [1e1,1e9],ystyle = 1
oplot,radii*dunit,LW_flux_rad_alpha,color = 120
oplot,radii*dunit,LW_flux_rad,color = 240
IF keyword_set(readarr) then oplot,rgas[where(ABS(g.z) lt 5e-6 AND g.dens gt 1)],lw[where(ABS(g.z) lt 5e-6 AND g.dens gt 1)],psym = 3
oplot,radii*dunit,total(LW_lum)/radii/dunit/radii/dunit/cmperkpc/cmperkpc,color = 50

oplot,rgas,lw,psym = 3
;file = 'h516.cosmo25cmb.1536g3HBWK.00504.halo.1'
readarr,dir+file+".lw",h,lw_prime,/ascii
lw = lw_prime[0:h.ngas - 1]*1d30/cmperkpc^2/dunit^2
oplot,rgas,lw,psym = 3,color = 30


IF KEYWORD_SET(outplot) then begin
    device,/close
    set_plot,'x'
ENDIF
!P.MULTI = 0
window,4,xsize = 800,ysize = 392

s_SD = weighted_histogram(rstar,input = s.mass*munit,locations = locations,min = 0,max = 2.0*maxrad*dunit,nbins = 100)
s_SD_num = histogram(rstar,locations = locations,min = 0,max = 2.0*maxrad*dunit,nbins = 100)
g_SD = weighted_histogram(rgas,input = g.mass*munit,locations = locations,min = 0,max = 2.0*maxrad*dunit,nbins = 100)
area = locations*2.0*maxrad*dunit/100*2*!PI
s_SD = s_SD/area
g_SD = s_SD/area
plot,locations,g_SD,/ylog,xrange = [0.1,10],xtitle = 'Radius [kpc]',ytitle = 'Stellar/Gas Surface Density [Msol kpc^2]'
oplot,locations,s_SD,color = 240
legend,['Gas','Stars'],color = [0,240],linestyle = [0,0],/right

fitpar =  [0,0,0,0,0]
rmin = 0.1
rmax = 10
ind = where(locations gt rmin and locations lt rmax)
;dblexpfit, locations[ind], s_SD[ind], s_SD[ind]/sqrt(s_SD_num[ind]), fitpar, red_chisq=exp2chi, /no_guess
;stop
;Make a guess at the exponential profile
A=[5.0,1d11]
expchi=1e11
weights = FLTARR(N_ELEMENTS(ind))+1.0
expFit=CURVEFIT(locations[ind],s_SD[ind],weights,A,sigma,FUNCTION_NAME='expdisk',CHISQ=expchi)
expchi=expchi/(n_elements(ind) - n_elements(A) - 1.)
print,'EXP: ',A,expchi
;Make a guess at the de Vaucouleurs profile
A=[1.0,1.0]
dVchi=1e11
weights = FLTARR(N_ELEMENTS(ind))+1.0
deVaucFit=CURVEFIT(locations[ind],s_SD[ind],weights,A,sigma,FUNCTION_NAME='deVauc',CHISQ=dVchi)
dVchi=dVchi/(n_elements(ind) - n_elements(A) - 1.)
print,'DeVauc: ',A,dVchi

if(dVchi lt expchi AND dVchi lt exp2chi) then best = '   deVauc' $
else if (exp2chi lt expchi AND exp2chi lt dVchi) then best = '   2 Exp'$
else best='   1 Exp'

;r1 = findgen(100)/100.*(fitpar[0])+3.
;r2 = findgen(100)/100.*(rmax-fitpar[0]+3.) + fitpar[0]-3.
;oplot, r1, fitpar[1]*exp(-r1/fitpar[2]), line = 1
;oplot, r2, fitpar[3]*exp(-r2/fitpar[4]), line = 1
;oplot, [fitpar[0],fitpar[0]], [1e4,1e11], line = 3

;stop

;sbox = s[box]
;nphots = TOTAL(LW_lum[box])
;print,''
;print,strtrim(total(nphots),2),' photons s^-1; '
;print,strtrim(volume,2),' cm^-3; '
;print,strtrim(total(nphots)/volume*leg/c,2),'  photons cm^-3'
;print,strtrim(total(nphots)/volume*leg/c*0.6,2),'  photons cm^-3'

;window,3
;!P.MULTI = 0
;plot,s.x*dunit,s.y*dunit,psym = 3,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',xrange = [-8,8],yrange = [-6,6]
;oplot,s[where(age lt 1e8)].x*dunit,s[where(age lt 1e8)].y*dunit,psym = 2,color = 30
;oplot,s[box_age].x*dunit,s[box_age].y*dunit,psym = 2,color = 100

;window,2
;plot,s.x*dunit,s.z*dunit,psym = 3,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',xrange = [-8,8],yrange = [-6,6]
;oplot,s[where(age lt 1e8)].x*dunit,s[where(age lt 1e8)].z*dunit,psym = 2,color = 30
;oplot,s[box_age].x*dunit,s[box_age].z*dunit,psym = 2,color = 100
;stop

end

function LW_lum,age1
age = age1
a0 = -332.61118
a1 = 210.37616
a2 = -43.684908
a3 = 4.0379845
a4 = -0.14145357
age[where(age lt 1e6)] = 1e6
lum = 10d0^(a0 + $
           a1*alog10(age) + $
           a2*alog10(age)*alog10(age) + $ 
           a3*alog10(age)*alog10(age)*alog10(age) + $
           a4*alog10(age)*alog10(age)*alog10(age)*alog10(age))
return,lum
end

PRO deVauc,x,A,F,pder
;x is the radius, 
;A[0]= effective radius, 
;A[1] = surface brightness at the effective radius
F = A[1]*EXP(-7.67*((X/A[0])^0.25 - 1))
d1 = EXP(-7.67*((X/A[0])^0.25 - 1))
d2 = A[1]*EXP(-7.67*((X/A[0])^0.25 - 1))*(-1.9175)*x^0.25*A[0]^(-1.25)
pder =[[d1],[d2]]
END


PRO expdisk,x,A,F,pder
;x is the radius, 
;A[0]= effective radius, 
;A[1] = surface brightness at the effective radius
F = A[1]*EXP((-1.0)*X/A[0])
d1 = EXP((-1.0)*X/A[0])
d2 = A[1]*EXP(x/A[0])*x*A[0]^(-2.0)
pder =[[d1],[d2]]
END
