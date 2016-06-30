pro envelope_LW,outplot = outplot
;See starburst99_Lw for the calculations from starburst

IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps'
    loadct,39
    device,filename = '/astro/users/christensen/code/MolecH/LWphotonden_calc.eps',/color,bits_per_pixel = 8,/times,xsize = 6,ysize = 8,/inch
    !P.MULTI = [0, 2, 3]
    !P.CHARSIZE = 2             ;3
    !P.THICK=1.5                ;4
    !P.CHARTHICK=1.5            ;4
    !P.font=0
    !X.THICK=1.5                ;4
    !Y.THICK=1.5                ;4
ENDIF ELSE BEGIN
    set_plot,'x'
    loadct,39
    !P.MULTI = [0, 2, 3]
    !P.CHARSIZE = 2             ;3
    !P.THICK=1.5                ;4
    !P.CHARTHICK=1.5            ;4
    !P.font=0
    !X.THICK=1.5                ;4
    !Y.THICK=1.5                ;4
    window,0,xsize = 800,ysize = 1200
ENDELSE

c = 3d10 ;cm s^-1
sigmaSB = 5.670400d-5  ;erg cm^-2 s^-1 K^-4
nu = 2.99d15  ; s^-1
dnu = 5.88e14 ; s^-1
h = 6.626068d-27 ;cm^2 g s^-1
k = 1.3806503d-16 ;cm^2 g s^-2 K^-1 
L_sol = 3.839d33 ;  erg s^-1

pfile = '/astro/net/scratch1/abrooks/FABIO/h277.cosmo50cmb.1536g2bwK/h277.cosmo50cmb.1536g2bwK.param'
pfile = '/astro/net/scratch2/christensen/MolecH/11M/Cloud_1e5/11M_1.param'
unit = tipsyunits(pfile)
dunit = 50000.000 ;kpc per sys
munit =  1.8479300d16 ;msol per sys
tunit =  3.8781500d10; yrs per sys
imass = 64021.535; solarmasses
cmperkpc = 3.08568025d21

stellar_types = [30000,10000,7500,6000,5200,3700]
mass_types = M_from_T(stellar_types)

nel = 10000
Tmin = 1000
Tmax = 100000
T = findgen(nel)/nel*(Tmax - Tmin) + Tmin

;indow,0
;plot,T,I_lw,xtitle = 'T [K]',ytitle = 'I_LW [ergs Hz^-1 cm^-2 s^-1]',/xlog,/ylog,xrange = [6000,100000],xstyle = 1
;oplot,[stellar_types[0],stellar_types[0]],[1e-20,1e5]
;oplot,[stellar_types[1],stellar_types[1]],[1e-20,1e5]
;oplot,[stellar_types[2],stellar_types[2]],[1e-20,1e5]
;oplot,[stellar_types[3],stellar_types[3]],[1e-20,1e5]
;oplot,[stellar_types[4],stellar_types[4]],[1e-20,1e5]
;oplot,[stellar_types[5],stellar_types[5]],[1e-20,1e5]

;window,1
;plot,T,L_from_T(T),/ylog,/xlog,xtitle = 'T [K]',xrange = [6000,100000],ytitle = 'Mass [M_sol] / Luminosity [L_sol]',yrange = [1,1e10],xstyle = 1
;oplot,T,M_from_T(T),linestyle = 1
;stelar_types = alog10(stellar_types)
;oplot,[stellar_types[0],stellar_types[0]],[1e-20,1e10]
;oplot,[stellar_types[1],stellar_types[1]],[1e-20,1e10]
;oplot,[stellar_types[2],stellar_types[2]],[1e-20,1e10]
;oplot,[stellar_types[3],stellar_types[3]],[1e-20,1e10]
;oplot,[stellar_types[4],stellar_types[4]],[1e-20,1e10]
;oplot,[stellar_types[5],stellar_types[5]],[1e-20,1e10]

;************************* Checked Against Wikipedia Stellar Classifications -- Looks Reasonable **************
;window,1
plot,M_from_T(T),L_from_T(T),/ylog,/xlog,ytitle = 'Luminosity [L_sol]',yrange = [1,1e7],xstyle = 1,xtitle = 'Mass [M_sol]',xrange = [1,100],ystyle = 1
;oplot,M_from_T(T),T,linestyle = 1
stelar_types = alog10(mass_types)
oplot,[mass_types[0],mass_types[0]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[1],mass_types[1]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[2],mass_types[2]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[3],mass_types[3]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[4],mass_types[4]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[5],mass_types[5]],[1e-20,1e10],linestyle = 2


;************************* Checked Against C&O Wien's Law Derivation, p 467 -- Looks Reasonable ************
I_lw = BB(T)
I_lwp = BBP(T)
;window,0
;plot,M_from_T(T),I_lw,ytitle = 'I_LW [ergs Hz^-1 cm^-2 s^-1]',/xlog,/ylog,xtitle = 'Mass [M_sol]',xstyle = 1,xrange = [1,100]
plot,M_from_T(T),I_lwp,ytitle = 'I_LW [LW photons s^-1]',/xlog,/ylog,xtitle = 'Mass [M_sol]',xstyle = 1,xrange = [1,100],yrange=[1d37,1d50],ystyle = 1
oplot,[mass_types[0],mass_types[0]],[1e-20,1d50],linestyle = 2
oplot,[mass_types[1],mass_types[1]],[1e-20,1d50],linestyle = 2
oplot,[mass_types[2],mass_types[2]],[1e-20,1d50],linestyle = 2
oplot,[mass_types[3],mass_types[3]],[1e-20,1d50],linestyle = 2
oplot,[mass_types[4],mass_types[4]],[1e-20,1d50],linestyle = 2
oplot,[mass_types[5],mass_types[5]],[1e-20,1d50],linestyle = 2

;************************* Checked Mass Normalization and Compared Shape to Plots -- Looks Reasonable ******************
;window,3
plot,M_from_T(T),IMF(M_from_T(T),/nstar),/xlog,/ylog,xtitle = 'Mass [M_sol]',ytitle = 'dN/dM',xrange = [1,100],xstyle = 1
oplot,M_from_T(T),IMF(M_from_T(T),/ms,/nstar),linestyle = 1
legend,['Kroupa','Miller-Scalo'],linestyle = [0,1],/top,/right,charsize = 0.6
oplot,[mass_types[0],mass_types[0]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[1],mass_types[1]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[2],mass_types[2]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[3],mass_types[3]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[4],mass_types[4]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[5],mass_types[5]],[1e-20,1e10],linestyle = 2
;mass_imfcheck = 10^(findgen(1000)/1000*5 - 3)
;m_bins = mass_imfcheck[1:999] - mass_imfcheck[0:998]
;imfcheck = (IMF(mass_imfcheck[1:999]) + IMF(mass_imfcheck[0:998]))/2.0*m_bins
;print,TOTAL((IMF(mass_imfcheck[1:999]) + IMF(mass_imfcheck[0:998]))/2.0*m_bins)

;************************** Sanity Check.  Looks Reasonable ************************************
MS_age = MS_age(M_from_T(T),fit = M_age_fit)
;window,2
;mass = M_from_T(T)
;plot,mass,10^(M_age_fit[0] + M_age_fit[1]*alog10(mass) + M_age_fit[2]*alog10(mass)*alog10(mass)),linestyle = 1,/xlog,/ylog
plot,M_from_T(T),MS_age(M_from_T(T)),/xlog,/ylog,xtitle = 'Mass [M_sol]',xrange = [1,100],ytitle = 'MS Lifetime [yr]',yrange = [1e6,1e10],ystyle = 1
oplot,[mass_types[0],mass_types[0]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[1],mass_types[1]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[2],mass_types[2]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[3],mass_types[3]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[4],mass_types[4]],[1e-20,1e10],linestyle = 2
oplot,[mass_types[5],mass_types[5]],[1e-20,1e10],linestyle = 2

;window,4
;plot,M_from_T(T),IMF(M_from_T(T),/nstar)*BB(T),/xlog,/ylog,xtitle = 'Mass [M_sol]',xrange = [1,100],xstyle = 1,ytitle = 'I_LW per Stellar Mass [ergs Hz^-1 cm^-2 s^-1 M_sol^-1]',yrange = [1e-14,1e-6]
;oplot,M_from_T(T),IMF(M_from_T(T),/ms,/nstar)*BB(T),linestyle = 1
plot,M_from_T(T),IMF(M_from_T(T),/nstar)*I_lwp,/ylog,xtitle = 'Mass [M_sol]',xrange = [1,100],xstyle = 1,ytitle = 'I_LW per SSP [photons s^-1 Msol^-1]',yrange = [1d37,1d45],/xlog;,yrange = [1d36,1d46]
oplot,M_from_T(T),IMF(M_from_T(T),/ms,/nstar)*I_lwp,linestyle = 1
;oplot,M_from_T(T),L_from_T(T),linestyle = 1
;oplot,M_from_T(T),IMF(M_from_T(T),/nstar)*1000,linestyle = 1
oplot,[mass_types[0],mass_types[0]],[1e-20,1d56],linestyle = 2
oplot,[mass_types[1],mass_types[1]],[1e-20,1d56],linestyle = 2
oplot,[mass_types[2],mass_types[2]],[1e-20,1d56],linestyle = 2
oplot,[mass_types[3],mass_types[3]],[1e-20,1d56],linestyle = 2
oplot,[mass_types[4],mass_types[4]],[1e-20,1d56],linestyle = 2
oplot,[mass_types[5],mass_types[5]],[1e-20,1d56],linestyle = 2

min_log_age = 6.
max_log_age = 9.6
nel = 500.
age_bins = 10^((findgen(nel) + 1)/(nel)*(max_log_age - min_log_age) + min_log_age)
;age_bins = [10^min_log_age,age_bins]
age_rad = fltarr(nel)
age_rad_ms = fltarr(nel)
min_mass_ind = fltarr(nel)
mass_lower = 1

FOR i = 0, nel - 1 DO BEGIN 
;    temp = MIN(ABS(age_bins[i]     - MS_age(M_from_T(T))),mass_upper_t)
    temp = MIN(ABS(age_bins[i] - MS_age(M_from_T(T))),ind) 
    min_mass_ind[i] = ind
    mass_ind = WHERE(0.5 lt M_from_T(T) AND MS_age(M_from_T(T)) ge age_bins[i])
    dm = ABS((M_from_T(T[mass_ind]) - M_from_T(T[mass_ind - 1])))
    age_rad[i]    = TOTAL((IMF(M_from_T(T[mass_ind])    ,/nstar    )*BBP(T[mass_ind]) + $
                           IMF(M_from_T(T[mass_ind + 1]),/nstar    )*BBP(T[mass_ind + 1]))/2*dm/1d30  )
    age_rad_ms[i] = TOTAL((IMF(M_from_T(T[mass_ind])    ,/nstar,/ms)*BBP(T[mass_ind]) + $
                           IMF(M_from_T(T[mass_ind + 1]),/nstar,/ms)*BBP(T[mass_ind + 1]))/2*dm/1d30  )
ENDFOR

;oplot,M_from_T(T[min_mass_ind]),IMF(M_from_T(T[min_mass_ind]),/nstar)*I_lwp[min_mass_ind]*imass,psym = 1,color = 240


;window,5
plot,age_bins,age_rad*1d30,xtitle = 'Time [yr]',/xlog,ytitle = 'I_LW per Stellar Mass [photons s^-1 M_sol^-1]',/ylog,yrange = [1e46,1e51]
;oplot,age_bins,age_rad*1d30,psym = 1,color = 240
oplot,age_bins,age_rad_ms*1d30,linestyle = 1

file = '/astro/net/scratch1/abrooks/FABIO/h277.cosmo50cmb.1536g2bwdK/h277.cosmo50cmb.1536g2bwdK.00512/h277.cosmo50cmb.1536g2bwdK.00512.1.std'
file = '/astro/net/scratch2/christensen/MolecH/11M/Cloud_1e5/o11M_00300.00001'
rtipsy,file,h,g,d,s
inrad = 6.18499e-05
outrad = 8.82041 ;0.000131241
height = 2.61542 ;2.27223e-06/2.

rstar = SQRT(s.x*s.x + s.y*s.y)
age = (MAX(s.tform) - s.tform)*unit.timeunit ;3.8781500d10
age[where(age lt 1e6)] = 1e6
box_age = where( rstar gt inrad AND rstar lt outrad AND ABS(s.z) le height AND age lt 1e8)
box = where( rstar gt inrad AND rstar lt outrad AND ABS(s.z) le height AND age lt 2e9)
volume = !PI*(outrad^2 - inrad^2) * height*dunit^3*(cmperkpc)^3
leg = volume^(1/3.0)/2.0
stop

;smassden = N_ELEMENTS(box)*imass/volume
;photon_den = smassden*1d50
;print,photon_den,' LW photons cm^-3 s^-1'

sbox = s[box]
sboxage = age[box]
s_box_age = sboxage[sort(sboxage)]
s_box_age[where(s_box_age lt 1e6)] = 1e6
nphots = spline(age_bins,age_rad,s_box_age)
nphots_ms = spline(age_bins,age_rad_ms,s_box_age) 
;oplot,s_box_age,nphots*1d30,psym = 2,color = 60
print,strtrim(total(nphots*1d30*imass),2),' (',strtrim(total(nphots_ms*1d30*imass),2),') photons s^-1; '
print,strtrim(volume,2),' cm^-3; '
print,strtrim(total(nphots*1d30*imass)/volume*leg/c,2),' (',strtrim(total(nphots_ms*1d30*imass*leg/c)/volume,2),')  photons cm^-3'
print,strtrim(total(nphots*1d30*imass)/volume*leg/c*.6,2),' (',strtrim(total(nphots_ms*1d30*imass*leg/c)/volume*.6,2),')  photons cm^-3'

fit = POLY_FIT(alog10(age_bins[where(age_bins lt 1e9)]),alog10(age_rad[where(age_bins lt 1e9)]*1d30),4)
age_fit = 10^(findgen(100)*(9 - 6)/100.0 + 6.)
age_rad_fit = 10^(fit[0] + $
              fit[1]*alog10(age_fit) + $
              fit[2]*alog10(age_fit)*alog10(age_fit) + $ 
              fit[3]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit) + $
              fit[4]*alog10(age_fit)*alog10(age_fit)*alog10(age_fit)*alog10(age_fit))
oplot,age_fit,age_rad_fit,linestyle = 2,color = 240

print,'Fit:'
print,'log10(L_LW) = ',strtrim(fit[0],2),' + ',strtrim(fit[1],2),'*log10(t) + ',strtrim(fit[2],2),'*log10(t)^2 + ',strtrim(fit[3],2),'*log10(t)^3 + ',strtrim(fit[4],2),'*log10(t)^4'
!P.MULTI = 0
window,1
lum = 10^(fit[0] + $
              fit[1]*alog10(age) + $
              fit[2]*alog10(age)*alog10(age) + $ 
              fit[3]*alog10(age)*alog10(age)*alog10(age) + $
              fit[4]*alog10(age)*alog10(age)*alog10(age)*alog10(age))*unit.istarmass
readcol,file+'.lw',lum_code
lum_code_star = lum_code[N_ELEMENTS(lum_code) - h.nstar:N_ELEMENTS(lum_code)-1]*1d30
plot,alog10(lum),alog10(lum_code_star),psym = 3
stop

IF KEYWORD_SET(outplot) THEN device,/close
set_plot,'x'

window,3
!P.MULTI = 0
plot,s.x*dunit,s.y*dunit,psym = 3,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',xrange = [-8,8],yrange = [-6,6]
oplot,s[where(age lt 1e8)].x*dunit,s[where(age lt 1e8)].y*dunit,psym = 2,color = 30
oplot,s[box_age].x*dunit,s[box_age].y*dunit,psym = 2,color = 100

window,2
plot,s.x*dunit,s.z*dunit,psym = 3,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',xrange = [-8,8],yrange = [-6,6]
oplot,s[where(age lt 1e8)].x*dunit,s[where(age lt 1e8)].z*dunit,psym = 2,color = 30
oplot,s[box_age].x*dunit,s[box_age].z*dunit,psym = 2,color = 100


stop


end

function MS_age,massMS,fit = fit
;Data from Carrol and Ostlie, p 486, ref Iben 1967
Mass = [15,9,5,3,2.25,1.5,1.25,1]
age1 = [1.010e7,2.144e7,6.547e7,2.212e8,4.802e8,1.553e9,2.803e9,7e9]
age2 = [2.270e5,6.053e5,2.173e6,1.042e7,1.647e7,8.10e7,1.824e8,2e9]
age = age1 + age2

;plot,mass,age,/xlog,/ylog,xtitle = 'Mass [M_sol]',ytitle = 'Age [Yr]'
fit = poly_fit(alog10(mass),alog10(age),2)
;oplot,mass,10^(fit[0] + fit[1]*alog10(mass) +fit[2]*alog10(mass)*alog10(mass)),linestyle = 1

lower_ind = WHERE(massMS lt MAX(Mass,max_ind),complement = upper_ind)
ageMS = massMS
ageMS[lower_ind] = 10^(fit[0] + fit[1]*alog10(massMS[lower_ind]) + fit[2]*alog10(massMS[lower_ind])^2)

x0 = alog10(mass[max_ind])
y0 = alog10(age[max_ind])  ;(fit[0] + fit[1]*x0 + fit[2]*x0^2)
slope = fit[1] + 2.0*fit[2]*x0
ageMS[upper_ind] = 10^(y0 + slope*(alog10(massMS[upper_ind]) - x0))
return,ageMS
end


function BB,T
I = 0.4/(exp(143497/T))
return,I
end

function BBP,T
c = 3d10 ;cm s^-1
sigmaSB = 5.670400d-5  ;erg cm^-2 s^-1 K^-4
nu = 2.99d15  ; s^-1
dnu = 5.88e14 ; s^-1
h = 6.626068d-27 ;cm^2 g s^-1
k = 1.3806503d-16 ;cm^2 g s^-2 K^-1 
L_sol = 3.839d33 ;  erg s^-1


BBconst = 2*!PI*(nu^2/c^2)*dnu/sigmaSB
BBexp = h*nu/k
Ip = BBconst/(exp(BBexp/T) - 1)*L_from_T(T)*L_sol/T^4
return,Ip
end

function M_from_T,T
;Balona, 1994
log_Teff = alog10(T)
log_L = log_Teff
log_m = log_Teff

highT = WHERE(log_Teff ge 3.82)
lowT =  WHERE(log_Teff lt 3.82 AND log_Teff ge 3.75)
vlowT = WHERE(log_Teff lt 3.75)

IF (highT[0] ne -1) THEN log_L[highT] = 3.4261    - 6.1589*log_Teff[highT]  + 1.4149*(log_Teff[highT])^2
IF (lowT[0] ne -1) THEN log_L[lowT]  = -239.9504 + 118.5614*log_Teff[lowT] - 14.5631*(log_Teff[lowT])^2

IF (highT[0] ne -1) THEN log_m[highT] = 31.6409 + 2.6152*log_L[highT] - 17.3446*log_Teff[highT] $
      + 0.08195*log_L[highT]^2 - 0.6982*log_L[highT]*log_Teff[highT] $
      + 2.3879*log_Teff[highT]^2 

IF (lowT[0] ne -1) THEN log_m[lowT]  = 21.3153 + 5.6698*log_L[lowT]  - 12.9151*log_Teff[lowT] $
      + 0.1233*log_L[lowT]^2   - 1.4917*log_L[lowT]*log_Teff[lowT] $
      + 1.9283*log_Teff[lowT]^2
IF (vlowT[0] ne -1) THEN log_m[vlowT] = 0

return,10.0^log_m
end


function L_from_T,T
;Balona, 1994
log_Teff = alog10(T)
log_L = log_Teff

highT = WHERE(log_Teff ge 3.82)
lowT =  WHERE(log_Teff lt 3.82 AND log_Teff ge 3.75)
vlowT = WHERE(log_Teff lt 3.75)

IF (highT[0] ne -1) THEN log_L[highT] = 3.4261    - 6.1589*log_Teff[highT]  + 1.4149*(log_Teff[highT])^2
IF (lowT[0] ne -1) THEN log_L[lowT]  = -239.9504 + 118.5614*log_Teff[lowT] - 14.5631*(log_Teff[lowT])^2
IF (vlowT[0] ne -1) THEN log_L[vlowT] = 0

return,10.0^log_L
end

function A_from_T,T
;Balona, 1994
log_Teff = alog10(T)
log_L = log_Teff
a = log_Teff

highT = WHERE(log_Teff ge 3.82 AND log_Teff le 4.64)
lowT =  WHERE(log_Teff lt 3.82 AND log_Teff ge 3.75)
vlowT = WHERE(log_Teff lt 3.75)
vhighT = WHERE(log_Teff gt 4.64)

log_L[highT] = 3.4261    - 6.1589*log_Teff[highT]  + 1.4149*(log_Teff[highT])^2
log_L[lowT]  = -239.9504 + 118.5614*log_Teff[lowT] - 14.5631*(log_Teff[lowT])^2
IF (vlowT[0] ne -1) THEN log_L[vlowT] = 0
log_L[vhighT] = 0

a[highT] = 10^(4.0013 - 0.6235*log_L[highT])*(21.5426 + 0.9811*log_L[highT] - 5.7133*log_Teff[highT])
a[lowT]  = 10^(4.1339 - 0.8206*log_L[lowT] )*(38.7918 + 1.2876*log_L[lowT]  - 10.2467*log_Teff[lowT])
a[vlowT] = 0
a[vhighT] = 0

return,10.0^log_L
end


function IMF, M, MS = MS,nstar = nstar
if KEYWORD_SET(MS) THEN imf = ms(M,nstar = nstar) ELSE imf = kroupa(M,nstar = nstar)
return,imf
end

FUNCTION kroupa,smass,nstar = nstar
;For stars with masses greater than 1 solar mass (parameters from A&A,
;315, 1996) -- normalized so that one solar mass of stars are formed
ind1 = WHERE(smass gt 1.0)
ind2 = WHERE(smass gt 0.5 AND smass le 1.0)

a = smass*0 + 0.3029*1.86606 ;for masses less than 0.08
b = smass*0 - 0.3
IF (ind2[0] ne -1) THEN BEGIN
    a[ind2] = 0.3029
    b[ind2] = -1.2
ENDIF
IF (ind1[0] ne -1) THEN BEGIN
    a[ind1] = 0.3029
    b[ind1] = -1.7
ENDIF
IF(KEYWORD_SET(nstar)) THEN b = b-1 ;so that it is no longer in terms of mass but in number of stars
nstars = a * smass^b
RETURN,nstars
END

FUNCTION ms,smass,nstar = nstar
;For stars with masses greater than 1 solar mass (parameters from A&A,
;315, 1996) -- normalized so that one solar mass of stars are formed
;Miller-Scalo IMF
ind1 = WHERE(smass gt 10.0)
ind2 = WHERE(smass gt 1.0 AND smass le 10.0)
c = 118.643 ;It sure looks like MS is scalled to 75 solar masses rather than 1
;c = 1.0

a = smass*0 + 42.0 ;for masses less than 0.08
b = smass*0 - 0.4
IF (ind2[0] ne -1) THEN BEGIN
    a[ind2] = 42.0
    b[ind2] = -1.5
ENDIF
IF (ind1[0] ne -1) THEN BEGIN
    a[ind1] = 240.0
    b[ind1] = -2.3
ENDIF
IF(KEYWORD_SET(nstar)) THEN b = b-1 ;so that it is no longer in terms of mass but in number of stars
nstars = a * smass^b
RETURN,nstars/c
END
