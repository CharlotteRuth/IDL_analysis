pro metaldistro
loadct,39
set_plot,'x'

file3 = '~/Scratch2/Sunrise/MW1lr_v3/mcrx.fits'
cam0_aux = 21
cam0_aux_data = mrdfits(file3,cam0_aux,h)
metallicity = cam0_aux_data[*,*,1]/cam0_aux_data[*,*,0]
temp = where(Finite(metallicity),complement = nan)
metallicity[nan] = 0

;peak = MAX(metallicity)
;theselevels = [peak/500.,peak/200.,peak/100.,peak/50.,peak/20.,peak/10.,peak/5.,peak/2.,peak]
;theselevels = [peak*0.001,peak*0.01,peak*0.1,peak*0.5,peak*0.7,peak*0.8,0.9*peak,peak]
;x=FINDGEN((SIZE(metallicity))[1])
;y=FINDGEN((SIZE(metallicity))[2])
;contour,metallicity,x,y,levels=theselevels,/fill
;,xstyle=1,ystyle=1,xthick=2,ythick=2,charsize=2.5,xmargin=[8,3],ymargin=[3,2],charthick=3
;stop
peak = 0.015
metallicity = metallicity[200:400,200:400]
x=FINDGEN((SIZE(metallicity))[1])+200
y=FINDGEN((SIZE(metallicity))[2])+200
theselevels = [peak*0.1,peak*0.15,peak*0.2,peak*0.25,peak*0.3,peak*0.35,peak*0.4,peak*0.45,peak*0.5,peak*0.55,peak*0.6,peak*0.65,peak*0.7,peak*0.75,peak*0.8,peak*0.85,0.9*peak,0.95*peak,peak]
window,0
contour,metallicity,x,y,levels=theselevels,/fill
;contour,metallicity,x,y,levels=theselevels
set_plot,'ps'
device,filename = 'metaldistro_MW1lr_v3.eps',/color,bits_per_pixel = 8,/times
contour,metallicity,x,y,levels = theselevels,/fill
device,/close
set_plot,'x'

rtipsy,'~/Scratch2/Sunrise/MW1lr_v3/MW1lr_001.std',h,g,d,s
unitlength_in_cm = 8.81609e25
unitlength_in_kpc = unitlength_in_cm/3.08568025e21
g.y = g.y*unitlength_in_kpc
g.x = g.x*unitlength_in_kpc
maxpos = MAX([g.y,g.x])
minpos = MIN([g.y,g.x])
g = g[where(g.x lt 25.0/2.0 AND g.x gt -25.0/2.0 AND g.y lt 25.0/2.0 AND g.y gt -25.0/2.0)]
range = maxpos - minpos
nbin = 600
binsize = range/nbin
x = findgen(nbin)*(range/nbin)+minpos
y = x
x = x[200:400]
y = y[200:400]
eps = d.eps[0]*unitlength_in_kpc/4.0*3.0
metallicity = fltarr(N_ELEMENTS(x),N_ELEMENTS(y))
for ctx = 0, N_ELEMENTS(x) - 1 DO BEGIN
    FOR cty = 0, N_ELEMENTS(y) - 1 DO BEGIN
        ind = where(g.x ge x[ctx]-eps AND g.x lt x[ctx]+binsize+eps AND g.y ge y[cty]-eps AND g.y lt y[cty]+binsize+eps)
       IF (ind[0] ne -1) THEN metallicity[ctx,cty] = TOTAL(g[ind].zmetal)/N_ELEMENTS(g[ind].zmetal) ELSE metallicity[ctx,cty] = 0
    ENDFOR
ENDFOR
window,1
x=FINDGEN((SIZE(metallicity))[1])+200
y=FINDGEN((SIZE(metallicity))[2])+200
contour,metallicity,x,y,levels=theselevels,/fill,ystyle = 1,xstyle = 1
set_plot,'ps'
device,filename = 'metaldistro_MW1lr_std.eps',/color,bits_per_pixel = 8,/times
contour,metallicity,x,y,levels = theselevels,/fill
device,/close
stop


;stop
;cam1_aux = 22
;cam1_aux_data = mrdfits(file3,cam1_aux,h)
;metallicity = cam1_aux_data[*,*,1]/cam1_aux_data[*,*,0]
;temp = where(Finite(metallicity),complement = nan)
;metallicity[nan] = 1
;peak = 0.015
;metallicity = metallicity[200:400,200:400]
;x=FINDGEN((SIZE(metallicity))[1])+200
;y=FINDGEN((SIZE(metallicity))[2])+200
;theselevels = [peak*0.1,peak*0.15,peak*0.2,peak*0.25,peak*0.3,peak*0.35,peak*0.4,peak*0.45,peak*0.5,peak*0.55,peak*0.6,peak*0.65,peak*0.7,peak*0.75,peak*0.8,peak*0.85,0.9*peak,0.95*peak,peak]
;contour,metallicity,x,y,levels=theselevels,/fill
end

;.r /astro/users/christensen/code/Dwarfs/prof
pro metal_profile
set_plot,'x'
msol = 3.17111799e15
kpc = 28571.0
rtipsy,'~/Scratch2/Sunrise/MW1lr_v2/MW1lr_001.std',h,g,d,s
g.mass = g.mass*msol
g.x = g.x*kpc
g.y = g.y*kpc
g.z = g.z*kpc
g_prof = prof(g, 'gas', h.time,rmax = 5.0,nbins = 20.0)
g_prof.rbins = g_prof.rbins;*kpc
g_prof.metal_den = g_prof.metal_den;/msol*kpc*kpc
g_prof.rho = g_prof.rho;*msol/kpc/kpc
set_plot,'ps'
device,filename = 'metalprofile_MW1lr.eps',/color,bits_per_pixel = 8,/times
!p.multi = [0,1,3,0,1]
multiplot
plot,g_prof.RBINS,g_prof.metal_den,ytitle = textoidl('\Sigma_{metal}')+' [M'+sunsymbol()+textoidl(' kpc^{-2}]'),/ylog,charsize = 1,title = 'Average Metal/Gas Profile, MW1lr'

multiplot
plot,g_prof.RBINS,g_prof.rho,ytitle = textoidl('\Sigma_{gas}')+' [M'+sunsymbol()+textoidl(' kpc^{-2}]'),/ylog,charsize = 1

multiplot
plot,g_prof.RBINS,g_prof.metal_den/g_prof.rho,xtitle = 'Radius [kpc]',ytitle = 'z',/ylog,charsize = 1
multiplot
multiplot,/reset
multiplot,/default
!p.multi = 0
device,/close

stop
end
