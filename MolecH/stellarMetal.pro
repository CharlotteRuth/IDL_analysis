;dir = '/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG/'
;filename = 'h603.cosmo50cmb.2304g5bwK.00512.halo.1.std'
;dKpcUnit =50000.
;dMsolUnit = 1.84793e16

;dir2 ='/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK.BUG/h603.cosmo50cmb.2304g2bwK.00512'
;filename2 = 'h603.cosmo50cmb.2304g2bwK.00512.halo.1.std'


;filename1 = '12M_hr.00800'
;dir1 = '/astro/net/nbody1/christensen/MolecH/MWHR/12M_hr'

pro stellarMetal,dir1,filename1,dir2,filename2,dKpcUnit,dMsolUnit
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790d+23
molec_weight = (0.76*1 + 0.24*4.0)
f_H = 0.764

set_plot,'x'
loadct,39
cd,dir1
rtipsy,filename1,h,g,d,s
s.x = s.x*dKpcUnit
s.y = s.y*dKpcUnit
s.z = s.z*dKpcUnit
s.mass = s.mass*dMsolUnit
g.x = g.x*dKpcUnit
g.y = g.y*dKpcUnit
g.z = g.z*dKpcUnit
g.mass = g.mass*dMsolUnit
cd,dir2
rtipsy,filename2,h2,g2,d2,s2
s2.x = s2.x*dKpcUnit
s2.y = s2.y*dKpcUnit
s2.z = s2.z*dKpcUnit
s2.mass = s2.mass*dMsolUnit
g2.x = g2.x*dKpcUnit
g2.y = g2.y*dKpcUnit
g2.z = g2.z*dKpcUnit
g2.mass = g2.mass*dMsolUnit

rmax = 40
rmin = 5
rclip = rmax
nbins = 200
fitpar =  [0,0,0,0,0]
fitpar2 =  [0,0,0,0,0]

cofm, h, g, d, s, /pot
s_prof = prof(s, 'star', h.time, nbins = nbins, rmax = rmax)
ind = where(s_prof.rbins gt rmin and s_prof.rbins lt rclip)
dblexpfit, s_prof.rbins[ind], s_prof.rho[ind], $
  s_prof.rho[ind]/sqrt(s_prof.num_in_bin[ind]), fitpar, red_chisq=exp2chi, /no_guess
cofm, h2, g2, d2, s2, /pot
s2_prof = prof(s2, 'star', h2.time, nbins = nbins, rmax = rmax)
ind = where(s2_prof.rbins gt rmin and s_prof.rbins lt rclip)
dblexpfit, s2_prof.rbins[ind], s2_prof.rho[ind], $
  s2_prof.rho[ind]/sqrt(s2_prof.num_in_bin[ind]), fitpar2, red_chisq=exp2chi, /no_guess

window,0
xtitle = 'Radius [kpc]'
ytitle = 'Stellar Density [M_sol pc^-2]'
plot, s_prof.rbins, s_prof.rho, /ylog, $
  ytitle = ytitle, title = title,$  
  ystyle=1, yrange=[1e2,1e11],$
  xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5
r1 = findgen(100)/100.*(fitpar[0])+3.
r2 = findgen(100)/100.*(rmax-fitpar[0]+3.) + fitpar[0]-3.
oplot, r1, fitpar[1]*exp(-r1/fitpar[2]), line = 1
oplot, r2, fitpar[3]*exp(-r2/fitpar[4]), line = 1
oplot, [fitpar[0],fitpar[0]], [1e4,1e11], line = 3
oplot,s2_prof.rbins, s2_prof.rho,color = 240
r1 = findgen(100)/100.*(fitpar2[0])+3.
r2 = findgen(100)/100.*(rmax-fitpar2[0]+3.) + fitpar2[0]-3.
oplot, r1, fitpar2[1]*exp(-r1/fitpar2[2]), line = 1,color = 240
oplot, r2, fitpar2[3]*exp(-r2/fitpar2[4]), line = 1,color = 240
oplot, [fitpar2[0],fitpar2[0]], [1e4,1e11], line = 3,color = 240

;**********************************************
maxdistance = 20.
nbins = 100.
bind = maxdistance/nbins
radius = findgen(nbins)/nbins*maxdistance + bind/2.0
sradii = sqrt(s.x*s.x + s.y*s.y + s.z*s.z)
gradii = sqrt(g.x*g.x + g.y*g.y + g.z*g.z)
ave_metal_s = radius
ave_metal_g = radius
sradii2 = sqrt(s2.x*s2.x + s2.y*s2.y + s2.z*s2.z)
gradii2 = sqrt(g2.x*g2.x + g2.y*g2.y + g2.z*g2.z)
ave_metal_s2 = radius
ave_metal_g2 = radius
FOR i = 0, nbins - 1 DO BEGIN
    ind = WHERE(sradii ge radius[i] - bind AND sradii lt radius[i] + bind)
    if (ind[0] ne -1) THEN ave_metal_s[i] = TOTAL(s[ind].metals*s[ind].mass)/TOTAL(s[ind].mass) else ave_metal_s[i] = 0

    ind = WHERE(gradii ge radius[i] - bind AND gradii lt radius[i] + bind)
    if (ind[0] ne -1) THEN ave_metal_g[i] = TOTAL(g[ind].zmetal*g[ind].mass)/TOTAL(g[ind].mass) else ave_metal_g[i] = 0

    ind = WHERE(sradii2 ge radius[i] - bind AND sradii2 lt radius[i] + bind)
    if (ind[0] ne -1) THEN ave_metal_s2[i] = TOTAL(s2[ind].metals*s2[ind].mass)/TOTAL(s2[ind].mass) else ave_metal_s2[i] = 0

    ind = WHERE(gradii2 ge radius[i] - bind AND gradii2 lt radius[i] + bind)
    if (ind[0] ne -1) THEN ave_metal_g2[i] = TOTAL(g2[ind].zmetal*g2[ind].mass)/TOTAL(g2[ind].mass) else ave_metal_g2[i] = 0
ENDFOR

ind_fit = where(radius le 20 AND radius ge 2)
fit1 = poly_fit(radius[ind_fit],alog10(ave_metal_s[ind_fit])  - alog10(0.0130215),1)
fit2 = poly_fit(radius[ind_fit],alog10(ave_metal_s2[ind_fit]) - alog10(0.0130215),1)

window,1
;set_plot,'ps'
;device,filename ='~/h603_metal.ps',/color,bits_per_pixel= 8,/times
plot,radius,alog10(ave_metal_s) - alog10(0.0130215),xtitle = 'Rg [kpc]',ytitle = '[Z/H]',yrange = [-1.5,0.5]
oplot,radius,alog10(ave_metal_g) - alog10(0.0130215),linestyle = 2
oplot,radius,fit1[0] + fit1[1]*radius,linestyle = 1

;oplot,radius,alog10(ave_metal_s2) - alog10(0.0130215),color = 240
;oplot,radius,alog10(ave_metal_g2) - alog10(0.0130215),linestyle = 2,color = 240
;oplot,radius,fit2[0] + fit2[1]*radius,linestyle = 1,color = 240
legend,['Stars','Gas'],linestyle = [0,2],/right
;device,/close

print,'Slope1: ',fit1[1]
print,'Slope2: ',fit2[1]
stop
END
