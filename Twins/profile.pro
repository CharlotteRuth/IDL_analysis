pro profile
infile = '/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512/h603.cosmo50cmb.2304g2bwK.00512.1.std'
outdir = '/astro/net/scratch1/christensen/Twins/h603.cosmo50cmb.2304g2bwK/'
galname = 'h603'
rmax = 25.
kpcunit = 50000.
msol = 1.84793e16
timeunit = 3.8785614e+10

;infile = '/astro/net/scratch1/abrooks/FABIO/h516.cosmo25cmb.2304g2bwK/h516.cosmo25cmb.2304g2bwK.00512/h516.cosmo25cmb.2304g2bwK.00512.1.std'
;outdir = '~/Scratch1/Twins/h516.cosmo25cmb.2304g2bwK/'
;galname = 'h516'
;rmax = 12.5
;kpcunit = 25000.
;msol = 2.310e15
;timeunit = 3.8785614e+10

loadct,39
set_plot,'ps'
;set_plot,'x'
!P.CHARSIZE = 1.25
!p.thick = 2.5
!X.style = 1.25
!Y.style = 1.25

SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
nbins = 50
rmin = 0

rtipsy,infile,h,g,d,s
;dbin = 4./kpcunit
;nbins = (MAX(SQRT(s.x*s.x + s.y*s.y + s.z*s.z)))/dbin
;stop
s_prof = prof(s, 'star', h.time, nbins = nbins,rmax = rmax/kpcunit)
g_prof = prof(g, 'gas', h.time, nbins = nbins,rmax = rmax/kpcunit)
d_prof = prof(d, 'dark', h.time, nbins = nbins,rmax = rmax/kpcunit)  
  
s_prof.rbins = s_prof.rbins*kpcunit
d_prof.rbins = d_prof.rbins*kpcunit
g_prof.rbins = g_prof.rbins*kpcunit
s_prof.rho = s_prof.rho*msol/kpcunit/kpcunit  ; Msol/kpc^2
d_prof.rho = d_prof.rho*msol/kpcunit/kpcunit*1e-6  
g_prof.rho = g_prof.rho*msol/kpcunit/kpcunit*1e-6 ; Msol/pc^2

fitpar =  [10,5e8,1.5,1e7,5]
dblexpfit, s_prof.rbins, s_prof.rho, $
  s_prof.rho/sqrt(s_prof.num_in_bin), fitpar
;set_plot,'x'
;window,2
device,filename = outdir +'profile.eps',/color,bits_per_pixel=8,/times
!p.multi = [0,1,3]
multiplot
plot,s_prof.rbins,s_prof.rho,/ylog,title = 'Surface Density Profiles of '+galname,ytitle = textoidl('\Sigma_{star}')+"[M"+sunsymbol()+" kpc!u-2!n]",xrange = [0,rmax],yrange = [MIN(s_prof.rho),max(s_prof.rho)]
r1 = findgen(100)/100.*(fitpar[0])
r2 = findgen(100)/100.*(rmax-fitpar[0]) + fitpar[0]
;oplot, r1, fitpar[1]*exp(-r1/fitpar[2]), line = 1,color = 240
;oplot, r2, fitpar[3]*exp(-r2/fitpar[4]), line = 1, color = 240
;oplot, [fitpar[0],fitpar[0]], [1,1e8], line = 3
multiplot
plot,g_prof.rbins,g_prof.rho,/ylog,ytitle = textoidl('\Sigma_{gas}')+" [M"+sunsymbol()+" pc!u-2!n]",xrange = [0,rmax],yrange = [MIN(g_prof.rho),max(g_prof.rho)]
; oplot,  [fitpar[0],fitpar[0]], [1e-5,1e5], line = 3
multiplot
plot,d_prof.rbins,d_prof.rho,/ylog,xtitle = 'Radius (kpc)',ytitle = textoidl('\Sigma_{dark}')+"[M"+sunsymbol()+" pc!u-2!n]",xrange = [0,rmax],yrange = [MIN(d_prof.rho),max(d_prof.rho)]
multiplot
device,/close


;set_plot,'x'
s_prof = prof(s, 'star', h.time, nbins = nbins,/sph,rmax = rmax/kpcunit)
g_prof = prof(g, 'gas', h.time, nbins = nbins,/sph,rmax = rmax/kpcunit)
d_prof = prof(d, 'dark', h.time, nbins = nbins,/sph,rmax = rmax/kpcunit)  
  
s_prof.rbins = s_prof.rbins*kpcunit
d_prof.rbins = d_prof.rbins*kpcunit
g_prof.rbins = g_prof.rbins*kpcunit
s_prof.rho = s_prof.rho*msol/kpcunit/kpcunit/kpcunit ;Msol/kpc^3
d_prof.rho = d_prof.rho*msol/kpcunit/kpcunit/kpcunit*1e-9 
g_prof.rho = g_prof.rho*msol/kpcunit/kpcunit/kpcunit*1e-9 

fitpar =  [10,1.3e7,2,1e6,5]
dblexpfit, s_prof.rbins, s_prof.rho, $
  s_prof.rho/sqrt(s_prof.num_in_bin), fitpar

device,filename = outdir+'profile_sph.eps',/color,bits_per_pixel=8,/times
;window,1
!p.multi = [0,1,3]
multiplot
plot,s_prof.rbins,s_prof.rho,/ylog,title = 'Density Profiles of '+galname+', Spherical Bins',ytitle = textoidl('\rho_{star}')+"[M"+sunsymbol()+" kpc!u-3!n]",xrange = [0,rmax];,yrange = [MIN(s_prof.rho),max(s_prof.rho)]
;rmax = MAX(s_prof.rbins)
r1 = findgen(100)/100.*(fitpar[0])
r2 = findgen(100)/100.*(rmax-fitpar[0]) + fitpar[0]
;oplot, r1, fitpar[1]*exp(-r1/fitpar[2]), line = 1,color = 240
;oplot, r2, fitpar[3]*exp(-r2/fitpar[4]), line = 1, color = 240
;oplot, [fitpar[0],fitpar[0]], [1,1e8], line = 3
multiplot
plot,g_prof.rbins,g_prof.rho,/ylog,ytitle = textoidl('\rho_{gas}')+" [M"+sunsymbol()+" pc!u-3!n]",xrange = [0,rmax];,yrange = [MIN(g_prof.rho),max(g_prof.rho)]
multiplot
plot,d_prof.rbins,d_prof.rho,/ylog,xtitle = 'Radius (kpc)',ytitle = textoidl('\rho_{dark}')+"[M"+sunsymbol()+" pc!u-3!n]",xrange = [0,rmax];,yrange = [MIN(d_prof.rho),max(d_prof.rho)]
multiplot
device,/close
;stop
END
