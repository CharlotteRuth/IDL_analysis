pro interpolate_mags
set_plot,'x'
sec_to_year = 365.*24.*3600.
timestep = 4e10
loadct,39

;Leo Girardi Isochrones ---
n_age_l = 42
n_metallicity_l = 15 
age_bin_l = 0.1
age_min_l = 6.0
age_vector_l = findgen(n_age_l)*age_bin_l + age_min_l
metallicity_vector_l = [0.0004,0.0005,0.0007,0.0010,0.0013,0.0017,0.0023,0.0030,0.0040,0.0054,0.0073,0.0098,0.0130,0.0234,0.03]
openr, 1, '/astro/net/scratch2/christensen/Database/procedures/data_leo/V_grid'
openr, 2, '/astro/net/scratch2/christensen/Database/procedures/data_leo/vr_grid'
openr, 3, '/astro/net/scratch2/christensen/Database/procedures/data_leo/ri_grid'
openr, 4, '/astro/net/scratch2/christensen/Database/procedures/data_leo/bv_grid'
openr, 5, '/astro/net/scratch2/christensen/Database/procedures/data_leo/vk_grid'
openr,6,'/astro/net/scratch2/christensen/Database/procedures/data_leo/ub_grid'

print,'LEO'
print,'---------------------'
print,'Range of Metallicities: '
print,MINMAX(metallicity_vector_l)
print,'Range of Ages:          '
print,MINMAX(age_vector_l)

V_grid_l =  fltarr(n_metallicity_l,n_age_l)
vr_grid_l = fltarr(n_metallicity_l,n_age_l)
ri_grid_l = fltarr(n_metallicity_l,n_age_l)
bv_grid_l = fltarr(n_metallicity_l,n_age_l)
vk_grid_l = fltarr(n_metallicity_l,n_age_l)
ub_grid_l = fltarr(n_metallicity_l,n_age_l)

readf, 1, V_grid_l
readf, 2, vr_grid_l
readf, 3, ri_grid_l
readf, 4, bv_grid_l
readf, 5, vk_grid_l
readf, 6, ub_grid_l
close, /all

;Starburst99 Isochrones ----
n_metallicity_s = 5
readcol,'/astro/net/scratch2/christensen/Database/procedures/data_star99/ages',age_vector_s
age_vector_s = ALOG10(age_vector_s)
age_min_s = MIN(age_vector_s)
age_max_s = MAX(age_vector_s)
n_age_s = N_ELEMENTS(age_vector_s)
metallicity_vector_s = [0.0004,0.004,0.008,0.02,0.05]

openr, 1, '/astro/net/scratch2/christensen/Database/procedures/data_star99/V_grid'
openr, 2, '/astro/net/scratch2/christensen/Database/procedures/data_star99/vr_grid'
openr, 3, '/astro/net/scratch2/christensen/Database/procedures/data_star99/ri_grid'
openr, 4, '/astro/net/scratch2/christensen/Database/procedures/data_star99/bv_grid'
openr, 5, '/astro/net/scratch2/christensen/Database/procedures/data_star99/vk_grid'
openr, 6, '/astro/net/scratch2/christensen/Database/procedures/data_star99/ub_grid'

print,''
print,'STARBURST99'
print,'---------------------'
print,'Range of Metallicities: '
print,MINMAX(metallicity_vector_s)
print,'Range of Ages:          '
print,MINMAX(age_vector_s)

V_grid_s =  fltarr(n_metallicity_s,n_age_s)
vr_grid_s = fltarr(n_metallicity_s,n_age_s)
ri_grid_s = fltarr(n_metallicity_s,n_age_s)
bv_grid_s = fltarr(n_metallicity_s,n_age_s)
vk_grid_s = fltarr(n_metallicity_s,n_age_s)
ub_grid_s = fltarr(n_metallicity_s,n_age_s)

readf, 1, V_grid_s
readf, 2, vr_grid_s
readf, 3, ri_grid_s
readf, 4, bv_grid_s
readf, 5, vk_grid_s
readf, 6, ub_grid_s
close, /all

all_metal_range=MAX([metallicity_vector_s,metallicity_vector_l])-MIN([metallicity_vector_s,metallicity_vector_l])
all_metal_min=MIN([metallicity_vector_s,metallicity_vector_l])
n_ind_m = 26
n_ind_a = 830   
m_low_l = MIN(metallicity_vector_l)
m_hi_l = MAX(metallicity_vector_l)
m_low_s = MIN(metallicity_vector_s)
m_hi_s = MAX(metallicity_vector_s)

a_low_s = MIN(age_vector_s)
a_hi_s = MAX(age_vector_s) 
a_low_l = MIN(age_vector_l)
a_hi_l = MAX(age_vector_l)    
m_spacing_l = (m_hi_l - m_low_l)/(n_ind_m-1)
m_vector_l = findgen(n_ind_m)*m_spacing_l + m_low_l
m_spacing_s = (m_hi_s - m_low_s)/(n_ind_m-1)
m_vector_s = findgen(n_ind_m)*m_spacing_s + m_low_s

a_spacing_l = (a_hi_l - a_low_l)/(n_ind_a-1)
a_vector_l = findgen(n_ind_a)*a_spacing_l + a_low_l
a_spacing_s = (a_hi_s - a_low_s)/(n_ind_a-1)
a_vector_s = findgen(n_ind_a)*a_spacing_s + a_low_s
V_grid2_l = TRI_SURF(V_grid_l, xvalues =  metallicity_vector_l, yvalues = age_vector_l,nx = n_ind_m, ny = n_ind_a)
V_grid2_s = TRI_SURF(V_grid_s, xvalues =  metallicity_vector_s, yvalues = age_vector_s,nx = n_ind_m, ny = n_ind_a)

window,2
plot,a_vector_l,V_grid2_l[*,0],xrange=[4,10.1],yrange=[-16,-6],xtitle='Log Age',ytitle='V Magnitude (6 Msol SSP)'
for i = 0, n_ind_m - 1 Do BEGIN
    oplot,a_vector_l,V_grid2_l[i,*]- 2.51*6.0,color=(m_vector_l[i]-all_metal_min)/all_metal_range*220.0+20
    oplot,a_vector_s,V_grid2_s[i,*],color=(m_vector_s[i]-all_metal_min)/all_metal_range*220.0+20
ENDFOR   
stop

plot,a_vector_l,V_grid2_l[*,0],xrange=[4,10.1],yrange=[-16,-8],xtitle='Log Age',ytitle='V Magnitude (6 Msol SSP)'
for i = 0, n_ind_m - 1 Do BEGIN
    oplot,a_vector_s,V_grid2_s[i,*],color=(m_vector_s[i]-all_metal_min)/all_metal_range*220.0+20
ENDFOR
for i = 0, 5 - 1 Do oplot,age_vector_s,V_grid_s[i,*]
 
stop
plot,a_vector_l,V_grid2_l[*,0],xrange=[6,8],yrange=[-16,-12],xtitle='Log Age',ytitle='V Magnitude (6 Msol SSP)'
for i = 0, n_ind_m - 1 Do BEGIN
    oplot,a_vector_s,V_grid2_s[i,*],color=(m_vector_s[i]-all_metal_min)/all_metal_range*220.0+20
ENDFOR
for i = 0, 5 - 1 Do oplot,age_vector_s,V_grid_s[i,*]

set_plot,'ps'
device,filename='interpolate_mags.eps',/color,bits_per_pixel=8
plot,a_vector_l,V_grid2_l[*,0],xrange=[4,10.1],yrange=[-16,-8],xtitle='Log Age',ytitle='V Magnitude (6 Msol SSP)'
for i = 0, n_ind_m - 1 Do BEGIN
    oplot,a_vector_s,V_grid2_s[i,*],color=(m_vector_s[i]-all_metal_min)/all_metal_range*220.0+20
ENDFOR
for i = 0, 5 - 1 Do oplot,age_vector_s,V_grid_s[i,*]
legend,['Starburst99, SSP at Constant z','Interpolation, z = 0.05','Interpolation, z = 0.0004'],color=[0,240,20]
device,/close
 
device,filename='interpolate_mags_zoom.eps',/color,bits_per_pixel=8
plot,a_vector_l,V_grid2_l[*,0],xrange=[6,8],yrange=[-16,-12],xtitle='Log Age',ytitle='V Magnitude (6 Msol SSP)'
for i = 0, n_ind_m - 1 Do BEGIN
    oplot,a_vector_s,V_grid2_s[i,*],color=(m_vector_s[i]-all_metal_min)/all_metal_range*220.0+20
ENDFOR
for i = 0, 5 - 1 Do oplot,age_vector_s,V_grid_s[i,*]
legend,['Starburst99, SSP at Constant z','Interpolation, z = 0.05','Interpolation, z = 0.0004'],color=[0,240,20] 
device,/close
END
