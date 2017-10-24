metals = star[1100000:1100001].metals
ages = log_age[1100000:1100001]
mags = Min_Curve_surf(V_grid, xvalues =  metallicity_vector, yvalues = age_vector, xpout = metals, ypout = ages,/tps)
plot,metallicity_vector,V_grid[*,37],yrange =[MIN(mags),MAX(mags)]
oplot,metals,mags,psym = 1
contour,V_grid,metallicity_vector,age_vector,NLEVELS = 10,ytitle = "Log(Age)",xtitle = "Metallicty"
oplot,star.metals,log_age,psym = 3

n_ind = 50
m_low = MIN(metallicity_vector)
m_hi = MAX(metallicity_vector)
a_low = MIN(age_vector)
a_hi = MAX(age_vector)
m_spacing = (m_hi - m_lo)/(n_ind - 1)
m_vector = findgen(n_ind)*m_spacing + m_low
a_spacing = (a_hi - a_low)/(n_ind - 1)
a_vector = findgen(n_ind)*a_spacing + a_low
mags_grid = TRI_SURF(V_grid, xvalues =  metallicity_vector, yvalues = age_vector,GS=[m_spacing, a_spacing], BOUNDS=[m_low,a_low,m_hi,a_hi])

contour,mags_grid,xvector,yvector,NLEVELS = 10,ytitle = "Log(Age)",xtitle = "Metallicty"
