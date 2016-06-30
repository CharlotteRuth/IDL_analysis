
PRO opt_depth,filename,kpcunit,massunit,massform,colors,psyms, outplot = outplot,overplot = overplot,guo = guo
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
molec_weight = (0.76*1 + 0.24*4.0)
;massunit = 1.84793e16
;kpcunit = 50000.
massform = massform*massunit
;htime = 1.0
diskheight = 2.5/2.0
rtipsy,filename + '.halo.1.std',h,g,d,s
rgas=SQRT(g.x^2+g.y^2)*kpcunit

maxrad = 20.0
nrad = 20
radii = findgen(nrad)*maxrad
surface_z = fltarr(nrad - 1)

FOR i = 0, nrad - 2 DO BEGIN
    indgas=WHERE(rgas LT radii[i + 1] AMD rgas gt radii[i] ,nindgas)
    surface_z[i] = TOTAL(g[indgas].mass*g[indgas].zmetal*massunit)/(!PI*(radii[i + 1]^2 - radii[i]^2))
END

plot,radii[0:nrad - 2] + maxrad/nrad,surface_z*0.4*5.86e-6

END
