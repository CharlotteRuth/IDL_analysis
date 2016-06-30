function integral,radius2,z
;d = distance from center of kernel
;z = distance along line of sight
;radius2 = distance between line of sight and center of kernel

d = sqrt(z*z + radius2)
ind = where(d LE 0, complement = ind0)
ind1 = where(d GT 2)
integral = fltarr(n_elements(z))
IF ind[0]  NE -1 THEN $
integral[ind ] = (1.0d0 - 1.5*radius2)*z[ind ] - 0.5*z[ind ]*z[ind ]*z[ind ] + 0.1875*z[ind ]*d[ind ]*d[ind ]*d[ind ] +         0.28125*radius2*z[ind ]*d[ind ] + 0.28125*radius2*radius2                *alog(z[ind ] + d[ind ])
IF ind0[0] NE -1 THEN $
integral[ind0] = (2.0d0 + 1.5*radius2)*z[ind0] + 0.5*z[ind0]*z[ind0]*z[ind0] - 0.0625*z[ind0]*d[ind0]*d[ind0]*d[ind0] - (1.5 + 0.09375*radius2)*z[ind0]*d[ind0] - (1.5*radius2 + 0.09375*radius2*radius2)*alog(z[ind0] + d[ind0])
IF ind1[0] NE -1 THEN integral[ind1] = 0
;stop
return,integral
END
