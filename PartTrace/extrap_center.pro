;CC, 12/12/11

;This function return the spline interpolated center of a galaxy given a
;time.
;It uses the output from read_center.pro, an array of
;(t,z,x_c,y_c,z_c) from the amiga.stats file 

FUNCTION extrap_center,time,center_stats

x = spline(center_stats[*,0],center_stats[*,2],time)
y = spline(center_stats[*,0],center_stats[*,3],time)
z = spline(center_stats[*,0],center_stats[*,4],time)

RETURN,[[x],[y],[z]]
END
