FUNCTION GBASIS,x1,x2,xline1,yline1,xline2,yline2
;Generates an othonormal basis function for which one of the basis
;vectors is the average slope of the two bounding isochrones near the
;tip
xline1 = reform(xline1)
yline1 = reform(yline1)
xline2 = reform(xline2)
yline2 = reform(yline2)

temp = MIN(ABS(xline1-x1),tip1)
temp = MIN(ABS(xline2-x2),tip2) 
x1 = ((xline1[tip1+1] - xline1[tip1-1]) + (xline2[tip2+1] - xline2[tip2-1]))/2
y1 = ((yline1[tip1+1] - yline1[tip1-1]) + (yline2[tip2+1] - yline2[tip2-1]))/2
mag = SQRT(x1*x1+y1*y1)
basis = [[x1/mag,y1/mag],[y1/mag,-1.0*x1/mag]]
RETURN,basis
END

FUNCTION GBASIS_tips,x1,x2,xline1,yline1,xline2,yline2
;Generates an othonormal basis function for which one of the basis
;vectors is the average slope of the line through the two tips
xline1 = reform(xline1)
yline1 = reform(yline1)
xline2 = reform(xline2)
yline2 = reform(yline2)

temp = MIN(ABS(xline1-x1),tip1)
temp = MIN(ABS(xline2-x2),tip2) 
;print, ((xline1[tip1+1] - xline1[tip1-1]) + (xline2[tip2+1] - xline2[tip2-1]))/2
;print,((yline1[tip1+1] - yline1[tip1-1]) + (yline2[tip2+1] - yline2[tip2-1]))/2
y1 = (xline1[tip1] - xline2[tip2])
x1 = (yline1[tip1] - yline2[tip2])*(-1.0)
;print,x1,y1
mag = SQRT(x1*x1+y1*y1)
basis = [[x1/mag,y1/mag],[y1/mag,-1.0*x1/mag]]
RETURN,basis
END


