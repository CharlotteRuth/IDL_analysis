FUNCTION interpol_color,x1,y1,y2,yline1,xline1,yline2,xline2
segment1 = where(yline1 LE y2 AND yline1 GE y1)
IF(segment1[0] EQ -1) THEN diff = MIN(ABS(yline1-y1),segment1)
segment1 = [MIN(segment1)-1,segment1,MAX(segment1)+1];This brackets off the points below and above the magnitdues we are interested in
segment2 = where(yline2 LE y2 AND yline2 GE y1)
IF(segment2[0] EQ -1) THEN diff = MIN(ABS(yline2-y1),segment2)
segment2 = [MIN(segment2)-1,segment2,MAX(segment2)+1]
xs1 = SPLINE(yline1[segment1],xline1[segment1],[y1,y2])
xs2 = SPLINE(yline2[segment2],xline2[segment2],[y1,y2])
IF(xs1[0] GT MAX(xline1[segment1]) OR xs1[0] LT MIN(xline1[segment1])) THEN BEGIN
    temp = MIN(ABS(yline1 - y1),nearest)
    xs1[0] = xline1[nearest]
ENDIF
IF(xs1[1] GT MAX(xline1[segment1]) OR xs1[1] LT MIN(xline1[segment1])) THEN BEGIN
    temp = MIN(ABS(yline1 - y2),nearest)
    xs1[1] = xline1[nearest]
ENDIF
IF(xs2[0] GT MAX(xline2[segment2]) OR xs2[0] LT MIN(xline2[segment2])) THEN BEGIN
    temp = MIN(ABS(yline2 - y1),nearest)
    xs2[0] = xline2[nearest]
ENDIF
IF(xs2[1] GT MAX(xline2[segment2]) OR xs2[1] LT MIN(xline2[segment2])) THEN BEGIN
    temp = MIN(ABS(yline2 - y2),nearest)
    xs2[1] = xline2[nearest]
ENDIF

weight1 = (x1-xs1[0])/(xs2[0]-xs1[0]) ;weight by how close to each boundery the point is
x2  = xs1[1] + weight1*(xs2[1]-xs1[1])
RETURN,x2
END
