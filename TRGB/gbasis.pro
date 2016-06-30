FUNCTION GBASIS,x1,x2,xline1,yline1,xline2,yline2
;Generates an orthonormal basis function for which one of the basis
;vectors is the average slope of the line through the two tips
xline1 = reform(xline1)
yline1 = reform(yline1)
xline2 = reform(xline2)
yline2 = reform(yline2)

temp = MIN(ABS(xline1-x1),tip1)
temp = MIN(ABS(xline2-x2),tip2) 
x_prime = (xline1[tip1] + xline2[tip2])/2
y_prime = (yline1[tip1] + yline2[tip2])/2
x =  xline2[tip2]- x_prime
y =  yline2[tip2]- y_prime
print,'X:       ',x,' Y:            ',y
;x_prime = (xline1[tip1] - xline2[tip2])
;y_prime = (yline1[tip1] - yline2[tip2])*(-1.0)
mag = SQRT(x*x+y*y)
print,mag
basis = [[x/mag,y/mag],[(-1.0)*y/mag,x/mag]]
print,basis
;stop
RETURN,basis
END

FUNCTION trans_rotate,x1,x2,xline1,yline1,xline2,yline2,points1,inverse = inverse,error = error
;This calls GBASIS to generate an orthonomatl basis function for which
;one of the basis vectors is the average slope of the line through the
;two tips
;It them moves the predicted tip to the origin, rotates such that the
;tip basis vector is horizontal and then moves the points back
;basis = GBASIS(x1,x2,xline1,yline1,xline2,yline2)
points = points1
xline1 = reform(xline1)
yline1 = reform(yline1)
xline2 = reform(xline2)
yline2 = reform(yline2)

temp = MIN(ABS(xline1-x1),tip1)
temp = MIN(ABS(xline2-x2),tip2) 
x_prime = (xline1[tip1] + xline2[tip2])/2
y_prime = (yline1[tip1] + yline2[tip2])/2
x =  xline2[tip2]- x_prime
y =  yline2[tip2]- y_prime
mag = SQRT(x*x+y*y)
;primes = [x_prime,y_prime]#basis
IF (KEYWORD_SET(inverse)) THEN BEGIN
    basis = [[x/mag,(-1.0)*y/mag],[y/mag,x/mag]]
    points[*,0] = points1[*,0] - x_prime;- primes[0]
    points[*,1] = points1[*,1] - y_prime;- primes[1]
    points = points#basis ;Rotate so that the line connecting the tips is the x-axes
    points[*,0] = points[*,0] + x_prime ;Translate such that the point between the two tips is the origin
    points[*,1] = points[*,1] + y_prime  
    IF(KEYWORD_SET(error)) THEN points = points1#basis ;No need to translate
ENDIF ELSE BEGIN
    basis = [[x/mag,y/mag],[(-1.0)*y/mag,x/mag]]    
    points[*,0] = points1[*,0] - x_prime ;Translate such that the point between the two tips is the origin
    points[*,1] = points1[*,1] - y_prime
    points = points#basis ;Rotate so that the line connecting the tips is the x-axes
    points[*,0] = points[*,0] + x_prime;+ primes[0] ;Translate back
    points[*,1] = points[*,1] + y_prime; primes[1]
    IF(KEYWORD_SET(error)) THEN points = points1#basis
ENDELSE
return,points
end




