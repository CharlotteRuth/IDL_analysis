function above,x,y,pointx,pointy
;This function will check if (pointx,pointy) lies above/right of the
;line defined by x,y
y_cutoff = (where(y lt pointy))[0] ;First time the y value is less than the x point
if (y_cutoff eq -1) then y_cutoff = (where(x gt pointx))[0]
if (y_cutoff lt 1) THEN RETURN,1
y2 = y[y_cutoff]
y1 = y[y_cutoff - 1]
x2 = x[y_cutoff]
x1 = x[y_cutoff - 1]
slope = (y2 - y1)/(x2 - x1)
if (pointy le (slope*(pointx - x1) + y1)) then return,1 
return,0
end
