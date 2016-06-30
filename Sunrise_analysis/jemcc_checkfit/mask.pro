;This program will print to the screen a list of the x and y
;coordinates of all pixels within r of the given center coordinate. 


pro mask,x_0,y_0,r

xmin=x_0-r
xmax=x_0+r
ymin=y_0-r
ymax=y_0+r
dist=0

openw, lun, 'mask.asc', /get_lun

for x=xmin,xmax,1 do begin
    for y=ymin,ymax,1 do begin
        dist=(((x-x_0)^2)+((y-y_0)^2))^.5
        if (dist le r) then begin
             printf, lun, x,y 
        endif
    endfor
endfor
close, lun

end



