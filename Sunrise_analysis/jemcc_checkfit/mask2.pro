;This program will print to screen/file a list of the x and y
;coordinates of all pixels that fall within the annulus defined by
;r_min and r_max

;The point of this is to obscure a central bar in a galaxy so that
;GALFIT can fit the bulge and disk without pissing itself

;Appropriate values for r_min and r_max may be found using e.g. checkfit


pro mask2, x_0, y_0, r_min, r_max

if (n_params() eq 0) then begin
 print,' Usage: mask, x center, y center, inner radius [pix], outer radius [pix]'
 stop
endif


xmin=x_0-r_max
xmax=x_0+r_max
ymin=y_0-r_max
ymax=y_0+r_max

openw, lun1, 'mask.asc', /get_lun;,/append

for x=xmin,xmax,1 do begin
    for y=ymin,ymax,1 do begin
        dist=(((x-x_0)^2.)+((y-y_0)^2.))^0.5
        if ( (dist le r_max) and (dist ge r_min) and (dist ne 0) )then begin
             printf, lun1, x,y 
        endif
    endfor
endfor

FREE_LUN, lun1

;readcol, 'mask.asc',xcol,ycol

;plot,xcol,ycol,xrange=[150,350],yrange=[150,350],psym=3,xtit='X',ytit='Y',/isotropic

end



