FUNCTION normalize,vector
;Returns the norlaized vecotr
mag = SQRT(TOTAL(vector * vector))
RETURN,vector/mag
END


PRO align_gals, list 
;align,'halo_list.dat'
msol = 2.362e5
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c

set_plot, 'ps'
IF(keyword_set(filename) EQ 0) THEN filename = 'align.ps'

;device, filename=filename, /encapsulated, ysize = 1.5*row, xsize = 2.5*column, /inches 
set_plot,'x'
loadct,39
;multiplot

;close,2
;openw,2,'fits.dat'

readcol, list, format = 'A60',files
FOR i = 0, n_elements(files)-1 DO BEGIN
    print,files[i]
    rtipsy, files[i], h, g, d, s
    points = [[g.x],[g.y],[g.z]]
    FOR ct = 0, 5 DO BEGIN
        print,'ITERATION: ',ct
        window,1
        plot,points[*,0],points[*,1],psym = 3,xtitle = 'X',ytitle = 'Y'
        result = POLY_FIT(points[*,0],points[*,1],1,chisq = chiz)
        slope_y_x = ABS(1.0/result[1])
        print,result
        x = findgen(100)/100.*(MAX(points[*,0]) - MIN(points[*,0])) + MIN(points[*,0])
        oplot,x,result[1]*x + result[0], color = 220
        
        window,2
        plot,points[*,2],points[*,1],psym = 3,xtitle = 'Z',ytitle = 'Y'   
        result = POLY_FIT(points[*,2],points[*,1],1,chisq = chix)
        slope_y_z = ABS(1.0/result[1])
        print,result
        x = findgen(100)/100.*(MAX(points[*,2]) - MIN(points[*,2])) + MIN(points[*,2])
        oplot,x,result[1]*x + result[0], color = 220
        
        window,3
        plot,points[*,0],points[*,2],psym = 3,xtitle = 'X',ytitle = 'Z'
        result = POLY_FIT(points[*,0],points[*,2],1,chisq = chiy)
        slope_z_x = ABS(1.0/result[1])
        print,result
        x = findgen(100)/100.*(MAX(points[*,0]) - MIN(points[*,0])) + MIN(points[*,0])
        oplot,x,result[1]*x + result[0],color = 220

        IF(chiz lt chiy AND chiz lt chix) THEN BEGIN
            perp = [1.0/slope_z_x,slope_y_z, 1.0]
            perp = normalize(perp)

            v1 = [1.0,0,0]
            v1 = v1 - TOTAL(v1*perp)*perp
            v1 = normalize(v1)
            
            v2 = [0,1.0,0] 
            v2 = v2 - TOTAL(v2*perp)*perp - TOTAL(v2*v1)*v1
            v2 = normalize(v2)

            rot_matrix = TRANSPOSE([[v1],[v2],[perp]])
        ENDIF ELSE BEGIN
            IF(chiy lt chix AND chiy lt chiz) THEN BEGIN

            ENDIF ELSE BEGIN
            ENDELSE
        ENDELSE


        print,rot_matrix
        points = points#rot_matrix
        stop
        ENDFOR
ENDFOR
;stop
;close,2
;device, /close

;multiplot, /reset
;multiplot, /default
;!p.multi = 0
;set_plot, 'x'

END
