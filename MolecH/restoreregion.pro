
FUNCTION restoreregion, xvec, yvec,   $
                       POLYFILE=polyfile, XPOLYVEC=xpolyvec, YPOLYVEC=ypolyvec, $
                       OUTSIDE = outside

;+
; ROUTINE:       restoreregion
;
; USAGE:         indx = restoreregion(xvec,yvec,POLYFILE='poly.dat')
;                indx = restoreregion(xvec,yvec,XPOLYVEC=xvec,YPOLYVEC=yvec)
;
; PURPOSE:       Recalculate the indices of the points that fall 
;                within a previously selected region
;
; INPUT: 
;   xvec         vector of quantity x
;   yvec         vector of quantity y
;
; OPTIONAL KEYWORD INPUT:
;
;   POLYFILE     file containing x-y coordinates of defined polygon
;
;                OR
;
;   XPOLYVEC     x-values of polygon vertices
;   YPOLYVEC     y-values of polygon vertices
;
; OUTPUT:        returns indices of points that fall within selected region
;                
; OPTIONAL KEYWORD OUPUT:
;
;   OUTSIDE      assigned to indices of those data points that fall outside of region
;
; PROCEDURE:     Assumes user has previously defined a polygon in x vs y
;                (possibly by running 'selectregion.pro').  
;
; EXAMPLE:       assume that a structure called "data" exists, with tags a, b, x, y
;
;          indx = restoreregion(data.x, data.y, POLYFILE='poly.dat')
;          trimmeddata = data[indx]             ; Data is trimmed to within region.
;          plot, trimmeddata.a, trimmeddata.b   ; a vs b properties of selected pts
;
;                OR
;
;          indx = restoreregion(data.x, data.y, XPOLYVEC=xvec, YPOLYVEC=yvec)  
;          trimmeddata = data[indx]             ; Data is trimmed to within region.
;          plot, trimmeddata.a, trimmeddata.b   ; a vs b properties of selected pts
;
; AUTHOR:        J. J. Dalcanton 2006
;                Basically a wrapper for 'inside'
;                      
;-


ON_ERROR,2                      ;Return to caller if an error occurs

IF (keyword_set(POLYFILE)) THEN BEGIN
    readcol, polyfile, xpolyvec, ypolyvec, FORMAT='F,F', /silent
    print, 'Read ', n_elements(xpolyvec), ' vertices from ', polyfile
ENDIF ELSE BEGIN
    IF (keyword_set(XPOLYVEC) AND keyword_set(YPOLYVEC)) THEN BEGIN
        print, 'Using ', n_elements(xpolyvec), ' vertices from user-supplied vectors'
    ENDIF ELSE BEGIN
        print, 'Must define polynomial vertices using POLYFILE or POLYVEC keywords'
    ENDELSE
ENDELSE

; return indices of points that are within the defined polygon
indx = inside(xvec, yvec, xpolyvec, ypolyvec)

; if requested, also return the indices that are _not_ inside the polygon.
IF (keyword_set(OUTSIDE)) THEN BEGIN
    junk = 0 * indgen(n_elements(xvec))
    junk[indx] = 1
    OUTSIDE = where(junk NE 1)
ENDIF

RETURN, indx

END


