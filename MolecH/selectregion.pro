
FUNCTION selectregion_tipsy, xvec, yvec,   $
                       XRANGEVEC=xrangevec, YRANGEVEC=yrangevec, PSYM=psym$
                       POLYFILE=polyfile, XPOLYVEC=xpolyvec, YPOLYVEC=ypolyvec, $
                       OUTSIDE = outside, color = color

;+
; ROUTINE:       selectregion
;
; USAGE:         indx = selectregion(xvec,yvec)
;                indx = selectregion(xvec,yvec,POLYFILE='poly.dat',xrangevec=[-1,4],
;                                    yrangevec=[22,28], outside=outsideindx)
;
; PURPOSE:       interactively select a region on a scatter plot, and return
;                the indices of the points that fall within the selected region
;
; INPUT: 
;   xvec         vector of quantity x, ploted on x axis of scatter plot
;   yvec         vector of quantity y, ploted on y axis of scatter plot
;
; 
; OPTIONAL KEYWORD INPUT:
;
;   xrangevec    vector indicating range in X over which to plot data
;   yrangevec    vector indicating range in Y over which to plot data
;   psym         integer number to override default PSYM=3
;
; OUTPUT:        returns indices of points that fall within selected region
;                
; OPTIONAL KEYWORD OUPUT:
;
;   POLYFILE     saves vertices of selected polygon to POLYFILE
;   XPOLYVEC     assigned to x-values of polygon vertices, if non-zero
;   YPOLYVEC     assigned to y-values of polygon vertices, if non-zero
;   OUTSIDE      assigned to indices of those data points that fall outside of region
;
; PROCEDURE:     Function begins by producing a scatter plot of the X and Y data.
;                User presses the left mouse button to select successive vertices
;                of the polygon.  The middle mouse button deletes the previous
;                point.  The right button ends the user input, and then closes the
;                polygon.
;
; EXAMPLE:       assume that a structure called "data" exists, with tags a, b, x, y
;
;          indx = selectregion(data.x, data.y)  ; User selects pts on x vs y plot.
;          trimmeddata = data[indx]             ; Data is trimmed to within region.
;          plot, trimmeddata.a, trimmeddata.b   ; a vs b properties of selected pts
;
;          indx = selectregion(data.x, data.y, outside=outindx)  
;          outsidedata = data[outindx]          ; Data is trimmed to outside of region.
;          plot, outsidedata.a, outsidedata.b   ; a vs b properties of nonselected pts
;
; AUTHOR:        J. J. Dalcanton 2006
;                Based partially on DEFROI      
;                      
;-


ON_ERROR,2                      ;Return to caller if an error occurs

; load rainbow+white color table

loadct,39

; redefine psym if needed

IF (keyword_set(PSYM)) THEN BEGIN
    psymnum = psym
ENDIF ELSE BEGIN
    psymnum = 3
ENDELSE

; plot the data points

IF (keyword_set(XRANGEVEC)) THEN BEGIN
    IF (keyword_set(YRANGEVEC)) THEN BEGIN
        plot, xvec, yvec, psym=psymnum, /ynozero, xrange=xrangevec, yrange=yrangevec, $
          /xstyle, /ystyle, /nodata
        oplot, xvec, yvec, color=255, psym=psymnum
    ENDIF ELSE BEGIN
        plot, xvec, yvec, psym=psymnum, /ynozero, xrange=xrangevec, /xstyle, /nodata
        oplot, xvec, yvec, color=255, psym=psymnum
    ENDELSE
ENDIF ELSE BEGIN
    IF (keyword_set(YRANGEVEC)) THEN BEGIN
        plot, xvec, yvec, psym=psymnum, /ynozero, yrange=yrangevec, /ystyle, /nodata
        oplot, xvec, yvec, color=255, psym=psymnum
    ENDIF ELSE BEGIN
        plot, xvec, yvec, psym=psymnum, /ynozero, /nodata
        oplot, xvec, yvec, color=255, psym=psymnum
    ENDELSE
ENDELSE

IF KEYWORD_SET(color) THEN BEGIN
    FOR i = 0, N_ELEMENTS(xvec) DO BEGIN
        oplot,[xvec[i],xvec[i]],[yvec[i],yvec[i]],psym = 3,color = color[i]
    ENDFOR
ENDIF

; give directions

print,'Press the left mouse button on begin selecting region'
print,'Press the middle mouse button to delete points'
print,'Press the right mouse button to Quit'

; set up vectors to hold ROI
maxpts = 200
xpoly = -666.0 + 0.0*findgen(maxpts)
ypoly = -666.0 + 0.0*findgen(maxpts)
npts = 0

;set color for drawing region
colorval = 250

cursor, x, y, /DATA, /DOWN

REPEAT BEGIN

    IF (!MOUSE.button EQ 1) THEN BEGIN
        
        print, npts, ' pts selected: X=',x,'  Y=',y
        xpoly[npts]=x
        ypoly[npts]=y

        IF (npts GT 0) THEN BEGIN
            plots,xpoly[npts-1:npts],ypoly[npts-1:npts], color=colorval
        ENDIF

	; reset mouse button value
        !MOUSE.button = 1

        ; increase counter
        npts = npts + 1
	
        print,'LEFT=select   MIDDLE=delete    RIGHT=quit'

    ENDIF

    ; middle button for normal mouse (2) or some microsoft compatible mice
    IF ((!MOUSE.button eq 2) or (!MOUSE.button eq 5)) THEN BEGIN
        
        print, 'Deleting last point'
        npts = npts - 1
        xpoly[npts]=0
        ypoly[npts]=0

        IF (npts GT 0) THEN BEGIN
            plots,xpoly[npts-2:npts-1],ypoly[npts-2:npts-1], color=colorval
        ENDIF

	; reset mouse button value
        !MOUSE.button = 1

        print,'LEFT=select   MIDDLE=delete    RIGHT=quit & close polygon'

    ENDIF

    cursor, x, y, /DATA, /DOWN

ENDREP UNTIL !MOUSE.button EQ 4

; close the polygon

xpoly[npts] = xpoly[0]
ypoly[npts] = ypoly[0]
plots,xpoly[npts-1:npts],ypoly[npts-1:npts],color=colorval

; trim the vector

xpoly=xpoly[0:npts]
ypoly=ypoly[0:npts]

; output the final corners of the polygon
FOR i = 0,npts DO BEGIN
    print, xpoly[i], ypoly[i]
ENDFOR

; Save polygon info, if needed
IF keyword_set(POLYFILE) THEN BEGIN
    print,'Saving polygons to ',POLYFILE
    n = n_elements(xpoly)
    OPENW, 1, POLYFILE
    FOR i=0,n-1 DO BEGIN
        printf, 1, xpoly[i], ypoly[i]
    ENDFOR
    CLOSE, 1
ENDIF

; pass polygon vertices back, if requested.
IF (keyword_set(XPOLYVEC)) THEN xpolyvec=xpoly
IF (keyword_set(YPOLYVEC)) THEN ypolyvec=ypoly

; calculate indices of points that are within the defined polygon
indx = inside(xvec, yvec, xpoly, ypoly)
oplot, xvec[indx], yvec[indx], psym=psymnum, color=colorval

; if requested, also return the indices that are _not_ inside the polygon.
IF (keyword_set(OUTSIDE)) THEN BEGIN
    junk = 0 * indgen(n_elements(xvec))
    junk[indx] = 1
    OUTSIDE = where(junk NE 1)
ENDIF

RETURN, indx

END


