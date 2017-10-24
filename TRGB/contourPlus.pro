PRO CONTOUR_PLUS, x,y,EXPORT=export,XBINSIZE=xbinsize,YBINSIZE=ybinsize,XMIN=xmin,YMIN=ymin,XMAX=xmax,YMAX=ymax,LEVELS=levels,THRESHOLD=threshold,REVERSE=reverse,LOG=LOG, NOLINES=NOLINES,NLEVELS=nlevels,nofill = nofill,_EXTRA=extra

; PROGRAM: CONTOUR_PLUS
;
; WRITTEN: Anil Seth, 2005
;
; PURPOSE: Created for use in color-magnitude diagrams, contour_plus is
; a program which makes a contour plot above a certain data density
; threshold and plot individual data points below that threshold.
;
; ARGUMENTS:
;           X --> data values for the xaxis (e.g. V-I)
;           Y --> data values for the yaxis (e.g. V)
;
; OPTIONAL KEYWORDS:
;           XBINSIZE --> size of bins in X direction
;           YBINSIZE --> size of bins in Y direction
;           XMIN,XMAX,YMIN,YMAX --> boundaries of the plot/binning
;           LEVELS --> density levels (points/bin) at which to draw contours
;             e.g. levels=[75,100,150,200,250,350,500,750,1000,1500,2000]
;           THRESHOLD --> density level above which to not draw points,
;              should be higher than lowest contour
;           REVERSE --> reverse the yaxis when plotting
;           _EXTRA --> applied to both PLOT and CONTOUR commands, can
;             include keywords such as C_COLOR, XMARGIN, etc.
;   MODIFICATIONS: Uses USERSYM rather than PSYM to avoid PDFing problems
;
; OPTIONAL OUTPUTS:
;           EXPORT --> RETURNS HISTOGRAM VALUES

IF NOT (KEYWORD_SET(xbinsize)) THEN xbinsize=.1
IF NOT (KEYWORD_SET(ybinsize)) THEN ybinsize=.1
IF NOT (KEYWORD_SET(threshold)) THEN threshold=100
IF NOT (KEYWORD_SET(xmin)) THEN xmin=MIN(x)
IF NOT (KEYWORD_SET(ymin)) THEN ymin=MIN(y)
IF NOT (KEYWORD_SET(xmax)) THEN xmax=MAX(x)
IF NOT (KEYWORD_SET(ymax)) THEN ymax=MAX(y)

IF (KEYWORD_SET(REVERSE)) THEN yrange=[ymax,ymin] ELSE yrange=[ymin,ymax]

nxbins=FLOOR((xmax-xmin) / xbinsize) + 1L
nybins=FLOOR((ymax-ymin) / ybinsize) + 1L
xbins=(FINDGEN(nxbins)*xbinsize)+xmin+xbinsize/2.
ybins=(FINDGEN(nybins)*ybinsize)+ymin+ybinsize/2.

hist=HIST_2D(x,y,bin1=xbinsize,bin2=ybinsize,min1=xmin,min2=ymin,max1=xmax,max2=ymax,_EXTRA=extra)
export=hist
PLOTSYM, 0, .1, /FILL

IF (MAX(hist) LT threshold*1.5) THEN BEGIN
    plot,x,y,xrange=[xmin,xmax],yrange=yrange,psym=8,_EXTRA=extra
;   COMMENTED TO GET USERSYM RATHER THAN PSYM=3
;    plot,x,y,xrange=[xmin,xmax],yrange=yrange,psym=3,_EXTRA=extra
ENDIF ELSE BEGIN
;now create array for indices of stars that fall outside of contours
    ind=[0l]
    FOR i=0l,N_ELEMENTS(x)-1 DO BEGIN
        mindistx=MIN(ABS(xbins-x[i]),xminpos)
        mindisty=MIN(ABS(ybins-y[i]),yminpos)
        IF (hist[xminpos,yminpos] LT threshold) THEN ind=[ind,i]
    ENDFOR
    ind=ind[1:*]
    plot,x[ind],y[ind],xrange=[xmin,xmax],yrange=yrange,psym=8,_EXTRA=extra
;   COMMENTED TO GET USERSYM RATHER THAN PSYM=3
;   plot,x,y,xrange=[xmin,xmax],yrange=yrange,psym=3,_EXTRA=extra

    IF NOT (KEYWORD_SET(nlevels)) THEN nlevels=MAX(hist)/threshold
    IF NOT (KEYWORD_SET(levels)) THEN levels=threshold*FINDGEN(nlevels)+threshold
    IF (KEYWORD_SET(log)) THEN levels=10.^(ALOG10(MAX(hist))-FINDGEN(nlevels))
    IF (KEYWORD_SET(nolines)) THEN contour,hist,xbins,ybins,yrange=yrange,levels=levels,/overplot,/FILL,_EXTRA=extra ,color=!P.BACKGROUND ELSE BEGIN
        IF (KEYWORD_SET(nofill)) THEN contour,hist,xbins,ybins,yrange=yrange,levels=levels,_EXTRA=extra,/overplot ELSE BEGIN
        contour,hist,xbins,ybins,yrange=yrange,levels=levels,_EXTRA=extra,/overplot,/FILL,color=!P.BACKGROUND ;,C_COLOR=REPLICATE(!P.BACKGROUND,nlevels)
        contour,hist,xbins,ybins,yrange=yrange,levels=levels,/overplot,color=!P.BACKGROUND
        ENDELSE
    ENDELSE

ENDELSE


END
