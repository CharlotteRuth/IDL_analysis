PRO CONTOUR_PLUS, x,y,weight = weight,XBINSIZE=xbinsize,YBINSIZE=ybinsize, $
                  XMIN=xmin,YMIN=ymin,XMAX=xmax,YMAX=ymax,LEVELS=levels,$
                  THRESHOLD=threshold,REVERSE=reverse,NLEVELS=nlevels, LOGLEVEL=loglevel,$
                  LENSCALE = LENSCALE,C_COLOR = C_COLOR, verbose = verbose,_EXTRA=EXTRA

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


IF NOT (KEYWORD_SET(xbinsize)) THEN xbinsize=.1
IF NOT (KEYWORD_SET(ybinsize)) THEN ybinsize=.1
IF NOT (KEYWORD_SET(threshold)) THEN threshold=100
IF NOT (KEYWORD_SET(xmin)) THEN xmin=MIN(x)
IF NOT (KEYWORD_SET(ymin)) THEN ymin=MIN(y)
IF NOT (KEYWORD_SET(xmax)) THEN xmax=MAX(x) + 0.1
IF NOT (KEYWORD_SET(ymax)) THEN ymax=MAX(y) - 0.1
IF NOT (KEYWORD_SET(weight)) THEN weigth = FLTARR(N_ELEMENTS(x)) + 1.0

IF (KEYWORD_SET(REVERSE)) THEN yrange=[ymax,ymin] ELSE yrange=[ymin,ymax]

nxbins=FLOOR((xmax-xmin) / xbinsize)
nybins=FLOOR((ymax-ymin) / ybinsize)
;xbins=(FINDGEN(nxbins)*xbinsize)+xmin+xbinsize/2.
;ybins=(FINDGEN(nybins)*ybinsize)+ymin+ybinsize/2.
xbins=FINDGEN(nxbins + 1)*(xmax - xmin)/nxbins + xmin
ybins=FINDGEN(nybins + 1)*(ymax - ymin)/nybins + ymin
xbins_plot=(xbins[0:nxbins - 1] + xbins[1:nxbins])/2
ybins_plot=(ybins[0:nybins - 1] + ybins[1:nybins])/2

hist=HIST_2D_weighted(x,y,nxbins = nxbins, nybins = nybins, weight = weight,bin1=xbinsize,bin2=ybinsize,min1=xmin,min2=ymin,max1=xmax,max2=ymax)

IF KEYWORD_SET(verbose) THEN BEGIN
    window,1
    yh = histogram(hist[where(hist ne 0)],locations = xh,nbins = 100)
    plot,xh,yh,psym = 10
    window,0
ENDIF

IF (KEYWORD_SET(LOGLEVEL)) THEN BEGIN
    hist = ALOG10(hist)
    threshold = ALOG10(threshold)
    highest = MAX(hist)
ENDIF ELSE highest = MAX(hist)

IF NOT (KEYWORD_SET(nlevels)) THEN BEGIN
    IF (KEYWORD_SET(levels)) THEN nlevels  = N_ELEMENTS(levels) ELSE nlevels = highest/threshold
ENDIF
IF NOT (KEYWORD_SET(levels)) THEN levels=FINDGEN(nlevels)*(highest - threshold)/nlevels + threshold;threshold*FINDGEN(nlevels)+threshold
IF NOT (KEYWORD_SET(C_COLOR)) THEN C_COLOR = (INDGEN(nlevels)+1)*255/nlevels
;ind = where(hist eq 0)
;hist[ind] = threshold - (highest - threshold)/nlevels
basecolor = C_COLOR[0]/2

PLOTSYM, 0, .2, /FILL
IF (MAX(hist) LT threshold*1.5) THEN BEGIN
    plot,x,y,xrange=[xmin,xmax],yrange=yrange,psym=8,_EXTRA=EXTRA
ENDIF ELSE BEGIN
;now create array for indices of stars that fall outside of contours
;    ind=[0l]
;    FOR i=0l,N_ELEMENTS(x)-1 DO BEGIN
;        mindistx=MIN(ABS(xbins-x[i]),xminpos)
;        mindisty=MIN(ABS(ybins-y[i]),yminpos)
;        IF xminpos ge 0 and xminpos lt N_ELEMENTS(xbins) AND yminpos ge 0 and yminpos lt N_ELEMENTS(ybins) THEN ind=[ind,i] ELSE $
;          IF (hist[xminpos,yminpos] LT threshold) THEN ind=[ind,i]
 ;   ENDFOR
 ;   ind=ind[1:*]
;    plot,x[ind],y[ind],xrange=[xmin,xmax],yrange=yrange,psym=8,_EXTRA=extra
    plot,x,y,xrange=[xmin,xmax],yrange=yrange,psym=8,_EXTRA=EXTRA

;    plot,x,y,xrange=[xmin,xmax],yrange=yrange,psym=8,_EXTRA=EXTRA
;    oplot,x,y,psym=8,color = basecolor
 
;    contour,hist,xbins_plot,ybins_plot,levels=levels,/FILL,C_COLOR = C_COLOR,/overplot,_EXTRA=EXTRA,closed = 1,/cell_fill ;,min_value = threshold
    contour,hist,xbins_plot,ybins_plot,yrange=yrange,levels=levels,/overplot,/FILL,_EXTRA=extra,C_COLOR=REPLICATE(!P.BACKGROUND,nlevels)
    contour,hist,xbins_plot,ybins_plot,yrange=yrange,levels=levels,/overplot
stop
;    contour,hist,xbins[0:N_ELEMENTS(xbins) - 2],ybins[0:N_ELEMENTS(ybins) - 2],levels=levels,/FILL,C_COLOR = C_COLOR,/overplot,_EXTRA=EXTRA
;yrange=yrange,levels=levels,/overplot,/FILL,_EXTRA=extra,C_COLOR = (INDGEN(nlevels)+1)*240/nlevels ;,C_COLOR=REPLICATE(!P.BACKGROUND,nlevels)
;     contour,hist,xbins[0:N_ELEMENTS(xbins) - 2],ybins[0:N_ELEMENTS(ybins) - 2],levels=levels,/overplot;yrange=yrange,levels=levels,/overplot
;,color=!P.BACKGROUND
ENDELSE
END
