PRO plot_colorscale, x, y, scale, min = min, max = max, overplot = overplot,_EXTRA = _extra

IF NOT keyword_set(min) THEN min = MIN(scale)
IF NOT keyword_set(max) THEN max = MAX(scale)
ind = where(scale GE min AND scale LE max)
;IF (where(scale LT min))[0] NE -1  THEN scale[where(scale LT min)] = min
;IF (where(scale GT max))[0] NE -1  THEN scale[where(scale GT max)] = max
colors = (scale - min)/(max - min)*254

IF NOT keyword_set(overplot) THEN plot,x,y,/nodata,_extra = _extra, psym = 3
IF ind[0] NE -1 THEN $
  FOR i = 0L, n_elements(ind) - 1 DO oplot,[x[ind[i]],x[ind[i]]],[y[ind[i]],y[ind[i]]],color = colors[ind[i]],psym = 3 

END
