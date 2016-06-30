;
;

function hist_2d_weighted, x, y, xbins = xbins, ybins = ybins, weight = weight,$
                           nxbins = nxbins, nybins = nybins, $
                           min1 = min1, max1 = max1, $
                           min2 = min2, max2 = max2, $
                           bin1 = bin1, bin2 = bin2

IF NOT KEYWORD_SET(weight) THEN weight = findgen(N_ELEMENTS(x)) + 1.0

IF NOT KEYWORD_SET(nxbins) THEN BEGIN
    IF (KEYWORD_SET(min1) AND KEYWORD_SET(max1) AND KEYWORD_SET(bin1)) THEN nxbins = (max1 - min1)/bin1 ELSE nxbins = 10.
ENDIF
IF NOT KEYWORD_SET(nybins) THEN BEGIN
    IF (KEYWORD_SET(min2) AND KEYWORD_SET(max2) AND KEYWORD_SET(bin2)) THEN nybins = (max2 - min2)/bin2 ELSE nybins = 10.
ENDIF

IF NOT KEYWORD_SET(min1) THEN xmin = 0. ELSE xmin = min1
IF NOT KEYWORD_SET(max1) THEN xmax = 0. ELSE xmax = max1
IF NOT KEYWORD_SET(min2) THEN ymin = 0. ELSE ymin = min2
IF NOT KEYWORD_SET(max2) THEN ymax = 0. ELSE ymax = max2

IF (KEYWORD_SET(min1) AND KEYWORD_SET(max1) AND KEYWORD_SET(nxbins)) THEN BEGIN
    dx = (xmax - xmin)/nxbins
ENDIF ELSE BEGIN
    IF NOT KEYWORD_SET(bin1) THEN dx = 10. ELSE dx = bin1
ENDELSE
IF (KEYWORD_SET(min2) AND KEYWORD_SET(max2) AND KEYWORD_SET(nybins)) THEN BEGIN
    dy = (ymax - ymin)/nybins
ENDIF ELSE BEGIN
    IF NOT KEYWORD_SET(bin2) THEN dy = 10. ELSE dy = bin2
ENDELSE

IF NOT KEYWORD_SET(xbins) THEN xarray = findgen(nxbins + 1)*dx + xmin ELSE xarray = xbins
IF NOT KEYWORD_SET(ybins) THEN yarray = findgen(nybins + 1)*dy + ymin ELSE yarray = ybins
nxbins = N_ELEMENTS(xarray) - 1
nybins = N_ELEMENTS(yarray) - 1
binned = fltarr(nxbins,nybins)

FOR x_ct = 0, nxbins-1 DO BEGIN
    FOR y_ct = 0, nybins-1 DO BEGIN
        ind = where((x ge xarray(x_ct) AND x lt xarray(x_ct + 1)) AND (y ge yarray(y_ct) AND y lt yarray(y_ct + 1)))
        IF (ind[0] ne -1) THEN BEGIN
            binned[x_ct,y_ct] = TOTAL(weight[ind])/N_ELEMENTS(ind)
        ENDIF ELSE BEGIN
            binned[x_ct,y_ct] = 0
        ENDELSE
    ENDFOR
ENDFOR
RETURN,binned
END
