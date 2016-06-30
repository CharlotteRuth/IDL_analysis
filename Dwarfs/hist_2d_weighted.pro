;
;

function hist_2d_weighted, x, y, weight = weight,$
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

IF NOT KEYWORD_SET(min1) THEN xmin = 50. ELSE xmin = min1
IF NOT KEYWORD_SET(max1) THEN xmax = 50. ELSE xmax = max1
IF NOT KEYWORD_SET(min2) THEN ymin = 50. ELSE ymin = min2
IF NOT KEYWORD_SET(max2) THEN ymax = 50. ELSE ymax = max2

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

xarray = findgen(nxbins + 1)*dx + xmin
yarray = findgen(nybins + 1)*dy + ymin
binned = fltarr(nxbins,nybins)
area = binned+dx*dy
FOR x_ct = 0, nxbins-1 DO BEGIN
    FOR y_ct = 0, nybins-1 DO BEGIN
        ind = where((x ge xarray(x_ct) AND x lt xarray(x_ct + 1)) AND (y ge yarray(y_ct) AND y lt yarray(y_ct + 1)))
        IF (ind[0] ne -1) THEN BEGIN
            binned[x_ct,y_ct] = TOTAL(weight[ind])/area[x_ct,y_ct]
        ENDIF ELSE BEGIN
            binned[x_ct,y_ct] = 0
        ENDELSE
    ENDFOR
ENDFOR
RETURN,binned
END
