function weighted_histogram,y,weight = weight,binsize=binsize,min = min,max = max,locations = locations,nbins = nbins, cum = cum
;if not keyword_set(min) AND min ne 0 then min = MIN(y)
;if not keyword_set(max) AND max ne 0 then max = MAX(y)
if not keyword_set(min) then min = MIN(y)
if not keyword_set(max) then max = MAX(y)
if not keyword_set(weight) then weight = fltarr(N_ELEMENTS(y)) + 1.0

if not keyword_set(binsize) then BEGIN
    if not keyword_set(nbins) then nbins = 100
    binsize = (max - min)/double(nbins)
endif else nbins = round((max - min)/double(binsize))
x = findgen(nbins + 1)*binsize + min

hist = fltarr(N_ELEMENTS(x)-1)
for i = 0L,N_ELEMENTS(hist)-1 do begin
    ind = where(y gt x[i] AND y lt x[i+1])
    if ind[0] ne -1 then hist[i] = TOTAL(weight[ind]) else hist[i] = 0
 endfor
IF keyword_set(cum) THEN BEGIN
    cumarr = fltarr(nbins,nbins)
    ii = lindgen((long(nbins))*(long(nbins)))
    cumarr(where(ii MOD nbins LE ii/nbins)) = 1
    hist = reform(hist # cumarr)
    ii = 0
    cumarr = 0
ENDIF
locations = (x[0:N_ELEMENTS(x) - 2] + x[1:N_ELEMENTS(x) - 1])/2
return,hist
end
