pro histogramp, array,binsize = binsize, input = input, locations = locations, max = max, min = min, nan = nan, nbins = nbins, omax = omax, omin = omin, l64 = L64,  REVERSE_INDICES=reverse_indices,weight = weight,overplot = overplot,normalize = normalize,cum = cum,_EXTRA=extra

IF NOT KEYWORD_SET(min) THEN min = MIN(array[where(FINITE(array))])
IF NOT KEYWORD_SET(max) THEN max = MAX(array[where(FINITE(array))])  
IF NOT KEYWORD_SET(nbins) AND max - min lt 1 THEN nbins = 100

IF KEYWORD_SET(weight) THEN $
   y = weighted_histogram(array,weight = weight,binsize = binsize,min = min,max = max,locations = locations,nbins = nbins) ELSE $
   y =          histogram(array,                binsize = binsize,min = min,max = max,locations = locations,nbins = nbins,nan = nan, omax = omax, omin = omin, L64 = L64,  REVERSE_INDICES=reverse_indices)

nbins = n_elements(y)
IF KEYWORD_SET(cum) THEN BEGIN
    cumarr=fltarr(nbins,nbins)
    ii=lindgen(nbins*nbins)
    cumarr(where(ii mod nbins le ii/nbins)) = 1
;    sfh0cum=reform(sfh0 # cumarr)
    y=reform(y # cumarr)
    ii=0
    cumarr=0
ENDIF

IF KEYWORD_SET(normalize) THEN BEGIN
   IF normalize EQ 1 THEN BEGIN
      IF KEYWORD_SET(weight) THEN norm = max(y) ELSE norm = N_ELEMENTS(array)
      y = y/FLOAT(norm)
    ENDIF ELSE y = y/normalize
ENDIF

IF KEYWORD_SET(overplot) THEN oplot,locations,y,psym = 10,_EXTRA=extra $
ELSE plot,locations,y,psym = 10,_EXTRA=extra

END
