FUNCTION z_to_t,z,uniqz = uniqz;,uniqt,z
IF NOT keyword_set(uniqz) THEN uniqz = z[uniq(z,sort(z))]
ageUniverse = wmap3_lookback(1000)
uniqt = (ageUniverse - wmap3_lookback(uniqz))/1e9
t = fltarr(n_elements(z))
FOR i = 0, n_elements(uniqz) - 1 DO t[where(z EQ uniqz[i])] = uniqt[i] 


;sortind = sort(uniqz)
;uniqz = uniqz[sortind]
;uniqt = uniqt[sortind]

;uniqz_fit = z[uniq(z,sort(z))]
;uniqt_fit = spline(uniqz,uniqt,uniqz_fit)
;match2,z,uniqz_fit,ind1,ind2
;t = fltarr(n_elements(z))
;t = uniqt_fit[ind1]

RETURN,t
END
