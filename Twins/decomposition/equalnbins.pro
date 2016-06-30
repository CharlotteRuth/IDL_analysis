;
;
;
; function to return bins that will contain equal numbers of input
; particles
;
;



function equalnbins, p, a, n


nbins = floor(n_elements(p)/n)

print, 'nbins = ', nbins

bins = fltarr(nbins)

sort = sort(a)

sorta = a[sort]

for j=2L,nbins do bins[j-1]=sorta[j*n-1]

return, bins

end

