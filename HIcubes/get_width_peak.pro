function get_width_peak, vaxis, spectrum, width, v1=v1, v2=v2, $
                         peak1=peak1, peak2=peak2, $
                         widthval1=widthval1, widthval2=widthval2, $
                         doplot=doplot, gal=gal, incl=incl

;vave = total(vaxis * spectrum) / total(spectrum)
;vaxis = vaxis - vave

; find the cumulative flux:
cumul = total(spectrum, /cumulative)
cumul = float(cumul) / max(cumul)

; figure out first side - within 20% of the cumul. flux.
ind = where(cumul lt 0.2)
peak1 = max(spectrum[ind], ind2)
widthval1 = peak1 * width

ind = where(vaxis lt vaxis[ind2])

blah = min(abs(spectrum[ind] - widthval1), i1)

if ((spectrum[i1] gt widthval1) and (i1 gt 0))then begin

    m = (spectrum[i1] - spectrum[i1-1]) / (vaxis[i1] - vaxis[i1-1])
    b = spectrum[i1] - m * vaxis[i1]
    v1 = (widthval1 - b) / m

endif else if ((spectrum[i1] lt widthval1) and $
               (i1 lt n_elements(vaxis)-1)) then begin

    m = (spectrum[i1+1] - spectrum[i1]) / (vaxis[i1+1] - vaxis[i1])
    b = spectrum[i1] - m * vaxis[i1]
    v1 = (widthval1 - b) / m

endif else v1 = vaxis[i1]



; figure out second side
ind = where(cumul gt 0.8, complement=icomp)
peak2 = max(spectrum[ind], ind2)
ind2 = ind2 + n_elements(icomp)
widthval2 = peak2 * width

ind = where(vaxis gt vaxis[ind2], complement=i2comp)


blah = min(abs(spectrum[ind] - widthval2), i2)
i2 = i2 + n_elements(i2comp)

if ((spectrum[i2] gt widthval2) and (i2 gt 0)) then begin

    m = (spectrum[i2+1] - spectrum[i2]) / (vaxis[i2+1] - vaxis[i2])
    b = spectrum[i2] - m * vaxis[i2]
    v2 = (widthval2 - b) / m

endif else if ((spectrum[i2] lt widthval2) and $
               (i2 lt n_elements(vaxis)-1)) then begin

    m = (spectrum[i2-1] - spectrum[i2]) / (vaxis[i2-1] - vaxis[i2])
    b = spectrum[i2] - m * vaxis[i2]
    v2 = (widthval2 - b) / m

endif else v2 = vaxis[i2]


; make sure that there's no some weird bumps and wiggles affecting
; the spectrum for both v1 and v2 - that is, v1 must be within dv of
; vaxis[i1] and the same for v2/vaxis[i2]
dv = vaxis[1] - vaxis[0]
if abs(vaxis[i1] - v1) gt dv then begin
    v1 = vaxis[i1]
endif
if abs(vaxis[i2] - v2) gt dv then begin
    v2 = vaxis[i2]
endif

w = abs(v2 - v1)



return,w

end
