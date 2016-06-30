function mag,array
IF (size(array))[0] eq 1 THEN BEGIN
    mag = 0
    FOR i = 0, (size(array))[1] - 1 DO mag = mag + array[i]*array[i]
    return,sqrt(mag)
ENDIF ELSE BEGIN
    mag = fltarr((size(array))[1])
    FOR i = 0, (size(array))[2] - 1 DO mag = mag + array[*,i]*array[*,i]
    return,sqrt(mag)
ENDELSE
end
