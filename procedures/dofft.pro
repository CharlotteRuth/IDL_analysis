function dofft,a

nea=n_elements(a)
i=findgen(nea)
; subtract off mean to elminate i=0 power
a=a-mean(a)

; Hanning window (see Numerical Methods)
; goes to 0 at beginning and end, 1 in the middle
w=0.5*(1.0-cos(2.*!pi*i/(nea-1.)))
a=a*w

; Padding to get more information about high frequencies
p=fltarr(3.*nea)
a=[a,p]

RETURN,fft(a)
END
