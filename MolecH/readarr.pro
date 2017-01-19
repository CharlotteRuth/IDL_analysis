;
;
;
;
; reads a tipsy array file
;
;
;
;



pro readarr, file, header, arr, part = part, type = type,ascii = ascii

IF KEYWORD_SET(ascii) THEN openr, 1, file ELSE openr, 1, file, /xdr

temp = 0.0

if(keyword_set(type) eq 0) then type = 'float'

if type eq 'float' then $
   arr = fltarr(header.n) $
else $
   arr = lonarr(header.n)

dummy = 0L

IF KEYWORD_SET(ascii) THEN readf,1,dummy ELSE readu,1,dummy

;byteorder, dummy

;print, dummy

temp = 't'
for i=0L, header.n-1 do begin
    IF KEYWORD_SET(ascii) THEN readf,1,temp,format = '(A)' ELSE readu,1,temp
    IF type EQ 'float' THEN $
      arr[i] = float(temp) $
    ELSE $
      arr[i] = long(temp)    
;    IF i gt 16856900 THEN BEGIN
;        print,arr[i],' ',temp
;    ENDIF
endfor


;if type eq 'float' then $
;  byteorder, arr, /xdrtof $
;else $
;  byteorder, arr

arrbk = arr
;stop

if(keyword_set(part)) then begin
    case part of
        'gas': ind = lindgen(header.ngas)
        'dark': ind = lindgen(header.ndark)+long(header.ngas)
        'star': ind = lindgen(header.nstar)+long(header.ngas)+long(header.ndark)
    endcase


    arr = arr[ind]
endif
;stop
close, 1

end


