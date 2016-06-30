pro read_tipsy_arr, file, header, arr, part = part, type = type
openr, 1, file
if(keyword_set(type) eq 0) then type = 'float'

if type eq 'float' then $
   arr = fltarr(header.n) $
else $
   arr = lonarr(header.n)

readf,1,dummy
readf,1,arr
if(keyword_set(part)) then begin
    case part of
        'gas': arr = arr[0:header.ngas - 1]
;findgen(header.ngas)
        'dark': arr = arr[header.ngas:header.ngas+header.ndark - 1]
;findgen(header.ndark)+header.ngas
        'star': arr = arr[header.ngas+header.ndark:N_ELEMENTS(arr) - 1]
;findgen(header.nstar)+header.ngas+header.ndark
    endcase


 ;   arr = arr[ind]
endif

close, 1

end
