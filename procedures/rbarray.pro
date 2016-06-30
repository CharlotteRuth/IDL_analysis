function rbarray,file,TYPE=type, TIME = time,VERBOSE = verbose,SWAP=swap,MAXSIZE=maxsize
;;; RTIPSY:  Tipsy reader for IDL
;;; Author:  James Wadsley
;;; 
if (N_PARAMS() eq 0) then begin
  print, "rbarray.pro  Reads single element binary array  files detecting the format: "
  print, "big endian, little endian, padded (standard) or non-padded header "
  print
  print, "Usage: "
  print, "        rbarray( filename [,TIME=time] [,/VERBOSE])"
  print
  print, "Input parameters: "
  print, "  filename  filename string"
  print, "  time      desired output time (optional)"
  print, "  /VERBOSE  print messages (optional)"
  print, "Return values:"
  print, "  The array"
  print, "Please read rbarray.pro for the structure definitions"
  print
  print, "Example: "
  print, "  array = rbarray( '/home/wadsley/usr5/mihos/mihos.std')"
  print, "  print, n_elements(array)"
  return,-1
endif

if ( keyword_set(maxsize) eq 0 ) then maxsize = 100000000L
if ( keyword_set(type) eq 0 ) then type='float'
;;; Note: IDL structures are never paddded 
arrsize= 0
dummy = arrsize

close,1
openr,1,file

Loop:  

;readu,1,arrsize
; Check how long long integers are (32 or 64 bit)
if ( arrsize eq 0 ) then begin
  point_lun,1,0
  ;arrsize=LONG64(0)
  arrsize=0L
  readu,1,arrsize
  verylong=1
endif else readu,1,dummy

; swap_endian if /swap is set.  Kind of clever trick to recognize when you
; want to do this 
endianswap = 0
if ( keyword_set(swap) OR ( arrsize lt 0 OR arrsize gt MAXSIZE ) ) then begin
  endianswap = 1
  arrsize=swap_endian(arrsize)
  if (keyword_set(verbose)) then print,"SWAP_ENDIAN"
endif

; create dummy array that will hold the contents of the file
CASE type OF
'float': blah=fltarr(arrsize)
'double': blah=dblarr(arrsize)
ELSE: blah=LONARR(arrsize)
ENDCASE

; read the file
readu,1,blah

if (endianswap eq 1) then blah=swap_endian(blah)

;;; Loop over output times if requested
if (keyword_set(time)) then begin
  if (abs(time-header.time) gt 1e-3) then begin
    on_ioerror, ReadError
    goto, Loop
  endif
endif

close,1
return,blah

ReadError:
print,"RTIPSY ERROR: Output time not found ",time
on_ioerror,NULL

close,1
return, -1

end
