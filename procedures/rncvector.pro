pro rncvector,file, header, blah, TIME = time,VERBOSE = verbose,MAXSIZE=maxsize
;;; RNCARRAY:  NChilada Array reader for IDL
;;; Author:  Greg Stinson
;;; 
if (N_PARAMS() eq 0) then begin
  print, "rncvector.pro  Reads single element binary array  files: "
  print
  print, "Usage: "
  print, "        rncvector( filename [,TIME=time] [,/VERBOSE])"
  print
  print, "Input parameters: "
  print, "  filename  filename string"
  print, "  time      desired output time (optional)"
  print, "  /VERBOSE  print messages (optional)"
  print, "Return values:"
  print, "  The array"
  print, "Please read rncvector.pro for the structure definitions"
  print
  print, "Example: "
  print, "  array = rncvector( 'mw_105k.00100/star/pos')"
  print, "  print, n_elements(array)"
  return
endif

if ( keyword_set(maxsize) eq 0 ) then maxsize = 100000000L
;;; Note: IDL structures are never paddded 
header = { magic:0L, time:double(0.0), n:LONG64(0), ndim:0L, code:0L, min:fltarr(3), max:fltarr(3) }

close,1
openr,1,file

Loop:  

readu,1,header
endianswap = 0
if (header.ndim lt 1 or header.ndim gt 3) then begin
  endianswap = 1
  header = swap_endian(header)
  if (keyword_set(verbose)) then print,"SWAP_ENDIAN"
endif

; create dummy array that will hold the contents of the file
CASE header.code OF
 3: blah=intarr(header.ndim,header.n)
 5: blah=lonarr(header.ndim,header.n)
 9: blah=fltarr(header.ndim,header.n)
 10: blah=dblarr(header.ndim,header.n)
 ELSE: print,"Couldn't read type code."
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
return

ReadError:
print,"RTIPSY ERROR: Output time not found ",time
on_ioerror,NULL

close,1
return

end
