pro rncarray,file,header,blah, TIME = time,VERBOSE = verbose,MAXSIZE=maxsize
;;; RNCARRAY:  NChilada Array reader for IDL
;;; Author:  Greg Stinson
;;; 
if (N_PARAMS() eq 0) then begin
  print, "rncarray.pro  Reads single element binary array  files: "
  print
  print, "Usage: "
  print, "        rncarray( filename [,TIME=time] [,/VERBOSE])"
  print
  print, "Input parameters: "
  print, "  filename  filename string"
  print, "  time      desired output time (optional)"
  print, "  /VERBOSE  print messages (optional)"
  print, "Return values:"
  print, "  The array"
  print, "Please read rncarray.pro for the structure definitions"
  print
  print, "Example: "
  print, "  array = rncarray( 'mw_105k.00100/star/timeform')"
  print, "  print, n_elements(array)"
  return
endif

if ( keyword_set(maxsize) eq 0 ) then maxsize = 100000000L
;;; Note: IDL structures are never paddded 
header = { magic:0L, time:double(0.0), n:long64(0), ndim:0L, code:0L, min:float(0.0), max:float(0.0) }

close,1
openr,1,file

readu,1,header
endianswap = 0
if (header.magic ne 1062053) then begin
  endianswap = 1
  header = swap_endian(header)
  if (keyword_set(verbose)) then print,"SWAP_ENDIAN"
endif

if (header.code ne 9) then begin
CASE header.code OF
 3: header = { magic:0L, time:double(0.0), n:long64(0), ndim:0L, code:0L, min:0, max:0 }
 5: header = { magic:0L, time:double(0.0), n:long64(0), ndim:0L, code:0L, min:0L, max:0L }
 6: header = { magic:0L, time:double(0.0), n:long64(0), ndim:0L, code:0L, min:0L, max:0L }
 10: header = { magic:0L, time:double(0.0), n:long64(0), ndim:0L, code:0L, min:double(0.0), max:double(0.0) }
 ELSE: print,"Couldn't read type code."
ENDCASE
point_lun,1,0
readu,1,header
if (endianswap EQ 1) then header=swap_endian(header)
ENDIF
; create dummy array that will hold the contents of the file
CASE header.code OF
 3: blah=intarr(header.n)
 5: blah=lonarr(header.n)
 6: blah=lonarr(header.n)
 9: blah=fltarr(header.n)
 10: blah=dblarr(header.n)
 ELSE: print,"Couldn't read type code."
ENDCASE

; read the file
readu,1,blah

if (endianswap eq 1) then blah=swap_endian(blah)

close,1

end
