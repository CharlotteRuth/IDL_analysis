pro rheader,file,header,TIME = time,VERBOSE = verbose
;;; RTIPSY:  Tipsy reader for IDL
;;; Author:  James Wadsley
;;; 
if (N_PARAMS() eq 0) then begin
  print, "rtipsy.pro  Reads tipsy files detecting the format: "
  print, "big endian, little endian, padded (standard) or non-padded header "
  print
  print, "Usage: "
  print, "        rtipsy, filename ,header [,g] [,d] [,s] [,TIME=time] [,/VERBOSE]"
  print
  print, "Input parameters: "
  print, "  filename  filename string"
  print, "  time      desired output time (optional)"
  print, "  /VERBOSE  print messages (optional)"
  print, "Return values:"
  print, "  header    tipsy header struct"
  print, "  g,d,s     gas, dark and star structures"
  print, "Please read rtipsy.pro for the structure definitions"
  print
  print, "Example: "
  print, "  rtipsy, '/home/wadsley/usr5/mihos/mihos.std',h,g,d"
  print, "  print, h.ndark"
  print, "  plot, d.x, d.y, psym=3"
  return
endif

;;; Note: IDL structures are never paddded 
header = { time:double(0.0), n:0L, ndim:0L, ngas:0L, ndark:0L, nstar:0L }

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

if (keyword_set(verbose)) then print,"Read time,n,ngas,nstar,ndark: ",header.time,header.n,header.ngas,header.ndark,header.nstar

fs = fstat(1)
;;; Explicitly pad header if required 
;if (fs.size eq 32UL+header.ngas*48+header.ndark*36+header.nstar*44) then begin
if (fs.size eq 32LL+header.ngas*48LL+header.ndark*36LL+header.nstar*44LL) then begin
  dummy = 1L 
  readu,1,dummy
;endif else if (fs.size ne 28UL+header.ngas*48+header.ndark*36+header.nstar*44) then begin  
endif else if (fs.size ne 28LL+header.ngas*48LL+header.ndark*36LL+header.nstar*44LL) then begin  
  print, "RTIPSY ERROR: Header and file size inconsistent"
  print, "Estimates: Header bytes:  28 or 32 (either is OK)"
  print, "     ngas: ",header.ngas," bytes:",48LL*header.ngas
  print, "    ndark: ",header.ndark," bytes:",36LL*header.ndark
  print, "    nstar: ",header.nstar," bytes:",44LL*header.nstar
;  print, "Actual File bytes:",fs.size,"  not one of:",32UL+header.ngas*48+header.ndark*36+header.nstar*44,28UL+header.ngas*48+header.ndark*36+header.nstar*44
  print, "Actual File bytes:",fs.size,"  not one of:",32LL+header.ngas*48LL+header.ndark*36LL+header.nstar*44LL,28LL+header.ngas*48LL+header.ndark*36LL+header.nstar*44LL
  close,1
  return
endif

close,1
return

ReadError:
print,"RTIPSY ERROR: Output time not found ",time
on_ioerror,NULL

close,1
return

end
