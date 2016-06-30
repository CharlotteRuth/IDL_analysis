; like read_ascii_array only reads a long integer rather than float
FUNCTION read_lon_array,filename, FLOAT=float

  openr,lun,filename,/get_lun
  readf,lun,arrsize
  readcol,filename,arrsize,format='L',numline=1, /silent
  if keyword_set(float) then array = fltarr(arrsize) else array = lonarr(arrsize)
  readf,lun,array
  close,lun
  free_lun,lun
  return, array

END
