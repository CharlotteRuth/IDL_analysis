FUNCTION read_mark_array,filename

  common ntot,ngas,nstar
  openr,lun,filename,/get_lun
  readf,lun,arrsize
  ntot=arrsize
  readf,lun,ngas
  readf,lun,nstar
  array = fltarr(arrsize)
  readf,lun,array
  close,lun
  free_lun,lun
  return, array

END
