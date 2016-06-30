pro rebin,filename1,filename2
  readcol,filename1,x,y,err
  sizu = floor ( x(size(x,/n_elements) - 1) - x(0) )
  u = fltarr(sizu)
  u = indgen(sizu)
  u = u + floor(x(0)) + .5
  newy = interpol(y,x,u)
  newerr = interpol(err,x,u)
  arr = fltarr(3,sizu)
  arr(0,*) = u
  arr(1,*) = newy
  arr(2,*) = newerr
  wascii, arr, filename2
end
