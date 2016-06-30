pro rcsm,filename1,filename2,width
  readcol,filename1,x,y,z
  w=smooth(y,width)
  arr = fltarr(3,size(x,/n_elements))
  arr(0,*) = x
  arr(1,*) = w
  arr(2,*) = z
  wascii, arr, filename2
end
