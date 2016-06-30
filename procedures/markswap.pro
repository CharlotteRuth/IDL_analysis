pro markswap,markfile,oldbase,newbase

  rdfloat,markfile,m,skip=1,/long
  rdfloat,oldbase+'.iord',oi,skip=1,/long
  ;rdfloat,oldbase+'.igasorder',oig,skip=1,/long
  rdfloat,newbase+'.iord',ni,skip=1,/long
  rtipsy,newbase,h,/justhead

  nm = lonarr(n_elements(m))
  iords = oi[m-1] ; need to subtract 1 since iord file is in fortran
		  ; array format

  ;fs = MIN( where( oig NE 0 ))
  ;soig = oig[fs:n_elements(oig)-1]

  nni = n_elements(ni)
  nnimo = nni -1
  i=0L
  j=0L
  while(j LT nnimo) do begin
    while (iords[i] NE ni[j] AND j LT nnimo) do j=j+1
    nm[i]=j+1
    i=i+1
  endwhile

  onm = nm[where(nm NE 0)]
  get_lun,lun
  openw,lun,oldbase+'.'+newbase+'.disk.mrk'
  printf,lun,string(h.n)+string(h.ngas)+string(h.nstar)
  for i=0L,n_elements(onm)-1 do printf,lun,strtrim(nm[i],2)
  close,lun
END
