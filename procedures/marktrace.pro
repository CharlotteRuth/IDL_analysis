pro marktrace,oldbase,newbase

  rdfloat,oldbase+'.mrk',marked,skip=1
  openr,lun,oldbase+'.mrk',/get_lun
  readf,lun,numotot
  readf,lun,numogas
  readf,lun,numostars
  close,lun
  numdark=2806399
  oio = read_ascii_array(oldbase+'.iord')
  nio = read_ascii_array(newbase+'.iord')
  nigo = read_ascii_array(newbase+'.igasorder')
  openw,lun,newbase+'.mrk',/get_lun
  numntot=strtrim(n_elements(nio),2)
  numnstars=n_elements(where(nigo NE 0))
  numngas = numntot-numnstars -numdark
  printf,lun,strtrim(numntot,2)+' '+string(numngas,format='(I)')+' '+strtrim(numnstars,2)
  iords = oio[marked]
  j=0L
  offset=0L
  for i=0L,n_elements(nio)-1 do begin
    if(iords[j] EQ nio[i]) then begin
      printf,lun,strtrim(i,2)
      j=j+1
    endif
    if(nio[i] NE i+offset) then begin
      if(iords[j] EQ i+offset) then begin
        endnewmark=[endnewmark,where(iords[j] EQ nigo[numngas+numdark+numostars:numntot])]
        j=j+1
      endif
      offset = offset+1
    endif
    if(j EQ n_elements(iords)-1) then break
  endfor
  if(n_elements(endnewmark) GT 0) then begin
    sortedend=endnewmark[sort(endnewmark)]
    for i=0,n_elements(endnewmark)-1 do printf,lun,strtrim(sortedend[i],2)
  endif

END
