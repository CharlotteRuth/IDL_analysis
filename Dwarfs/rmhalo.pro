PRO rmHalo,filename

;This is taken from Robert Marcus's nodark IDl procedure.  It will
;read in a tipsy file of a galaxy, remove the dark matter halo and
;output a new tipsy file

  res = ['1E2R','1E3R','1E4R','1E5R','5E1R']
  mass = ['10M','11M','12M','13M','7M','8M','9M']
  for i=0,n_elements(res)-1 do begin
    for j=0,n_elements(mass)-1 do begin
     infile=res[i]+'/'+mass[j]+'/'+mass[j]+'.std'
     outfile=res[i]+'/'+mass[j]+'/'+mass[j]+'_nodark.std'
     rtipsy,infile,h,g,d,s
     d=''
     h.n=h.n-h.ndark
     h.ndark=0
     h.time=0.2
  ;   s.eps=3.25e-6
     g.h=6.50396e-6
     wtipsy,outfile,h,g,d,s,/standard
     close,1
     close,2
    endfor
  endfor
end
