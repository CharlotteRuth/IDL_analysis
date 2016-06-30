; This code uses the output of the .star.history trace to find the progenitor halo that 
; contained most of the z=0 stars

pro haloid,simname,haloid

files = file_search(simname+'.0*/*grp*' + haloid + '*.star.history.fits')
;command =  "ls "+dir+simname+".0*/*grp*" + haloid + "*.star.history.fits | cut -d '.' -f1-7"
;command = "ls "simname+".0*/*.star.history.fits | grep star.history.fits | sed 's/.star.history.fits//g'"
;files = file_search(simname+'.0*/*.star.grp2.history.fits')
;command = "ls "+simname+".0*/*.star.grp2.history.fits | grep star.grp2.history.fits | sed 's/.star.grp2.history.fits//g'"
;spawn, command, filebase

filebase = files
endpos = strpos(files[0],'.grp')
FOR i = 0, n_elements(files) - 1 DO filebase[i] = strmid(files[i],0,endpos)
;startpos = strpos(files[0],'/')
;FOR i = 0, n_elements(files) - 1 DO filebase[i] = strmid(files[i],0,startpos)

openw, lun, simname+'.grp' + haloid + '.haloid.dat', /get_lun
;printf, lun, format='(A75, 2x, I4)', simname+'.00512/'+simname+'.00512', 1
;openw, lun, simname+'.haloid.grp2.dat', /get_lun
;printf, lun, format='(A75, 2x, I4)', simname+'.00128/'+simname+'.00128', 1
for i=n_elements(files)-1,0,-1 do begin
  test = mrdfits(files[i],1)
  if i eq n_elements(files)-1 then begin
   ind1 = where(test.haloid ne 0 and test.temp eq 0.0) 
   res = histogram(test[ind1].haloid, locations=x)
  endif else begin
   ind1 = where(test[ind2].haloid ne 0 and test[ind2].temp eq 0.0)
   res = histogram(test[ind2(ind1)].haloid, locations=x)
  endelse
   maj = x(where(res eq max(res)))
   printf, lun, format='(A75, 2x, I4)', filebase[i], maj
   ind2 = where(test.haloid eq maj[0])
  ;maj1 = x(where(res eq max(res)))
  ;maj2 = x(where(res eq max(res(where(x ne maj1[0])))))
endfor
close, lun
free_lun, lun

print, "Check and edit the file ", simname + '.grp' + haloid +'.haloid.dat', " before using it in other programs."

end

