; This code uses the output of the .star.history trace to find the progenitor halo that 
; contained most of the z=0 stars

PRO haloid, filebase, finalid = finalid
IF NOT keyword_set(finalid) THEN finalid = '1'

files = file_search(filebase+'.0*/*grp' + finalid + '.star.history.fits')
;command =  "ls "+dir+filebase+".0*/*grp*" + finalid + "*.star.history.fits | cut -d '.' -f1-7"
;command = "ls "filebase+".0*/*.star.history.fits | grep star.history.fits | sed 's/.star.history.fits//g'"
;files = file_search(filebase+'.0*/*.star.grp2.history.fits')
;command = "ls "+filebase+".0*/*.star.grp2.history.fits | grep star.grp2.history.fits | sed 's/.star.grp2.history.fits//g'"
;spawn, command, filename
filename = files
endpos = strpos(files[0],'.grp')
FOR i = 0, n_elements(files) - 1 DO filename[i] = strmid(files[i],0,endpos)
;startpos = strpos(files[0],'/')
;FOR i = 0, n_elements(files) - 1 DO filename[i] = strmid(files[i],0,startpos)

openw, lun, filebase+'.grp' + finalid + '.haloid.dat', /get_lun
;printf, lun, format='(A75, 2x, I4)', filebase+'.00512/'+filebase+'.00512', 1
;openw, lun, filebase+'.haloid.grp2.dat', /get_lun
;printf, lun, format='(A75, 2x, I4)', filebase+'.00128/'+filebase+'.00128', 1
for i=n_elements(files)-1,0,-1 do begin
   test = mrdfits(files[i],1)
   if i eq n_elements(files)-1 then begin
      ind1 = where(test.haloid ne 0 and test.temp eq 0.0) 
      IF ind1[0] EQ -1 THEN return
      IF n_elements(ind1 EQ 1) THEN BEGIN
          res = [1]
          x = [0]
      ENDIF ELSE  res = histogram(test[ind1].haloid, locations=x)
   endif else begin
      ind1 = where(test[ind2].haloid ne 0 and test[ind2].temp eq 0.0)
      IF n_elements(ind1) GT 1 THEN res = histogram(test[ind2(ind1)].haloid, locations=x) ELSE BEGIN
         res = [1]
         IF n_elements(ind1) EQ 1 THEN x = [test[ind2(ind1)].haloid] ELSE x = [0]
      ENDELSE
   endelse
   maj = x(where(res eq max(res)))
   printf, lun, format='(A84, 2x, I4)', filename[i], maj
   ind2 = where(test.haloid eq maj[0])
  ;maj1 = x(where(res eq max(res)))
  ;maj2 = x(where(res eq max(res(where(x ne maj1[0])))))
endfor
close, lun
free_lun, lun

print, "Check and edit the file ", filebase + '.grp' + finalid +'.haloid.dat', " before using it in other programs."

end

