;tfile = 'h799.cosmo25cmb.3072g14HBWK.00512'
;rtipsy,tfile,h,g,d,s
;array = [g.mass,d.mass,s.mass]
;writetipsyarray,tfile+'.mass',array

PRO writetipsyarray,filename,array,type = type
IF NOT keyword_set(type) THEN type = "float"
openw,lunhi,filename,/get_lun
printf,lunhi,n_elements(array)
for j=0LL,n_elements(array)-1 do begin
    IF type EQ "int" THEN printf,lunhi,array[j],format = '(I)'  ELSE printf,lunhi,array[j],format = '(g)'
endfor
close,lunhi
END
;16756355
