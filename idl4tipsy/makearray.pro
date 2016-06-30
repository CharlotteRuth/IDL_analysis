;filename = 'h516.cosmo25cmb.3072g14HBWK.00492'
PRO makearray, filename, halo = halo

;ext = ['OxMassFrac', 'FeMassFrac', 'HeI', 'HeII']
;ext = ['eCoolH2']
;ext = ['coolontime']
ext = ['decomp']

IF NOT keyword_set(halo) THEN halo = '1'
rtipsy, filename + '.halo.' + halo + '.std', hshort,g,d,s,/justhead
readarr,filename + '.halo.' + halo + '.iord',hshort,iordshort,/ascii,type = 'long'
rtipsy, filename,h,g,d,s,/justhead
readarr,filename + '.iord',h,iord,/ascii,type = 'long'
;iordshort = iordshort
;iord = long(iord)

match,iord,iordshort,ind,temp
print,min(ind[sort(ind)] EQ ind)
print,min(temp[sort(temp)] EQ temp)
print,hshort.n EQ n_elements(ind)
;stop

FOR i = 0, n_elements(ext) - 1 DO BEGIN
    readarr,filename + '.' + ext[i],h,extarr,/ascii;,part = 'gas'
;    array = fltarr(hshort.ngas)
;    array[0:hshort.ngas - 1] = extarr[ind]
    array = extarr[ind]

    fileout = filename + '.halo.' + halo + '.' + ext[i]
    openw,lunhi,fileout,/get_lun
    printf,lunhi,hshort.n
    FOR j = 0L, hshort.n - 1 DO $
      IF (ext[i] eq 'decomp' OR ext[i] eq 'iord' ) THEN printf,lunhi,array[j],format = '(I)'  ELSE printf,lunhi,array[j],format = '(g)'
    close,lunhi
ENDFOR

END


PRO makearray2, filename, halo = halo

ext = ['OxMassFrac', 'FeMassFrac', 'HeI', 'HeII','HI','H2']

IF NOT keyword_set(halo) THEN halo = '1'
rtipsy, filename + '.' + halo + '.std',hshort,g,d,s,/justhead
readarr,filename + '.' + halo + '.iord',hshort,iordshort,/ascii,type = 'long'
rtipsy, filename,h,g,d,s,/justhead
readarr,filename + '.iord',h,iord,/ascii,type = 'long'
;iordshort = iordshort
;iord = long(iord)

match,iord,iordshort,ind,temp
print,min(ind[sort(ind)] EQ ind)
print,min(temp[sort(temp)] EQ temp)
print,hshort.n EQ n_elements(ind)
stop
FOR i = 0, n_elements(ext) - 1 DO BEGIN
    readarr,filename + '.' + ext[i],h,extarr,/ascii;,part = 'gas'
    array = extarr[ind]

    fileout = filename + '.' + halo + '.' + ext[i]
    openw,lunhi,fileout,/get_lun
    printf,lunhi,hshort.n
    FOR j = 0L, hshort.n - 1 DO $
      IF (ext[i] eq 'decomp' OR ext[i] eq 'iord' ) THEN printf,lunhi,array[j],format = '(I)'  ELSE printf,lunhi,array[j],format = '(g)'
    close,lunhi
ENDFOR

END
