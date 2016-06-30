PRO igasord_from_starlog, tipsyfile, iordfile, starlogfile, molecularh = molecularh
rtipsy,tipsyfile,h,g,d,s,/justhead
read_tipsy_arr,iordfile,h,iord,type = 'long'
stariord = iord[h.ngas+h.ndark:h.n - 1]
starlog = rstarlog(starlogfile)
igasorder = fltarr(h.n)
igasorder[h.ngas+h.ndark:h.n - 1] = starlog.iordergas

openw,lunhi,tipsyfile + '.igasorder',/get_lun
printf,lunhi,h.n
FOR j=0L,h.n-1 DO printf,lunhi,igasorder[j],format = '(I)' 
close,lunhi

END
