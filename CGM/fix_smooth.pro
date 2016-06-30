;Fixes the smoothing length of the file using hsmooth
PRO fix_smooth,filename,filename2 = filename2
IF NOT keyword_set(filename2) THEN filename2 = filename
rtipsy,filename + '.std',h,g,d,s
read_tipsy_arr,filename2 + '.hsm',h,smooth,part = 'gas'
;readarr,filename2 + '.smoothlength',h,smooth,/ascii,part = 'gas'
g.h = smooth
wtipsy,filename + '.std',h,g,d,s,/standard
END
