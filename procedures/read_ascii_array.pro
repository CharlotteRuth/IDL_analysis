;******************************************************************************
;This functions reads in ascii arrays output by the std_array_convert
;program.  The default output is float, but /DOUBLE or /LONG will output
;arrays of the appropriate type.
;******************************************************************************
FUNCTION read_ascii_array, filename,DOUBLE=double,LONG=long,int = int
numlines = 0L
openr, lun, filename,/get_lun
readf, lun, numlines,format='(I)'
;stop
;numlines=LONG(numlines)
;stop
;Get the array into the proper format
IF KEYWORD_SET(double) THEN BEGIN
    array=DBLARR(numlines) 
ENDIF ELSE BEGIN
    IF KEYWORD_SET(long) THEN BEGIN
        array=LONARR(numlines)
    ENDIF ELSE BEGIN
        IF KEYWORD_SET(INT) THEN BEGIN
            array=INTARR(numlines)
        ENDIF ELSE BEGIN
            array=FLTARR(numlines)
        ENDELSE
    ENDELSE
ENDELSE
;stop
readf, lun, array
close, lun
free_lun, lun

return, array

END

