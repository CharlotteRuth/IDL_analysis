PRO paperps, FILENAME=filename, charsize=charsize,thick=thick,noencaps=noencaps,_EXTRA=extra

!P.FONT=1
if (keyword_set(charsize)) then !P.CHARSIZE=charsize ELSE !P.CHARSIZE=1.3
if (keyword_set(thick)) then begin
  !P.THICK=thick 
  !P.CHARTHICK=thick
  !X.THICK=thick
  !Y.THICK=thick
ENDIF ELSE BEGIN 
  !P.THICK=4
  !P.CHARTHICK=4
  !X.THICK=4
  !Y.THICK=4
ENDELSE

set_plot,'PS'
if(keyword_set(noencaps) EQ 1) then device,filename=filename,landscape=0,SET_FONT='Times',/tt,encapsulated=0,_EXTRA=extra else $
device,filename=filename,/encapsulated,landscape=0,SET_FONT='Times',/tt,_EXTRA=extra

END
