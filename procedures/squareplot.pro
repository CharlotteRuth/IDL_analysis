pro squareplot,filename=filename,close=close
;plots postscripts without having to set all of that crap.

;if opening, type psplot,filename='filename.ps'
;if closing, type psplot,/close

!P.THICK=1.5;4
!P.CHARTHICK=1.2;4
!X.THICK=1.5;4
!Y.THICK=1.5;4
!p.charsize=1.2;1.8
!p.font=0


if keyword_set(close) then begin
 device,/close
 set_plot,'x'

 !P.THICK=1
 !P.CHARTHICK=1
 !X.THICK=1
 !Y.THICK=1
 !p.charsize=1
 !p.font=-1

 return
endif else begin

 set_plot,'ps'
 
device,/color,bits_per_pixel=8,filename=filename,/times,ysize=5,xsize=5,/inch
 return

endelse

END
