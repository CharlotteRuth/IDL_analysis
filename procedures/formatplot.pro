
PRO formatplot, outplot = outplot,thick = thick
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps'
    device, /helvetica
    IF KEYWORD_SET(THICK) THEN BEGIN
        !P.CHARTHICK=6
        !X.THICK=6
        !Y.THICK=6
        !p.charsize=1.5
        !x.charsize=1.25    
        !y.charsize=1.25
        !P.thick = 4
        !X.MARGIN = [14,3]
        !Y.MARGIN = [6,2]
    ENDIF ELSE BEGIN
        !P.CHARTHICK=4
        !X.THICK=4
        !Y.THICK=4
        !p.charsize=1.0
        !x.charsize=1.2         ;2.25
        !y.charsize=1.2         ;2.25
        !P.thick = 2
        !X.MARGIN = [12,3]
        !Y.MARGIN = [6,2]
    ENDELSE
    !p.font=0 
ENDIF ELSE BEGIN
    device, decomposed=0
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
ENDELSE
!p.multi = 0

end
