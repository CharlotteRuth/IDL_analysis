pro fill_between, xx, y1,y2 , _extra=e, n_big=n_big
;+
; NAME:  fill_between
;
;
;
; PURPOSE:  shade in between two curves
;
;
;
; CATEGORY:  plotting
;
;
;
; CALLING SEQUENCE:
;
;
;
; INPUTS:  x, y1,y2
;           ;if all else fails, just switch y1 and y2
;
;
; OPTIONAL INPUTS:
;              Use standard contour keywords to define how the region
; is to be filled.  Good ones include:
;            /fill :  if you just want a solid color, use this
;            c_color:  set 0 to 255 for a nice color
;            c_orientation :If filling with lines, set the angle to
; lines make
;            c_linestyle:  linestyle
;            c_thick:  how thick to make the lines
;            /cell_fill:  Looks like /fill is a better choice.
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:
;
;
;
; OPTIONAL OUTPUTS:
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;-


;generate a 2-d array
x=xx
  if not keyword_set(n_big) then n_big=1000. ;keyword to make sure the
                                ;contour isn't too pixelated
z=fltarr(n_elements(x), n_big);n_elements(y))

;generate a y-array
spread=max(y2)-min(y1)
y=findgen(n_big)/n_big*spread+min(y1)

;gx=x(where(y1 lt y2))
;gy=y(where(y gt y1 and y lt y2))

;z(gx,gy)=1

for i=0L,n_elements(x)-1 do begin
;  if y1(i) lt y2(i) then begin
    if max(where(y ge y1(i) and y le y2(i))) ne -1 then begin
      z(i,(where(y ge y1(i) and y le y2(i))))=1
    endif
    if max(where(y ge y2(i) and y le y1(i))) ne -1 then begin
      z(i,(where(y ge y2(i) and y le y1(i))))=1
    endif
;  endif
endfor





;stop

contour, z, x, y, _extra=e

;stop
end
