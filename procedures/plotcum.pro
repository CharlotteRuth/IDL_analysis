pro plotcum,array,binsize=sbin,OVERPLOT=overplot,_EXTRA=_extra

  plothist,array,x,h,bin=sbin,/noplot
  tot=0.0
  y=fltarr(n_elements(h))
  for i=n_elements(h)-1,0,-1 do begin
    tot=h[i]+tot
    y[i]=tot
  endfor

  if keyword_set(overplot) then oplot,x,y,_EXTRA=_extra else begin
    plot,x,y,_EXTRA=_extra
  endelse
end
