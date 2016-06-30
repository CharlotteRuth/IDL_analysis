pro plottempdist,gas=g,file=FILE,MASSUNIT=massunit,OVERPLOT=overplot,XLABEL=xlabel,YLABEL=ylabel,_EXTRA=_extra

  if(keyword_set(massunit) EQ 0) THEN massunit=13.6e16
  if(keyword_set(file)) THEN rtipsy,file,h,g,d,s
  if(keyword_set(gas) EQ 0 AND keyword_set(file) EQ 0) THEN BEGIN
    print,'Need some input:   plottempdist,[gas=g | file=filename]'
    return
  ENDIF
  sortind = reverse(sort(g.tempg))
  sorttemps = g[sortind].tempg
  sum = fltarr(n_elements(sorttemps))
  totalmass=0.0
  for i=1L,n_elements(sorttemps)-1 DO BEGIN
    sum[i] = totalmass +massunit*g[sortind[i]].mass
    totalmass= totalmass +massunit*g[sortind[i]].mass
  ENDFOR
  if (keyword_set(overplot)) then oplot,sorttemps,sum/totalmass,_EXTRA=_extra ELSE $
    if(keyword_set(xlabel)) then plot,sorttemps,sum/totalmass,xrange=[1e3,1e7],/xlog,/ylog,xtit='                                       Temperature [K]',_EXTRA=_extra ELSE $
     if(keyword_set(ylabel)) then plot,sorttemps,sum/totalmass,xrange=[1e3,1e7],/xlog,/ylog,ytit='dM(>T)/M!ltot!n                              ',_EXTRA=_extra ELSE $
      plot,sorttemps,sum/totalmass,xrange=[1e3,1e7],/xlog,/ylog,_EXTRA=_extra

END
