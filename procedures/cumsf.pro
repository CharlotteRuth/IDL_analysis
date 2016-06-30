pro cumsf,s,massunit=massunit,timeunit=timeunit,overplot=overplot,normalized=normalized,_extra=_extra

  tfs = sort(s.tform)
  if (keyword_set(massunit) eq 0) then massunit = 2.325e5
  if (keyword_set(timeunit) eq 0) then timeunit = 1e9
  os = s[tfs]
  netfs = n_elements(tfs)

  smass = fltarr(netfs)
  smass[0]=os[0].mass*massunit
  for i=1L,netfs-1 do begin
    smass[i]=smass[i-1]+os[i].mass*massunit
  endfor

  if(keyword_set(normalized)) then smass = smass/(total(s.mass)*massunit)
  tform = os.tform*timeunit/1e9

  if (keyword_set(overplot)) then $
    oplot,tform,smass,_extra=_extra $
  else if (keyword_set(normalized)) then $
    plot,tform,smass,xtit="Formation time [Gyr]",ytit="Normalized Cumulative Mass [M"+sunsymbol()+"]",_extra=_extra $
  else $
    plot,tform,smass,xtit="Formation time [Gyr]",ytit="Cumulative Mass [M"+sunsymbol()+"]",_extra=_extra

end
