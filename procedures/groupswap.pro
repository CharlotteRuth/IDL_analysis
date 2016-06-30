pro groupswap,oldbase,newbase,og,oi,ig

  if (keyword_set(og) EQ 0) then og = read_ascii_array (oldbase+'.grp',/long)
  if (keyword_set(oi) EQ 0) then oi = read_ascii_array (oldbase+'.iord',/long)
  ni = read_ascii_array(newbase+'.iord',/long)

  get_lun,lun
  openw,lun,oldbase+'.'+newbase+'.grp'

  printf,lun,strtrim(n_elements(ni),2)
  offset = 0L
  offsetdelt = 0L
  noimo=n_elements(oi-1)
  for i=0L,n_elements(ni)-1 do begin
  if ( (i mod 1e5) EQ 0) then print, i
   if(i GT noimo) then printf,lun,nsg[i-noimo] else begin
    if(oi[i-offset] NE ni[i]) then begin
      ; What follows is all the stuff to do the first time through
      ; once we figure out if we're going from late output to early or
      ; early to late
      if(offsetdelt EQ 0) then begin
        if (oi[i+1] EQ ni[i]) then begin 
          offsetdelt = -1
          if (keyword_set(ig) EQ 0) then ig = read_ascii_array(newbase+'.igasorder',/long)
        endif else begin 
          offsetdelt = 1
          if (keyword_set(ig) EQ 0) then ig = read_ascii_array(oldbase+'.igasorder',/long)
        endelse
        fs = MIN( where( ig NE 0 ))
        soig = ig[fs:n_elements(ig)-1]
        nsg = lonarr(n_elements(soig))
      endif

      ; Do real stuff
      if ( offsetdelt GT 0 ) then begin
        ta = og[fs + where(ni[i] EQ soig)]
        printf, lun,strtrim(ta[0],2)
      endif else nsg[where(oi[i] EQ soig)]=og[i-offset]
      ;printf, lun,strtrim(0,2)
      offset = offset+offsetdelt
    endif else printf,lun,strtrim(og[i-offset],2)
   endelse
  endfor
  print,"Offset:"+string(offset)
  print,"Offset Delta:"+string(offsetdelt)
  close,lun
END
