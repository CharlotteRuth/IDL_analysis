pro temphist
!ORDER=1
  ;prefix="HSBc17s0.045V45."
  prefix="cosdwf."
  iOutputInt=1
  for i=1,256 do begin
    if (i*iOutputInt LT 10) then begin
      infile = prefix + '0000' + strcompress(string(i*iOutputInt),/remove_all)
    endif else begin
      if (i*iOutputInt LT 100) then infile = prefix+'000' + strcompress(string(i*iOutputInt),/remove_all) else begin
        if (i*iOutputInt LT 1000) then infile = prefix+'00' + strcompress(string(i*iOutputInt),/remove_all)
      endelse
    endelse
    if (i LT 10) then begin
      ppmout = "temp00"+strcompress(string(i),/remove_all)+".ppm"
    endif else if (i LT 100) then begin
      ppmout = "temp0"+strcompress(string(i),/remove_all)+".ppm"
    endif else ppmout = "temp"+strcompress(string(i),/remove_all)+".ppm"
    rtipsy,infile,h,g,d,s

    plot,1e-29*g.dens,g.tempg,psym=3,xtit="Density (g/cc)",ytit="Temperature (K)",/xlog,/ylog,xrange=[1e-34,1e-22],yrange=[1e3,1e7]
    img=tvrd()
    write_ppm,ppmout,img
  endfor
!ORDER=0

end
