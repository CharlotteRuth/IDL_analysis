pro trhoprof,start=start,finish=finish

  restore,'/net/grads-1/stinson/trace/tipsyprofile.sav'
if (keyword_set(start)) then begin
  a = read_ascii(start+'.prof',template=temp)
  c = read_ascii(start+'.gas.prof',template=temp)
endif else begin
  a = read_ascii('MW100.prof',template=temp)
  c = read_ascii('MW100.gas.prof',template=temp)
endelse
if (keyword_set(finish)) then begin
  b = read_ascii(finish+'.prof',template=temp)
  d = read_ascii(finish+'.gas.prof',template=temp)
endif else begin
  ;b = read_ascii('adiabatic/MW.00300.prof',template=temp)
  ;d = read_ascii('adiabatic/MW.00300.gas.prof',template=temp)
  b = read_ascii('nosf/MW.00170.prof',template=temp)
  d = read_ascii('nosf/MW.00170.gas.prof',template=temp)
endelse
  !p.multi = [0,2,1]
  plot,c.r,c.T,/xlog,/ylog,xtit='radius [kpc]',ytit='Temperature [K]'
  oplot,d.r,d.T,linestyle=2

  plot,a.r,a.rho,/xlog,/ylog,yrange=[1e-3,1e4],xtit='radius [kpc]',ytit='Density'
  oplot,b.r,b.rho,linestyle=2
  oplot,c.r,c.rho,line=1
  oplot,d.r,d.rho,linestyle=3
END
