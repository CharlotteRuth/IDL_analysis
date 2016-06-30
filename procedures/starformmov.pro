pro starformmov,tipsyfile
  paramfile = strmid(tipsyfile,0,strlen(tipsyfile)-5) + "param"
  dKpcUnit = 1e5
  a = timeunits(paramfile)
  dtUnit = a[0]
  dt  = 1e6
  twidth = 1e7
  rformfile=tipsyfile+".rform"
  rdfloat,rformfile,rform,skip=1
  rtipsy,tipsyfile,h,g,d,s
  ns = h.nstar
  ntot = h.n
  ngd= h.ngas + h.ndark

  tform = s.tform * dtUnit

  ind = where(tform gt 0)
  mint = min(tform[ind])
  maxt = max(tform[ind])
  
  r = replicate({x:0.,y:0.,z:0.},ns)
  r.x = rform[ngd:ntot-1]
  r.y = rform[ntot+ngd:2*ntot-1]
  r.z = rform[2*ntot+ngd:3*ntot-1]

!ORDER=1
 ;mpegid = mpeg_open([400,400]) 
  for t=mint,maxt,dt do begin
    i = fix(t/dt + 0.1)
    in = where( (tform gt t-twidth) AND (tform lt t+twidth),nin )
    if (nin gt 0) then begin
      plot, r[in].x*dKpcUnit, r[in].y*dKpcUnit,xrange=[-20,20],yrange=[-20,20],psym=4,/isotropic,xtit='X [kpc]',ytit = 'Y [kpc]'
      if (i LT 10) then begin
        ppmout = "sf00"+strcompress(string(i),/remove_all)+".ppm"
      endif else if (i LT 100) then begin
        ppmout = "sf0"+strcompress(string(i),/remove_all)+".ppm"
      endif else ppmout = "sf"+strcompress(string(i),/remove_all)+".ppm"
      img = tvrd()
      write_ppm,ppmout,img

 ;     mpeg_put,mpegid,window=0,frame=fix(t/dt + 0.1)
    endif
  endfor
!ORDER=0  
 ;mpeg_save,mpegid,filename='starform.mpg'
 ;mpeg_close,mpegid
end
