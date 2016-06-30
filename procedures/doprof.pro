;paperps,file='starprof.eps',charsize=2.0

dirs=['30','25','20','15']
for i=0,n_elements(dirs)-1 do begin
  file=dirs[i]+'kms/'+dirs[i]+'kms.01350'
print,file
  rtipsy,file,h,g,d,s
;Use minimum potential as center
  xcen=s[where(s.phi EQ min(s.phi))].x
  ycen=s[where(s.phi EQ min(s.phi))].y
  zcen=s[where(s.phi EQ min(s.phi))].z
  s.x=s.x-xcen
  s.y=s.y-ycen
  s.z=s.z-zcen

  if (keyword_set(xz)) then Rsqr = (s.x^2. + s.z^2.) $
  else if (keyword_set(yz)) then Rsqr = (s.y^2. + s.z^2.) $
  else Rsqr = (s.x^2. + s.y^2.)
  R = sqrt(Rsqr)

  prof=profile(R,s.mass,nbin=50,/xz)
  minsig=min(prof.sigma[where(prof.sigma GT 0)])
  maxsig=max(prof.sigma)
  if(i GT 0) then oplot,prof.r,prof.sigma,psym=i+1 $
  else plot,prof.r,prof.sigma,psym=i+1,/ylog,yrange=[1e-5,maxsig]
    ;xrange=[0,10]
  fit=mpfitprof(prof.sigma,prof.r)
  oplot,prof.r,mpexp(prof.r,fit.fit),line=i
print,"Max r",max(prof.r)
  ;oplot,prof.r,mpexp(prof.r,fit.twofit),line=2
endfor
  legend,dirs+' km/s',psym=indgen(n_elements(dirs))+1,line=indgen(n_elements(dirs)),/right
;starprof,'15kms/15kms.01000',xmargin=[7.5,0.8],ymargin=[3.2,0.4]
;ppsclose

end
