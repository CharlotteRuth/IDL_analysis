
loadct,39
!p.multi = [0,2,2]
dirs = ['15','20','25','30']
ned = n_elements(dirs)
nouts=27
scalelength = fltarr(ned,nouts)

paperps,file='timeprofs.eps',/color,charsize=1.6
for j=0,n_elements(dirs)-1 do begin
for i=1,nouts do begin
  t=string(i*50,format='(I05)')
  file=dirs[j]+'kms/'+dirs[j]+'kms.'+t
print,file
  rtipsy,file,h,g,d,s
;Use minimum potential as center
  xcen=s[where(s.phi EQ min(s.phi))].x
  ycen=s[where(s.phi EQ min(s.phi))].y
  zcen=s[where(s.phi EQ min(s.phi))].z
  s.x=s.x-xcen
  s.y=s.y-ycen
  s.z=s.z-zcen
  Rsqr = (s.x^2. + s.y^2.)
  R = sqrt(Rsqr)

  prof=profile(R,s.mass,nbin=50,/xz)
  minsig=min(prof.sigma)
  maxsig=max(prof.sigma)
  if(i GT 1) then oplot,prof.r,prof.sigma,psym=1,color=i*254/nouts $
  else plot,prof.r,prof.sigma,psym=1,/ylog,yrange=[minsig,50*maxsig], $
    xtit='R [kpc]',ytit='Density [M'+sunsymbol()+' kpc!u-2!n]', $
    xmargin=[9.0,0.3],ymargin=[3.2,0.4]
  fit=mpfitprof(prof.sigma,prof.r)
  oplot,prof.r,mpexp(prof.r,fit.fit),color=i*254/nouts

  scalelength[j,i-1]=fit.fit[1]
endfor
endfor
ppsclose

loadct,39
paperps,file='timescalelengths.eps',/color,charsize=1.6
  plot,0.5*(findgen(nouts)+1),scalelength[0,*],psym=-1, $
     xtit='time [Gyr]',ytit='Scalelength [kpc]', $
     xmargin=[5.5,0.9],ymargin=[3.0,0.3],yrange=[0.3,1.],ystyle=1
for j=1,ned-1 do begin
  oplot,0.5*(findgen(nouts)+1),scalelength[j,*],psym=-(j+1),color=254*j/ned
endfor
legend,dirs+' km s!u-1!n',psym=-(indgen(ned)+1),color=254*indgen(ned)/ned
ppsclose
END
