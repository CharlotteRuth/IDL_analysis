pro f
 !p.multi=[0,2,3]
 for i=0,4 do begin
  j=strcompress(string(i),/remove_all)
  readcol,"p1runge"+j+".dat",t,x,v
  readcol,"p1leap"+j+".dat",t2,x2,v2
  real=2^(.5)*cos(t-!pi/4)
  plot,t2,x2-real,title="Stepsize ="+string(1.6*(0.5^i)),xtit="t",ytit="error",linestyle=1,xrange=[0,7]
  oplot,t,x-real,linestyle=2
 endfor
  legend,['R-K','Leapfrog'],linestyle=[2,1],/left
end
