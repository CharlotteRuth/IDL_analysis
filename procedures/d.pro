pro d
 !p.multi=[0,2,3]
 for i=0,4 do begin
  j=strcompress(string(i),/remove_all)
  readcol,"p1runge"+j+".dat",t,x,v
  readcol,"p1leap"+j+".dat",t2,x2,v2
  realt=(7.0/100.0)*findgen(100)
  plot,x2,v2,title="Stepsize ="+string(1.6*(0.5^i)),xtit="x",ytit="v",linestyle=1
  oplot,x,v,linestyle=2
  oplot,2^(0.5)*cos(realt),-2^(0.5)*sin(realt)
 endfor
  legend,['anal','R-K','Leapfrog'],linestyle=[0,2,1],position=[2.0,0.75]
end
