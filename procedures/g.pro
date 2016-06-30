pro g
 !p.multi=[0,2,2]
  readcol,"p2runge5.dat",t,x,y,vx,vy,E
  readcol,"p2leap5.dat",t2,x2,y2,vx2,vy2,E2
  readcol,"p2arunge5.dat",t3,x3,y3,vx3,vy3,E3
  readcol,"p2aleap5.dat",t4,x4,y4,vx4,vy4,E4
  in1 = min((x^2+y^2)^(0.5))
  out1 = max((x^2+y^2)^(0.5))
  print,"in="+string(in1)
  print,"out="+string(out1)
  oplot,x2,y2,linestyle=1
  plot,x3,y3,title="Stepsize = .05",xtit="x",ytit="y",linestyle=2
  oplot,x4,y4,linestyle=1
  plot,t,E,title="Energy behavior",xtit="time",ytit="Energy",linestyle=2,yrange=[-0.208,-0.207]
  oplot,t2,E2,linestyle=1
  plot,t3,E3,title="Energy behavior",xtit="time",ytit="Energy",linestyle=2,yrange=[-0.50001,-0.499999]
  oplot,t4,E4,linestyle=1
  legend,['R-K','Leapfrog'],linestyle=[2,1],/bottom
end
