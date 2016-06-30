pro b,file1,file2
  readcol,file1,ft,fnstar,fngas,fmstar,fmgas
  readcol,file2,ot,onstar,ongas,omstar,omgas
  ;readcol,"plottable",t,nstar,ngas,mstar,mgas
  dMSolUnit=1.35985e17
 ; msolstar=mstar*dMSolUnit
  fmsolstar=fmstar*dMSolUnit
  omsolstar=omstar*dMSolUnit
 ; msolgas=mgas*dMSolUnit
  fmsolgas=fmgas*dMSolUnit
  omsolgas=omgas*dMSolUnit
 ; msoltot=msolstar+msolgas
  fmsoltot=fmsolstar+fmsolgas
  omsoltot=omsolstar+omsolgas
  !p.multi=[0,2,2]
  plot,ft,fnstar,xtitle="expansion factor",ytitle="# of stars",linestyle=2
  oplot,ot,onstar,linestyle=3
 ; oplot,t,nstar
  legend,['new','fabio','old'],linestyle=[0,2,3]
  plot,ft,fngas,xtitle="expansion factor",ytitle="# of SPHs",yrange=[1e5,1.5e5],linestyle=2
 ; oplot,t,ngas
  oplot,ot,ongas,linestyle=3
  plot,ft,fmsolstar,xtitle="expansion factor",ytitle="Solar masses of stars",linestyle=2
 ; oplot,t,msolstar
  oplot,ot,omsolstar,linestyle=3
  plot,ft,fmsoltot,xtitle="expansion factor",ytitle="Solar masses of gas + stars",yrange=[6.3e10,6.4e10],linestyle=2
 ; oplot,t,msoltot
  oplot,ot,omsoltot,linestyle=3
end
