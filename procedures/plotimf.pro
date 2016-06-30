pro plotimf

  mass=findgen(400)/10.0
  ;Miller Scalo
  a=[0.3547,-0.4,0.1 ]
  b=[0.3547,-1.5,1.0 ]
  c=[2.0272,-2.3,10.0]
  nms=imfn(mass,a,b,c) ;Fraction of stars below each mass, Millo Scalo
  nsnms = fltarr(n_elements(nms))
  basemsn = max(nms[where(mass EQ 8.0)])
  nsnms[where(mass GT 8.0)]=nms[where(mass GT 8.0)]-basemsn
  ;Kroupa
  a=[0.56523,-0.3,0.08]
  b=[0.3029, -1.2,0.5 ]
  c=[0.3029, -1.7,1.0 ]
  nkr=imfn(mass,a,b,c) ;Fraction of stars below each mass, Kroupa
  nsnkr = fltarr(n_elements(nkr))
  basekrn = max(nkr[where(mass EQ 8.0)])
  nsnkr[where(mass GT 8.0)]=nkr[where(mass GT 8.0)]-basekrn
!p.multi=[0,2,2]
  plot,mass,nms,/xlog,xrange=[1e-1,40.0],xstyle=1
  oplot,mass,nkr,linestyle=2
  legend,['Kroupa','Miller-Scalo'],linestyle=[2,0],/bottom,/right
  plot,mass,1e4*nsnms,/xlog,xrange=[1e-1,40.0],xstyle=1,xtit="Star mass (M_sol)",ytit="SNII/1e4 M_sol "
  oplot,mass,1e4*nsnkr,linestyle=2
  ;plot,mass,fems,xtit="Star mass (M_sol)",ytit="Iron"
  ;oplot,mass,fekr
  ;plot,mass,oms,xtit="Star mass (M_sol)",ytit="Oxygen"
  ;oplot,mass,okr

END
