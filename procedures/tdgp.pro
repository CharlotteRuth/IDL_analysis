
restore,'../trace/tipsyprofile.sav'
a=read_ascii('15kms.01000.allprof',template=temp)
g=read_ascii('15kms.01000.gasprof',template=temp)
paperps,file='tipsdensprofs.eps'
  plot,a.r,a.rho,/xlog,/ylog,xrange=[0.01,300], $
    yrange=[min(a.rho[where(a.rho GT 0)]),max(a.rho)], $
    xmargin=[6.5,1.2]
  oplot,g.r,g.rho,line=2
ppsclose
!p.multi=[0,1,1]

END
