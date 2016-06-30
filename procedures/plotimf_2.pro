function dMSCumNumber,mass, mmax, a, b, c
;mmax = 100
n = fltarr(N_ELEMENTS(mass))

ind = where(mass lt mmax)
n[ind] = 0

ind = where(mass ge c[2], COMPLEMENT = nind)
IF (ind[0] ne -1) THEN n[ind]  = c[0]/c[1]*(mmax^c[1] - mass[ind]^c[1])
;n[nind] = c[0]/c[1]*(mmax^c[1] - c[2]^c[1]) 

ind = where(mass ge b[2] AND mass lt c[2], COMPLEMENT = nind)
IF (ind[0] ne -1) THEN n[ind]  = c[0]/c[1]*(mmax^c[1] - c[2]^c[1]) $
        + b[0]/b[1]*(c[2]^b[1] - mass[ind]^b[1])
;n[nind] = n[nind] + b[0]/b[1]*(c[2]^b[1] - b[2]^b[1])    

ind = where(mass ge a[2] AND mass lt b[2], COMPLEMENT = nind) 
IF (ind[0] ne -1) THEN n[ind]  = c[0]/c[1]*(mmax^c[1] - c[2]^c[1]) $
        + b[0]/b[1]*(c[2]^b[1] - b[2]^b[1]) $ 
        + a[0]/a[1]*(b[2]^a[1] - mass[ind]^a[1])
;n[nind] = n[nind] + a[0]/a[1]*(b[2]^a[1] - a[2]^a[1])

return,n                   ;
END


function imfn,mass,a,b,c

n = fltarr(N_ELEMENTS(mass))

ind = where(mass LT b[2], COMPLEMENT = nind)
IF (ind[0] ne -1) THEN n[ind]  = a[0]/(a[1]+1.0)*(mass[ind]^(a[1]+1.0)-a[2]^(a[1]+1.0))
;n[nind] = a[0]/(a[1]+1.0)*(b[2]^(a[1]+1.0)     -a[2]^(a[1]+1.0))

ind = where((mass LT c[2]) AND (mass GE b[2]), COMPLEMENT = nind)
IF (ind[0] ne -1) THEN n[ind]  = a[0]/(a[1]+1.0)*(b[2]     ^(a[1]+1.0)-a[2]^(a[1]+1.0)) $
        + b[0]/(b[1]+1.0)*(mass[ind]^(b[1]+1.0)-b[2]^(b[1]+1.0))
;n[nind] = n[nind] + b[0]/(b[1]+1.0)*(c[2]^(b[1]+1.0)     -b[2]^(b[1]+1.0))

ind = where(mass GE c[2], COMPLEMENT = nind)
IF (ind[0] ne -1) THEN n[ind]  = a[0]/(a[1]+1.0)*(b[2]     ^(a[1]+1.0)-a[2]^(a[1]+1.0)) $
        + b[0]/(b[1]+1.0)*(c[2]     ^(b[1]+1.0)-b[2]^(b[1]+1.0)) $
        + c[0]/(c[1]+1.0)*(mass[ind]^(c[1]+1.0)-c[2]^(c[1]+1.0))

return,n

END

;Cumulative mass of stars with mass greater than "mass".

FUNCTION dMSCumMass,mass, mmax, a, b, c
;mmax = 100
n = fltarr(N_ELEMENTS(mass))

ind = where(mass lt mmax)
n[ind] = 0

ind = where(mass ge c[2], COMPLEMENT = nind)
IF (ind[0] ne -1) THEN n[ind]  = c[0]/(c[1] + 1)*(mmax^(c[1] + 1) - mass[ind]^(c[1] + 1))
;n[nind] = c[0]/(c[1] + 1)*(mmax^(c[1] + 1) - c[2]^(c[1] + 1)) 

ind = where(mass ge b[2] AND mass lt c[2], COMPLEMENT = nind)
IF (ind[0] ne -1) THEN n[ind]  = c[0]/(c[1] + 1)*(mmax^(c[1] + 1) - c[2]^(c[1] + 1)) $ 
        + b[0]/(b[1] + 1)*(c[2]^(b[1] + 1) - mass[ind]^(b[1] + 1))
;n[nind] = n[nind] + b[0]/(b[1] + 1)*(c[2]^(b[1] + 1) - b[2]^(b[1] + 1))    

ind = where(mass ge a[2] AND mass lt b[2], COMPLEMENT = nind) 
IF (ind[0] ne -1) THEN n[ind]  = c[0]/(c[1] + 1)*(mmax^(c[1] + 1) - c[2]^(c[1] + 1)) $    
        + b[0]/(b[1] + 1)*(c[2]^(b[1] + 1) - b[2]^(b[1] + 1)) $
        + a[0]/(a[1] + 1)*(b[2]^(a[1] + 1) - mass[ind]^(a[1] + 1))
;n[nind] = n[nind] + a[0]/(a[1] + 1)*(b[2]^(a[1] + 1) - a[2]^(a[1] + 1))

return,n ;
END

pro plotimf_2
  mmax = 100.0
  mass=findgen(400*mmax)/mmax
  ;Miller Scalo
  ams=[0.3547,-0.4,0.1 ]
  bms=[0.3547,-1.5,1.0 ]
  cms=[2.0272,-2.3,10.0]
  nms = dMSCumNumber(mass,mmax,ams,bms,cms) ;Millo Scalo
  mms = dMSCumMass(mass,mmax,ams,bms,cms)
  ms = imfn(mass,ams,bms,cms)
  nsnms = fltarr(n_elements(nms))
  basemsn = max(nms[where(mass EQ 8.0)])
  nsnms[where(mass GT 8.0)]=nms[where(mass GT 8.0)]-basemsn
  ;Kroupa
  akr=[0.56523,-0.3,0.08]
  bkr=[0.3029, -1.2,0.5 ]
  ckr=[0.3029, -1.7,1.0 ]
  nkr = dMSCumNumber(mass,mmax,akr,bkr,ckr) ;Kroupa
  mkr = dMSCumMass(mass,mmax,akr,bkr,ckr)
  kr = imfn(mass,akr,bkr,ckr)
  nsnkr = fltarr(n_elements(nkr))
  basekrn = max(nkr[where(mass EQ 8.0)])
  nsnkr[where(mass GT 8.0)]=nkr[where(mass GT 8.0)]-basekrn

  window,0
  plot,mass,nms/MAX(nms),/xlog,xrange=[1e-1,40.0],xstyle=1,xtitle = 'Stellar Mass (M_sol)',ytitle = 'Num Above a Mass',/ylog,yrange=[1e-4,1]
  oplot,mass,nkr/MAX(nkr),linestyle=2
  legend,['Kroupa','Miller-Scalo'],linestyle=[2,0],/top,/right
;  stop

  window,1
  plot,mass,mms,/xlog,xrange=[1e-1,40.0],xstyle=1,xtitle = 'Stellar Mass (M_sol)',ytitle = 'Mass Above a Mass',yrange = [0,1]
  oplot,mass,mkr,linestyle=2
  legend,['Kroupa','Miller-Scalo'],linestyle=[2,0],/top,/right
;  stop

  window,2
  plot,mass,ms,/xlog,xrange=[1e-1,40.0],xstyle=1,xtitle = 'Stellar Mass (M_sol)',ytitle = 'Mass Below a Mass',yrange = [0,1]
  oplot,mass,kr,linestyle=2
  legend,['Kroupa','Miller-Scalo'],linestyle=[2,0],/bottom,/right
  

  print,'Number of type II SN: ',dMSCumNumber([8],40,ams,bms,cms),dMSCumNumber([8],40,akr,bkr,ckr),dMSCumNumber([8],40,ams,bms,cms)/dMSCumNumber([8],40,akr,bkr,ckr)
  print,'Mass of type II SN:   ',dMSCumMass([8],40,ams,bms,cms),dMSCumMass([8],40,akr,bkr,ckr),dMSCumMass([8],40,ams,bms,cms)/dMSCumMass([8],40,akr,bkr,ckr)
  print,'Number of type Ia SN: ',dMSCumNumber([3],16,ams,bms,cms),dMSCumNumber([3],16,akr,bkr,ckr),dMSCumNumber([3],16,ams,bms,cms)/dMSCumNumber([3],16,akr,bkr,ckr)
  print,'Mass of type Ia SN:   ',dMSCumMass([3],16,ams,bms,cms),dMSCumMass([3],16,akr,bkr,ckr),dMSCumMass([3],16,ams,bms,cms)/dMSCumMass([3],16,akr,bkr,ckr)
  stop  

  

  plot,mass,1e4*nsnms,/xlog,xrange=[1e-1,40.0],xstyle=1,xtit="Star mass (M_sol)",ytit="SNII/1e4 M_sol "
  oplot,mass,1e4*nsnkr,linestyle=2
  ;plot,mass,fems,xtit="Star mass (M_sol)",ytit="Iron"
  ;oplot,mass,fekr
  ;plot,mass,oms,xtit="Star mass (M_sol)",ytit="Oxygen"
  ;oplot,mass,okr

END
