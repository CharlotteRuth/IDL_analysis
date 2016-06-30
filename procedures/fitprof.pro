function fitprof,rho,Rarray

  ; Find r=0
  mt = min(where(rho EQ max(rho)))
  dhguess = (Rarray[mt+10]-Rarray[mt])/alog(rho[mt]/rho[mt+10])

; Initial guess
  a=[max(rho),dhguess]
print,Rarray[mt]
print,a
  minrho = min(rho[where(rho GT 0)])
  maxR = max(Rarray[where(rho GT minrho)])
  weights=fltarr(n_elements(Rarray));+1.0
  weights[mt:min(where(rho LT max(rho)/1e4 AND rho GT 0))] = 1.0
;for i=0,nbin-1 do print,Rarray[i],weights[i]
  fit=curvefit(Rarray,rho,weights,a,offset,function_name='expfit')
  print,"Inner Scale length:",a[1]

  ;outrho = rho - fit
  ;mt=min(where(outrho EQ max(outrho)))
  ;hhguess = (Rarray[mt+3]-Rarray[mt])/alog(outrho[mt]/outrho[mt+3])
  ;outa=[max(outrho),dhguess]
  twoa=[max(rho),dhguess,max(rho-fit),a[1]*2.]
  weights=fltarr(n_elements(Rarray));+1.0
  weights[where(rho GT max(rho)/1e6)]=1.0
  twofit=curvefit(Rarray,rho,weights,twoa,offset,function_name='twoexpfit')
  print,"Outer Scale length:",twoa[3]

  return,{fit:fit,a:a,twofit:twofit,twoa:twoa}
end
