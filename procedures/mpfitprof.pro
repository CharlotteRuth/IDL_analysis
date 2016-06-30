function mpfitprof,rho,Rarray

  good=where(rho GT 0 AND rho LT max(rho)/2.)
  rhog = rho[good]
  rg = Rarray[good]
  ; Find r=0
  mt = min(where(rg LT max(rg)/2.))
  neg = n_elements(good)
  dhguess = (rg[neg-1]-rg[mt])/alog(rhog[mt]/rhog[neg-1])
; Initial guess
  a=[max(rhog),dhguess]
print,Rarray[mt]
print,a
  minrho = min(rho[good])
  maxR = max(Rarray[where(rho GT minrho)])
  fit=mpfitfun('mpexp',Rarray[good],rho[good],sqrt(rho[good]),a,perror=perror,bestnorm=bestnorm,/quiet)
  print,"Inner Scale length:",fit[1]

  ;twoa=[max(rho),a[1],max(rho-fit),a[1]*2.]
  ;twofit=mpfitfun('mptwoexp',Rarray[good],rho[good],sqrt(rho[good]),twoa,perror=perror,bestnorm=bestnorm)
  ;print,"Outer Scale length:",twofit[3]

  return,{fit:fit};,twofit:twofit}
end
