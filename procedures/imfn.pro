function imfn,mass,a,b,c

  n=mass

  ind = where(mass LT b[2])
  n[ind]=a[0]/(a[1]+1.0)*(mass[ind]^(a[1]+1.0)-a[2]^(a[1]+1.0))
  ind = where((mass LT c[2]) AND (mass GE b[2]))
  n[ind]=a[0]/(a[1]+1.0)*(b[2]^(a[1]+1.0)-a[2]^(a[1]+1.0))+ $
    b[0]/(b[1]+1.0)*(mass[ind]^(b[1]+1.0)-b[2]^(b[1]+1.0))
  ind = where(mass GE c[2])
  n[ind]=a[0]/(a[1]+1.0)*(b[2]^(a[1]+1.0)-a[2]^(a[1]+1.0))+ $
    b[0]/(b[1]+1.0)*(c[2]^(b[1]+1.0)-b[2]^(b[1]+1.0)) + $
    c[0]/(c[1]+1.0)*(mass[ind]^(c[1]+1.0)-c[2]^(c[1]+1.0))

  return,n

END
