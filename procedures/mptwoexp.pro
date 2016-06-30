FUNCTION mptwoexp,t,a

  f = a[0]*exp(-t/a[1])+a[2]*exp(-t/a[3])
  return,f

END
