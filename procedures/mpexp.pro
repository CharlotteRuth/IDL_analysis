function mpexp,x,a

  f = a[0]*exp(-x/a[1])
  return,f
END
