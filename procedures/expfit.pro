PRO expfit,t,a,f,pder

  f = a[0]*exp(-t/a[1])

  if n_params() GE 4 THEN $
    pder = [[exp(-t/a[1])], [(a[0]/a[1]^2.)*t*exp(-t/a[1])] ]

END
