PRO sandagetau,t,a,f,pder

  btau=1.d/(1.0-exp(-(a[0]^2.0)/(2.d*a[1]^2.0)))
  expt=exp(-((t)^2.0)/(2.d*(a[1]^2.0)))
  expa=exp((a[0]^2.0)/(2.d*(a[1]^2.0)))
  mult=expt*expa
  f = btau*((t)/(a[1]^2.0))*expt

  if n_params() GE 4 THEN $
    pder = [[-a[0]*(t)*mult*btau*btau/(a[1]^(4.d))], [btau*btau*mult*(t)*(a[0]*a[0]-(expa-1.d)*(2.d*a[1]*a[1]-(t)*(t)))/(a[1]^(5.d))] ]

END
