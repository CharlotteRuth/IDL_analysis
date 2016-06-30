PRO NFW
;check out the NFW profile.  

r=FINDGEN(100000)*1.e-8
c=15.d0
m200=8.24e-6
r200=.00214
cx = r*(c/r200)
dSoft = 5.0e-5
eps = c*dSoft/r200
;m_const=1.
g=1.
m_const = m200/( (1./3.)*eps*eps/(1+eps)/(1+eps) - 1/(1+eps) + alog((1+c)/(1+eps)) + 1/(1+c) )
a=FINDGEN(N_ELEMENTS(r))			
  ind=WHERE(cx LT eps,comp=comp)
  a[ind]=G*M_const*(1./3.)/(eps*(1+eps)*(1+eps))*(c*c*c)/(r200*r200*r200)


  a[comp] = G*M_const*((1./3.)*eps*eps/(1+eps)/(1+eps)- 1/(1+eps) + ALOG((1+cx[comp])/(1+eps)) + 1/(1+cx[comp]))/(r[comp]*r[comp]*r[comp]);


x=ALOG10((1+cx[comp])/(1+eps)) + 1/(1+cx[comp])

;  a[comp]=G*M_const*( (1./3.)*eps*eps/(1+eps)/(1+eps) - 1/(1+eps) + alog10((1+cx[comp])/(1+eps)) + 1/(1+cx[comp]) ) /(r[comp]*r[comp]*r[comp])
  a=-a*r ; multiply by r to get force
LOADCT,0
plot,r,ABS(a),/xlog,/ylog,xrange=[1e-5,1e-3],yrange=[.9,80]
LOADCT,13
  
    oplot,r[ind],ABS(a[ind]),COLOR=100,psym=3
    oplot,r[comp],ABS(a[comp]),COLOR=200,psym=3
STOP
END

