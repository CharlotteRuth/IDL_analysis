pro fitfunc,x,a,y,pder
;This is the function to which I fit the isochrone
y = a[0] + a[1]*x + a[2]*x^2 + a[3]*x^3 + a[4]*x^4 + a[5]*x^5
x=reform(x)
;If the procedure is called with four parameters, calculate the  
;partial derivatives.  
  IF (N_PARAMS() GE 4) THEN $  
    pder = [[replicate(1.0, N_ELEMENTS(X))],[x],[x^2],[x^3],[x^4],[x^5]] 
end
