pro sd,filename
   readcol,filename,r,num,den,m,vc,vr,sr,vth,sth,j,jth,jph
   n = n_elements(m)
   SD = fltarr(n-1)
   for i=0,n-2 do begin
     SD[i] = (m[i+1] - m[i])/(!PI*(r[i+1]^2.0 - r[i]^2.0))
   endfor
   plot,r,SD,xtit="radius",ytit="surface density",/ylog
   j = 5
   SD0 = SD[j]
   while ( SD0/2.71828 LT SD[j] ) do j=j+1
   print,"scale length ="+string(r[j]-r[5])+" r["+string(j)+"]="+string(r[j])+" r[5]="+string(r[5])
end
