pro rotcur,filename
  readcol,filename,r,num,den,m,vc,vr,sr,vth,sth,j,jth,jphi
  plot,r*1e5,vc*2418,xtit="radius (kpc)",ytit="velocity (km/s)",tit="Rotation Curve"
stop

end
