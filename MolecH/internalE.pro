

PRO internalE,filename,Y,dens,Pres,u_code
;  zsolar = 0.0130215
  CL_Eerg_gm_degK = 8.2494e7
  unit_conv = 1e10
  gamma = 5./3.

  rtipsy,filename,h,g,d,s
  readarr,filename + '.HI',h,HI,part = 'gas',/ascii
  IF file_test(filename + '.H2') THEN readarr,filename + '.H2',h,H2,part = 'gas',/ascii ELSE h2 = HI*0
  readarr,filename + '.HeI',h,HeI,part = 'gas',/ascii
  readarr,filename + '.HeII',h,HeII,part = 'gas',/ascii
  Y_H = HI*0
  Y_He = HI*0
 ; g.zmetal = g.zmetal/zsolar

  lowz = where(g.zmetal LE 0.1, complement = highz)
  Y_He[lowz] = (0.236 + 2.1*g[lowz].zmetal)/4.0
  IF highz[0] NE -1 THEN Y_He[highz] =  (-0.446*(g[highz].zmetal - 0.1)/0.9 + 0.446)/4.0
  Y_H = 1.0 - Y_He*4.0 - g.zmetal ; 
  HII = Y_H - HI - 2*H2
  HeIII = Y_He - HeI - HeII
  Y = Y_H - H2 + HII + Y_He + HeII + 2.0*HeIII

  u = 3./2.*Y*g.tempg*CL_Eerg_gm_degK
  u_code = u/unit_conv
  Pres = u*g.dens*(gamma - 1)

  dens = g.dens
END

;device,/color,filename = 'u_dens.eps'
;plot,dens2,u2,/ylog,/xlog,psym = 3,xrange = [1e-6,1e6],yrange = [1e-1,1e5]
