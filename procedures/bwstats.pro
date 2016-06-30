pro bwstats,tmax=tmax

  dMSolUnit=13.5985e16
  dpcUnit=1e8
  InitStarMass21k=3.6e-12*dMSolUnit
  InitStarMass105k=7.2e-13*dMSolUnit
  InitStarMass525k=1.44e-13*dMSolUnit

  !p.multi=[0,2,2]
  restore,'bwtemplate.sav'
  o21=read_ascii('21k.bwout',template=temp)
  o105=read_ascii('105k.bwout',template=temp)
  o525=read_ascii('525k.bwout',template=temp)
  nsn21 = o21.esn/1.8d54/InitStarMass21k
  nsn105 = o105.esn/1.8d54/InitStarMass105k
  nsn525 = o525.esn/1.8d54/InitStarMass525k
  ;nsn21 = o21.nsn/InitStarMass21k/1.8d3
  ;nsn105 = o105.nsn/InitStarMass105k/1.8e3
  ;nsn525 = o525.nsn/InitStarMass525k/1.8e3
  r21 = o21.r*1e8*(InitStarMass525k/InitStarMass21k)^(1./3.)
  r105 = o105.r*1e8*(InitStarMass525k/InitStarMass105k)^(1./3.)
  r525 = o525.r*1e8*(InitStarMass525k/InitStarMass525k)^(1./3.)

paperps,filename='bwstats.eps',charsize=1.3
  print,'NSN'
  mom=moment(nsn21)
  print,'21k:  ',mom[0], ' +/-',sqrt(mom[1])
  bin21=sqrt(mom[1])/5
  mom=moment(nsn105)
  print,'105k:  ',mom[0], ' +/-',sqrt(mom[1])
  bin105=sqrt(mom[1])/5
  mom=moment(nsn525)
  print,'525k:  ',mom[0], ' +/-',sqrt(mom[1])
  bin525=sqrt(mom[1])/5
  plothist,nsn21*1e8,peak=100,xtit='SN Energy Eff (10!u8!n E!lSN!u*!n/m!lstar!nc!u2!n)',bin=bin21*1e8,ytit='N',xmin=0,xstyle=1,yrange=[0,110],ystyle=1,line=1,xmargin=[6,0.5],ymargin=[3.7,0.2]
  plothist,nsn105*1e8,peak=100,/over,linestyle=2,bin=bin105*1e8
  plothist,nsn525*1e8,peak=100,/over,bin=bin525*1e8

  print,'number'
  mom=moment(o21.n)
  print,'21k:  ',mom[0], ' +/-',sqrt(mom[1])
  mom=moment(o105.n)
  print,'105k:  ',mom[0], ' +/-',sqrt(mom[1])
  mom=moment(o525.n)
  print,'525k:  ',mom[0], ' +/-',sqrt(mom[1])
  plothist,o21.n,peak=100,xtit='Particles with Cooling Disabled',ytit='N',xmin=0,xstyle=1,yrange=[0,110],ystyle=1,line=1,xmargin=[6,0.8],ymargin=[3.7,0.2]
  plothist,o105.n,peak=100,/over,linestyle=2
  plothist,o525.n,peak=100,/over

  print,'Blast Radius (pc)'
  mom=moment(r21)
  print,'21k:  ',mom[0], ' +/-',sqrt(mom[1])
  bin21=sqrt(mom[1])/5
  mom=moment(r105)
  print,'105k:  ',mom[0], ' +/-',sqrt(mom[1])
  bin105=sqrt(mom[1])/5
  mom=moment(r525)
  print,'525k:  ',mom[0], ' +/-',sqrt(mom[1])
  bin525=sqrt(mom[1])/5
  plothist,r21,peak=100,bin=bin21,xtit='Blast Radius*(M!lstar,hr!u1/3!n / M!lstar!u1/3!n) [pc]',ytit='N',xmin=0,xstyle=1,yrange=[0,110],ystyle=1,line=1,xmargin=[6,0.5],ymargin=[3.7,0.2]
  plothist,r105,peak=100,/over,linestyle=2,bin=bin105
  plothist,r525,peak=100,/over,bin=bin525

  print,'Shutoff time (yrs)'
  mom=moment(4e10*o21.t)
  print,'21k:  ',mom[0], ' +/-',sqrt(mom[1])
  mom=moment(4e10*o105.t)
  print,'105k:  ',mom[0], ' +/-',sqrt(mom[1])
  mom=moment(4e10*o525.t)
  print,'525k:  ',mom[0], ' +/-',sqrt(mom[1])
  tbin=max(4e10*o525.t/100)
  tmaxdef=max(4e10*o525.t)
  if (keyword_set(tmax)) then plothist,4e10*o21.t,peak=100,bin=tbin,xtit='Cooling Shutoff Time (yrs)',xrange=[0,tmax],ytit='N',yrange=[0,110],ystyle=1,line=1,xmargin=[6,0.8],ymargin=[3.7,0.2] $
else plothist,4e10*o21.t,peak=100,bin=tbin,xtit='Cooling Shutoff Time (yrs)',/xlog,xrange=[5e6,1.1e8],xstyle=1,ytit='N',yrange=[0,110],ystyle=1,line=1,xmargin=[6,0.8],ymargin=[3.7,0.2]
  plothist,4e10*o105.t,peak=100,bin=tbin,/over,linestyle=2
  plothist,4e10*o525.t,peak=100,bin=tbin,/over

device,/close
end
