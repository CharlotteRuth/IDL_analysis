pro sfrMaster,tipsyfile,dMsolUnit,dKpcUnit

timeunit=3.872d3*SQRT((dKpcUnit*3.0856805d21)^3/(dMsolUnit*1.98892d33))/31556926.0

components = tipsyfile+'.cmp'

loadct,39
set_plot,'ps'
device,/portrait,filename=tipsyfile+'SFH.ps', /color

!P.thick = 1.5
!X.Charsize = 1.25
!Y.Charsize = 1.25
!P.CHARSIZE = 1.25

;;Component Index Numbers in the .cmp file

diskcomp = 1
halocomp = 2
bulgecomp = 3
thickcomp = 4
psbulgecomp = 5
othercomp = 6

rtipsy,tipsyfile,h,g,d,s
readcol,components,data
star = data[h.n-h.nstar+1:h.n ]

diskind = where(star eq 1 or star eq 4)
diskstars = s[diskind]

bulgeind = where(star eq 3 or star eq 2 or star eq 5)
bulgestars = s[bulgeind]

satind = where(star eq 6)
satstars = s[satind]

;;Plot the total SFH

sfr,s,massu = dMsolUnit,timeunit = timeunit,OVERPLOT = 0,title = 'SFH for h258',xtitle = 'Time [Gyrs]', ytitle = 'Average SFR [yr^-1]',linestyle = 0
sfr,s,massu = dMsolUnit,timeunit = timeunit,OVERPLOT = 1,linestyle = 0,color=250

;;Plot the Disk Stars (disk + thickdisk).
sfr,diskstars,massu=dMsolUnit,timeunit = timeunit,OVERPLOT=1,linestyle = 0,color=150

;;Plot the Bulge Stars (bulge + psbulge + halo)
sfr,bulgestars,massu=dMsolUnit,timeunit = timeunit,OVERPLOT = 1,linestyle = 0,color=80

;;Plot the Satellites
sfr,satstars,massu=dMsolUnit,timeunit = timeunit,OVERPLOT = 1,linestyle = 0,color=200

legend,['Total','Disk + Thickdisk','Bulge + Psuedo Bulge + Halo','Other'],linestyle = [0,0,0,0],color = [250,150,80,200],charsize=1


device,/close
set_plot,'x'

stop

END
