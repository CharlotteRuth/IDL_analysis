pro phaseplot, gid

  basefile='MW.00512'
  outputs=['MW.00010','MW.00018','MW.00044','MW.00058','MW.00081','MW.00123','MW.00128','MW.00160','MW.00184','MW.00186','MW.00192','MW.00200','MW.00208','MW.00216','MW.00219','MW.00224','MW.00232','MW.00240','MW.00248','MW.00256','MW.00262','MW.00264','MW.00272','MW.00280','MW.00288','MW.00296','MW.00304','MW.00312','MW.00320','MW.00321','MW.00328','MW.00336','MW.00344','MW.00352','MW.00360','MW.00368','MW.00376','MW.00384','MW.00392','MW.00400','MW.00401','MW.00408','MW.00416','MW.00424','MW.00432','MW.00440','MW.00448','MW.00452','MW.00456','MW.00464','MW.00472','MW.00480','MW.00488','MW.00496','MW.00504','MW.00512'];

!p.thick=4
!p.charsize=1.5

  for i=0,n_elements(outputs)-1 do begin
    print,'Processing:  '+outputs[i]
    rtipsy,'../'+outputs[i],h,g,d,s
    grp= read_ascii_array('../'+basefile+'.'+outputs[i]+'.grp',/long)
    ggrp = grp[0:h.ngas-1]
    ind = where(ggrp EQ gid)
    plot,1e-29*g[ind].dens,g[ind].tempg,psym=1,/xlog,/ylog,xtit="Density (g/cc)", ytit="Temperature (K)",xrange=[1e-31,1e-20],yrange=[1e3,1e7]
    legend,['z = '+1./h.time -1.0],box=0
    write_png,'phasemov.'+string(gid)+'.'+outputs[i]+'.png',tvrd()
  endfor 

END
