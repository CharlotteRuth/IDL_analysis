pro phasemarkplot

outputs=['MW1.1024g1.00016','MW1.1024g1.00018','MW1.1024g1.00032','MW1.1024g1.00035','MW1.1024g1.00048','MW1.1024g1.00058','MW1.1024g1.00064','MW1.1024g1.00080','MW1.1024g1.00081','MW1.1024g1.00096','MW1.1024g1.00112','MW1.1024g1.00128','MW1.1024g1.00144','MW1.1024g1.00208']

!p.thick=4
!p.charsize=1.5

  for i=0,n_elements(outputs)-1 do begin
    print,'Processing:  '+outputs[i]
    rtipsy,outputs[i],h,g,d,s
    rdfloat,outputs[i]+'.mrk',marks,skip=1
    ind = marks[where(marks LT h.ngas-1)]
    plot,1e-29*g[ind].dens,g[ind].tempg,psym=1,/xlog,/ylog,xtit="Density (g/cc)", ytit="Temperature (K)",xrange=[1e-31,1e-20],yrange=[10,1e7]
    legend,['z = '+strtrim(string(1./h.time -1.0),2)],box=0
    write_png,'gidmark.953/'+outputs[i]+'.png',tvrd()
  endfor 

END
