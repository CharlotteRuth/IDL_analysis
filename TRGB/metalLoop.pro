PRO metalLoop,galaxy,index,data,points
rgb_v = data.MAG1_ACS[index]
rgb_i = data.MAG2_ACS[index]
rgb_verr = data.MAG1_STD[index]
rgb_ierr = data.MAG2_STD[index]
rgb_vi = rgb_v - rgb_i
minC = MIN(rgb_vi)
;minC = 0.5
maxC = MAX(rgb_vi)
;maxC = 1.9
dC = (maxC - minC)/points
print,"minC: ",minC," maxC: ",maxC," dC: ",dC
results = fltarr(points,3)

FOR i=0.0, points-1 DO BEGIN
    subindex = where((rgb_vi LT (dC*(i+1)+minC)) AND (rgb_vi GE (dC*i+minC)))
    plot,rgb_vi[subindex],rgb_v[subindex],psym=3,xrange=[-1,4],yrange=[28,22]
    print,"i: ",i
    results[i,2] = 1
    results[i,1]=trgb(galaxy,rgb_v[subindex],rgb_verr[subindex],rgb_i[subindex],rgb_ierr[subindex])
    results[i,0]=dc*(i+.5)+minC
    print,results[i,2]
    stop
END

fit = poly_fit(results[*,0],results[*,1],1,MEASURE_ERRORS = results[*,2],SIGMA = SIGMA)
print,'Poly Fit: ',fit
print,'Poly Fit Errors: ',sigma

fit_linex = findgen(20)*((maxC-MinC)/20) + minC
print,fit_linex
fit_liney = fit[0] + fit[1]*fit_linex


set_plot,'x'
;device,filename=galaxy+'_metalDependence.eps',/color, bits_per_pixel=8
plot,rgb_v-rgb_i,rgb_v,XTITLE='V-I',YTITLE='V',TITLE='Metal Dependancy',psym = 3,yrange = [30,22]
oplot,results[*,0],results[*,1],psym=4,color = 254
oplot,fit_linex,fit_liney,color = 254
;device,/close

;infile = "galaxy+'/'+galaxy +'_trgb_vi.dat'
;openr,1,outfile
;readf,1,data

;plot,data(0,*),
    
end    
