pro fourier
;This graphs the (star fromation rate)/(gass mass) over time for galaxies with
;different resolution.  One graph per galactic mass is produced.

r=['5E1R', '1E2R', '1E3R', '1E4R', '1E5R']
m=['9M', '10M', '11M','12M', '13M']

dataread = dblarr(2,7)
integrate=fltarr(5)
binsize=1.e8
timeunit=1e9
massunit=2.325e5
root = "/astro/net/scratch1/christensen/DwarfResearch/ResMassTests";
;set_plot,'ps'
set_plot,'x'
cd,root

print,"Begining Graph"
for mct=0,4 do begin
;    device,filename=m[mct]+'_fourierStars.eps',/color, bits_per_pixel=8
    for rct=3,4 do begin
        cd,r[rct]    
        cd,m[mct]
        print,r[rct],"/",m[mct]
        input = "amps.stars"
        close,1
        openr,1,input
        readf,1,dataread
        data = [TRANSPOSE(dataread[*,6]),TRANSPOSE(dataread[*,0:5])]
        print,data
        if (rct eq 3) then plot,data[*,0],data[*,1],title = "Fourier Analysis of Stars "+trim(m[mct]),linestyle = 0,yrange=[0,0.3]  else oplot,data[*,0],data[*,1],linestyle = rct-3
        cd,'../..'
    endfor
legend,['10^5 ','10^4 '],linestyle=[0,1],/right
;    device,/close
    stop
endfor

for mct=0,4 do begin
;    device,filename=m[mct]+'_fourierGas.eps',/color, bits_per_pixel=8
    for rct=1,4 do begin
        cd,r[rct]    
        cd,m[mct]
        print,r[rct],"/",m[mct]
        input = "amps.gas"
        close,1
        openr,1,input
        readf,1,dataread
        data = [TRANSPOSE(dataread[*,6]),TRANSPOSE(dataread[*,0:5])]
        if (rct eq 1) then plot,data[*,0],data[*,1],title = "Fourier Analysis of Gas "+trim(m[mct]),linestyle = 0,yrange=[0,1] else oplot,data[*,0],data[*,1],linestyle = rct-1
        cd,'../..'
    endfor
legend,['10^5','10^4','1000','100'],linestyle=[0,1,2,3],/right
;    device,/close
    stop
endfor



END
