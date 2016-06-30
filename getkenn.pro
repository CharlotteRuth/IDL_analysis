function getkenn

dirs=[['1E4R/10M/o10M','1E5R/11M/o11M'],['1E4R/11M/o11M','1E5R/12M/o12M'],['1E4R/12M/o12M','1E5R/13M/o13M']]
 ;the number of particles per mass remains the same
m=2.362e5
t=1e9
sfr=fltarr(2,3);Starformation rate for each sim
gas=fltarr(2,3);Gas density for each sim

loadct,39 ;load color table
for i=0,1 do begin
    for j=0,2 do begin
;iterate through different scale simulations
        file=dirs[i,j]+'.00300'
print,file
        if (i eq 0) AND (j eq 0) then begin
            sfr[i,j]=0
            gas[i,j]=2.83297
        endif else begin
            sfrgas=sigsfrgas(file,mass=m,time=t)
            sfr[i,j]=sfrgas.sfr
            gas[i,j]=sfrgas.gas
        endelse
        print,sfr,' ',gas
    endfor
endfor

return,{sfr:sfr,gas:gas}

END
