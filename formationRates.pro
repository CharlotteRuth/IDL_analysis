pro formationRates,directory

cd,directory

data0 = PTR_NEW(FLTARR(3,5000))
data1 = PTR_NEW(FLTARR(3,5000))
data2 = PTR_NEW(FLTARR(3,5000))
data3 = PTR_NEW(FLTARR(3,5000))
;data4 = PTR_NEW(FLTARR(3,5000))

a=['dESN0','dESN.001','dESN.005','dESN.01']
;'dESN.1']
sizea = size(a)
lengtha = sizea[1] - 1
set_plot, 'ps'
device,filename = directory + 'starForm.ps'
loadct,12
FOR i = 0,lengtha DO BEGIN
    data = [data0,data1,data2,data3]
    file = a[i]+'/StarFormation.txt'
    OPENR,1,file
    j = 0
    WHILE (NOT EOF(1)) DO BEGIN
        READF, 1, s, m, t
        k = j/10
        (*data[i])[0,k] = s + (*data[i])[0,k]
        (*data[i])[1,k] = m + (*data[i])[1,k]
        (*data[i])[2,k] = t
 ;       IF k lt 10 then print,k,s,m,t,', ',(*data[i])[2,k]
        j = j + 1
    ENDWHILE
;    print,size(*data[i])
;    print,(*data[i])[0:2,0:10]
    (*data[i]) = (*data[i])[0:2,0:k]
;    print,(*data[i])[0:2,1590:k]
;    print,size(*data[i])
;    print,(*data[i])[0:2,0:10]
;   IF i eq 0 then print,data[i,2,*],xrange=[0.200000,0.239]
    IF i eq 0 then plot,(*data[0])[2,*],(*data[0])[0,*],xrange=[0.2,0.239], $
      xtitle = 'Mass of Stars Formed',ytitle = 'dTime',title = directory
   ;   color = 130
    IF i ne 0 then oplot,(*data[i])[2,*],(*data[i])[0,*], linestyle = i
 ;     color = 110 - 20*i
    CLOSE,1
ENDFOR
legend,a,linestyle = [0,1,2,3], colors=[130,90,70,50,30],/right
device,/close
cd,'..'
END
