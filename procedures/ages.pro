dirs=['15','20','25','30']
m=2.325e5
t=1e9
ned=n_elements(dirs)
xmargin=[7.5,0.8]
ymargin=[3.0,0.2]
!p.multi=[0,2,2]
paperps,file='ages.eps',charsize=1.5
for j=0,ned-1 do begin
    file=dirs[j]+'kms/'+dirs[j]+'kms.01000'
print,file
    rtipsy,file,h,g,d,s
    sfr,s,sfh=sfh,maxt=maxt,massu=m,time=t,tarray=tarray,/age,/xlog, $
          xmargin=xmargin, ymargin=ymargin, xrange=[-0.5,9.5],xstyle=1
endfor
ppsclose
    
END

