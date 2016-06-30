
dir=['15','20','25','30']
;dir=['nofb']
loadct,39
for i=0,n_elements(dir)-1 do begin
 paperps,file=dir[i]+'kms.rhoevol.eps',/color
  file=dir[i]+'kms/'+dir[i]+'kms.00050'
  gasprof,file,/all,tit=dir[i]+' km/s',xmargin=[6.5,1.5]
 for j=1,20 do begin
  file=dir[i]+'kms/'+dir[i]+'kms.0'+string(j*50,format='(I04)')
print,file
  gasprof,file,/all,/over,color=(255-254*j/20)
 endfor
 ppsclose
endfor
!p.multi=[0,1,1]

END
