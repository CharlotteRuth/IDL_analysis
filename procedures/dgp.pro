
dir=['15','20','25']
lowsfr=['0600','0350','0800']
hisfr=['1000','0250','0200']
for i=0,n_elements(dir)-1 do begin
!p.multi=[0,2,1]
 paperps,file=dir[i]+'kms.rhoprof.eps'
  file=dir[i]+'kms/'+dir[i]+'kms.0'+lowsfr[i]
  lotime = lowsfr[i]/100
  gasprof,file,/all,tit=string(lotime)+' Gyr (after low SFR)',xmargin=[6.5,1.5]
  gasprof,file,/over,line=2
  file=dir[i]+'kms/'+dir[i]+'kms.0'+hisfr[i]
  hitime = hisfr[i]/100
  gasprof,file,/all,tit=string(hitime)+' Gyr (after hi SFR)',xmargin=[6.5,1.5]
  gasprof,file,/over,line=2
 ppsclose
endfor
!p.multi=[0,1,1]

END
