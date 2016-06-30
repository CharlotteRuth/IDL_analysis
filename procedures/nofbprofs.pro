
loadct,39
 paperps,file='nofb.rhoevol.eps',/color
  file='15kms/15kms.std'
  gasprof,file,/all,tit='No FB',xmargin=[6.5,1.5]
 for j=1,20 do begin
  file='15kms/nofb/15kms.0'+string(j*50,format='(I04)')
print,file
  gasprof,file,/all,/over,color=(255-254*j/20)
 endfor
 ppsclose
!p.multi=[0,1,1]

END
