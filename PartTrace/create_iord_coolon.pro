PRO create_iord_coolon,file,prevfile,nextfile,starlog
file = 'h603.cosmo50cmb.3072g14HBWK.00476'
prevfile = '../h603.cosmo50cmb.3072g14HBWK.00472.dir/h603.cosmo50cmb.3072g14HBWK.00472'
nextfile = '../h603.cosmo50cmb.3072g14HBWK.00480.dir/h603.cosmo50cmb.3072g14HBWK.00480'
starlog = '../../h603.cosmo50cmb.3072g14HBWK'

rtipsy,file,h,g,d,s
rtipsy,prevfile,h1,g1,d1,s1
rtipsy,nextfile,h2,g2,d2,s2

readarr,prevfile + '.iord',h1,i1,/ascii,type = 'long'
i1gas = i1[0:h1.ngas - 1]
i1dark = i1[h1.ngas:h1.ngas + h1.ndark - 1]
i1star = i1[h1.ngas + h1.ndark:h1.ngas + h1.ndark + h1.nstar - 1]
readarr,nextfile + '.iord',h2,i2,/ascii,type = 'long'
i2gas = i2[0:h2.ngas - 1]
i2dark = i2[h2.ngas:h2.ngas + h2.ndark - 1]
i2star = i2[h2.ngas + h2.ndark:h2.ngas + h2.ndark + h2.nstar - 1]
sl = rstarlog(starlog,/molecularH)

idark = i1dark
;sfind = where(sl.timeform 
istar = i2star[0:h.nstar - 1]
match2,i1gas,i2gas,ind1,ind2
dispind = where(ind1 EQ -1)
dispi = i1gas(dispind)

readarr,nextfile + '.igasorder',h2,i2gasord,part = 'star',/ascii,type = 'long'
i2starnew = i2star[h.nstar:h2.nstar - 1]
i2gasnew = i2gasord[h.nstar:h2.nstar - 1]
match2,dispi,i2gasnew,ind1,ind2
dispiearly = dispi[where(ind1 EQ -1)] ;Gas gone be the step in question
;dispilate = dispi[where(ind1 NE -1)]
match2,dispiearly,i1gas,ind1,ind2
igas = i1gas[where(ind2 EQ -1)]
iord = [h.n,igas,idark,istar]
writecol,file + '.iord',iord,format='(I)'
END
