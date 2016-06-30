PRO bulgeSFH,filename

rtipsy, filename,h,g,d,s
readarr,filename+'.decomp',h,deco,part = 'star',/ascii
;readarr,filename+'.massform',h,massform,part = 'star',/ascii
;readarr,filename+'.iord',h,iord,part = 'star',/ascii,type = 'long'
;data = rstarlog('../../h986.cosmo50cmb.3072g14HBWK.starlog',/molecularH)
;match,data.iorderstar,iord,ind,temp
;Bulge: 3
;PsuedoBulge: 5
bulge = s[where(deco EQ 3 OR deco EQ 5,complement = otherind)]
other = s[otherind]

END
