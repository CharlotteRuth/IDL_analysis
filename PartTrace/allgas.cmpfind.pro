allgasiord = mrdfits('../alyson/grp1.allgas.iord.fits',0)
iord = read_ascii_array('MW1.1024g1bwK.00512.iord')
igasord = read_ascii_array('MW1.1024g1bwK.00512.igasorder')
rtipsy, 'MW1.1024g1bwK.00512', h,g,d,s
igasord = igasord[(h.ngas+h.ndark):(h.n-1)]
gasind = binfind(iord, allgasiord)
gasind = gasind(where(gasind NE -1))
starind2 = binfind(allgasiord, igasord)
starind = where(starind2 NE -1)
ind = [gasind, starind]
cmp = read_ascii_array('MW1.1024g1bwK.00512.cmp')
cmpind = cmp[ind]
cmpstar = cmp[starind]
cmpgas = cmp[gasind]
starind3 = starind2(where(starind2 NE -1))
iordgas = iord[gasind]
gasind2 = binfind(allgasiord, iordgas)

;Too complicated.  Trying something else.
stop

end

