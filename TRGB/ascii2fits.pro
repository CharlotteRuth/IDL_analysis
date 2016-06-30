pro ascii2fits,asciifile,templatefile

restore,templatefile
data = my_read_ascii( asciifile, template=template )
mwrfits, data, asciifile+'.fits', /create

end
