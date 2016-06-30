function makeheader, cubehead



sxaddpar, header, 'CTYPE1', 'KPC-LIN'
sxaddpar, header, 'CUNIT1', 'KPC'
sxaddpar, header, 'CRVAL1', cubehead.crval1
sxaddpar, header, 'CDELT1', cubehead.cdelt1
sxaddpar, header, 'CRPIX1', cubehead.crpix1

sxaddpar, header, 'CTYPE2', 'KPC-LIN'
sxaddpar, header, 'CUNIT2', 'KPC'
sxaddpar, header, 'CRVAL2', cubehead.crval2
sxaddpar, header, 'CDELT2', cubehead.cdelt2
sxaddpar, header, 'CRPIX2', cubehead.crpix2

return,header

end
