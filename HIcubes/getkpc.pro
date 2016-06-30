pro getkpc, header, xaxis=xaxis, yaxis=yaxis, vaxis=vaxis

xaxis = dindgen(header.naxis1) * header.cdelt1 + header.crval1
yaxis = dindgen(header.naxis2) * header.cdelt2 + header.crval2
vaxis = dindgen(header.naxis3) * header.cdelt3 + header.crval3

end
