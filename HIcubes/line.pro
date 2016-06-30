function line, x, y
if n_elements(x) ne 2 then begin
    print,"params = line(x, y)"
    print,"x and y must be two element vectors"
    return,0
endif
if n_elements(y) ne 2 then begin
    print,"params = line(x, y)"
    print,"x and y must be two element vectors"
    return,0
endif

b = float(y[0] - y[1]) / (x[0] - x[1])
a = y[1] - b * x[1]

return, [a, b]
end
