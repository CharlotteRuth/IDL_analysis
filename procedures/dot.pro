
FUNCTION dot, a, b 
compile_opt strictarr 
IF n_elements(a) NE n_elements(b) THEN $ 
message, /ioerror, 'DOT(): the two vectors must have same number of elements' 
ans = total(a * b,1) 
return, ans 
END 
