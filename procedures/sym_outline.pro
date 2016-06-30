;this function takes a symbol code for a filled symbol and returns the
;symbol that is the unfilled version

FUNCTION sym_outline,sym

outline = sym
CASE sym OF
   14: outline = 4 ; Diamond
   15: outline = 6 ; Square
   16: outline = 9 ; Circle
   17: outline = 5 ; Upward triangle
   18: outline = 11 ; Downward triangle
   19: outline = 12 ; Right-facing triangle
   20: outline = 13 ; Left-facing triangle
   22: outline = 21 ; Hourglass  
   24: outline = 23 ; Bowtie
   26: outline = 25 ; Standing bar
   28: outline = 27 ; Laying bar
   34: outline = 33 ; Big cross
   38: outline = 37 ; Upper half circle
   40: outline = 39 ; Lower half circle
   42: outline = 42 ; Left half circle
   44: outline = 43 ; Right half circle
   46: outline = 45 ; Star
ELSE: outline = sym
ENDCASE

RETURN,outline
END
