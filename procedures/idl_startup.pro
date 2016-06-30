
;--------------------idl libraries---------------
;------------mycodes-------------

DEFINE_KEY, ESC=string([27B, 79B, 68B]), 'ALT_L', /BACK_CHAR
DEFINE_KEY, ESC=string([27B, 79B, 67B]), 'ALT_R', /FORWARD_CHAR
DEFINE_KEY, ESC=string([27B, 79B, 65B]), 'ALT_UP', /PREVIOUS_LINE
DEFINE_KEY, ESC=string([27B, 79B, 66B]), 'ALT_DN', /NEXT_LINE

  ;Set backing store to be handled by IDL
  ;(so window redraws when uncovered)

device, retain = 2

  ;Set IDL to use the true color (24 bit) visual class
  ;instead of direct color, which it will default to, if available.

;device, true_color = 24

  ;Have IDL colors behave like it's using a pseudo color, 8 bit display.
  ;You can set decomposed = 1 to go back to true color behavior.
  ;(mainly for use with LOADCT, etc.)

device, decomposed = 0









