; docformat = 'rst'
;
; NAME:
;   SetDecomposedState
;
; PURPOSE:
;   Provides a device-independent way to set the color decomposition state of the
;   current graphics device. Devices that do not have a DECOMPOSED keyword to the
;   DEVICE command are assumed to be in indexed color mode always.
;
;******************************************************************************************;
;                                                                                          ;
;  Copyright (c) 2010, by Fanning Software Consulting, Inc. All rights reserved.           ;
;                                                                                          ;
;  Redistribution and use in source and binary forms, with or without                      ;
;  modification, are permitted provided that the following conditions are met:             ;
;                                                                                          ;
;      * Redistributions of source code must retain the above copyright                    ;
;        notice, this list of conditions and the following disclaimer.                     ;
;      * Redistributions in binary form must reproduce the above copyright                 ;
;        notice, this list of conditions and the following disclaimer in the               ;
;        documentation and/or other materials provided with the distribution.              ;
;      * Neither the name of Fanning Software Consulting, Inc. nor the names of its        ;
;        contributors may be used to endorse or promote products derived from this         ;
;        software without specific prior written permission.                               ;
;                                                                                          ;
;  THIS SOFTWARE IS PROVIDED BY FANNING SOFTWARE CONSULTING, INC. ''AS IS'' AND ANY        ;
;  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES    ;
;  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT     ;
;  SHALL FANNING SOFTWARE CONSULTING, INC. BE LIABLE FOR ANY DIRECT, INDIRECT,             ;
;  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED    ;
;  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;         ;
;  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND             ;
;  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT              ;
;  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS           ;
;  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                            ;
;******************************************************************************************;
;
;+
; :Description:
;   Provides a device-independent way to set the color decomposition state of the
;   current graphics device. Devices that do not have a DECOMPOSED keyword to the
;   DEVICE command are assumed to be in indexed color mode always.
;
; :Categories:
;    Graphics, Utilities
;    
; :Params:
;    state: in, required, type=integer, default=0
;         Set to 1 to set the current graphics device to decomposed color. Set to 0 
;         to set the current graphics device to indexed color. Devices lacking a 
;         DECOMPOSED keyword are assumed to be perpetually in indexed color mode.
;       
; :Keywords:
;     CurrentState: out, optional, type=integer
;         The current decomposition state of the current graphics device when the
;         program is called. A 1 indicates decomposed color. A 0 indicates indexed 
;         color.
;          
; :Examples:
;       IDL> SetDecomposedState, 0, CurrentState=mode
;       
; :Author:
;       FANNING SOFTWARE CONSULTING::
;           David W. Fanning 
;           1645 Sheely Drive
;           Fort Collins, CO 80526 USA
;           Phone: 970-221-0438
;           E-mail: davidf@dfanning.com
;           Coyote's Guide to IDL Programming: http://www.dfanning.com
;
; :History:
;     Change History::
;        Written, 16 November 2010. DWF.
;        Changes to include SET_PIXEL_DEPTH in Z-graphics buffer. 19 Nov 2010. DWF.
;        Allow PostScript 7.0 to set the decomposition keyword. 12 Dec 2010. DWF.
;
; :Copyright:
;     Copyright (c) 2010, Fanning Software Consulting, Inc.
;-
PRO SetDecomposedState, state, CURRENTSTATE=currentState

    Compile_Opt idl2
    
    ; Set to indexed color if you are not told differently.
    state = Keyword_Set(state)
    
    ; Need to know the current state?
    IF Arg_Present(currentState) THEN currentState = DecomposedColor()
    
    ; Set the decomposition state, if you are able. Otherwise assume
    ; indexed color mode.
    CASE StrUpCase(!D.Name) OF
    
       'PS': BEGIN ; PostScript
             IF Float(!Version.Release) GE 7.0 THEN Device, Decomposed=state
             END
             
        'Z': BEGIN ; Z-Graphics Buffer
             IF Float(!Version.Release) GE 6.4 THEN BEGIN
                CASE state OF
                    0: Device, Decomposed=state, Set_Pixel_Depth=8
                    1: Device, Decomposed=state, Set_Pixel_Depth=24
                 ENDCASE
             ENDIF
             END
             
        'MAC': BEGIN
             IF (Float(!Version.Release) GE 5.2) THEN Device, Decomposed=state
             END
             
         'X': Device, Decomposed=state
         
         'WIN': Device, Decomposed=state
         
         ELSE: Message, 'Unrecognized device. Assuming indexed color state.', /INFORMATIONAL
         
    ENDCASE
    
END    
