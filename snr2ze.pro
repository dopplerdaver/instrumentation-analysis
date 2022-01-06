function snr2ze, snr, hdr, txn_power, raw=raw, rg_offs=rg_offs
; Calculate the Reflectivity Ze from SNR (linear).
; 
; Assuming Raighley conditions and no attenuation the dBZ reflectivity
; represents the integral over the (D^6)-weighted droplet
; distribution, where D is the Droplet diameter. 
;
;     Z = int_over_D( n(D)*D^6 ) 
;
; For a small Droplet size range deltaD the droplet distribution is defined by:
;
;     deltaD * n(D) = Number of drops per m^3 with Diameter +/- deltaD/2
;
; For further information see the documentation.
;
; Matthias Bauer
;

; parameters of radar calibration
    calparams = {         $
        Pt0    : 30000. , $  ; W Puls Leistung
        tau0   : 200e-9 , $  ; s
        prf    : 5000.  , $  ; Hz
        h0     : 5000.  , $  ; m
        c0     : 0.00724  $  ; -21.4
    }
    ; see the MIRA36 documentation section 2.6 "radar equ. of meteorol. radar"
    ; 6.12.2005: nachdem ich die Doku überarbeitet habe, komme ich auf c0=-26.1.
    ;     Der wesentliche Fehler war dass ich die Verluste im Sende-Pfad vergessen
    ;     hatte. Damit die Bilder http://metekgmbh.dyndns.org/fzk/pngs/2005/*/*.z1.png
    ;     konsisten bleiben lasse ich c0=-26.9
    ; 28.11.2006 fuers DWD radar sollten wir uns mal einigen: 
    ;     -21.7 hab ich berechen, 
    ;     -22.4 nehme ich seit es viewpds gibt, 
    ;     -21.9 hat Uli bei MRR-Vergleichen rausgefunden 
    ;     -21.4 fuer NASA etwas realistischere Beruecksichtigung der Verluste in der Antenne
    ;     

    if n_elements(raw) ne 1 then raw=0
    if n_elements(rg_offs) ne 1 then rg_offs=0

    ; pulse lenght in s
    tau0 = calparams.tau0        ; as used in calculation of calparams.c0
    tau	 = hdr.pulse_width       ; actually selected

    ; mean transmitted power in W
    Ptav0 = calparams.Pt0 * calparams.tau0 * calparams.prf  ; as used in calculation of calparams.c0
    Ptav  = txn_power * 27.8/34.0 / 100.           ; actually indicated by the thermistor
                                                   ; 27.8/34 NASA entsprechend den Messungen vom 4.11.2008

    ;reference height
    h0    = calparams.h0
    if raw eq 0 then h = (findgen(hdr.mom_zrg)+rg_offs) * hdr.delta_h + hdr.mom_h_min $
                else h = (findgen(hdr.raw_zrg)+rg_offs) * hdr.delta_h + hdr.raw_h_min

    Ze = snr * calparams.c0 * Ptav0 * tau0 * (h/h0)^2 / (Ptav*tau) 

    return, Ze
end

