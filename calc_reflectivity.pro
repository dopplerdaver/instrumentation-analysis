Pro calc_reflectivity, _hdr, _mom, _calparams, rfl_cop, rfl_crp, rfl_min

; Calculates the reflectivity factor Z in DBZ which is equal 10log(Z[mm6/m3])
;
; Parameters:
;
;	_hdr	: structure contains header information
;	_mom	: structure contains moment valuesa and technical information
;_calparams : structure contains radar constants and the corresponding values of tau, Pt
;
; For further information see documentation and user manual.
;
; UG, 18.02.2004



; pulse lenght
tau0 	= _calparams.tau0 ; in m
tau		= _hdr.pdr *1e-9 ; in m

;receiver bandwidth
B0		= 1./tau0
B		= 1./tau

;transmitted peak power
Pt		= *_mom.tpow/100  ; mean transmitted power in W
Pt0		= _calparams.Pt0 ; Peak power in W

Pt 		= Pt/(_hdr.prf * tau ); peak power (mean Pt / (Pulse repetition frequency * pulse lenght)

; receiver gain (only different from zero when receiver operates with attenuation)
Gr		= 1
Gr0		= 1

; radar constant for given values
rc0		= _calparams.c0

;reference height
h0		= _calparams.h0
h		= _mom.height


n_dwells = n_elements(pt)
n_height = n_elements(h)

rfl_cop = replicate(!values.f_nan, n_dwells, n_height)
rfl_crp = replicate(!values.f_nan, n_dwells, n_height)
rfl_min = replicate(!values.f_nan, n_dwells, n_height)


snr_cop 		= *_mom.snr_cop
snr_crp 		= *_mom.snr_crp


rc		= rc0 - 10*alog10(tau/tau0) - 10*alog10(B0/B) - 10*alog10(Gr/Gr0) - 10*alog10(Pt/Pt0)

; calculation of the smallest detectable snr (due to signal processing); there is to add an offset value _hdr.of1
snr_min = 5/(_hdr.sft * sqrt(_hdr.avc)) + _hdr.of1


for i = 0, n_dwells-1 do begin

  	rfl_cop[i,*] 	= rc[i] + 20*alog10(h/h0) + 10*alog10(snr_cop[i,*])
	rfl_crp[i,*] 	= rc[i] + 20*alog10(h/h0) + 10*alog10(snr_crp[i,*])
	rfl_min[i,*]	= rc[i] + 20*alog10(h/h0) + 10*alog10(snr_min)

endfor


end