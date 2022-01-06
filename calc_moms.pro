function calc_moms,spc,xx,noise_snr
; calculate the three spectral moments from the selected spectral points spc which are at the given positions xx
    pwx = total(spc)
    if pwx gt 0.0 then begin
        vlx = total(spc*xx)/pwx
        vvx = total(spc*xx*xx)/pwx
        if n_elements(spc) gt 1 then begin
            return, [pwx/noise_snr, -vlx,sqrt(vvx - vlx*vlx)]
        endif else begin 
            return, [pwx/noise_snr, -vlx, 0.0] 
        endelse
    endif else return,[pwx,-9999.,-9999.]
end
