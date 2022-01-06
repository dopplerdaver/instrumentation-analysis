pro remove_clutter,spc,noise,det_thres,pp
    ; remove clutter from a spectrum which is distrubuted symmetrically around zero velocity (index 0 in spc)
    ; spc = Spectrum of one range gate, zero velocity at index 0
    ; noise = Noise level
    ; det_thres = detection threshold
    ; pp = Structure containing the following Parameters:
    ; Parameters tested with Data from September + October 2005:
    ; pp.max_sholder_dbc=-0.20  ; theoretically the shoulders should be -6 dbc.
    ;                  Non clutter signal may cover the the clutter. The typically one of the shoulders is higer.
    ;                  In this case no clutter removal is done.
    ; pp.fftwin=[0.75, 0.3, 0.16, 0.065, 0.023, 0.008, 0.0025] ; decaying factor for the neighbors
    n = n_elements(spc)
    nop = n_elements(pp.fftwin)
    pc = spc[0]   ; power of the spectral component at zero Doppler shift
    signif_neg = 2*noise - det_thres
    pc9n = pc - noise 
    ;plot,db(spc),yra=[-7,30]
    if pc ge det_thres then begin   ; just for saving cpu time
        max_sholders = spc[1] > spc[n-1]  
        if (max_sholders-noise) / pc9n lt linear(pp.max_sholder_dbc) then begin
            spc[0] = noise 
            for i=0,nop-1 do begin
                pcs = pc9n * pp.fftwin[i]
                spc[i+1] = (spc[i+1] - pcs) > signif_neg
                spc[n-i-1] = (spc[n-i-1] - pcs) > signif_neg
            endfor
        endif
    endif
    ;oplot,db(spc),col=1
    ;stop
end
