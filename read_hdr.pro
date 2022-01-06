Function hdr2phys,ppar
; calculation the tags of the physically meaningful structure hdr
    c = 2.998e8
    xmt = 35.5e9
    prf= 2500. * (ppar.prf + 1)
    pulse_width = 100. * (ppar.pdr + 1L) * 1e-9
    delta_h = 0.5 * c * pulse_width
    mom_h_min   = ppar.ihp * delta_h
    mom_zrg = ppar.chg-2                     ; -2 is applied as chg includes the two noise 
                                             ; estimates which are appended after the higest range gate
    ; mom_height = [findgen(mom_zrg) * delta_h + mom_h_min,15000.,15000.]

    raw_h_min = (ppar.raw_gate1 )*delta_h
    ; für Daten vor dem 10.7.2004 (nasa) und vor dem 1.Sept (dwd)
    ;raw_h_min = (ppar.raw_gate1 + ppar.ihp )*delta_h 

    raw_zrg = ppar.raw_gate2 - ppar.raw_gate1+1
    ; raw_height = [findgen(raw_zrg) * delta_h + raw_h_min,15000.,15000.]

    prf = 2500. * (ppar.prf + 1)
    nfft= 128L * 2L ^ ppar.sft
    ave = float(ppar.avc) * float(nfft) / prf

    hdr={ $
        name:string(ppar.name),$        ; Original filename
        sn:'dwd',$                      ; device
        prf: prf, $                     ; Pulse repetition frequency
        pulse_width: pulse_width, $     ; puls diameter as choosen in the magnetron mode
        nfft: nfft, $                   ; number of FFT points
        avc: ppar.avc, $                ; number of spectral averages
        ave: ave, $                     ; averaging time
        c: c, $                         ; speed of light
        xmt: xmt, $                     ; transmitting frequency. This is a constant and not from the BITE data!
        delta_h: delta_h , $            ; range spacing
        mom_h_min: mom_h_min, $         ; height of the lowes range gate of the moment data
        mom_zrg: mom_zrg, $             ; number of range gates of the raw data
        raw_h_min: raw_h_min, $         ; height of the lowes range gate of the raw (or spectrum) data
        raw_zrg: raw_zrg, $             ; number of range gates of the raw data
        do_raw: ppar.raw, $             ; 0: no raw or spectra saving; 1: spectra; 2: raw data
        vuar:  c * prf / (2.0*xmt), $   ; diameter of the range of velocities which can be alocated without ambiguity 
                                    $   ; (Velocity Un-Ambigiuty Range)
        vny: c * prf / (4.0*xmt), $     ; nyquist velocity = vuar/2
        la: c / xmt, $                  ; wave length
        pol: ppar.pol  $                ; cross polarized data on/off
    }
    return,hdr
end

Function read_hdr, lun, hdr
; read the header which appears only once at the begining of each pds-file
common initialize,init_extract,init_header_counter
init_extract=1 ; forces recalculation of init_extract in read1pds_well
init_header_counter=1 ; be silent if header_indes dose not match

    ppar    = { name      : bytarr(32),$
            time         : bytarr(32),$
            oper         : bytarr(64),$
            place        : bytarr(128),$
            descr        : bytarr(256), $
            prf          : 0l, $
            pdr          : 0l, $
            sft          : 0l, $
            avc          : 0l, $
            ihp          : 0l, $
            chg          : 0l, $
            pol          : 0l, $
            att          : 0l, $
            tx           : 0l, $
            adcgain0     : 0., $
            adcgain1     : 0., $
            wnd          : 0l, $
            pos          : 0l, $
            add          : 0l, $
            len          : 0l, $
            cal          : 0l, $
            nos          : 0l, $
            of0          : 0l, $
            of1          : 0l, $
            swt          : 0l, $

            sum          : 0l, $
            osc          : 0l, $
            tst          : 0l, $
            cor          : 0l, $
            ofs          : 0l, $
            hdb          : 0l, $

            gainbblock    : 0., $

            calibrpower_m : 0., $
            calibrsnr_m   : 0., $
            calibrpower_s : 0., $
            calibrsnr_s   : 0., $

            raw_gate1     : 0l, $
            raw_gate2     : 0l, $
            raw           : 0l, $
            prc           : 0l $
    }
    dummy=bytarr(372)
    readu, lun, ppar, dummy

    hdr = hdr2phys(ppar)

    return, ppar
end
