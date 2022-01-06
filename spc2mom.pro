pro spc2mom,arguments
; 
; Reads the spectra from pds files and estimates 3 spectral moments of the nmulti most 
; dominant peaks. If nmulti=1 then it estimates the 3 spectral moments of the complete 
; spectra as the DSP software does.
;
; 
; unix> idl spc2mom infile=data/20040602_1000.pds nmulti=9
;
; IDL > spc2mom, 'infile=data/20040602_1000.pds'
; 
; infile      File to read the pds data from. Multible files with contain successive 
;             (in time) data may be given as array of strings.
;
; outfile     File for writing the output data.
;
; nmulti      number of most dominant peaks from which moments should be saved.
; 
;
; Matthias Bauer-Pfundstein bauer@metek.de
;
common startup,startup_done
if n_elements(startup_done) ne 1 then begin
@startup
endif

; Arguments can be provided by two methods:
; - If spc2mom is called by an interactive IDL session the first (and only) 
; string parameter arguments which contains a space separated list is split to
; the strarr argv.
; - Aditionally strings contained in the environment variables ARG1, ARG2, ..., 
; ARG9, ARG10, ... are prepended to argv. The script idl packs its calling arguments
; these environment variables. This is done as IDL can not extract the shell 
; variables like $1, $2,...
argv = env2argv(arguments)

kp = { $
    infile:'', $
    outfile:'' $
}
argv2stru,argv,kp

; Default values for the Ploting parameters. These Parameters can be changed by the command "set"
pp={ $
nmulti:1,$
of0:-1,$     
hilde_method:3,$ 
ofdiv:9.,$
ndiv:16,$
peak_trh:0.8 ,$ 
gap4dom:7, $
nnave:16 $
}
argv2stru,argv,pp

; online descripion for the above parameters
pph={ $
nmulti: "Number of  peaks for which the spectra should be saved" , $
of0:      "ofs0, parameter for detection threshold", $
hilde_method: "enumerate: 0 take noise from noise gate, 1 hildebrand-sekhon, 2 hildebrand-fischer-henemut, 3 cheap div",$
ofdiv: "cheap hildebrand, similar to ofs0 for determination of detection threshold. divisions from minimum to this thresold are regarded as noise",$
ndiv: "cheap hildebrand, number of parts to which the spectrum is divided",$
peak_trh: "Threshold for peak separation" ,$ 
gap4dom: "Peaks are separated in any case if more than gap4dom values are below detection thresold between them",$
nnave: "number of averages for the noise values, if hilde_method=0 is selected" $
}

; open the pds file and read the fileheader
get_lun,lun & openr,lun,kp.infile
get_lun,outlun & openw,outlun,kp.outfile
; read obsolete header which is in the first 1024 bytes of each file
ppar = read_fhdr1024(lun,hdr)
write_fhdr1024,outlun, ppar, pp.nmulti

; read the header which appears in the data each time any parameter was changed and at the begining
ppar = read_ppar_frame( lun, hdr)
write_ppar_frame,outlun, ppar, hdr.len_ppar
zchan=2 ; number of channels
nmulti=pp.nmulti
maxchg = 512L
nmom=3
zrg = -1L
lig3=lindgen(nmom)
nfft=0  ; nfft of the spectra in momory. If hdr.nfft differs then a new buffer for noisspectra is created.
while not eof(lun) do begin
    ; read one dwell of spectra
    while read_mainframe_header(lun, hdr) eq 0 do begin
        ppar = read_ppar_frame( lun, hdr)
        write_ppar_frame,outlun, ppar, hdr.len_ppar
    endwhile
    if zrg ne hdr.raw_zrg then begin
        ; create some arrays
        zrg = hdr.raw_zrg
        moments = fltarr(zchan, maxchg, nmulti, nmom)
    endif
    if hdr.nfft ne nfft then begin
        nfft=hdr.nfft
        noise_spectras = fltarr(nfft, zchan, pp.nnave)
        noise_values = fltarr(zchan, pp.nnave)
        inave=0
        nnave=0
        ix2ord = (lindgen(nfft) + nfft/2 ) mod nfft
        dx = hdr.VUAR/nfft
        x = cmod(dx * findgen(nfft), hdr.VUAR,0)
        ed_laprf = 1./(hdr.la * hdr.prf)
        xord=x[ix2ord]*ed_laprf
    endif
    if n_elements(p_srvinfo) then heap_free,p_srvinfo
    if n_elements(p_data) then heap_free,p_data
    ; read the spectra here:
    read1pds_dwell, lun, hdr, ppar, p_srvinfo, p_data, extract=['SP1','SP2']
    if n_elements( (*(*p_data).sp1)[*,0] ) ne nfft then begin
        if eof(lun) then goto,ende
        print,'corrupted data'
        stop
    endif
    noise_spectras[0,0,inave] = (*(*p_data).sp1)[*,hdr.raw_zrg]
    noise_spectras[0,1,inave] = (*(*p_data).sp2)[*,hdr.raw_zrg]
    noise_values[0,inave] = total(noise_spectras[*,0,inave])/float(nfft)
    noise_values[1,inave] = total(noise_spectras[*,1,inave])/float(nfft)
    inave = ++inave mod pp.nnave
    if nnave lt pp.nnave then ++nnave
    if n_elements(sp_tagnr) ne 2 then begin
        sp_tagnr = [(where(tag_names(*p_data) eq 'SP1'))[0],(where(tag_names(*p_data) eq 'SP2'))[0]]
        cal_tagnr = [(where(tag_names(*p_srvinfo) eq 'CPW1'))[0],(where(tag_names(*p_srvinfo) eq 'CPW2'))[0]]
        ; cal = (*p_srvinfo).(cal_tagnr[pp.chan])
    endif

    if pp.of0 gt 0 then of0 = pp.of0 else of0 = ppar.of0

    ; do the moment estimation for each range gate and each channel
    for ichan=0,zchan-1 do begin
        noise = total(noise_values[ichan,0:nnave-1])/float(nnave)
        if pp.hilde_method eq 0 then begin
            det_thres_from_noisegate = noise * ( 1 + of0/10. / sqrt(ppar.avc))
        endif else begin
            det_thres_fac = ( 1 + of0/10. / sqrt(ppar.avc))
        endelse
        noise_snr = noise*float(nfft)
        for irg = 0, zrg-1 do begin
            spc = (*(*p_data).(sp_tagnr[ichan]))[ix2ord,irg]
            if pp.hilde_method eq 0 then begin
                det_thres = det_thres_from_noisegate > max(spc) / linear(ppar.SWT)
            endif else begin
                noise = hilde(spc,method=pp.hilde_method, ofdiv=pp.ofdiv, ndiv=pp.ndiv, n_spcave=ppar.avc)
                det_thres = noise * det_thres_fac
            endelse
            if nmulti eq 1 then begin
                ; estimate moments of the complete spectrum as done by the dsp software
                ix_above_thres = where(spc ge det_thres, nix_above)
                if nix_above gt 0 then begin
                    ; noise removal and thresholding
                    spc_above_9noise = spc[ix_above_thres] - noise
                    xx = xord[ix_above_thres]
                    ; moment estimation
                    moments[ichan, irg, 0, lig3] = calc_moms(spc_above_9noise,xx,noise_snr)
                endif else begin
                    moments[ichan, irg, 0, lig3] = [0.0, -9999.0, -9999.0]
                endelse
            endif else begin
                ; use canion algorithm for identifying the most dominant peaks    
                spcord = spc[ix2ord] 
                spc_trh = (spcord - noise) * (spcord gt det_thres)
                peak_ifrqs = identify_peaks_by_canion( spc_trh, peak_trh_fac = pp.peak_trh/sqrt(ppar.avc), gap4dom=pp.gap4dom)
                npeaks=n_elements(peak_ifrqs[0,*]) 
                if peak_ifrqs[0,0] eq -1 then npeaks=0
                for ipeak=0,npeaks-1 do begin
                    ix4spc = ab2ix4spc(peak_ifrqs[0,ipeak],peak_ifrqs[1,ipeak],nfft)
                    spc_of_peak = spc_trh[ix4spc]
                    xx_of_peak = xord[ix4spc]
                    moments[ichan, irg, ipeak, lig3] = calc_moms(spc_of_peak, xx_of_peak,noise_snr)
                endfor
                for ipeak=npeaks, nmulti-1 do begin
                    moments[ichan, irg, ipeak, lig3] = [0.0, -9999.0, -9999.0]
                endfor
            endelse
        endfor
    endfor
    ; write the moments to file
    len_mainframe = 8L + hdr.len_srvi + 8L + 0L + 3L * (8L + maxchg * 4 * 2) + 8L 
    writeu,outlun,'FZKF',len_mainframe
    writeu,outlun,'SRVI',hdr.len_srvi,*p_srvinfo
    writeu,outlun,'FFTD',0L
    writeu,outlun,'SNRD', maxchg * 4 * 2 * nmulti & writeu,outlun,moments[*, *, *, 0]
    writeu,outlun,'VELD', maxchg * 4 * 2 * nmulti & writeu,outlun,moments[*, *, *, 1]
    writeu,outlun,'RMSD', maxchg * 4 * 2 * nmulti & writeu,outlun,moments[*, *, *, 2]
    writeu,outlun,'EXPD',0L

endwhile

ende:
close,lun
close,outlun
if n_elements(p_srvinfo) then heap_free,p_srvinfo
if n_elements(p_data) then heap_free,p_data
end
