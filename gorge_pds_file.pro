Pro gorge_pds_file, infile, ppar, hdr, parr_srvinfo, parr_data, extract=extract, bdate_begin=bdate_begin, $
    bdate_end=bdate_end,max_dwells=max_dwells,progress_bar=progress_bar,max_mem=max_mem,lun=lun
;
; Reads a complete data file (*.pds) produced by the cloud radar MIRA36 . 
;
; infile        string Filename of file to be read in.
; lun           dummy parameter which is needed in the gorge_pds_file routine of the fzk radar and later
; ppar          header structure with appears onc at the begining of the pds file
; hdr            structure with human readable parameters calculated from ppar
; parr_srvinfo  array of pointers which point to the serverinfo structure of each dwell
; parr_data     array of pointers which point to the profiles of radar data structure of each dwell
; extract       optional array of strings which indicates which profiles should be read and
;               which should be skipped. it should be a subset of
;                      extrac=['p1','v1','r1','p2','v2','r2','sp1','sp2'].
;               If not given all available profiles are read.
; bdate_begin   see hour_end
; hour_end      only dwells within the the given time period are saved in the moment and spectra structure
;
; After reading the file header by the function ppar=read_hdr(lun) swallow_pds_file performes one call of 
;     read1pds_dwell, lun, hdr, ppar, p_srvinfo, p_data, extract=extract
; for each dwell (averaging time).
;
; Examples for obtaining Data from parr_srvinfo  and  parr_data:
;
; IDL > itime=100 & irg=30                                                     ; dwell and range gate number used below
; 
; IDL > utime = lonarr(n_elements(PARR_SRVINFO))                               ; read the time stamps to utime
; IDL > for itime=0, n_elements(PARR_SRVINFO)-1 do utime[itime] = (*parr_srvinfo[itime]).TIME_T
;
; IDL > plot, hdr.mom_height, db((*parr_data[itime]).P1),yrange=[-40,70],yst=1 ; plot SNR-Profile of dwell number itime 
;   
; IDL > hight_profile_of_cross_spectra = (*(*parr_data[itime]).SP2) 
; IDL > plot,dB(hight_profile_of_cross_spectra[*,irg])                         ; plot one spectum
;
; !!! See also un_gorge.pro which is a convenient tool for extracitng slices of parr_data.
; 
; If the pointer arrays parr_srvinfo and parr_data already contain data then the new data is appended.
; The structures ppar and hdr will be over-written nevertheless.
; 
; !!! If PARR_SRVINFO and PARR_DATA should be reused without appending then free heap by:
;
; IDL > heap_free,PARR_SRVINFO & heap_free,PARR_DATA
; 
; Description of structures which are provided by gorge_pds_file:
;
;   ppar is the header structure as described in the pds file documentation:
; 
;   hdr={ $
;       prf: 2500. * (ppar.prf + 1)     ; Pulse repetition frequency
;       pulse_width: (ppar.pdr+1)*1e-7  ; puls diameter as choosen in the magnetron mode
;       nfft: 128 * (ppar.sft + 1)      ; number of FFT points
;       c: c                            ; speed of light
;       xmt: xmt                        ; transmitting frequency. This is a constant and not from the BITE data!
;       delta_h: 0.5 * c * pulse_width  ; range spacing
;       mom_h_min: ppar.ihp * delta_h   ; height of the lowes range gate of the moment data
;       mom_zrg: ppar.chg-2             ; number of range gates of the raw data
;       mom_height: findgen(mom_zrg)*delta_h+mom_h_min ; Vector with range values for each range gate of mom data
;       raw_h_min:(ppar.raw_gate1+ppar.ihp)*delta_h  ; height of the lowes range gate of the raw (or spectrum) data
;       raw_zrg: ppar.raw_gate2-ppar.raw_gate1+1  ; number of range gates of the raw data
;       raw_height:findgen(raw_zrg)*delta_h+raw_h_min ; Vector with range values for each range gate of raw data
;       do_raw: ppar.raw, $             ; 0: no raw or spectra saving; 1: spectra; 2: raw data
;       vuar:  c * prf / (2.0*xmt), $   ; diameter of the range of velocities which can be alocated without ambiguity
;                                   $   ; (Velocity Un-Ambigiuty Range)
;       vny: c * prf / (4.0*xmt), $     ; nyquist velocity = vuar/2
;       la: c / xmt, $                  ; wave length
;       pol: ppar.pol  $                ; cross polarized data on/off
;   }
;
;   parr_srvinfo is an array of pointers. Each elements points to the following structure which contains the
;   _server_info structure which is decribed in the pds file description
;   *parr_srvinfo[itime] = { 
;                time_t :   unix time
;                tpow   :   Power measured by thermistor from the BITE data
;                npw1   :   integral power of the spectrum at the noise gate * 356515840.   co-polarized channel 
;                npw2   :   integral power of the spectrum at the noise gate * 356515840.   cross-polarized channel
;                           The noise gate is at a height where no signal is expected.
;                           npwx is an average of the most recent 16 spectra. 
;                           All spectra (including the calib. noise spectra) are normalized to npwx/(356515840*nfft)   
;                cpw1   :   integral power of the spectrum of the calib. noise source , co-polarized channel
;                cpw2   :   integral power of the spectrum of the calib. noise source , cross-polarized channel
;                           This signal is introduced at the end of every pulse cycle by pin switches.
;                ps_err :   Error flag of the power supply unit.
;                te_err :   Error flags of the transmitter. 
;                tt_err :   Flags which shows due to which threshold the transmitter was switched off
;                rc_err :   Error flags of the receiver
;   }
;
;   parr_data is an array of pointers. Each elements points to the following structure which contains the 
;   data of one dwell:
;    *parr_data[itime] = { 
;       p1  : p1      ; SNR profile               co-polarized signal
;       v1  : v1      ; Velocity profile          co-polarized signal
;       r1  : r1      ; RMS (Peak width) profile  co-polarized signal
;       p2  : p2      ; SNR profile               cross-polarized signal
;       v2  : v2      ; Velocity profile          cross-polarized signal
;       r2  : r2      ; RMS (Peak width) profile  cross-polarized signal
;       sp1 : sp1     ; pointer to profile of spectra *sp1 = spc[ifft,irg]  co-polarized signal
;       sp2 : sp2     ; pointer to profile of spectra *sp1 = spc[ifft,irg]  cross-polarized signal
;   }
;
;
;
; U. Görsdorf, Matthias Bauer-Pfundstein bauer@metek.de
get_lun,lun & openr,lun,infile

ppar = read_hdr(lun,hdr)

default_value,max_dwells,7000L
default_value,max_mem,500000000L
default_value,progress_bar,0

itime = n_elements(parr_srvinfo)
if itime gt 0 then itime = itime * ptr_valid(parr_srvinfo[0])
if itime gt 0 and n_elements(parr_data) eq itime then begin
    parr_srvinfo = [parr_srvinfo, ptrarr(max_dwells)]
    parr_data    = [parr_data,    ptrarr(max_dwells)]
    smax_dwells = n_elements(parr_srvinfo)
endif else begin
    parr_srvinfo = ptrarr(max_dwells)
    parr_data = ptrarr(max_dwells)
    itime=0
    smax_dwells=max_dwells
endelse

skip=0
ocur_ptr = (fstat(lun)).cur_ptr
bytes_dwell=0L
while not eof(lun) do begin   ; dwells
    if memory(/current) ge max_mem then begin
        printf,-2,'max_mem = ',max_mem,' is too small'
        break
    endif
    if itime ge smax_dwells then begin
        printf,-2,'max_dwells too small'
        break
    endif
    fstatx = fstat(lun)
    if fstatx.cur_ptr ne ocur_ptr then begin
        bytes_dwell = fstatx.cur_ptr - ocur_ptr
        ocur_ptr = fstatx.cur_ptr
    endif
    if fstatx.size-fstatx.cur_ptr lt bytes_dwell then break
    if progress_bar and bytes_dwell ne 0 then progress_bar,fstatx.size/bytes_dwell,fstatx.cur_ptr/bytes_dwell
    read1pds_dwell, lun, hdr, ppar, p_srvinfo, p_data, extract=extract, skip=skip
    unixtime_of_dwell = (*p_srvinfo).time_t
    bdate_of_dwell = unixtime2bdate(unixtime_of_dwell,/utc)  
    if n_elements(unixtime_begin) eq 0 then begin
         unixtime_begin = unixtime_of_dwell
         if KEYWORD_SET(bdate_begin) then begin
            if bdate_begin ne 'default' then begin
                bdate_begin = complete_bdate(bdate_begin,bdate_of_dwell)
                unixtime_begin = bdate2unixtime(bdate_begin)
            endif else unixtime_begin = unixtime_of_dwell
         endif

         unixtime_end = unixtime_of_dwell + 1000000L
         if KEYWORD_SET(bdate_end) then begin
            if bdate_end ne 'default' then begin
                bdate_end = complete_bdate(bdate_end,bdate_of_dwell)
                unixtime_end = bdate2unixtime(bdate_end)
            endif
         endif
    endif

    if unixtime_of_dwell ge unixtime_begin and unixtime_of_dwell le unixtime_end then begin
        if skip then skip=0 else begin
            parr_srvinfo[itime] = p_srvinfo
            parr_data[itime]    = p_data
            itime = itime + 1
        endelse
    endif  else begin
        skip=1
    endelse
endwhile
if itime gt 0 then begin
    parr_srvinfo=parr_srvinfo[0:itime-1]
    parr_data=parr_data[0:itime-1]
endif else begin
    printf,-2,'gorged no data from '+infile
    parr_srvinfo=parr_srvinfo[0:0]
    parr_data=parr_data[0:0]
endelse

free_lun,lun & lun=-1000
end
