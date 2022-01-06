Function read_srvinfo, lun
; DWD version
; read the header which precedes each dwell. It should contain the important log parameters which
; may change from dwell to dwell
common read_srvinfo,frame_cnt
common initialize,init_extract,init_header_counter
    srvinfo = { frame_cnt  : 0ul, $
                 time_t : 0ul, $
                 tpow   : 0.,  $
                 npw1   : 0.,  $
                 npw2   : 0.,  $
                 cpw1   : 0.,  $
                 cpw2   : 0.,  $
                 ps_err : 0ul, $
                 te_err : 0ul, $
                 tt_err : 0ul, $
                 rc_err : 0ul  $
    }
    ens = n_elements(srvinfo)
    readu, lun, srvinfo
    if n_elements(srvinfo) ne ens then return,-1
    if init_header_counter then begin
        frame_cnt = srvinfo.frame_cnt - 1
        init_header_counter = 0
    endif
;    if srvinfo.frame_cnt ne frame_cnt+1 then $
;        printf,-2,'Skipped ',srvinfo.frame_cnt-frame_cnt-1,' dwells at dwell count',srvinfo.frame_cnt,' '
     frame_cnt=srvinfo.frame_cnt
    return, srvinfo
end

