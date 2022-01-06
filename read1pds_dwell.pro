Function read_cop_header, lun
; read in the rather the header which precedes the profiles of moments from the cop or the crp channel
; ? This header may be useful to skip the whole cop or the whole crp section. But siz is 8 bytes too big!!
    sgnt = 'FCOP'
    if ~ find_sgnt(lun,sgnt) then print,-2,'Signature FCOP not foud'
    siz = 0l & readu,lun,siz
    _cophdr = { sig     : sgnt, $
                siz   : siz}
    return, _cophdr
end

Function read_crp_header, lun
; read in the rather the header which precedes the profiles of moments from the cop or the crp channel
; ? This header may be useful to skip the whole cop or the whole crp section. But siz is 8 bytes too big!!
    sgnt = 'FCRP'
    if ~ find_sgnt(lun,sgnt) then print,-2,'Signature FCRP not foud'
    siz = 0l & readu,lun,siz
    _cophdr = { sig     : sgnt, $
                siz   : siz}
    return, _cophdr
end


Function read_mom_section, lun, taget_sgnt, do_read
; read one mom_profile on the moment data (i.e. cop_snr)
    if ~ find_sgnt(lun,taget_sgnt) then printf,'-2','Signature "'+taget_sgnt+'" expected but not found'    
    siz = 0l & readu,lun,siz
    
    if do_read then begin
        zrg_pz   = siz / 4     ; _pz means Plus Zwei <=> +2
        if zrg_pz le 0 then mom_profile = !values.f_nan else begin
            mom_profile = fltarr(zrg_pz,/nozero)
            readu, lun, mom_profile
        endelse
    endif else begin
        ;fstatx=fstat(lun)
        ;point_lun, lun, fstatx.cur_ptr + siz
        seek,lun,siz
        mom_profile = !values.f_nan
    endelse
    return, mom_profile
end


Function read_fft_raw_section, lun,hdr, do_read
; read one profile of I-Q data or spectra
    if ~ find_sgnt(lun,'FFFT') then print,-2,'Signature "FFFT" not found'
    siz = 0l & readu,lun,siz

    raw_zrg_pz = hdr.raw_zrg+2
    expected_siz = hdr.nfft*raw_zrg_pz*4*hdr.do_raw
    if expected_siz ne siz then printf,-2,'raw data has a size of ',siz,$
        ' expected',expected_siz

    if hdr.do_raw eq 0 then raw_pofile = !values.f_nan else begin
        if do_read then begin
            if hdr.do_raw eq 1 then begin ;spectra
                raw_pofile = ptr_new(fltarr(hdr.nfft,raw_zrg_pz,/nozero) )
                readu, lun, *raw_pofile
            endif

            if hdr.do_raw eq 2 then begin ; I and Q values
                raw_pofile = ptr_new(replicate({i: 0, q: 0}, siz/4))
                readu, lun, *raw_pofile
            endif
        endif else begin
            ;fstatx=fstat(lun)
            ;point_lun, lun, fstatx.cur_ptr+siz
            seek,lun,siz
            raw_pofile = ptr_new(!values.f_nan)
        endelse
    endelse 

    return, raw_pofile
end


Pro read1pds_dwell, lun, hdr, ppar, p_srvinfo, p_data, extract=extract, skip=skip
; read one dwell from pds-files
; Parameter
;   lun         logical unit of the file which has opened before
;   hdr          header information, obtained before by read_hdr(lun)
;   p_srvinfo   pointer to struct with dwell header containing information from the log data
;   p_data      pointer to struct with returned data
;   extract     optional array of strings which indicates which profiles should be read and 
;               which should be skipped. it should be a subset of
;                      extrac=['p1','v1','r1','p2','v2','r2','sp1','sp2'].
;               If not given all available profiles are read.
;   skip        skip data from this dwell, this cannot be achieved by setting extract=['']
;               as extract is evaluated only at the first call of read1pds_dwell after reading 
;               the file header.
;
;  See the source of gorge_pds_file.pro for information about the structures
;
; U. Görsdorf, Matthias Bauer-Pfundstein bauer@metek.de
;
    if n_elements(skip) le 0 then skip=0
    dwd_me = hdr.pol eq 1 ? -1 : 1 

    common initialize,init_extract,init_header_counter
    common read1pds_dwell,extr
    if n_elements(init_extract) eq 0 then begin
        printf,-2,'read_hdr has to be called before read1pds_dwell'
        init_extract = 0
    endif
    if init_extract then begin
        extract_all=['p1','v1','r1','p2','v2','r2','sp1','sp2']
        if n_elements(extract) le 0 then extract=extract_all
        if extract[0] eq 'default' then extract=extract_all
        extract = strupcase(extract)
        extr = { p1: max(extract eq 'P1' ), $
                 v1: max(extract eq 'V1' ), $
                 r1: max(extract eq 'R1' ), $
                 p2: max(extract eq 'P2' ), $
                 v2: max(extract eq 'V2' ), $
                 r2: max(extract eq 'R2' ), $
                 sp1:max(extract eq 'SP1'), $
                 sp2:max(extract eq 'SP2')  $
        }
    endif

    srvi_sgnt='FHDR' 
    if ~ find_sgnt(lun,srvi_sgnt) then begin
        print,'Signature '+srvi_sgnt+' not found'
        goto,ioerr
    endif
    siz=0L & readu,lun,siz
    p_srvinfo = ptr_new(read_srvinfo(lun))

    ; reading header and data of co-polarizised channel (cop)
    _cophdr = read_cop_header(lun)
        ;print, '_sig ', _cophdr.sig,' _siz ', _cophdr.siz

    sp1 = read_fft_raw_section(lun, hdr, extr.sp1 and skip eq 0)

    p1 =                             read_mom_section(lun, 'FINT', extr.p1 and skip eq 0)
    v1 = dwd_me * hdr.la * hdr.prf * read_mom_section(lun, 'FVEL', extr.v1 and skip eq 0)
    r1 =          hdr.la * hdr.prf * read_mom_section(lun, 'FRMS', extr.r1 and skip eq 0)

    if hdr.pol eq 1 then begin
         ; reading header and data of cross-polarizised channel (crp)
         _crphdr = read_crp_header(lun)
             ;print, '_sig ', _crphdr.sig,' _siz ', _crphdr.siz

         sp2 = read_fft_raw_section(lun, hdr, extr.sp2 and skip eq 0)

         p2 =                     read_mom_section(lun, 'FINT', extr.p2 and skip eq 0)
         v2 = dwd_me * hdr.la * hdr.prf * read_mom_section(lun, 'FVEL', extr.v2 and skip eq 0)
         r2 =       hdr.la * hdr.prf * read_mom_section(lun, 'FRMS', extr.r2 and skip eq 0)
    endif else begin
        sp2 = !values.f_nan
        p2 = !values.f_nan
        v2 = !values.f_nan
        r2 = !values.f_nan
    endelse

    p_data = ptr_new({ $
        p1  : p1  ,   $
        v1  : v1  ,   $
        r1  : r1  ,   $
        p2  : p2  ,   $
        v2  : v2  ,   $
        r2  : r2  ,   $
        sp1 : sp1 ,   $
        sp2 : sp2     $
    })
ioerr:
end

