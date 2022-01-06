pro  limit_itr2contiguous,hours,max_ave,pp
; soubroutine of viewpds 
; which limits the time range (pp.itr) to a piece with contiguous data
; (no gap in the time sequence)
    if pp.itr[1]-pp.itr[0] gt 2 then begin 
        delta_hour = hours[pp.itr[0]+1:pp.itr[1]] - hours[pp.itr[0]:pp.itr[1]-1]
        ix_skip = where(delta_hour*36d2 gt max_ave*pp.tolerable_gap, nix_skip)
        if nix_skip gt 0 then begin ; an un-tolerable gap in the data was detected
            print,'Number of gaps in the time sequence:',nix_skip
            print,'Depending on takeat_gap the time range will be set to the first or last piece of contiguous data'
            print,'tolerable_gap=',pp.tolerable_gap,'  takeat_gap=', pp.takeat_gap 
            case pp.takeat_gap of
                'left': pp.itr[1]=pp.itr[0]+ix_skip[0]  
                'right': pp.itr[0]=pp.itr[0]+ix_skip[nix_skip-1]+1
                else:
            endcase
        endif
    endif
end

pro viewpds,arguments
; Tool for visualizing profiles of moments and Spectra of the MIRA-36 cloud radar
; provided as pdf files which are saved by the radarserver. 
;
; viewpds can be called from the IDL prompt or from a unix shell:
; 
; unix> idl viewpds infile=data/20040602_1000.pds:data/20040602_1200.pds \
;                extract=p1:p2:sp1 bdate_begin=DT120000 bdate_end=180000
;
; IDL > viewpds, 'infile=data/20040602_1000.pds extract=p1:p2:sp1 bdate_begin=110000'
; 
; infile      File to read the pds data from. Multible files with contain successive 
;             (in time) data may be given as array of strings.
;
; extract     An array of tag names which specify which profiles should be extracted.
;             If omitted all profiles which are in the pds file are extracted.
;             If the file is biger than the system memory it is useful to extract
;             only a subset of the profiles. The following 'tags' are allowed: 
;             p1, p2, v1, v2, r1, r2, sp1, sp2.
;
; bdate_begin 
; bdate_end   Specify a date-time range from which data should be extracted.
;             Both must be provided as String in the form yymmddHHMMSS. If leading
;             chars are omitted they are replaced to match the date-time of the
;             first dwell contained in the pds-file.
;
; max_dwells  limits the number of dwells which is loaded per file.
; max_mem     viewpds will stopp goging data if this much memory is used.
;
; After launching the pds files are gorged and a cmd line interface is launched.
; Different visualizations can be created by the following commands:
;
;-----------------------------------------------------------------------------------------------------------------
;     View PDS > sp [sp1|sp2]
;
; Plot the profile of spectra at dwell number itime as image plot. The profile of velocities from the moment data 
; is oplotted into the image plot. The dwellnumber itime can be selected by one of the following commands:
;
;     View PDS > + [increment]
;     View PDS > - [decrement]
;     View PDS > ++                # go to the begining of the next file if multible files where loaded
;
; Note multible commands seperated by " ; " (blancs are necessary) can be given in one line:
;
;     View PDS > + ; sp  return return return ...
;
; Like this profiles can be stepped. This works because blanc command lines cause repeating 
; of the preveous (non-blanc) command line.
;
; If multible files where loaded then ++ (--) seeks to the begining of the next 
; (preveous) file.
;
;-----------------------------------------------------------------------------------------------------------------
;     View PDS > thc [p1|p2|v1|v2|s1|s2|ldr]
;
; Plot a Time Height Cross section of the profiles idicated by the first argument.
; The default value of the z-range of each type of profiles is stored in the
; parameter <f>range where <f> is the first letter of the profiles name. For 
; example the default value of the z-range for the ldr plots cam be changed by
;
;    View PDS > set lrange -30 -1 ; thc ldr
;
;-----------------------------------------------------------------------------------------------------------------
;    View PDS > set <parameter> <valu1> [<value2> ...]
;
; Change the default value of a plotting parameter. If <value1> is omited the
; current value is plotted. All the plotting parameters which 
; are available can be shown by
;
;    View PDS > set help
;-----------------------------------------------------------------------------------------------------------------
;    View PDS > tr <begin_hour> <end_hour>
;
; limits the time range of the thc plots. <begin_hour> <end_hour> are given as decimal numbers:
; 17.5 means 17:30. The hours ar countet from 00:00 of the first day contained in the data
; which is loaded.
;-----------------------------------------------------------------------------------------------------------------
;    View PDS > postscript <filename>
;
; Directs the output of the following plot commands to the given postscript filename.
; for flushing and closing the postscript file the command
;
;    View PDS > close
;
; must be given. After that the plot out put is directed to the screen again.
;----------------------------------------------------------------------------------------------------------------
;    View PDS > png <filename>
;
; exports the screen to png-format
;-----------------------------------------------------------------------------------------------------------------
;    View PDS > info
;
; gives some information from the file header
;-----------------------------------------------------------------------------------------------------------------
;    View PDS > x
;
; exits ViewPDS
;----------------------------------------------------------------------------------------------------------------
; If no command is recognized then the command is executed by the IDL interpreter.
;
; Matthias Bauer-Pfundstein bauer@metek.de
;
common startup,startup_done
if n_elements(startup_done) ne 1 then begin
@startup
endif

; Arguments can be provided by two methods:
; - If viewpds is called by an interactive IDL session the first (and only) 
; string parameter arguments which contains a space separated list is split to
; the strarr argv.
; - Aditionally strings contained in the environment variables ARG1, ARG2, ..., 
; ARG9, ARG10, ... are prepended to argv. The script idl packs its calling arguments
; these environment variables. This is done as IDL can not extract the shell 
; variables like $1, $2,...
argv = env2argv(arguments)

kp = { $
    infile:'',$
    extract:'default', $
    bdate_begin:'default', $
    bdate_end:'default', $
    max_mem:500000000L, $
    progress_bar:1, $
    max_dwells:9000L $
}
argv2stru,argv,kp

infiles = strsplit(kp.infile,':',/extract)

if n_elements(infiles) eq 1 and infiles[0] eq '' then begin
    infiles = DIALOG_PICKFILE(/read,/must_exist,/multiple_files, path = path, filter='*.pds')
endif

extracts = strsplit(kp.extract,':',/extract)

; Default values for the Ploting parameters. These Parameters can be changed by the command "set"
pp={ $
hr:[240.,13000.],$          ; height range for plotting
itr:[0,100L], $             ; range of dwell numbers which are shown in thc-plots
badcol :15, $               ; grey, color for plotting values which are not plausible.
srange:[-15,75.], $         ; zrange for spectra
prange:[-25,70.], $         ; zrange for snr profiles
zrange:[-60,60.], $         ; zrange for dbz
vrange:[!values.f_nan,!values.f_nan], $  ; zrange for velocities profiles
rrange:[0.0,4.], $          ; zrange for    Peak with (RMS) profiles
lrange:[-35.,0.], $         ; zrange for    LDR profiles
min_txn_pow:1800., $        ; thc plots only the profiles where srvinfo.tpow (transmitted avg power * 100) grater than min_txn_pow
min_npw:40.,$               ; thc plots only the profiles where srvinfo.npw1 > min_npw
rec_err_good_bits:'0000'X,$ ; thc plots only the profiles where rec_err_good_bits in tsrvinfo.re_err are set
rec_err_bad_bits:'0000'X, $ ; thc plots only the profiles where rec_err_bad_bits in tsrvinfo.re_err are not set
xtu:'Minutes',$             ; parameter to control the number of ticks on the time axis. 'minutes'|'time'|'hours'..
tolerable_gap:6d0, $        ; a gap in the time sequence will be detected when difference between the time stamps of two dwells is greater than tolerable_gap*median dwelltime
takeat_gap:'right', $      ; describes how itr should be changed when a gap is detected
charsize:1.1, $  
thick:1.0, $
b4tit:'',$
colt:'smooth' $
}
argv2stru,argv,pp

startcolor=30
ocolt='none'
td=1000.

!p.multi=[0,1,1]
if !d.name eq "X" then begin
    if !d.window ge 0 then wdelete
    screen_size=get_screen_size()
    scale=screen_size/[1280.,1024]
    window,/free,retain=2,colors=maxcolor,title="viewpds",xpos=90*scale[0],ypos=250*scale[1], $
    xsize=1180*scale[0],ysize=800*scale[1]
endif

zfiles = n_elements(infiles)
last_zt = 0L
for ifile=0,zfiles-1 do begin
   infile=infiles[ifile]
   gorge_pds_file, infile, _hdr, hdr, parr_srvinfo, parr_data, extract=extracts, bdate_begin=kp.bdate_begin, $
        bdate_end=kp.bdate_end,max_dwells=kp.max_dwells,max_mem=kp.max_mem,progress_bar=kp.progress_bar
   print
   if ifile eq 0 then begin
       hdrs=replicate(hdr,zfiles)
       _hdrs=replicate(_hdr,zfiles)
       ztimes = lonarr(zfiles)
   endif
   ztimes[ifile] = n_elements(parr_srvinfo) * ptr_valid(parr_srvinfo[0])
   hdrs[ifile] = hdr
   _hdrs[ifile] = _hdr
   print,infile,': ',unixtime2bdate((*parr_srvinfo[last_zt]).TIME_T,nice=2,utc=1), $
               ' - ',unixtime2bdate((*parr_srvinfo[ztimes[ifile]-1]).TIME_T,nice=2,utc=1)
   last_zt = ztimes[ifile]
endfor
max_ave = max(hdrs[*].ave,min=min_ave)
if max_ave ne min_ave then print,'Files with different averaging time have been loaded. The time axis of thc plots will be wrong. Use ++, -- !'

ztime = n_elements(parr_srvinfo) * ptr_valid(parr_srvinfo[0])
if ztime eq 0 then message,'No Data gorged'
unix_times = lonarr(ztime) 
for itime=0, ztime-1 do unix_times[itime] = (*parr_srvinfo[itime]).TIME_T
hours = unix_times/3600d0
first_day_0000 = long(hours[0]/24d0)*24d0
hours = hours - first_day_0000
pp.itr = [0,ztime-1] ; range of dwell numbers which are shown in thc-plots
juldate = unixtime2jul(unix_times)

txn_power = fltarr(ztime)
for itime=0, ztime-1 do txn_power[itime] = (*parr_srvinfo[itime]).tpow

txn_err = intarr(ztime)
for itime=0, ztime-1 do txn_err[itime] = (*parr_srvinfo[itime]).te_err

rec_err = intarr(ztime)
for itime=0, ztime-1 do rec_err[itime] = (*parr_srvinfo[itime]).rc_err

npw1 = fltarr(ztime)
for itime=0, ztime-1 do npw1[itime] = (*parr_srvinfo[itime]).npw1

cpw1 = fltarr(ztime)
for itime=0, ztime-1 do cpw1[itime] = (*parr_srvinfo[itime]).cpw1

ztime=n_elements(parr_data)
itime=0
ipcs=0
last_draw_arg='thc'
while 1 do begin
    if ocolt ne pp.colt then begin
        loadcolor,startcolor,colorfile=pp.colt,zcolor=zcolor 
        ocolt=pp.colt
    endif
    
    arg = myprompt('View PDS >')
artarg:
    cmd = arg(0) &  narg=n_elements(arg)
    case 1 of
        cmd eq 'redraw' or cmd eq 'r': begin
            arg = last_draw_arg
            goto,artarg
        end
        cmd eq 'sp': begin
            last_draw_arg = arg
            if narg ge 2 then tag_sp=arg(1) else tag_sp='sp1'
            ; plot a profile of spectra and oplot the profile of velocities from the moment data 
            hdr = hdrs[first(where(itime lt ztimes))] ; pick the header of the file from which itime was loaded
            dsp2vel = hdr.la*hdr.prf
            nffth = hdr.nfft/2 
            ; toord = [lindgen(nffth)+nffth,lindgen(nffth)] ; DWD
            toord = hdr.nfft-1 - [lindgen(nffth)+nffth,lindgen(nffth)] ; idx-vector for centering the spectra to zero velocity
            raw_height = findgen(hdr.raw_zrg) * hdr.delta_h + hdr.raw_h_min
            raw_ihr=[first(where(raw_height ge pp.hr[0])),last(where(raw_height le pp.hr[1]))]
            if min(raw_ihr) lt 0 or raw_ihr[0] ge raw_ihr[1] then begin
                printf,-2,'height range pp.hr',pp.hr,' should be between',first(raw_height),last(raw_height)
                break
            endif
            spcs_plt = un_gorge(hdr,parr_data,itime=itime,irg=raw_ihr,tag=tag_sp)
            if not finite(spcs_plt[0,0]) then break
            mimage, db(spcs_plt[toord,*]), [-0.5,.5]*hdr.vuar, raw_height[raw_ihr]/td , pp.srange, $
                zcolor, startcolor,pp.badcol,charsize=pp.charsize,thick=pp.thick, $
                title='Profile of Spectra ' + unixtime2bdate(unix_times[itime], nice=2,utc=1), $
                ztitle='Intensity',xtitle='Velocity [m/s]',$
                heitit='Height AGL [km]'
            tag_velo = 'v'+strmid(tag_sp,2,1)
            vr_plt = un_gorge(hdr,parr_data,itime=itime,tag=tag_velo)
            if not finite(vr_plt[0]) then break
            mom_height = [findgen(hdr.mom_zrg) * hdr.delta_h + hdr.mom_h_min,15000.,15000.]
            ;oplot, vr_plt , mom_height / td + 999.*(abs(vr_plt) gt 100.),color=0, max=900.,$
            ;    thick = pp.thick
            oplot, vr_plt , mom_height / td + 999.*(abs(vr_plt) gt 100.),color=0, max=900.,$
                thick = pp.thick, psym=4
        end
        cmd eq 'thc': begin
            last_draw_arg = arg
            limit_itr2contiguous,hours,max_ave,pp
            ; time height cross section
            if narg ge 2 then tag_mom=arg(1) else tag_mom='p1'
            hdr = hdrs[first(where(pp.itr[0] lt ztimes))] ; pick the header of the file from which itime was loaded
            mom_height = findgen(hdr.mom_zrg) * hdr.delta_h + hdr.mom_h_min
            mom_ihr=[first(where(mom_height ge pp.hr[0])),last(where(mom_height le pp.hr[1]))]
            case strupcase(strmid(tag_mom,0,1)) of 
              'P': begin
                ztitle='SNR [dB]'
                p_tmp = un_gorge(hdr,parr_data,itime=pp.itr,irg=mom_ihr,tag=tag_mom)
                if finite(p_tmp[0,0]) and n_elements(p_tmp) gt 1 then prof=db(p_tmp) else prof = p_tmp
                zrange = pp.prange
              end
              'Z': begin
                ztitle='dBZ [dB]'
                tag_mom = 'P'+strmid(tag_mom,1)
                p_tmp = un_gorge(hdr,parr_data,itime=pp.itr,irg=mom_ihr,tag=tag_mom)
                if finite(p_tmp[0,0]) and n_elements(p_tmp) gt 1 then begin
                    prof = db(p_tmp) 
                    good = prof ge pp.prange[0]
                    for iti=pp.itr[0],pp.itr[1] do begin
                        iti_prof = iti - pp.itr[0]
                        prof[0,iti_prof] = dbz(prof[*,iti_prof], hdr, txn_power[iti]) $
                            - 1000. * (good[*,iti_prof] eq 0)
                    endfor
                endif else prof = p_tmp
                zrange = pp.zrange
              end
              'V': begin
                ztitle='Velocity [m/s]'
                prof =  un_gorge(hdr,parr_data,itime=pp.itr,irg=mom_ihr,tag=tag_mom) 
                if finite(pp.vrange[0]) then zrange=pp.vrange else zrange=[-0.5,.5]*hdr.vuar
              end
              'R': begin
                ztitle='Peak Width [m/s]'
                prof =  un_gorge(hdr,parr_data,itime=pp.itr,irg=mom_ihr,tag=tag_mom) 
                zrange = pp.rrange 
              end
              'L': begin
                if strupcase(tag_mom) eq 'LDR' then begin
                  ztitle='LDR [dB]'
                  p1 = db(un_gorge(hdr,parr_data,itime=pp.itr,irg=mom_ihr,tag='p1'))
                  p2 = db(un_gorge(hdr,parr_data,itime=pp.itr,irg=mom_ihr,tag='p2'))
                  good = p1 ge pp.prange[0] and p2 ge pp.prange[0]
                  if n_elements(p1) gt 1 and n_elements(p2) gt 1 then begin
                      prof = (p2 - p1) * good + 100.* (good eq 0)
                  endif else prof = !values.f_nan
                  zrange = pp.lrange
                endif else prof = !values.f_nan
              end
              else: begin
                  prof = !values.f_nan
                  zrange=[0,1]
              end
            endcase

            if finite(prof[0,0]) then begin
                good_times = txn_power[pp.itr[0]:pp.itr[1]] ge pp.min_txn_pow $
                    and ( rec_err[pp.itr[0]:pp.itr[1]] and pp.rec_err_bad_bits) eq 0 $
                    and ( rec_err[pp.itr[0]:pp.itr[1]] and pp.rec_err_good_bits) eq pp.rec_err_good_bits $
                    and npw1[pp.itr[0]:pp.itr[1]] ge pp.min_npw
                zrg = n_elements(prof[*,0])                          ; should be mom_ihr[1]-mom_ihr[0]+1
                prof = transpose(   prof +  replicate(1.,zrg) # ( (good_times ne 1) * 999999.)    ) 
                sxtf=!x.TICKFORMAT & !x.TICKFORMAT='LABEL_DATE'
                sxtu=!X.TICKUNITS[0] & !X.TICKUNITS[0]=pp.xtu
                dummy = LABEL_DATE(DATE_FORMAT='%h:%i')
                mimage, prof, juldate[pp.itr], mom_height[mom_ihr]/td , zrange, $
                    zcolor, startcolor,pp.badcol,charsize=pp.charsize,thick=pp.thick, $
                    title = pp.b4tit + unixtime2bdate(unix_times[pp.itr[0]], nice=2,utc=1) $ 
                        +  ' - '     + unixtime2bdate(unix_times[pp.itr[1]], nice=2,utc=1), $
                    ztitle=ztitle,xtitle='Time UT',$
                    heitit='Height AGL [km]'
                !x.TICKFORMAT = sxtf & !X.TICKUNITS[0]=sxtu
            endif else printf,-2,'No data for plotting ' + tag_mom
        end
        cmd eq '+' or cmd eq '-': begin 
            if narg eq 2 then stride = long(arg[1]) else stride = 1L 
            if cmd eq '-' then stride *= -1
            itime = ((itime + stride) < (ztime-1) ) > 0
        end
        cmd eq '++' or cmd eq '--': begin
           if cmd eq '++' then stride=1L else stride=-1L
           if  ipcs eq 0 then pp.itr = [0,ztimes[ipcs]-1] $
                         else pp.itr = [ztimes[ipcs-1],ztimes[ipcs]-1]
           ipcs=(ipcs+stride+zfiles) mod zfiles
           if cmd eq '++' then itime = pp.itr[0] else itime = pp.itr[1]
        end
        cmd eq 'tr': begin
            case narg of
              1: pp.itr = [0,ztime-1]
              2: begin
                if strupcase(arg[1]) ne 'NAN' then begin
                    if double(arg[1]) ge 0d0 then begin
                        dumy=min(abs(hours-double(arg[1])),itr0)
                        pp.itr[0] = itr0
                    endif else begin
                        dumy=min(abs(hours-hours[ztime-1]-double(arg[1])),itr0)
                        pp.itr = [itr0,ztime-1]
                    endelse
                endif
              end
              3: begin
                if strupcase(arg[1]) ne 'NAN' then begin
                    dumy=min(abs(hours-double(arg[1])),itr0)
                    pp.itr[0] = itr0
                endif
                if strupcase(arg[2]) ne 'NAN' then begin
                    dumy=min(abs(hours-double(arg[2])),itr1)
                    pp.itr[1] = itr1
                endif
                if pp.itr[0] ge  pp.itr[1] then begin
                    printf,-2,'Timerange crossed over. Resetting time range'
                    pp.itr = [0,ztime-1]
                endif
              end
              else: printf,-2,'wrong number of arguments'
            endcase
        end
        cmd eq 'postscript' or cmd eq 'ps': begin
            if narg eq 2 then psfile = arg[1] else psfile="viewpds.ps"
            set_plot,'PS'
            device,/color,bits_per_pixel=8,filename=psfile ; ,yoffset= 27,xoffset=0.7
        end
        cmd eq 'png': begin
            if narg ge 2 then png_fn = arg[1] else begin
                png_fn = infile + '_' + tag_mom + '.png'
            endelse
            screen_dump = TVRD(true=1)
            WRITE_PNG, png_fn, screen_dump
        end
        cmd eq 'close':begin
            device,/close
            set_plot,'X'
        end
        cmd eq 'set': begin
            if narg ge 2 then set_tag_name = arg(1) else set_tag_name='HELP'
            if strupcase(set_tag_name) eq 'HELP' then begin
                pp_tags = tag_names(pp)
                for i=0,n_elements(pp_tags)-1 do print,pp_tags[i]+' :',pp.(i)
            endif else begin
                set_tag_number = (where(strupcase(set_tag_name) eq tag_names(pp)))[0]
                if set_tag_number ge 0 then begin
                    if narg eq 2 then print,set_tag_name+' :',pp.(set_tag_number) else begin
                        for i=0,(narg-3) < (n_elements(pp.(set_tag_number))-1) do begin
                            if strupcase(arg[i+2]) ne 'SKIP' then pp.(set_tag_number)[i] = arg[i+2]
                        endfor
                    endelse
                endif else printf,-2,'unknown variable name'
            endelse
        end
        cmd eq 'info': begin
            help,hdrs[first(where(itime lt ztimes))],/struct
            if narg eq 2 then help,_hdrs[first(where(itime lt ztimes))],/struct
        end
        cmd eq 'stop': if lmgr(/vm) then goto,ende else stop
        cmd eq 'x': goto,ende
        else: begin
            setcmd = ''
            for i=0 , narg - 1 do setcmd = setcmd + arg(i) + ' '
            if lmgr(/vm) then printf,-2,'VM license can not execute '+setcmd $
            else dummi=execute( setcmd )
        end
    endcase
endwhile

ende:
if !d.name eq "X" then wdelete
; heap_free,parr_data
; heap_free,parr_srvinfo
end
