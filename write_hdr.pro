pro write_hdr,lun, ppar, nmulti
; write the 1024 byte pds file header
    ppar.oper[57:63] = byte('MULTI'+string(nmulti,format='(i2.2)')) 
    writeu,lun,ppar
    writeu,lun,replicate(0B,372)
end
