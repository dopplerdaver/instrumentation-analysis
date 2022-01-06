function whoami,iam
    common magic4chs,magic4chs,FzkPdsFmt
    if n_elements(iam) ne 1 || iam eq 'nobody' then iam='JUEL'
case iam of
    'DWD': begin
        FzkPdsFmt=0b
        magic4chs='FHDR'
    end
    'NASA':begin
        FzkPdsFmt=0b
        magic4chs='FHDR'
    end
    'FZK':begin
        FzkPdsFmt=1b
         magic4chs='FZKF'
    end
    'MPI':begin ; MPI2 is meant as completely other software is used for MPI1
        FzkPdsFmt=1b
        magic4chs='FZKF'
    end
    'HALO':begin
        FzkPdsFmt=1b
         magic4chs='HALO'
    end
    'NUIG':begin
        FzkPdsFmt=1b
        magic4chs='NUIG'
    end
    'CNRP':begin
        FzkPdsFmt=1b
        magic4chs='CNRP'
    end
    'JUEL':begin
        FzkPdsFmt=1b
        magic4chs='0709'
    end
    'IFT':begin
        FzkPdsFmt=1b
        magic4chs='NMRA'
    end
    else:begin
       print,'I dont know who Iam: '+kp.iam 
       exit
    end
endcase
    return,iam
end
