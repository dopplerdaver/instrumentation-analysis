OWLeS_slwsonde_launches.m

This code gives parameters for cropping the csv files from the OWLeS field campaign launches

1.10.2014 - dserke

% line numbers of data file header
file_hdr = 1:17;  

cell_hdr = flight,    date, location,           start_line, stop_line;
cell_l1  = ow001, 1/9/2014, SUNY Oswego Campus, 1225,       1487;
cell_l2  = ow002, 1/9/2014, SUNY Oswego Campus, 408,        747;
cell_l3  = ow003, 1/9/2014, SUNY Oswego Campus, 81,         end;
cell_l4  = ow004, 1/9/2014, SUNY Oswego Campus, 189,        end;

% setup code to only accept and plot freq vs. hgt if freq is >39 & < 45 
