% MEDIAN_2D_FILTER Median filter 2D 
%
% median_2D_filter.c -- Routine to compute the median filter of a matrix.
%
% Usage: [B] = median_2D_filter(A,imed,m_window,n_window)
%
% Inputs:
%       A               2-D matrix of data to be median filtered
%       m_window        # of rows in the filter window, an odd number > 0
%       n_window        # of columns in the filter window, an odd number > 0
%       imed            method used to calculate the median 
%                       1=quicksort, 2=heapsort, 3=select_kth_smallest
%                       These are all Numerical Recipes methods with   
%                       modifications for double array, no nrutil.h 
%                       dependence, and adjusted the indexing to begin 
%                       at 0 rather than 1.
%
% Outputs:
%        B               2-D matrix of filtered data
%
%
%
%  Example:
%          If A = [  0.1934    0.1509    0.8537    0.8216
%	             0.6822    0.6979    0.5936    0.6449
%		     0.3028    0.3784    0.4966    0.8180
%		     0.5417    0.8600    0.8998    0.6602 ]
%          
%          then median_2D_filter(A,1,3,3) is 
%	  
%	             0.4378    0.6379    0.6714    0.7333
%		     0.3406    0.4966    0.6449    0.7314
%		     0.6119    0.5936    0.6602    0.6526
%		     0.4600    0.5191    0.7391    0.7391
%  
%	  
% median_2D_filter is a MEX-file for MATLAB.
% Type "mex median_2D_filter.c" in Matlab to create the 
% executable file "median_2D_filter.mexglx".

% Version 1 of 1-May-2002 Beth Chorbajian @ NCAR.
%
% help(median_2D_filter)
%
%/**************************************************************************/
%/**************************************************************************/
%
% Copyright (c) 1999-2000 UCAR
% University Corporation for Atmospheric Research(UCAR)
% National Center for Atmospheric Research(NCAR)
% Research Applications Program(RAP)
% P.O.Box 3000, Boulder, Colorado, 80307-3000, USA
% All rights reserved. Licenced use only.
% Date: 2002/05/01
%



























