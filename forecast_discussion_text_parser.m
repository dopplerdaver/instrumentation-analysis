%-------------------------------------------------------------------
% Name:    forecast_discussion_text_parser.m
%
% Purpose:
%
% Inputs:  1. 
%          2. 
%
% Usage:   matlab -nodesktop
%          > forecast_discussion_text_parser 
%
% Created: 12.11.2013 dserke
%
% Coding tasks list:
% Task Status
% 1.   0      
% 2.   x     
%
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% hardcoded variables
%---------------------------------------------------------------------
% input station
station = 'CLE';

%---------------------------------------------------------------------
% initialize fields to be passed on to NIRSS volumetric testbed 
%---------------------------------------------------------------------
% assume can interpolate icing hazard to whole radar domain unless
%  determined otherwise by this program
interp_haz_to_radardomain_now = 1;
interp_haz_to_radardomain_fut = 1;

%---------------------------------------------------------------------
% read data from url into text string variable
%---------------------------------------------------------------------
text = urlread(['http://forecast.weather.gov/product.php?site=NWS&issuedby=' station '&product=AFD&format=TXT&version=1&glossary=1']);

%---------------------------------------------------------------------
% define start and end indices of desired info blocks in text discussion 
%---------------------------------------------------------------------
ind_start_synopsis = findstr('.SYNOPSIS...',text);
ind_end_synopsis   = findstr('&&',text);
ind_start_nearterm = findstr('.NEAR TERM /',text);
ind_end_nearterm   = findstr('&&',text);
ind_start_aviation = findstr('.AVIATION /',text);
ind_end_aviation   = findstr('&&',text);

%---------------------------------------------------------------------
% define desired info blocks in text discussion 
%---------------------------------------------------------------------
str_synopsis_block = text(ind_start_synopsis+12:ind_end_synopsis(1)-1);
str_nearterm_block = text(ind_start_nearterm+11:ind_end_synopsis(2)-1);
str_aviation_block = text(ind_start_aviation+10:ind_end_aviation(end-2)-1);
%time_valid_until   = str_nearterm_block(start:end);

%---------------------------------------------------------------------
% determine the yyyy mm dd hh mm of the discusion in variable 'text'  
%---------------------------------------------------------------------

%---------------------------------------------------------------------
% search info blocks for keywords and possible associated info  
%---------------------------------------------------------------------
time_date_current    = text(ind_start_synopsis-29:ind_start_synopsis-1);  
ind_time_valid_until = find(str_nearterm_block == '/');
time_valid_until     = str_nearterm_block(ind_time_valid_until(1):ind_time_valid_until(2));

keyword = {'FRONT'};

for j=1:length(keyword)

  ind_front = findstr(keyword{j},str_synopsis_block);
  count     = 1;

  for i=1:length(ind_front)
    if str_synopsis_block(ind_front(i)+5:ind_front(i)+8) == '</a>'
      % ...then read rest of the sentence  
      ind_end_sent                  = findstr('.',str_synopsis_block);
      ind_end_keyword_sent          = find(ind_end_sent > ind_front(i));
      ind_end_keyword_sent          = ind_end_sent(min(ind_end_keyword_sent));
      str_sent_after_keyword{count} = str_synopsis_block(ind_front(i)+8:ind_end_keyword_sent);
      count                         = count + 1;
    end 
  end

end      % end of for j=1:length(keyword)

% try to discern where the front is and when it will affect the area of interest





