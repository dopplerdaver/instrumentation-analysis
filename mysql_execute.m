function data = mysql_execute(query,varargin)
% mysql_execute: Get result from mysql query (or other command)
%  TO work a mysql client must be available.  Note that this is probably not a 
%  good option if the query returns large results, especially if you have non-numbers
%  (dates,strings) coming out of mysql.  If you have all numbers being returned
%  you should use the assume_all_numbers optimization (much faster).  
%
%  WARNING: This uses the mysql client to get the data so it must be installed.  Also,
%  because it is using the client, the output is just a string.  It is possible to confuse 
%  it if there are tabs in the data (as it will be interpreted as a field seperator).
%
%  WARNING: if the output type is a float (probably only get this if the table column
%  is a float and the query does not manipulate it, the precision is only 6 digits - 
%  when mysql sends the data (less than mysql saves at).  To get around this you can
%  use the matlab round function, or add 0.0 (which will convert to double).  
%
%  WARNING: CTRL-C may not work when trying to stop execution of this function
%  if matlab is waiting for the mysql client to finish.  You may have to kill
%  the mysql client (via unix kill command), to get your session back.  This is
%  on the TODO list to resolve.  
%
%  usage data = mysql_execute(query,'param1',val1,...)
%    query is a string containing the mysql statement (any valid mysql statement will
%    work).
%
%    data is either a string (for output_type 1) or a struct.  If no records are
%    returned then data will always be an empty string.
%
%    where params can be: (default is shown
%
%    args = ''; Arguments to pass to the mysql client.  e.g. '--host=delphi9 -u insitu_data 
%               --password=XXXXX gturb'  
%               Note that any password given should not be visible from 'ps' and such.  
%               It will just be in your script or session.  To avoid any security issues
%               you should set up a ~/.my.cnf file with the username/password or if you
%               need multiple connections, some other file that you then tell the mysql
%               client to use via the --defaults-extra-file=# option.  Either config file
%               would have something like:
%                 [client]
%                 user=insitu_data
%                 password=XXXXX
%               The format assumed for any translations (see output_type) other than 
%               output_type = 0, are tab delimted with column header line.  So if the 
%               output_type is something other than 0, you should not use options like
%               -t, --table, -v, -N, or --skip-column-names -H, --html or -X, --xml.
%    quote = ''''; When calling the mysql client, a unix shell is used and the query string 
%                is enclosed in the quote specified by this variable.  So if you need
%                quotes is the query string, you need to use the opposite one.  So given
%                the default, you should always use " in your query.  You probably don't
%                want to mess with this.
%    output_type = 2; 0: raw string.  No translation into matlab is done so any mysql
%                        options are probably fine (see args)
%                     1: matlab struct where the fieldnames are derived from the column 
%                        names.  The mysql options that will not crash this function are
%                        limited (See args).  Note that the column names MUST be valid
%                        matlab fieldnames from structures.  Use the AS command in your
%                        mysql select statement to control column names.
%                     2: matlab struct with 3 fields: header (cell array of strings),
%                        data (cell array of data), and ref (an inline function to assist in 
%                        pulling the data out. e.g. data.ref(data,'measurement_time') 
%                        will return the data from the column corresponding to 
%                        measurement_time.) The mysql options that will not crash this 
%                        function are limited (See args)
%    assume_all_numbers = 0; This is only used if output_type is not 0.  If this is set
%                 to 1, then the conversions will occur assuming that all data are numbers.
%                 This can speed up the conversion drastically, depending on the size of the
%                 query results.  Note that NULL's are OK.  They are automatically converted
%                 to NaN's.
%    autoconvert_dates = 0;  This is only used if output_type is not 0 and assume_all_numbers
%                 is 0.  This function will attempt to convert the dates into matlab 
%                 datenum's (type 'help datenum', also see autoconvert_date_formats)
%    autoconvert_date_formats = {'yyyy-mm-dd HH:MM:SS','yyyy-mm-dd','HH:MM:SS'}; Only used if
%                 autoconvert_dates occurs (see autoconvert_dates).  This is a cell array of 
%                 strings of date formats (see datestr) to convert the mysql date into 
%                 a matlab datenum.
%    use_temp_file = 0; You can instead force this function to use a temporary file for the
%                 mysql output instead of capturing the STDOUT from the unix call.  This can
%                 be useful if calling unix gives you extra output
%    temp_dir =  '/tmp'; Specify the directory to put the temp file in.  Make sure that your
%                your temp partition does not run out of space!  Only used if  use_temp_file=1
%    echo_query = 0; turn to 1 to print the final query to screen
%    echo_unix = 0; turn to 1 to print the final unix command to screen
%
%
%   examples:   d1 = mysql_execute('select latitude,longitude,altitude,unix_timestamp(measurement_time) as MT,peak_edr from insitu_orig limit 3','args','--host=delphi9 --table -u insitu_data --password=XXXXX  insitu','output_type',0)
%               d2 = mysql_execute('select latitude,longitude,altitude,unix_timestamp(measurement_time) as MT,peak_edr from ude_usa limit 3','args','--host=delphi9 -u insitu_data --password=XXXXX  gturb','output_type',1,'assume_all_numbers',1)
%               d3 = mysql_execute('select latitude,longitude,altitude,measurement_time,peak_edr from ude_usa limit 3','args','--host=delphi9 -u insitu_data --password=INSITU gturb','output_type',1,'autoconvert_dates',1)

sql = sql_core(query,'query_call',@(query,args) sprintf('echo %s | mysql %s | sed -re ''s/(\t|^)NULL($|\t)/\\1NaN\\2/g;s/(\t|^)NULL($|\t)/\\1NaN\\2/g;s/\t/,/g''',query,args),varargin{:});

data = sql.run();
return
end










% $$$ args = '';
% $$$ quote = '''';
% $$$ output_type = 2;
% $$$ assume_all_numbers = 0;
% $$$ guess_format_lines = 0;
% $$$ autoconvert_dates = 0;
% $$$ autoconvert_date_formats = {'yyyy-mm-dd HH:MM:SS','yyyy-mm-dd','HH:MM:SS'};
% $$$ provide_format = '';
% $$$ use_temp_file = 0;
% $$$ temp_dir = '/tmp';
% $$$ 
% $$$ echo_query = 0;
% $$$ echo_unix = 0;
% $$$ 
% $$$ paramparse(varargin);
% $$$ 
% $$$ if guess_format_lines
% $$$   assume_all_numbers = 0;
% $$$   limit = {guess_format_lines};
% $$$ else
% $$$   limit = {};
% $$$ end;
% $$$ 
% $$$ if ~isempty(provide_format)
% $$$   assume_all_numbers = all(provide_format(2:2:end)=='f');
% $$$ end
% $$$ 
% $$$ ref = inline('this.data{find(strcmp(this.header,fld))}','this','fld');
% $$$ 
% $$$ str = sprintf('echo %s%s%s | mysql %s | sed -re ''s/(\t|^)NULL($|\t)/\\1NaN\\2/g;s/(\t|^)NULL($|\t)/\\1NaN\\2/g;s/\t/,/g'' ',quote,query,quote,args);
% $$$ 
% $$$ if use_temp_file
% $$$   tmpfile = tempname(tempdir);
% $$$   str = [str '> ' tmpfile];
% $$$ end
% $$$ 
% $$$ if echo_query
% $$$   fprintf('%s\n',query);
% $$$ end
% $$$ 
% $$$ if echo_unix
% $$$   fprintf('%s\n',str);
% $$$ end
% $$$ 
% $$$ try
% $$$   [r,o] = unix(str);
% $$$ catch
% $$$   r = 1;
% $$$   o = lasterr;
% $$$ end
% $$$ 
% $$$ if r
% $$$   if use_temp_file
% $$$     try
% $$$       delete(tmpfile);
% $$$     end
% $$$   end
% $$$   disp(r)
% $$$   error(o);
% $$$ end
% $$$ 
% $$$ if use_temp_file
% $$$   fid = fopen(tmpfile,'rt');
% $$$   o = fscanf(fid,'%c');
% $$$   %delete(tmpfile);
% $$$ end
% $$$   
% $$$   
% $$$ if output_type == 0
% $$$   data = o;
% $$$   return
% $$$ end
% $$$ 
% $$$ if isempty(o)
% $$$   data = '';
% $$$   return
% $$$ end
% $$$ 
% $$$ ind = find(sprintf('\n')==o(1:min(10000,length(o))));
% $$$ if length(ind) == 0 
% $$$   o = sprintf('%s\n',o);
% $$$   ind = length(o);
% $$$ end;
% $$$ 
% $$$ ind = ind(1);
% $$$ header = char2delim(o(1:ind-1),',','IgnoreRepeated','on','Reverse','tocell').';
% $$$ 
% $$$ o = o(ind+1:end);
% $$$ params = varstruct({'assume_all_numbers','provide_format','limit','autoconvert_dates','autoconvert_date_formats'});
% $$$ 
% $$$ [d,string] = parse_data(o,header,params);
% $$$  
% $$$ if guess_format_lines
% $$$   params.limit = {};
% $$$   params.provide_format = string;
% $$$   d = parse_data(o,header,params);
% $$$ end
% $$$ 
% $$$ if output_type==1
% $$$   % look_for dups
% $$$   try
% $$$     header = handle_dups(header);
% $$$     data = cell2struct(d,header,2);
% $$$     return
% $$$   catch
% $$$     warning(sprintf('Could not convert to struct.  Falling back to output_type 2.  REASON: %s',lasterr));
% $$$     output_type = 2;
% $$$   end;
% $$$ end 
% $$$ 
% $$$ if output_type==2
% $$$   data.header = header;
% $$$   data.data = d;
% $$$   data.ref = ref;
% $$$   return
% $$$ end  
% $$$ 
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% $$$ function [d,string] = parse_data(o,header,p)
% $$$ if p.assume_all_numbers
% $$$   string = repmat('%f',1,length(header));
% $$$ elseif ~isempty(p.provide_format)
% $$$   string = p.provide_format;
% $$$ else
% $$$   string = repmat('%q',1,length(header));
% $$$ end
% $$$ d = textscan(o,string,p.limit{:},'delimiter',',','emptyvalue',NaN,'ReturnOnError',0);
% $$$ 
% $$$ if ~p.assume_all_numbers & isempty(p.provide_format)
% $$$   for l = 1:length(d)
% $$$     try
% $$$       if ~isnumeric(d{l})
% $$$         res = str2double(d{l});
% $$$         if any(isnan(res))
% $$$           % if any NaN's show up, check to see if they exactly correspond to the 'NaN' string
% $$$           if all(isnan(res)==strcmp('NaN',d{l}))
% $$$             % if they do, then replace the d{l} with the doubles
% $$$             % otherwise leave the d{l} alone, since there must be strings
% $$$             d{l} = res;
% $$$             string(2*l) = 'f';
% $$$           elseif p.autoconvert_dates
% $$$              for k = 1:length(p.autoconvert_date_formats)
% $$$                try
% $$$                  d{l} = datenum(d{l},p.autoconvert_date_formats{k});
% $$$                  break;
% $$$                end;
% $$$              end;
% $$$           end
% $$$         else
% $$$           % if there were no NaN's then replace with double version
% $$$           d{l} = res;
% $$$           string(2*l) = 'f';
% $$$         end
% $$$       end
% $$$     end;
% $$$   end;
% $$$ end
% $$$ 
% $$$ %%%%%%%%%%%%%%%%%5
% $$$ %%% function handle_dups
% $$$ function h = handle_dups(h)
% $$$ 
% $$$ % find the unique field names 'u' and find the mapping from h to 'u'
% $$$ [u,~,h_inds] = unique(h);
% $$$ % now, count how many each u has
% $$$ ct = hist(h_inds,1:length(u));
% $$$ 
% $$$ % find the problems
% $$$ nonunique = find(ct>1);
% $$$ % loop over the problem fields
% $$$ for ll = 1:length(nonunique)
% $$$   % find the inds in h that correspond to nonunique(ll)
% $$$   inds = find(strcmp(u{nonunique(ll)},h));
% $$$   for kk = 2:length(inds)
% $$$     new_fld = sprintf('%s_%02i',h{inds(kk)},kk);
% $$$     if ~strcmp(new_fld,h)
% $$$       warning(sprintf('Renaming "%s"[%i] to "%s"',h{inds(kk)},kk,new_fld));
% $$$       h{inds(kk)} = new_fld;
% $$$     else
% $$$       error('Duplicate field.  Attempt to rename "%s" failed because there is already "%s".',h{inds(kk)},new_fld);
% $$$     end
% $$$   end
% $$$ end
% $$$ return