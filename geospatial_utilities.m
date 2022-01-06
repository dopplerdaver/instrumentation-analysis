
%------------------------------------------------------------------------------%
%
% Author: Steve Mueller (August-December 2006)
%
% Purpose:
%  "Master switch" for geospatial utilties.
%
% Input:
%  subfunction_name - String representing specific processing function.
%  Currently, valid values are:
%     'get_distance_on_sphere'
%     'get_local_earth_radius'
%     'get_solid_angle'
%
% Output:
%  Returns the processed value(s) requested.
%
% Usage:
%  geospatial_utilities('subfunction_name', arguments);
%   Examples ---
%      distance = geospatial_utilities('get_distance_on_sphere', {59.0 -134.0 59.0 -133.0});
%      dist_km  = geospatial_utilities('get_distance_on_sphere', {59.0 -134.0 59.0 -133.0})*(1/1000);
%      earth_radius = geospatial_utilities('get_local_earth_radius', 58.0);
%      solid_angle = geospatial_utilities('get_solid_angle', {(pi/180)*58.0 (pi/180)*(-134.0) (pi/180)*59.0 (pi/180)*(-133.0)});
%
% Notes:
%  The somewhat awkward manner in which these geospatial utilities
%   must be accessed is a consequence of the fact that Matlab does not
%   provide a convenient means of bundling multiple functions into a
%   single file. The only other alternative would require a separate
%   file for each function.
%
%
%------------------------------------------------------------------------------%
function varargout=geospatial_utilities(subfunction_name, varargin)
   varargout=cell(1, max(1, nargout));
   [varargout{:}]=feval(subfunction_name, varargin);
   return
% End of read_data()


%----------------------------------------------------------------------%
%
% Name: get_distance_on_sphere()
%
% Author: Steve Mueller (August 2005)
%
% Purpose:
%  Determines distance along the surface of a sphere based on lat-lon pairs.
%
% Input:
%  Cell array of lat-lon pairs (e.g., {lat1_degrees lon1_degrees lat2_degrees lon2_degrees}).
%     lat1_degrees - Float for 1st latitude  in degrees.
%     lon1_degrees - Float for 1st longitude in degrees.
%     lat2_degrees - Float for 2nd latitude  in degrees.
%     lon2_degrees - Float for 2nd longitude in degrees.
%
% Output:
%  Distance in meters.
%
% Comments:
%    Refer to Geodynamics by Turcotte and Schubert (pg. 206) for details.
%
% Usage example:
%   distance = geospatial_utilities('get_distance_on_sphere', {59.0 -134.0 59.0 -133.0});
%   dist_km = geospatial_utilities('get_distance_on_sphere', {59.0 -134.0 59.0 -133.0})*(1/1000);
%
%----------------------------------------------------------------------%
function distance=get_distance_on_sphere(cell_array_input)
   FLATTENING = 'FALSE';
   EQUATORIAL_RADIUS = 6.378139E06;

   % Determine solid angle between lat-lon pairs.
   lat1_radians = cell_array_input{1}{1}*(pi/180.0);
   lon1_radians = cell_array_input{1}{2}*(pi/180.0);
   lat2_radians = cell_array_input{1}{3}*(pi/180.0);
   lon2_radians = cell_array_input{1}{4}*(pi/180.0);
   solid_angle = get_solid_angle(lat1_radians, lon1_radians, lat2_radians, lon2_radians);

   if(strcmp(FLATTENING, 'TRUE'))
      lat_avg_degrees = (cell_array_input{1}{1} + cell_array_input{1}{3})/2.0;
      earth_radius = get_local_earth_radius(lat_avg_degrees);
   else
      earth_radius = EQUATORIAL_RADIUS;
   end

   distance = solid_angle*earth_radius;
% End of get_distance_on_sphere() method.


%----------------------------------------------------------------------%
%
% Name: get_local_earth_radius()
%
% Author: Steve Mueller (August 2005)
%
% Purpose:
%    Calculates Earth radius (m) as a function of latitude by compensating
%    for rotational oblation.
%
% Input:
%    lat_deg - Float latitude in degrees.
%
% Output:
%    Local Earth radius in meters.
%
% Comments:
%    Refer to Geodynamics by Turcotte and Schubert (pg. 206) for details.
%
% Usage example:
%  earth_radius = geospatial_utilities('get_local_earth_radius', 58.0);
%
%----------------------------------------------------------------------%
function local_earth_radius =get_local_earth_radius(lat_degrees)
   GEOFLATTENING = (1.0/298.256);
   EQUATORIAL_RADIUS = 6.378139E06;

   lat_radians = lat_degrees{1}*(pi/180.0);

   oblate_factor = 1.0/sqrt(1 + ((2*GEOFLATTENING + GEOFLATTENING^2)/(1 - GEOFLATTENING)^2)*sin(lat_radians)^2);
   local_earth_radius = EQUATORIAL_RADIUS*oblate_factor;
% End of get_local_earth_radius() method.


%----------------------------------------------------------------------%
%
% Name: get_solid_angle()
%
% Author: Steve Mueller (August 2005)
%
% Purpose:
%    Calculates solid angle between lat-lon pairs.
%
% Input:
%    Cell array with lat-lon pairs ({lat1_radians lon1_radians lat2_radians lon2_radians}).
%       lat1_radians - Float first  latitude  in radians.
%       lon1_radians - Float first  longitude in radians.
%       lat2_radians - Float second latitude  in radians.
%       lon2_radians - Float second longitude in radians.
%
% Output:
%    Float angle in radians.
%
% Comments:
%    Input values must be RADIANS and the output is radians.
%
% Usage example:
%  solid_angle = geospatial_utilities('get_solid_angle', {(pi/180)*58.0 (pi/180)*(-134.0) (pi/180)*59.0 (pi/180)*(-133.0)});
%
%----------------------------------------------------------------------%
function solid_angle = get_solid_angle(input_cell_array);
   lat1_radians = input_cell_array{1}{1};
   lon1_radians = input_cell_array{1}{2};
   lat2_radians = input_cell_array{1}{3};
   lon2_radians = input_cell_array{1}{4};
   solid_angle = acos(((sin(lat1_radians))*(sin(lat2_radians)) + (cos(lat1_radians)*cos(lat2_radians))*(cos(lon1_radians-lon2_radians))));
% End of get_solid_angle() method
