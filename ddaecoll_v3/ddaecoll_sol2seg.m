function prob = ddaecoll_sol2seg(prob, oid, varargin)
%COLL_SOL2SEG   Append 'coll' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% PROB     = COLL_SOL2SEG(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB - Continuation problem structure.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_sol2seg.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'ddaecoll'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run  = str.get;
if ischar(str.peek)
  soid = str.get;
else
  soid = oid;
end
lab = str.get;

[sol data] = ddaecoll_read_solution(soid, run, lab);  % Extract solution and toolbox data from disk
data       = ddaecoll_get_settings(prob, tbid, data); % Get toolbox settings
data       = ddaecoll_init_data(data, sol.x, sol.y, sol.p);  % Build toolbox data
sol        = ddaecoll_init_sol(data, sol.t, sol.x, sol.y, sol.p);  % Build initial solution guess
prob       = ddaecoll_construct_seg(prob, tbid, data, sol); % Append continuation problem

end
