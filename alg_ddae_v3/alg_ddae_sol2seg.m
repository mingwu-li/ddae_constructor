function prob = alg_ddae_sol2seg(prob, oid, varargin)
%COLL_SOL2SEG   Append 'coll' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% PROB     = COLL_SOL2SEG(PROB, OID, VARARGIN)
% VARARGIN = { GIDX RUN [SOID] LAB }
%
% PROB - Continuation problem structure.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_sol2seg.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'alg_ddae'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
gidx = str.get;
run  = str.get;
if ischar(str.peek)
  soid = str.get;
else
  soid = oid;
end
lab  = str.get;

fid               = sprintf('gseg.%s',gidx);
soid              = coco_get_id(soid, 'alg_ddae');
soid              = coco_get_id(soid, fid);
data              = coco_read_solution(soid, run, lab);
fbid              = coco_get_id(oid, 'ddaecoll');        
[fdata, uidx]     = coco_get_func_data(prob, fbid, 'data', 'uidx');
data.ddaecoll_seg = fdata; % get updated data
data.ddae_uidx    = uidx;
data              = alg_ddae_init_data(data);                 % Build toolbox data
prob              = alg_ddae_construct_seg(prob, tbid, data); % Append continuation problem

end
