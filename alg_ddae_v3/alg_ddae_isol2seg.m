function prob = alg_ddae_isol2seg(prob, oid, varargin)
%COLL_ISOL2SEG   Append 'coll' instance constructed from initial data.
%
% Parse input sequence to construct toolbox data and initial solution guess
% and use this to construct an instance of 'coll'.
%
% PROB     = COLL_ISOL2SEG(PROB, OID, VARARGIN)
% VARARGIN = { @F [(@DFDX | '[]') (@DFDy | '[]') [(@DFDP | '[]')]] T0 X0 [PNAMES] P0 }
% VARARGIN = {gidx, gdim, {@s^2,@s^3,...,@s^nd-1},{g1,g2,...,g_nd}}
% g1 = {@gfunc,(@tau^x | '[]'),(@tau^y | '[]'),(@tau^t | '[]'),(@phi |'[]'),dg,dtaux,dtauy,dtaut}
% dg = {@dgdt,@dgdx,@dgdy,@dgdxs,@dgdys,@dgdphi,@dgdp}
% dtaux = {dtauxDT0,dtauxDT,dtauxDP}
% dtauy = {dtauyDT0,dtauyDT,dtauyDP}
% dtaut = {dtautDT0,dtautDT,dtautDP}

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_isol2seg.m 2839 2015-03-05 17:09:01Z fschild $

fbid              = coco_get_id(oid, 'ddaecoll'); % Create toolbox instance identifier
[fdata, uidx]     = coco_get_func_data(prob, fbid, 'data', 'uidx');
data.ddaecoll_seg = fdata;
data.ddae_uidx    = uidx;

tbid = coco_get_id(oid, 'alg_ddae'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
data.gidx    = str.get;
data.gdim    = str.get;
sj           = str.get;
nd           = numel(sj);
data.nseg    = nd+1;
data.sj      = cell(nd+2,1);
data.sj{1}   = @sinit;
data.sj{end} = @send;
for k=1:nd
    data.sj{k+1} = sj{k};
end

data.g = str.get;





% if is_empty_or_func(str.peek)
%   data.taux = str.get;
%   if is_empty_or_func(str.peek)
%      data.taut = str.get;
%      if is_empty_or_func(str.peek)
%          data.phi = str.get;
%      end
%   end
% end

% data.ghan = str.get;
% data.dgdthan  = [];
% data.dgdxhan  = [];
% data.dgdyhan  = [];
% data.dgdphan  = [];
% data.dgdtdthan  = [];
% data.dgdtdxhan  = [];
% data.dgdtdyhan  = [];
% data.dgdtdphan  = [];
% data.dgdxdxhan  = [];
% data.dgdxdyhan  = [];
% data.dgdxdphan  = [];
% data.dgdydyhan  = [];
% data.dgdydphan  = [];
% data.dgdpdphan  = [];
% 
% if is_empty_or_func(str.peek)
%   data.dgdthan = str.get;
%   if is_empty_or_func(str.peek)
%     data.dgdxhan = str.get;
%     if is_empty_or_func(str.peek)
%         data.dgdyhan = str.get; 
%         if is_empty_or_func(str.peek)
%             data.dgdphan = str.get;
%             if is_empty_or_func(str.peek)
%                 data.dgdtdthan = str.get;
%                 if is_empty_or_func(str.peek)
%                     data.dgdtdxhan = str.get;
%                     if is_empty_or_func(str.peek)
%                         data.dgdtdyhan = str.get;
%                         if is_empty_or_func(str.peek)
%                             data.dgdtdphan = str.get;
%                             if is_empty_or_func(str.peek)
%                                 data.dgdxdxhan = str.get;
%                                 if is_empty_or_func(str.peek)
%                                     data.dgdxdyhan = str.get;
%                                     if is_empty_or_func(str.peek)
%                                         data.dgdxdphan = str.get;
%                                         if is_empty_or_func(str.peek)
%                                             data.dgdydyhan = str.get;
%                                             if is_empty_or_func(str.peek)
%                                                 data.dgdydphan = str.get;
%                                                 if is_empty_or_func(str.peek)
%                                                     data.dgdpdphan = str.get;
%                                                 end
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%   end
% end

data = alg_ddae_init_data(data);        % Build toolbox data
prob = alg_ddae_construct_seg(prob, tbid, data); % Append continuation problem

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end


function y = sinit(T0,T,p)
y = 0;
end


function y = send(T0,T,p)
y = 1;
end
