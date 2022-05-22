% Copyright (c) 2012, Adam Cohen
% All rights reserved.
% 
% Redistribution in source or binary forms, with or without modification,
% is not permitted
%        
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function out = spikefind(dat, thresh);
% Function out = spikefind(dat, thresh);
% for each time that dat goes above thresh, returns the index of the local
% maximum in dat
% DRH and AEC 24 Feb. 2011

spikeon = find(dat(2:end) > thresh & dat(1:end-1) < thresh);
spikeoff = find(dat(2:end) < thresh & dat(1:end-1) > thresh);


if isempty(spikeon) | isempty (spikeoff)
    out = [];
    return
end;

if spikeoff(1) < spikeon(1);
    spikeoff(1) = [];
end;
if spikeon(end) > spikeoff(end);
    spikeon(end) = [];
end;

for j = 1:length(spikeon);
    [y, indx] = max(dat(spikeon(j):spikeoff(j)));
    out(j) = indx + spikeon(j)-1;
end;

