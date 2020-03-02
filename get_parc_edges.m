function [ edge_verts ] = get_parc_edges(vert_weights,nbrs)
% find some edges yo
%
% INPUTS
%
% vert_weights:         vector (length == number of verticies) of label
%                       weights
% nbrs:                 matrix (size1 == number of verticies) of neighbors
%                       per vertex
%
% OUTPUTS
%
% edge_verts            index vector of the edges 
%
% Josh Faskowitz
% Indiana University
% Computational Cognitive Neurosciene Lab
% See LICENSE file for license
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do a little input checking

% if ~exist('nbrs_size','var') || isempty(nbrs_size)
%     % most local
%     nbrs_size = 1 ;
% end
% 
% if nbrs_size < 1
%    error('invalid nbrs_size') 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize seed for potential random numbers
rng(4321)

% copy weights_to_dilate
edge_verts = zeros(size(vert_weights,1),1) ;

% get rid of any possible zero neighbors
nbrs(nbrs==0) = -1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first loop to find verts to fill, 
% we want verts that border a label area

for idx = 1:size(edge_verts,1)

    fill_nbrs = nbrs(idx,:) ;
    % trim any of the potential -1
    fill_nbrs = fill_nbrs(fill_nbrs > 0) ;
    fill_nbrs_vals = vert_weights(fill_nbrs) ; 
    % if this lil set of vals is not all equal
    if sum(fill_nbrs_vals(1) ~= fill_nbrs_vals(:)) > 0
        edge_verts(idx) = 1 ;
    end
end