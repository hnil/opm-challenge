function G = computeBoundingBoxes(G, varargin)
% Add cell and/or face bounding boxes.
%
% SYNOPSIS:
%   G = computeBoundingBoxes(G, varargin)
%
% DESCRIPTION:
%   This routine adds fields 'bbox' to cells and/or faces repersenting a 
%   'minimal' bounding box in the givven coordinate system. Useful for grid 
%   searching. 
%
% REQUIRED PARAMETERS:
%   G  - The grid.
%
%  OPTIONAL PARAMETERS:
%   'cells'  : compute bbox for all cells (default true)
%   'faces'  : compute bbox for all faces (default true)
%
% RETURNS:
%   Updated grid G
%
% SEE ALSO:
% findFacesCloseToSegment

opt = struct('faces',    true, ...
             'cells',    false);        
opt = merge_options(opt, varargin{:});

% bbox for cells
if opt.cells
    cno  = cellNodes(G);
    nn   = accumarray(cno(:,1), ones(size(cno,1), 1));
    npos = cumsum([1;nn]);
    % max and min ix after sorting
    minIx = npos(1:end-1);
    maxIx = npos(2:end)-1;
    bbox = nan(G.cells.num, G.griddim);
    for d = 1:G.griddim
        tmp = sortrows([cno(:,1), G.nodes.coords(cno(:,3),d)]);
        bbox(:,d) = tmp(maxIx, 2) - tmp(minIx, 2);
    end
    G.cells.bbox = bbox;
end

% bbox for faces
if opt.faces
    npos = G.faces.nodePos;
    fno  = rldecode((1:G.faces.num)', diff(npos));
    % max and min ix after sorting
    minIx = npos(1:end-1);
    maxIx = npos(2:end)-1;
    bbox = nan(G.faces.num, G.griddim);
    for d = 1:G.griddim
        tmp = sortrows([fno(:,1), G.nodes.coords(G.faces.nodes, d)]);
        bbox(:,d) = tmp(maxIx, 2) - tmp(minIx, 2);
    end
    G.faces.bbox = bbox;
end
end
    