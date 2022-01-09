function [Xmesh,Ymesh,Xvec,Yvec] = get_XY_grid( X_list, Y_list, X_step, Y_step )

X = round(X_list,3); 
Y = round(Y_list,3);

% Xsort = sort(unique(X));
% Ysort = sort(unique(Y));
% 
% Xstep = mode(abs(Xsort(2:end)-Xsort(1:end-1)));
% Ystep = mode(abs(Ysort(2:end)-Ysort(1:end-1)));

Xvec = min(X) : X_step : max(X);
Yvec = min(Y) : Y_step : max(Y); min(Y)

if Xvec(end) < max(X), Xvec(end+1) = Xvec(end)+X_step; end
if Yvec(end) < max(Y), Yvec(end+1) = Yvec(end)+Y_step; end

[Xmesh,Ymesh] = meshgrid( Xvec,Yvec );

end
