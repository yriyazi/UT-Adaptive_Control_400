if numel(A)~=numel(A_m)
    error('numel(A)~=numel(A_m)')
end

if numel(B)~=numel(B_m)
    error('numel(A)~=numel(A_m)')
end

% if (conv(B_plus,B_minus))~=B
%     error('B_plus*B_minus~=B')
% end
% 
% if conv(B_prim,B_minus)~=B_m/DC_gain
%     error('B_prim*B_minus~=B_m')
% end