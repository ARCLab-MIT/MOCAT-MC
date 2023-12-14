function LinkZoom(fig)
%
% LinkZoom(fig)
%
%   links zoom scale of the 2 sublots in figure.
%
if ~exist('fig','var')
    fig=gcf;
end
figure(fig)

% Get handles for the axes on this figure whose Tag is blank
axes_handles = findobj(fig,'Type','axes','Tag','');
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
current_axis = find(axes_handles == gca);
other_axis = find(axes_handles ~= gca);

% Determine original limits for all axes_handles.
n_axes = length(axes_handles);
orig_limits = zeros(n_axes,4);
for  n = 1:n_axes
    set(fig,'CurrentAxes',axes_handles(n));
    orig_lims(n,1:4) = zoom(fig,'getlimits');
end
%compute normalized x,y zoom limits for current axis.
xzoom_norm = [(xlim(1)-orig_lims(current_axis,1))/(orig_lims(current_axis,2)-orig_lims(current_axis,1))...
    (xlim(2)-orig_lims(current_axis,1))/(orig_lims(current_axis,2)-orig_lims(current_axis,1))];
yzoom_norm = [(ylim(1)-orig_lims(current_axis,1))/(orig_lims(current_axis,2)-orig_lims(current_axis,1))...
    (xlim(2)-orig_lims(current_axis,3))/(orig_lims(current_axis,4)-orig_lims(current_axis,3))];

%compute x,y, limits for other axis from normalized zoom limits.
xlim_other = orig_lims(other_axis,1)+xzoom_norm.*(orig_lims(other_axis,2)-orig_lims(other_axis,1));
ylim_other = orig_lims(other_axis,3)+yzoom_norm.*(orig_lims(other_axis,4)-orig_lims(other_zxis,3));

set(axes_handles(other_axis),'Xlim',xlim_other,'YLim',ylim_other);

return