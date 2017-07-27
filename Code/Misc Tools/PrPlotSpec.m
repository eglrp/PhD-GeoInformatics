function PrPlotSpec(data)
    s = size(data);
    ll = getlablist(data);
    nl = getnlab(data);
    nll = unique(nl);
    figure;
    colours = get(gca, 'ColorOrder');
    lines = {'-','--',':','.-'};
    markers = {'none';'x';'o';'^';'square';'>';'<'};
    for i = 1:length(nll)
        h_ = plot(+data(nl==nll(i), :)', 'Color', colours(1 + rem(i, size(colours, 1)), :), ...
            'LineStyle', lines{1 + rem(floor(i/size(colours, 1)), length(lines))}); %, 'Marker', markers{1 + rem(floor(i/size(colours, 1)), length(markers))});
        hold on;
        h(i) = h_(1);
    end
    if isstr(ll)
        legend(h, cellstr(ll));
    else
        legend(h, cellstr(num2str(ll(:))));
    end
    grid on;
return

