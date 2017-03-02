function [hmap] = hmap(obj)
    stats = obj.clusterStats{end}('centroids');
    sizes = obj.clusterStats{end}('size');
    data = [];
    for cluster = 1:size(stats)
       row = [];
       for channel = 1:size(obj.channels)
          row = [row, stats{cluster}(contractor.channels{channel})];
       end
       if cluster == 1
         data = [data;0, sizes(cluster), [row]];
       else
         data = [data;norm(row - data(1,3:end)), sizes(cluster), [row]];
       end
    end
    [~,I]=sort(data(1,:));
    data=data(:,I);
    hmap = HeatMap(data(:,3:end), 'RowLabels', data(:,2), 'ColumnLabels', obj.channels);
end