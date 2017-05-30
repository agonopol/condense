function sanky(assigments, fields, output)
    rows = {};
    colors = {};
    for cluster = unique(assigments(end,:))
        id = name(size(assigments,1), fields, cluster);
        branches  = drill(assigments, fields, size(assigments,1) - 1, size(assigments,1), cluster, id, colorassigments(assigments));
        for branch = branches
            rows{size(rows, 2) + 1} = sprintf('[ "%s", "%s", %d ]', ...
                                                branch{:}{1}, ...
                                                branch{:}{2}, ...
                                                branch{:}{3} );
            if isempty(colors)
               colors = {branch{:}{1}, sprintf('"%s"', branch{:}{4})};
            end
            if (not(ismember(colors(:,1), branch{:}{1})))
                colors = [colors; {branch{:}{1}, sprintf('"%s"', branch{:}{4})}];
            end
            if (not(ismember(colors(:,1), branch{:}{2})))
                colors = [colors; {branch{:}{2}, sprintf('"%s"', branch{:}{4})}];
            end
        end
    end
    
    write_sanky(strjoin(rows(:), ',\n '), ...
        output, ...
        strjoin(colors(:, 2), ',\n '));
end


function id = name(iteration, fields, cluster)
    if iteration == 1
       id = sprintf('%s', string(fields(cluster)));
    else
       id = sprintf('I%d/N%d', iteration, cluster);
    end
end

function colors = colorassigments(assigments)
    colors = cell(size(assigments));
    set = distinguishable_colors(length(unique(assigments(end, :))));
    for c = 1:size(colors, 2)
        colors(:,c) = {rgb2hex(set(assigments(end, c),:))};
    end
end

function tree = drill(assigments, fields, level, rlevel, rid, rname, colors)
    tree = {};
    clusters = unique(assigments(level, find(assigments(rlevel,:) == rid)));
    if size(clusters, 2) == 1
        if level == 1
            tree{size(tree, 2) + 1} = {rname, ...
                        name(level, fields, clusters), ...
                        1, ...
                        colors{1, clusters} };
        else
            branches = drill(assigments, fields, level - 1, rlevel, rid, rname, colors);
            for branch = branches
               tree{size(tree, 2) + 1} = branch{:};
           end
        end
    else
        for cluster = clusters
           to = name(level, fields, cluster);
           tree{size(tree, 2) + 1} = {rname, ...
                                      to, ... 
                                      sum(assigments(level, :) == cluster), ...
                                      colors{1, find(assigments(level,:) == cluster, 1)} };
           branches = drill(assigments, fields, level - 1, level, cluster, to, colors);
           for branch = branches
               tree{size(tree, 2) + 1} = branch{:};
           end
        end
    end
end


