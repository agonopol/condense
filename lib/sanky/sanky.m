function sanky(assigments, fields, output)
    rows = {};
    for cluster = unique(assigments(end,:))
        id = name(size(assigments,1), fields, cluster);
        branches = drill(assigments, fields, size(assigments,1) - 1, size(assigments,1), cluster, id);
        for branch = branches
            rows{size(rows, 2) + 1} = branch{:};
        end
    end
    st = ShinyTemplate();
    [path, ~, ~] = fileparts(mfilename('fullpath'));
    t = fileread(fullfile(path, '/sanky.tpl'));
    st.loadString(t);
    matrix = strjoin([rows{:}], ',\n ');
    c = struct('matrix', matrix);
    fid = fopen(output,'wt');
    fprintf(fid, st.render(c));
    fclose(fid);
end


function id = name(iteration, fields, cluster)
    if iteration == 1
       id = string(fields(cluster));
    else
       id = sprintf('I%d/N%d', iteration, cluster);
    end
end

function tree = drill(assigments, fields, level, rlevel, rid, rname)
    tree = {};
    clusters = unique(assigments(level, find(assigments(rlevel,:) == rid)));
    if size(clusters, 2) == 1
        if level == 1
            tree{size(tree, 2) + 1} = {sprintf('[ "%s", "%s", %d ]', rname, name(level, fields, clusters), 1)};
        else
            branches = drill(assigments, fields, level - 1, rlevel, rid, rname);
            for branch = branches
               tree{size(tree, 2) + 1} = branch{:};
           end
        end
    else
        for cluster = clusters
           to = name(level, fields, cluster);
           tree{size(tree, 2) + 1} = {sprintf('[ "%s", "%s", %d ]', ...
                                        rname, ...
                                        to, ... 
                                        sum(assigments(level, :) == cluster))};
           branches = drill(assigments, fields, level - 1, level, cluster, to);
           for branch = branches
               tree{size(tree, 2) + 1} = branch{:};
           end
        end
    end
end


