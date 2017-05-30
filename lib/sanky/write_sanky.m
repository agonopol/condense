function write_sanky(matrix, output, colors)
    st = ShinyTemplate();
    [path, ~, ~] = fileparts(mfilename('fullpath'));
    t = fileread(fullfile(path, '/sanky.tpl'));
    st.loadString(t);
    
    c = struct('matrix', matrix, 'colors', colors);
    fid = fopen(output,'wt');
    fprintf(fid, st.render(c));
    fclose(fid);
end