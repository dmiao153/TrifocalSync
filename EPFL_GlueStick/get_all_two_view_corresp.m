function two_view_corresp = get_all_two_view_corresp(dataset_name, n)
two_view_corresp_name = "./EPFL_GlueStick/"+dataset_name+"_results/";
two_view_corresp = cell(n,n);

for i = 0:n-1
    for j = 0:n-1
        if i < j
            cur_data = [];
            corresp_name = two_view_corresp_name + i + "_" + j + "_matches.txt";
            fid = fopen(corresp_name, 'r');
            for k= 1:4
                tline = fgetl(fid);
            end
            tline = fgetl(fid);

            while length(tline) > 1
                dat = str2double(split(tline, " "));
                cur_data = [cur_data; dat(2:end)'];
                tline = fgetl(fid);
                
            end
            
            
            two_view_corresp{i+1,j+1} = cur_data; 
            fclose(fid);
        end


    end
end


end
