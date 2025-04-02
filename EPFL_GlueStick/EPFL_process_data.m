%%%% -=============================================================- %%%%
%%%% This function calculates the line triplet correspondences from
%%%% pairwise matches from GlueStick
%%%% -=============================================================- %%%%

function [triplet_points, triplet_lines, num_cams] = EPFL_process_data(dataset_name, n)

num_cams = n;
filenames = dir("EPFL_GlueStick/" + dataset_name +"_results/*.txt");
suffix = "EPFL_GlueStick/" + dataset_name +"_results/";
point_matches = cell(num_cams, num_cams);
line_matches = cell(num_cams, num_cams);

for i = 1:length(filenames)
    cur_name = split(filenames(i).name, '_');
    im1_num = str2double(cur_name{1});
    im2_num = str2double(cur_name{2});
    [pmatches,lmatches] = ETH3D_read_results(suffix + filenames(i).name);
    line_matches{im1_num+1, im2_num+1} = lmatches;
    point_matches{im1_num+1,im2_num+1} = pmatches;
end

triplet_lines = cell(num_cams, num_cams, num_cams);
for i = 1:num_cams
    for j = 1:num_cams
        for k= 1:num_cams
            if i < j && j < k 
                cur_lmatches = line_matches{i,j};
                cur_lmatches2 = line_matches{j,k};
                if ~size(cur_lmatches,1) || ~size(cur_lmatches2,1)
                    continue
                end

                % Can further restrict lines that appear in all three pairs
                % cur_lmatches3 = line_matches{i,k};

                [ijjk,ijind,jkind] = intersect(cur_lmatches(:,5:8), cur_lmatches2(:,1:4), 'rows');
                % triplet_lines{i,j,k} = calc_lines_from_endpoints([cur_lmatches(ijind, :), cur_lmatches2(jkind,5:8)]);
                triplet_lines{i,j,k} = [cur_lmatches(ijind, :), cur_lmatches2(jkind,5:8)];
            end
        end
    end
end

triplet_points = ETH3D_process_pmatches(point_matches, num_cams);

end
%%%% -=============================================================- %%%%


%%%% -=============================================================- %%%%
%%%%     This function processes the pmatches from GlueStick          %%%%
%%%% -=============================================================- %%%%

function triplet_points = ETH3D_process_pmatches(pmatches, num_cams)
triplet_points = cell(num_cams,num_cams,num_cams);
for i = 1:num_cams
    for j = 1:num_cams
        for k = 1:num_cams
            if i < j && j < k
                cur_pmatches = pmatches{i,j};
                cur_pmatches2 = pmatches{j,k};
                if ~size(cur_pmatches,1) || ~size(cur_pmatches2,1)
                    continue;
                end
                [ijjk,ijind,jkind] = intersect(cur_pmatches(:,3:4), cur_pmatches2(:,1:2), 'rows');
                % triplet_lines{i,j,k} = calc_lines_from_endpoints([cur_lmatches(ijind, :), cur_lmatches2(jkind,5:8)]);
                triplet_points{i,j,k} = [cur_pmatches(ijind, :), cur_pmatches2(jkind,3:4)];

            end


        end
    end
end



end


%%%% -=============================================================- %%%%
%%%%     This function processes the results from GlueStick          %%%%
%%%% -=============================================================- %%%%
function [pmatches,lmatches] = ETH3D_read_results(filename)

fid = fopen(filename, 'r');
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);
tline = fgetl(fid);

pmatches = [];

while ischar(tline) && length(tline) > 0
    num_line = str2double(split(tline,' '));
    pmatches = [pmatches; num_line(2:end)'];
    tline = fgetl(fid);
end

tline = fgetl(fid);

lmatches = [];

while ischar(tline) && length(tline) > 0
    num_line = str2double(split(tline,' '));
    lmatches = [lmatches; num_line(2:end)'];
    tline = fgetl(fid);
end

fclose(fid);

end
%%%% -=============================================================- %%%%

%%%% -=============================================================- %%%%
%%%%     This function calcualtes lines from a N x 12 matrix and returns a 
%%%%     N x 9 matrix of the line from endpoints. 
%%%% -=============================================================- %%%%
function calc_lines = calc_lines_from_endpoints(endpoints)
    calc_lines = zeros(size(endpoints,1), 9);
    for i = 1:size(endpoints,1)
        hp1 = [endpoints(i,1:2),1];
        hp2 = [endpoints(i,3:4),1];
        line1 = cross(hp1, hp2);
        calc_lines(i,1:3) = line1/norm(line1);

        hp1 = [endpoints(i,5:6),1];
        hp2 = [endpoints(i,7:8),1];
        line2 = cross(hp1, hp2);
        calc_lines(i,4:6) = line2/norm(line2);

        hp1 = [endpoints(i,9:10),1];
        hp2 = [endpoints(i,11:12),1];
        line3 = cross(hp1, hp2);
        calc_lines(i,7:9) = line3/norm(line3);
    end
end
%%%% -=============================================================- %%%%
