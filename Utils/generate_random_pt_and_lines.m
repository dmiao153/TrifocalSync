function [real_pts, im_pts, im_lines, real_pvl] = generate_random_pt_and_lines(temp_P)
    rand_pt1 = 10*(rand(1,3) - 0.5);
    pt_real_1 = [rand_pt1, 1]';
    
    rand_pt2 = 10*(rand(1,3)-0.5);
    pt_real_2 = [rand_pt2, 1]';
    

    alpha = 10 * (rand(1)-0.5);

    rand_pt_3 = rand_pt1 + alpha*(rand_pt2-rand_pt1);
    pt_real_3 = [rand_pt_3, 1]';

    
    im1_pt = temp_P{1} * pt_real_3;
    im2_pt = temp_P{2} * pt_real_3;
    im3_pt = temp_P{3} * pt_real_3;
    
    im1_pt = im1_pt/im1_pt(3);
    im2_pt = im2_pt/im2_pt(3);
    im3_pt = im3_pt/im3_pt(3);
    
    
    pvl = vgg_line3d_pv_from_XY(pt_real_1, pt_real_2)';

    pvp_1 = vgg_line3d_Ppv(temp_P{1})';
    im1_line = pvp_1 * pvl;

    im1_line = im1_line/im1_line(3);
    % im1_line = im1_line/norm(im1_line);

    pvp_2 = vgg_line3d_Ppv(temp_P{2})';
    im2_line = pvp_2 * pvl;
    im2_line = im2_line/im2_line(3);
    % im2_line = im2_line/norm(im2_line);

    pvp_3 = vgg_line3d_Ppv(temp_P{3})';
    im3_line = pvp_3 * pvl;

    im3_line = im3_line/im3_line(3);
    % im3_line = im3_line/norm(im3_line);

    
    im_lines = [im1_line, im2_line, im3_line];
    real_pvl = pvl;
    real_pts = [pt_real_1, pt_real_2, pt_real_3];
    im_pts = [im1_pt, im2_pt, im3_pt];
    
end

