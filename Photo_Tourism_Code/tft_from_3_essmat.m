% [U,V,np,F] = generate_random_projection_cameras(3, false, false);
% 
% 
% F21 = vgg_F_from_P(np{1}, np{2});
% F31 = vgg_F_from_P(np{1}, np{3});
% F32 = vgg_F_from_P(np{2}, np{3});
% 
% [P1,P2,P3] = tft_from_3_essmat2(F21, F31, F32);
% T_from_P({P1,P2,P3})/frob(T_from_P({P1,P2,P3}))
% T_from_P({np{1},np{2},np{3}}) / frob(T_from_P({np{1},np{2},np{3}}))
% 
% (F31 / norm(F31)) - (vgg_F_from_P(P1, P3) / norm(vgg_F_from_P(P1, P3)))
% 
% (F32 / norm(F32)) - (vgg_F_from_P(P2, P3) / norm(vgg_F_from_P(P2, P3)))
% 
% (F21 / norm(F21)) - (vgg_F_from_P(P1, P2) / norm(vgg_F_from_P(P1, P2)))

% np{3}' * F31 * np{1} 
% 
% P1 = eye(3,4);
% P2 = vgg_P_from_F(F21);
% 
% fp1 = F31 * P1;
% fp2 = F32 * P2;
% A1 = order_equations(fp1);
% A2 = order_equations(fp2);
% 
% [u,s,v] = svd([A1;A2]);
% 
% P3 = reshape(v(:,end), [3,4]);
% 
% (F32 / norm(F32)) - (vgg_F_from_P(P2, P3) / norm(vgg_F_from_P(P2, P3)))



function [P1,P2,P3] = tft_from_3_essmat(F21, F31, F32)

P1 = eye(3,4);
P2 = vgg_P_from_F(F21);

fp1 = F31 * P1;
fp2 = F32 * P2;

A1 = order_equations(fp1);
A2 = order_equations(fp2);
[u,s,v] = svd([A1;A2]);

P3 = reshape(v(:,end), [3,4]);
end



function A = order_equations(fp)
%%% diagonal equations
A = zeros(10,12);
A(1,:)=[fp(1,1),fp(2,1),fp(3,1),0,0,0,0,0,0,0,0,0];
A(2,:)=[0,0,0,fp(1,2),fp(2,2),fp(3,2),0,0,0,0,0,0];
A(3,:)=[0,0,0,0,0,0,fp(1,3),fp(2,3),fp(3,3),0,0,0];
A(4,:)=[0,0,0,0,0,0,0,0,0,fp(1,4),fp(2,4),fp(3,4)];

%%% offdiagonal equations
A(5,:)=[fp(1,2),fp(2,2),fp(3,2),fp(1,1),fp(2,1),fp(3,1),0,0,0,0,0,0]; % 12
A(6,:)=[fp(1,3),fp(2,3),fp(3,3),0,0,0,fp(1,1),fp(2,1),fp(3,1),0,0,0]; % 13
A(7,:)=[fp(1,4),fp(2,4),fp(3,4),0,0,0,0,0,0,fp(1,1),fp(2,1),fp(3,1)]; % 14
A(8,:)=[0,0,0,fp(1,3),fp(2,3),fp(3,3),fp(1,2),fp(2,2),fp(3,2),0,0,0]; % 23
A(9,:)=[0,0,0,fp(1,4),fp(2,4),fp(3,4),0,0,0,fp(1,2),fp(2,2),fp(3,2)]; % 24
A(10,:)=[0,0,0,0,0,0,fp(1,4),fp(2,4),fp(3,4),fp(1,3),fp(2,3),fp(3,3)]; % 34

end

