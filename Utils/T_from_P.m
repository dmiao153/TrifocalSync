%%% T_from_P calculates the trifocal tensor using the determinant method
%%% described in Hartley and Zisserman in Chapter 17
%%  P: Cell array of three projection matrices
function T = T_from_P(P)

P1 = P{1};
P2 = P{2};
P3 = P{3};

T=zeros(3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
%             temp_P = zeros(2,4);
%             if i == 1
%                 temp_P = P1([2 3],:);
%             end
%             if i == 2
%                 temp_P = P1([1 3],:);
%             end
%             if i == 3
%                 temp_P = P1([1 2],:);
%             end
%             T(j,k,i)=(-1)^(i+1)*det([temp_P;P2(j,:);P3(k,:)]);

            T(j,k,i)=(-1)^(i+1)*det([P1([1:(i-1) (i+1):3],:);P2(j,:);P3(k,:)]);
        end
    end
end
 
% T=T/norm(T(:));

end








