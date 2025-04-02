function stacked_exterior_square = calculate_stacked_exterior_square(np,n)

    stacked_exterior_square = zeros(3*n, 6);
    
    for i=1:n
        stacked_exterior_square(3*(i-1)+1:3*i, :) = exterior_square(np{i});
    end
    
end

        
%     for i=1:n        
%         matrix = zeros(3,6);
%         camera = np{i};
%         for j = 1:3
%             if j==1
%                 m = 2;
%                 n = 3;
%             end
%             if j==2
%                 m = 1;
%                 n = 3;
%             end
%             if j==3
%                 m = 1;
%                 n = 2;
%             end
%             
%             matrix(j,1) = det([camera(2,1), camera(2,2); camera(3,1), camera(3,2)]);                
%             matrix(j,2) = det([camera(2,1), camera(2,3); camera(3,1), camera(3,3)]);
% 
% 
%                 
% 
%               
% 
%             
%         end
%         
%     end

% end