[leftB,~] = qr(randn(3,3));
leftB = (sign(det(leftB)))*leftB; %for safety, doesn't appear to be necessary 
rightB = randn(3,1);
rightB = rightB/norm(rightB,'fro');%wlog right column of ground-truth B has norm 1
B0 = [leftB rightB];
[leftC,~] = qr(randn(3,3));
leftC = (sign(det(leftC)))*leftC; %for safety, doesn't appear to be necessary 
C0 = [leftC randn(3,1)];
l0 = randn(5,1);
g = [l0(1) 0 0 0;
     0 l0(1) 0 0;
     0 0 l0(1) 0;
     l0(2) l0(3) l0(4) l0(5)];
 B = B0*inv(g);
 C = C0*inv(g);
 [B1,C1,l1]=calibrate(B,C);
 diff1 = norm(B1 - B0, 'fro')+norm(C1 - C0, 'fro');
 diff2 = norm(B1*diag([1 1 1 -1]) - B0, 'fro')+norm(C1*diag([1 1 1 -1]) - C0, 'fro');
 min(diff1,diff2) %should be 0... and it is :)
 
 %should test robustness to noise too... but daniel can try that

  