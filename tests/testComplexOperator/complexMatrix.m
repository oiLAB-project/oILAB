clc; clear; close all;


n=100
B=rand(n,n);
k=20;
B=(B+transpose(B))/2;
C=rand(n,n);
C=(C-transpose(C))/2;

fileID1 = fopen('realMatrix.txt','w');
fileID2 = fopen('complexMatrix.txt','w');

for j=1:1:size(B,2)
    for i=1:1:size(B,1)
        fprintf(fileID1,'%12.8f\n',B(i,j));
        fprintf(fileID2,'%12.8f\n',C(i,j));

    end
end

fclose(fileID1);
fclose(fileID2);

D=B+C*1i;
A=[B -C;
   C  B];


evalues_D= eigs(D,k);
evalues_A= eigs(A,2*k);


fileID = fopen('evalues.txt','w');
fprintf(fileID,'%12.8f\n',evalues_A);
fclose(fileID);



