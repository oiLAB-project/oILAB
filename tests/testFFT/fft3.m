clc; clear; close all;

m3D=rand(10,15,20);

fileID = fopen('m3D.txt','w');

for k=1:1:size(m3D,3)
    for j=1:1:size(m3D,2)
        for i=1:1:size(m3D,1)
            fprintf(fileID,'%12.8f\n',m3D(i,j,k));
        end
    end
end
fclose(fileID);

n3D=fftn(m3D);
fileID = fopen('fft_m3D.txt','w');

for k=1:1:size(m3D,3)
    for j=1:1:size(m3D,2)
        for i=1:1:size(m3D,1)
            fprintf(fileID,'%20.12f   %20.12f\n',real(n3D(i,j,k)),imag(n3D(i,j,k)));
        end
    end
end
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%  FFT 2D   %%%%%%%%%%%%%%%%%%

m2D=rand(15,20);

fileID = fopen('m2D.txt','w');


for j=1:1:size(m2D,2)
    for i=1:1:size(m2D,1)
        fprintf(fileID,'%12.8f\n',m2D(i,j));
    end
end
fclose(fileID);

n2D=fftn(m2D);
fileID = fopen('fft_m2D.txt','w');

for j=1:1:size(m2D,2)
    for i=1:1:size(m2D,1)
        fprintf(fileID,'%20.12f   %20.12f\n',real(n2D(i,j)),imag(n2D(i,j)));
    end
end

fclose(fileID);