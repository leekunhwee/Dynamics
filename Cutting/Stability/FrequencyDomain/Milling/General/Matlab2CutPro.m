% Generate the .frf file for modal analysis in CutPro  
function Matlab2CutPro(f,HX,HY)
CC1='Freq'; CC2='Hz'; CC3='Real';
CC4='m/N'; CC5='Imag'; CC6='Real';

fileID = fopen('FRFX.frf','w+'); % Creat/Clear the FRFX.frf file in write mode
AA=[f(1)';HX(1,1)';HX(1,2)'];
fprintf(fileID,'%14.6e\t%14.12e\t%14.12e\t',AA);
fclose(fileID);
fileID = fopen('FRFX.frf','a+'); % Append to the end 
fprintf(fileID,'%s\t',CC1);
fprintf(fileID,'%s\t',CC2);
fprintf(fileID,'%s\t',CC3);
fprintf(fileID,'%s\t',CC4);
fprintf(fileID,'%s\t',CC5);
fprintf(fileID,'%s\t\n',CC6);
fclose(fileID);
fileID = fopen('FRFX.frf','a+'); % Append to the end 
AA=[f(2:end)';HX(2:end,1)';HX(2:end,2)'];
fprintf(fileID,'%14.6e\t%14.12e\t%14.12e\t\n',AA);
fclose(fileID);

fileID = fopen('FRFY.frf','w+');
BB=[f(1)';HY(1,1)';HY(1,2)'];
fprintf(fileID,'%14.6e\t%14.12e\t%14.12e\t',BB);
fclose(fileID);
fileID = fopen('FRFY.frf','a+');
fprintf(fileID,'%s\t',CC1);
fprintf(fileID,'%s\t',CC2);
fprintf(fileID,'%s\t',CC3);
fprintf(fileID,'%s\t',CC4);
fprintf(fileID,'%s\t',CC5);
fprintf(fileID,'%s\t\n',CC6);
fclose(fileID);
fileID = fopen('FRFY.frf','a+'); 
BB=[f(2:end)';HY(2:end,1)';HY(2:end,2)'];
fprintf(fileID,'%14.6e\t%14.12e\t%14.12e\t\n',BB);
fclose(fileID);