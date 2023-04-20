

fileID = fopen ('p_T_test.txt','w');

k = 1;
while k < 500
[pfinal, error] = first(365);
      err = error;
      p11 = pfinal;
      k = k +1;
   if error < 5000
     fprintf(fileID,'%12s\n','error');  
     fprintf(fileID,'%12.8f\n',err); 
     fprintf(fileID,'%12s\n','parameters');
     fprintf(fileID,'%12.8f\n',p11);
     fprintf(fileID,'%12s\n','----------');
   end
end   

fclose(fileID);