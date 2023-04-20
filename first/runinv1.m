

fileID = fopen ('p_T_inv.txt','w');

k = 1;
while k < 500
[pfinal, error] = invfirst(365);
      err = error;
      p11 = pfinal;
      k = k +1;
   if error < 100
     fprintf(fileID,'%12s\n','error H = 20/abc');  
     fprintf(fileID,'%12.8f\n',err); 
     fprintf(fileID,'%12s\n','parameters');
     fprintf(fileID,'%12.8f\n',p11);
     fprintf(fileID,'%12s\n','----------');
   end
end   

fclose(fileID);