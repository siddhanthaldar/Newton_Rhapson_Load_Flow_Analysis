% Read no. of variables
prompt = 'No. of variable:';
var_num = input(prompt);

% Read bus_dat
fileID = fopen('gauss_var.txt','r');
formatSpec ='';
for i=1:var_num+1
    formatSpec = strcat(formatSpec,'%f');
    if i<=var_num
        formatSpec = strcat(formatSpec,' ');
    end
end


gauss_var = fscanf(fileID,formatSpec);
gauss_var =reshape(gauss_var,var_num+1,var_num);  %each column contains data for one bus
sol = zeros(var_num,1);

% Obtain augmented matrix
for i=1:var_num-1
    for j=i+1:var_num
       for k=var_num+1:-1:1
           gauss_var(k,j) = gauss_var(k,j) - gauss_var(k,i)*gauss_var(i,j)/gauss_var(i,i);
       end
    end
end

% Solve
for i=var_num:-1:1
    for j=var_num:-1:i+1
        sol(i,1) = sol(i,1) + sol(j,1)*gauss_var(j,i);
    end
    sol(i,1) = (gauss_var(var_num+1,i)-sol(i,1))/gauss_var(i,i);
end

