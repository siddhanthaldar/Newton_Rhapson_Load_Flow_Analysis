% Read no. of buses and machines
prompt = 'No. of buses:';
nbs = input(prompt);
prompt = 'No. of machines:';
nms = input(prompt);

% Read bus_dat
fileID = fopen('bus_dat.txt','r');
formatSpec = '%d %d %f %f %f %f %f %f';
bus_dat = fscanf(fileID,formatSpec);
bus_dat =reshape(bus_dat,8,nbs);  %each column contains data for one bus

% Read line_dat
fileID = fopen('line_dat.txt','r');
formatSpec = '%d %d %f %f %f %f';
line_dat = fscanf(fileID,formatSpec);
n_lines = length(line_dat)/6;
line_dat =reshape(line_dat,6,n_lines);  %each column contains data for one line

% Formulate Ybus
Ybus = zeros(nbs,nbs);
for x=1:nbs   %row
    for y=x:nbs
       if x==y
           for k=1:n_lines
               if (line_dat(1,k) == x)
                   z = line_dat(3,k)+1i*line_dat(4,k);
                   Ybus(x,y) = Ybus(x,y) + (1/(abs(line_dat(6,k)))^2)*(1/z + 1i*line_dat(5,k)/2);
               elseif (line_dat(2,k) == x)
                   z = line_dat(3,k)+1i*line_dat(4,k);
                   Ybus(x,y) = Ybus(x,y) + (1/z + 1i*line_dat(5,k)/2);
               end
           end
       else 
           for k=1:n_lines
               if (line_dat(1,k) == x) && (line_dat(2,k) == y)
                   z = line_dat(3,k)+1i*line_dat(4,k);
                   Ybus(x,y) = Ybus(x,y) + (-1/conj(line_dat(6,k)))*(1/z);
               elseif (line_dat(1,k) == y) && (line_dat(2,k) == x)
                   z = line_dat(3,k)+1i*line_dat(4,k);
                   Ybus(x,y) = Ybus(x,y) + (-1/conj(line_dat(6,k)))*(1/z);
               end
           end
           Ybus(y,x) = Ybus(x,y);
       end
    end
end

% Newton Rhapson Solver
% tolerance
tol = 0.001;

while 1
    % Mismatch vector calculation
    del_P = zeros(1,nbs-1);
    del_Q = zeros(1,nbs-nms);
    cnt_P = 0;
    cnt_Q = 0;
    P = zeros(1,nbs-1);
    Q = zeros(1,nbs-nms);

    for i=1:nbs
        if bus_dat(2,i) == 101
            cnt_P = cnt_P + 1;
            cnt_Q = cnt_Q + 1;
            p_spec = bus_dat(5,i) - bus_dat(7,i);
            q_spec = bus_dat(6,i) - bus_dat(8,i);
            for it=1:nbs
                ang = bus_dat(4,i)-bus_dat(4,it);
                re = real(Ybus(i,it));
                img = imag(Ybus(i,it)); 
                P(1,i) = P(1,i) + bus_dat(3,i)*bus_dat(3,it)*(re*cos(ang)+img*sin(ang));
                Q(1,i) = Q(1,i) + bus_dat(3,i)*bus_dat(3,it)*(re*sin(ang)-img*cos(ang));
            end
            del_P(cnt_P) = p_spec - P(1,i);
            del_Q(cnt_Q) = q_spec - Q(1,i);
        elseif bus_dat(2,i) == 102
            cnt_P = cnt_P + 1;
            p_spec = bus_dat(5,i) - bus_dat(7,i);
            for it=1:nbs
                ang = bus_dat(4,i)-bus_dat(4,it);
                re = real(Ybus(i,it));
                img = imag(Ybus(i,it)); 
                P(1,i) = P(1,i) + bus_dat(3,i)*bus_dat(3,it)*(re*cos(ang)+img*sin(ang));
            end

            del_P(cnt_P) = p_spec - P(1,i);
        end
    end

    mismatch = [del_P,del_Q];

    if sum(mismatch>0.1)==0
        break;
    end

    % Jacobian
    dim = 2*nbs-nms-1;  %size of Jacobian matrix
    J = zeros(dim,dim);

    % Jp,theta
    for x=1:nbs-1
        for y=1:nbs-1
            if x==y
                if x<=nbs-nms
                    J(x,y) = -Q(1,x)-bus_dat(3,x)*bus_dat(3,y)*imag(Ybus(x,y));
                else
                    for z=1:nbs
                        if z==x
                            continue
                        else
                            ang = bus_dat(4,x)-bus_dat(4,z);
                            re = real(Ybus(x,z));
                            img = imag(Ybus(x,z));
                            J(x,y) = J(x,y)- bus_dat(3,x)*bus_dat(3,z)*(re*sin(ang)-img*cos(ang)); 
                        end
                    end
                end
            else
                ang = bus_dat(4,x)-bus_dat(4,y);
                re = real(Ybus(x,y));
                img = imag(Ybus(x,y));
                J(x,y) = bus_dat(3,x)*bus_dat(3,y)*(re*sin(ang)-img*cos(ang));
            end
        end
    end

    % Jq,theta
    for x=1:nbs-nms%dim-(nbs-1)
        for y=1:nbs-1
            if x==y
                J(x+nbs-1,y) = P(1,x)- bus_dat(3,x)*bus_dat(3,y)*real(Ybus(x,y));
            else
                ang = bus_dat(4,x)-bus_dat(4,y);
                re = real(Ybus(x,y));
                img = imag(Ybus(x,y));
                J(x+nbs-1,y) = -1.0*bus_dat(3,x)*bus_dat(3,y)*(re*cos(ang)+img*sin(ang));
            end        
        end
    end

    % Jp,V
    for x=1:nbs-1
        for y=1:nbs-nms
            if x==y
                J(x,y+nbs-1) = P(1,x)+ bus_dat(3,x)*bus_dat(3,y)*real(Ybus(x,y));
            else
                ang = bus_dat(4,x)-bus_dat(4,y);
                re = real(Ybus(x,y));
                img = imag(Ybus(x,y));
                J(x,y+nbs-1) = bus_dat(3,x)*bus_dat(3,y)*(re*cos(ang)+img*sin(ang));
            end        
        end
    end

    % Jq,V
    for x=1:nbs-nms
        for y=1:nbs-nms
            if x==y
                J(x+nbs-1,y+nbs-1) = Q(1,x)- bus_dat(3,x)*bus_dat(3,y)*imag(Ybus(x,y));
            else
                ang = bus_dat(4,x)-bus_dat(4,y);
                re = real(Ybus(x,y));
                img = imag(Ybus(x,y));
                J(x+nbs-1,y+nbs-1) = bus_dat(3,x)*bus_dat(3,y)*(re*sin(ang)-img*cos(ang));
            end        
        end
    end

    % Gradients using inverse
    grad = J\mismatch';        

    % Update theta
    for x=1:nbs-1
        bus_dat(4,x) = bus_dat(4,x) + grad(x);
    end
    % Update V
    for x=1:nbs-nms
        bus_dat(3,x) = bus_dat(3,x) + grad(x+nbs-1)*bus_dat(3,x);  % since del_V/V calculated, multiple by V
    end

end












