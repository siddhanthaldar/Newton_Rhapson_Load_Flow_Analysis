% Read no. of buses and machines
prompt = 'No. of buses:';
nbs = input(prompt);
prompt = 'No. of machines:';
nms = input(prompt);

% Read bus_dat
fileID = fopen('test_bus.txt','r');
formatSpec = '%d %d %f %f %f %f %f %f';
bus_dat = fscanf(fileID,formatSpec);
bus_dat =reshape(bus_dat,8,nbs);  %each column contains data for one bus

% Read line_dat
fileID = fopen('test_line.txt','r');
formatSpec = '%d %d %f %f %f %f';
line_dat = fscanf(fileID,formatSpec);
n_lines = length(line_dat)/6;
line_dat =reshape(line_dat,6,n_lines);  %each column contains data for one line

load_bus = zeros(8,nbs-nms);
PV_bus = zeros(8,nms-1);
slack = zeros(8,1);
load_it = 0;
PV_it = 0;

% Arrange bus_data in order load_bus, PV bus, slack bus
for i=1:nbs
    if bus_dat(2,i) == 101
        load_it = load_it + 1;
        load_bus(:,load_it) = bus_dat(:,i);
    elseif bus_dat(2,i) == 102
        PV_it = PV_it + 1;
        PV_bus(:,PV_it) = bus_dat(:,i);
    elseif bus_dat(2,i) == 103
        slack(:,1) = bus_dat(:,i);
    end
end
load_bus
PV_bus
%sort load bus
load = zeros(nbs-nms,2);
bus = zeros(8,nbs-nms);
for i=1:nbs-nms
    load(i,1) = i;
    load(i,2) = load_bus(1,i);
end

load(:,2) = sort(load(:,2));
for i=1:nbs-nms
    bus(:,load(i,1)) = bus_dat(:,load(i,2));
end
load_bus = bus;


%sort PV bus
load1 = zeros(nms-1,2);
bus1 = zeros(8,nms-1);
for i=1:nms-1
    load1(i,1) = i;
    load1(i,2) = PV_bus(1,i);
end
load1(:,2) = sort(load1(:,2));
for i=1:nms-1
    bus1(:,load1(i,1)) = bus_dat(:,load1(i,2));
end
PV_bus = bus1;

bus_dat = [load_bus, PV_bus, slack];


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
tol = 0.0001;

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
    
    if abs(mismatch)<tol
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

    %{
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
    %}

    %%{
    % Gradients using Gaussian Elimination
    gauss_var = [J';mismatch];
    [d,var_num] = size(gauss_var);
    grad = zeros(var_num,1);

    % Obtain augmented matrix
    for x=1:var_num-1
        for y=x+1:var_num
           for z=var_num+1:-1:1
               gauss_var(z,y) = gauss_var(z,y) - gauss_var(z,x)*gauss_var(x,y)/gauss_var(x,x);
           end
        end
    end

    % Solve
    for x=var_num:-1:1
        for y=var_num:-1:x+1
            grad(x,1) = grad(x,1)+ grad(y,1)*gauss_var(y,x);
        end
        grad(x,1) = (gauss_var(var_num+1,x)-grad(x,1))/gauss_var(x,x);
    end

    % Update theta
    for x=1:nbs-1
        bus_dat(4,x) = bus_dat(4,x) + grad(x,1);
    end
    % Update V
    for x=1:nbs-nms
        bus_dat(3,x) = bus_dat(3,x) + grad(x+nbs-1,1)*bus_dat(3,x);  % since del_V/V calculated, multiple by V
    end

end

% Print bus voltage magnitude(in p.u) and angle after convergence
disp('   ')
disp('Bus Voltage Magnitude(in p.u.) : ')
disp(bus_dat(3,:))

disp('Angle :')
disp(bus_dat(4,:))

% Active and reactive power calculation at all buses
P = zeros(1,nbs);
Q = zeros(1,nbs);
for i=1:nbs
    for it=1:nbs
        ang = bus_dat(4,i)-bus_dat(4,it);
        re = real(Ybus(i,it));
        img = imag(Ybus(i,it)); 
        P(1,i) = P(1,i) + bus_dat(3,i)*bus_dat(3,it)*(re*cos(ang)+img*sin(ang));
        Q(1,i) = Q(1,i) + bus_dat(3,i)*bus_dat(3,it)*(re*sin(ang)-img*cos(ang));
    end
end

disp('Active power at all buses :')
disp(P)
disp('Reactive power at all buses :')
disp(Q)
disp('Total active power loss')
disp(sum(P))

% Reactive power generation at PV buses
react_gen = zeros(1,nms-1);
tot_react_pow_gen = 0;
for i=nbs-nms+1:nbs-1
    react_gen(1,i-(nbs-nms)) = Q(1,i) + bus_dat(8,i);
    tot_react_pow_gen = tot_react_pow_gen + react_gen(1,i-(nbs-nms));
end

disp('Reactive power generated at P-V buses :')
disp(react_gen)
disp('Total reactive power generated at P-V buses :')
disp(tot_react_pow_gen)

% Calculate complex power flow
complex_pow_flow = zeros(n_lines,3);
for i=1:n_lines
    Z = line_dat(3,i) + 1i*line_dat(4,i);
    A = 1 + Z*(1i*line_dat(5,i))/2;
    B = Z;
    Vs = bus_dat(3,line_dat(1,i))*exp(1i*bus_dat(4,line_dat(1,i)));%*180/pi);
    Vr = bus_dat(3,line_dat(2,i))*exp(1i*bus_dat(4,line_dat(2,i)));%*180/pi);
    Ir = (Vs - A*Vr)/B;
    
    complex_pow_flow(i,1) = line_dat(1,i); % Sending end
    complex_pow_flow(i,2) = line_dat(2,i); % Receiving end
    complex_pow_flow(i,3) = Vr*conj(Ir);
end

disp('Complex power flow in lines(from sending to receiving end) :')
disp(complex_pow_flow)












