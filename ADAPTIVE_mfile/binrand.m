%% EXER_Start UP
%%system identification exercises-first series


% -starat up with coding m-file, which could generate a vector of random
% signals,u, which could have three values of -1,0 and 1  u=binrand(T,N,Tmin,t0,dist)
% In which T is the time vector, N is the number of pulses, tmin is the
% minimum pulse width, to is considered as starting time, and dist as the
% random distribution, by which the signal varies between the three values
function u=binrand(T,N,tmin,t0,dist) %%Describing Binrand Function

u=zeros(1,length(T));
pulse_value=zeros(1,N);

% considering t0; 
 for i=1:1:t0
     u(1,i)=0;
 end
 
%%generating random number considering three values and dist form
for pulse=1:1:N
    switch dist
        case 'uniform'
            bni=unifrnd(-1.5,1.5,1,1);
        case 'normal'
            bni=normrnd(0,1,1,1);
    end
    if bni>0.5
        pulse_value(pulse)=1;
    elseif bni<-0.5
        pulse_value(pulse)=-1;
    else
        pulse_value(pulse)=0;
    end
end

%considering number of pulses the pulse value must be switched
%randint generates a random scalar that is either 0 or 1
for i=2:1:N
    if pulse_value(i)== pulse_value(i-1)
        switch pulse_value(i)
            case 1
                pulse_value(i)=randint-1;
            case 0
                pulse_value(i)=2*randint-1;
            case -1
                pulse_value(i)=randint;
        end
    end
end 

%out = randint(m,n,rg) generates an m-by-n integer matrix. If rg is zero,
%out is a zero matrix. Otherwise, the entries are uniformly distributed and independently chosen from the range
feas=true;
while feas
    length_pulse=randint(1,N-1,[t0,T(1,end)]);
    bni=[t0,sort(length_pulse),T(1,end)];
    feas=false;
    for j=1:1:N
        if bni(1,j+1)-bni(1,j)<tmin % checking pulse width and increasing them if they are less than tmin
            if j==N
                feas=true;
            end
            bni(1,j+1)=bni(1,j)+tmin;
        end
    end
    
end

for i=1:1:N
    
      u(1,bni(1,i):T(1,end))=pulse_value(i);  
end
u=u';
end
