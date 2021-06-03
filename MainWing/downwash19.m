function[X,Y,thetaij,deltagamma,w,downwash] = downwash19(x,y,gamma,theta,gammam,npartition,nd)

%winglets‚É‚æ‚é—U“±R—Í‚ÌŒ¸­—¦
lets = 0;

%À•W‚ğ—¼—ƒ‚Åİ’è
X = [-fliplr(x), 0, x];
Y = [fliplr(y), 0, y];

%‚ ‚é2‚Â‚Ì•ªŠ„“_‚ğŒ‹‚Ôü‚É‚’¼‚Èü‚Ìx²‚Æ‚ÌŠp“x(0 <= thetaij < pi)
thetaij = zeros(2*npartition+1,2*npartition+1);

for i = 1:2*npartition+1
     for j = 1:2*npartition+1
          if i == j 
              
              thetaij(i,j) = 0;
              
          elseif i + j == 2*npartition+2
              
              thetaij(i,j) = pi/2;
          
          else
              
              thetaij(i,j) = atan(-(X(1,i) - X(1,j))/(Y(1,i) - Y(1,j)));
          end    
     end
end

for i = 1:2*npartition+1
     for j = 1:2*npartition+1
         if thetaij(i,j) < 0
             thetaij(i,j) = thetaij(i,j) + pi;
         end
     end
end


%w(i,j):•ªŠ„“_j‚ª•ªŠ„“_i‚É—U“±‚·‚é—U“±‘¬“x
%ãŒü‚«—U“±‘¬“x‚ğ³

w = zeros(2*npartition+1,2*npartition+1);

%”n’û‰Q
deltagamma = zeros(1,npartition);

for i = 1:npartition
    if i == 1
        deltagamma(1,1) = gamma(1,1) - gammam;%•‰
        
    elseif i == npartition
        
        deltagamma(1,npartition) = -gamma(1,npartition-1)*(100 - lets)/100;%•‰
    
    else
        deltagamma(1,i) = gamma(1,i) - gamma(1,i-1);%•‰
    end
end

deltagamma = [fliplr(deltagamma) , 0 , deltagamma] ;

for i = npartition+2:2*npartition+1
    for j = 1:2*npartition+1
        
        if j <= npartition
             w(i,j) = deltagamma(1,j)*nd/(4*pi*sqrt((X(1,i)-X(1,j))^2 + (Y(1,i)-Y(1,j))^2));
        
        elseif j >= npartition+1 && j <= i-1
             w(i,j) = -deltagamma(1,j)*nd/(4*pi*sqrt((X(1,i)-X(1,j))^2 + (Y(1,i)-Y(1,j))^2));
        
        elseif j == i
             w(i,j) = 0;
             
        else 
             w(i,j) = deltagamma(1,j)*nd/(4*pi*sqrt((X(1,i)-X(1,j))^2 + (Y(1,i)-Y(1,j))^2));
        end
    end
end

we = fliplr(flipud(w));

for j = 1:2*npartition+1
    if j == npartition+1
        w(npartition+1,j) = 0;
    else
        w(npartition+1,j) = deltagamma(1,j)*nd/(4*pi*sqrt((X(1,j))^2 + (Y(1,j))^2));
    end
end

w = we + w;

%—U“±‘¬“x
downwash = zeros(1,npartition);
    
for i = 1:npartition
    
    downwash(1,i) = sum((cos(thetaij(i+npartition+1,:)-theta(1,i)-pi/2)).*w(i+npartition+1,:));
end


    
    
    
    