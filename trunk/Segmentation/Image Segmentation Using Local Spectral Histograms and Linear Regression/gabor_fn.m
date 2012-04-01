function gabout = gabor_fn(S,theta)
% S: std
% theta: orientation

phs=0;
gamma=1;
wins=ceil(S*2);
f=1/(wins+.5);
G=zeros(2*wins+1);

for x = -wins:wins
    for y = -wins:wins
        xPrime = x * cos(theta) + y * sin(theta);
        yPrime = y * cos(theta) - x * sin(theta);
        G(wins+x+1,wins+y+1) =1/(2*pi*S*S)*exp(-.5*((xPrime)^2+(yPrime*gamma)^2)/S^2)*cos(2*pi*f*xPrime+phs);
    end
end

gabout=G;