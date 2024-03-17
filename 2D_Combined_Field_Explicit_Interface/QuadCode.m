X = [0.1012865073235; 0.794269853531; 0.1012865073235; 0.4701420641051; 0.4701420641051; 0.0597158717898; 1/3];
Y = [1 1 2 6 4 4 7];
W = [0.1259391805448; 0.1259391805448; 0.1259391805448; 0.1323941527885; 0.1323941527885; 0.1323941527885; 0.225];

sum = zeros(6,6);

Result = zeros(6,6);

for i =1:length(X)

L = [1-X(i)-X(Y(i)) X(i) X(Y(i))];
phi = [L(1)*(2*L(1)-1) L(2)*(2*L(2)-1) L(3)*(2*L(3)-1) 4*L(2)*L(3) 4*L(3)*L(1) 4*L(3)*L(1)];
Result = W(i)*(phi'*phi);
sum = sum+ Result;

end

sum
