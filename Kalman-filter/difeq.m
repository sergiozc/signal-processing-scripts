function [dxp] = difeq(t, xp)
% Función que construye las ecuaciones diferenciales

global GM;
global Q; 

x = xp(1:4); % Las 4 primeras son las variables de estado
r = sqrt(x(1)^2 + x(2)^2); % Radio
dx(1) = x(3); % Derivada de la primera variable de estado es la tercera variable de estado y etc
dx(2) = x(4);
dx(3) = (-GM / r^3) * x(1);
dx(4) = (-GM / r^3) * x(2);

P = reshape(xp(5:20),[4,4]); % Hay que reorganizarlas en la matriz P con un reshape
% (Comprobar que la matriz P sea diagonal)

% Construimos la matriz de transición F
elemento31 = GM * (3*x(1)^2 - r^2) / r^5;
elemento32 = GM * (3*x(1) * y(1)) / r^5;
elemento41 = GM * (3*x(1) * y(1)) / r^5;
elemento42 = GM * (3*x(1)^2 - r^2) / r^5;
F = [0 0 1 0; 0 0 0 1; elemento31 elemento32 0 0; elemento41 elemento42 0 0];

% Derivada de la matriz P
dP = F * P + P * F' + Q; % Habiendo definido Q como global

dp_array = reshape(dP, [1, 16]);

% dxp tiene que ser un vector (formato de la función ode45)
dxp = [dx(1) dx(2) dx(4) dp_array];
end

