% Matthew Mangione
% AMATH 481
% Bose-Einstein condensate final project


% INPUT
tspan = [0:.5:4];
n = 16;
L = pi;
A = [-1 -1 -1];
B = -A;

% VAR SETUP
x = linspace(-L, L, n+1);
[X, Y, Z] = meshgrid(x(1:n), x(1:n), x(1:n));

% spacial linear multiplier
spacial = (A(1)*sin(X).^2 + B(1)) .* (A(2)*sin(Y).^2 + B(2)) .* (A(3)*sin(Z).^2 + B(3));

% fourier laplace multiplier
k = (pi/L) * [0:(n/2 -1) (-n/2):-1]; 
k(1) = 10^-6;
[KX, KY, KZ] = meshgrid(k, k, k);
K = -(KX.^2 + KY.^2 + KZ.^2);

% initial conditions (A) 
psi0 = cos(X) .* cos(Y) .* cos(Z);

% initial conditions (B)
%psi0 = sin(X) .* sin(Y) .* sin(Z);

psi0_f = fftn(psi0);
psi0_fvec = reshape(psi0_f, n^3, 1);

% solve in Fourier space for t = [0: .5 : 4]
[t, psi_fvec] = ode45(@(t, psi_fvec) gross_pita_rhs(t, ...
                psi_fvec, K, spacial, n), tspan, psi0_fvec);

% revert back to proper shape and time-domain
psi_f = reshape(psi_fvec(length(tspan),:), n, n, n);
psi = ifftn(psi_f);

% particle density
p = psi .* psi;

%}
for i= 1:length(tspan)
    
    % reshape vector, square it.
    psi_f = reshape(psi_fvec(i,:), n, n, n);
    psi = real(ifftn(psi_f)).^2;
    clf
    
    isosurface(X, Y, Z, psi, -2)

    % edge lines for clarity
    line([pi,pi], [-pi,pi], [-pi, -pi], 'LineWidth', 1, 'Color', 'k');
    line([-pi,pi], [pi,pi], [-pi, -pi], 'LineWidth', 1, 'Color', 'k');
    line([-pi,-pi], [-pi,-pi], [-pi, pi], 'LineWidth', 1, 'Color', 'k');
    line([pi,pi], [-pi,-pi], [-pi, pi], 'LineWidth', 1, 'Color', 'k');
    line([pi,pi], [pi, pi], [-pi, pi], 'LineWidth', 1, 'Color', 'k');
    line([-pi,pi], [pi, pi], [pi, pi], 'LineWidth', 1, 'Color', 'k');
    line([-pi,pi], [-pi,-pi], [pi, pi], 'LineWidth', 1, 'Color', 'k');
    line([pi, pi], [-pi,pi], [pi, pi], 'LineWidth', 1, 'Color', 'k');
    line([-pi,-pi], [-pi,pi], [pi, pi], 'LineWidth', 1, 'Color', 'k');
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    pause(.1)
end


% semi-spectral time step: Gross-Pitaevskii system
function psi_fvec = gross_pita_rhs(~, psi_fvec, K, spacial, n)

    % non-linear part (evaluated in normal space)
    psi_f = reshape(psi_fvec, n, n, n);
    psi = ifftn(psi_f);
    
    N = (-psi.* conj(psi) + spacial) .* psi;
    Nf = fftn(N);
    
    % linear part (evaluated in Fourier Space)
    Lf = 0.5 * K .* psi_f;
    
    % formulate psi_t
    psi_new_f = 1i* (Nf + Lf);
    psi_fvec = reshape(psi_new_f, n^3, 1);
        
end
