% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4] = solution()
 [consoleout, A1, A2, A3, A4] = evalc('student_solution(0)'); 
end

function [A1, A2, A3, A4] = student_solution(dummy_argument)
    
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
    
    % initial conditions (a) ----------------------------------
    psi0 = cos(X) .* cos(Y) .* cos(Z);
    psi0_f = fftn(psi0);
    psi0_fvec = reshape(psi0_f, n^3, 1);
    tic
    [t, psi_fvec] = ode45(@(t, psi_fvec) gross_pita_rhs(t, psi_fvec, K, spacial, n), tspan, psi0_fvec);

    A1 = real(psi_fvec);
    A2 = imag(psi_fvec);
    size(psi_fvec)
    % initial conditions (b) ---------------------------------
    %{
    psi0 = sin(X) .* sin(Y) .* sin(Z);
    psi0_f = fftn(psi0);
    psi0_fvec = reshape(psi0_f, n^3, 1);
    tic
    [t, psi_fvec] = ode45(@(t, psi_fvec) gross_pita_rhs(t, psi_fvec, K, spacial, n), tspan, psi0_fvec);

    A3 = real(psi_fvec);
    A4 = imag(psi_fvec);
    
    %
    psi_f = reshape(psi_fvec, n, n, n);
    psi = ifftn(psi_f);
    
    p = patch(isosurface(X, Y, Z, psi, -3));
    isonormals(x,y,z,v, p)
    p.FaceColor = 'red';
    p.EdgeColor = 'none';
    daspect([1 1 1])
    view(3)
    camlight; lighting phong
    %}
    
end

% semi-spectral time step: Gross-Pitaevskii system
function psi_fvec = gross_pita_rhs(t, psi_fvec, K, spacial, n)

    % non-linear part (evaluated in normal space)
    psi_f = reshape(psi_fvec, n, n, n);
    psi = ifftn(psi_f);
    
    N = (-psi.* conj(psi) + spacial) .* psi;
    Nf = fftn(N);
    
    % linear part (evaluated in Fourier Space)
    Lf = 0.5 * K .* psi_f;
    
    % formulate psi_t
    psi_new_f = -2i* (Nf + Lf);
    psi_fvec = reshape(psi_new_f, n^3, 1);
        
end