# Continuous-Piecewise-Quadratic-FEM
A MATLAB implementation to compute and visualize the convergence of the finite element method (FEM) with piecewise quadratic elements for the dirichlet boundary‑value problem: -u'' = pi sin(pi x) in (0,1)

# Description 
This repository contains a MATLAB script FEMPiecewiseQuadraticErrConvMain that solves the above BVP using P2 (quadratic) elements on uniform meshes of sizes. It computes the L²-norm and energy-norm errors via 3-point Gauss–Legendre quadrature, plots error convergence, and displays observed rates

# Code Overview 
The main script FEMPiecewiseQuadraticErrConvMain:
- Defines Gauss–Legendre nodes and weights for 3-point quadrature.
- Loops over mesh refinements h = 2^{-j}, j= 1,..,J_max. Builds node and midpoint coordinates for quadratic DOFs.
- Assembles global stiffness matrix and load vector via helper functions GStiffQuad and GLoadQuad.
- Applies zero Dirichlet boundary conditions, solves for interior DOFs.

Computes error norms on each element:
- Energy norm: integrates (u'-u'_h)^2 with quadrature.
- L² norm: integrates (u-u_h)^2 with quadrature.

Produces a log–log plot of both errors vs. h alongside reference slopes O(h^3) (L²) and O(h^2) (energy). Prints observed convergence rates for each norm. Finally, plots the FEM solution vs. exact solution on the finest mesh.

# Main Function 
```matlab
function FEMPiecewiseQuadraticErrConvMain()
    clear; clc;

    Jmax = 10;
    H    = 2.^-(1:Jmax);
    errL2 = zeros(Jmax,1);
    errE  = zeros(Jmax,1);

    % 3‑point Gauss‑Legendre on [-1,1]
    xi = [-sqrt(3/5); 0; sqrt(3/5)];
    wi = [5/9; 8/9; 5/9];

    for j = 1:Jmax
        h      = H(j);
        x      = (0:h:1)';
        N      = length(x)-1;

        % Build DOF coords: nodes and midpoints
        coords        = zeros(2*N+1,1);
        coords(1:2:end) = x;
        coords(2:2:end) = (x(1:end-1)+x(2:end))/2;

        % Assemble stiffness and load
        A = GStiffQuad(coords);
        F = GLoadQuad(coords);

        % Solve with Dirichlet BCs
        u      = zeros(size(coords));
        free   = 2:length(coords)-1;
        u(free) = A(free,free) \ F(free);

        % Compute energy norm
        errE(j)  = computeErrEnergy(coords,u,xi,wi);
        % Compute L2 norm
        errL2(j) = computeErrL2(coords,u,xi,wi);
    end

    % Plot convergence
    figure; hold on; box on;
    loglog(H, errL2, 'cy--','LineWidth',1.2);
    loglog(H, errE,  'r--','LineWidth',1.2);
    loglog(H, H.^3,  '-','LineWidth',1.2);
    loglog(H, H.^2,  '-','LineWidth',1.2);
    legend('L^2 norm','Energy norm','O(h^3)','O(h^2)','Location','SouthEast');
    xlabel('h'); ylabel('Error');
    title('P2 FEM Convergence for u(x)=sin(\pi x)');
    grid on;

    % Print rates
    displayQuadRates(H, errL2, errE);

    % Plot approximate vs exact on finest mesh
    figure;
    plot(coords, u, 'b--','LineWidth',1.2); hold on;
    plot(coords, sin(pi*coords), 'r-','LineWidth',1.2);
    legend('u_h (P2)','u_{exact}','Location','Best');
    xlabel('x'); ylabel('u(x)');
    title(sprintf('P2 FEM vs exact, h=%.4f', H(end)));
end
```

## Helper Functions
```matlab
# Stiffness Matrix
function A = GStiffQuad(coords)
    M = length(coords);
    N = (M-1)/2;
    A = zeros(M);
    for i = 1:N
        idx = [2*(i-1)+1,2*(i-1)+2,2*i+1];
        a   = coords(idx(1));
        b   = coords(idx(3));
        h   = b-a;
        Ae  = (1/h)*[7/3,-8/3,1/3; -8/3,16/3,-8/3; 1/3,-8/3,7/3];
        A(idx,idx) = A(idx,idx) + Ae;
    end
end

# Load Vector
function F = GLoadQuad(coords)
    M = length(coords);
    N = (M-1)/2;
    F = zeros(M,1);
    for i = 1:N
        idx = [2*(i-1)+1,2*(i-1)+2,2*i+1];
        a   = coords(idx(1));
        m   = coords(idx(2));
        b   = coords(idx(3));
        h   = b-a;
        fa  = pi^2*sin(pi*a);
        fm  = pi^2*sin(pi*m);
        fb  = pi^2*sin(pi*b);
        F(idx) = F(idx) + (h/6)*[fa;4*fm;fb];
    end
end

```


# Usage 
```matlab
% In MATLAB
FEMPiecewiseQuadraticErrConvMain;
```

## License
This project is licensed under the MIT License - see the LICENSE file for details.
```





