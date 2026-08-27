% Plot aggregate consumption C(B,N) across versions v1--v5.
%
% The household budget identity is
%
%   dB/dt = w(B,N) + r(B,N)*B - C(B,N).
%
% Since PLM_finegrid is the perceived law of motion h(B,N) for B,
% aggregate consumption is C(B,N) = w(B,N) + r(B,N)*B - h(B,N).

clearvars
close all
clc

script_dir = fileparts(mfilename('fullpath'));
model_dir = fileparts(script_dir);

versions = {'v1','v2','v3','v4','v5'};
required_variables = {'PLM_finegrid','BB_grid_2D','NN_grid_2D', ...
                      'alpha','Zeta','delta','sigma2'};

results = struct('version',versions, ...
                 'B',cell(size(versions)), ...
                 'N',cell(size(versions)), ...
                 'h',cell(size(versions)), ...
                 'r',cell(size(versions)), ...
                 'w',cell(size(versions)), ...
                 'C',cell(size(versions)));

C_min = Inf;
C_max = -Inf;

for i_version = 1:numel(versions)
    version = versions{i_version};
    workspace_file = fullfile(model_dir,version,'z_FinalWorkspace.mat');

    if ~isfile(workspace_file)
        error('d2_C:WorkspaceNotFound', ...
              'Workspace not found: %s',workspace_file);
    end

    workspace = load(workspace_file,required_variables{:});
    missing_variables = setdiff(required_variables,fieldnames(workspace));
    if ~isempty(missing_variables)
        error('d2_C:VariablesNotFound', ...
              'Missing variables in %s: %s',workspace_file, ...
              strjoin(missing_variables,', '));
    end

    B_grid = workspace.BB_grid_2D;
    N_grid = workspace.NN_grid_2D;
    h_grid = workspace.PLM_finegrid;

    if ~isequal(size(B_grid),size(N_grid),size(h_grid))
        error('d2_C:InconsistentGridSizes', ...
              ['BB_grid_2D, NN_grid_2D, and PLM_finegrid must have ', ...
               'the same size in %s.'],workspace_file);
    end
    if any(~isfinite(B_grid(:))) || any(~isfinite(N_grid(:))) || ...
            any(N_grid(:) <= 0) || any(B_grid(:) + N_grid(:) <= 0)
        error('d2_C:InvalidGrid', ...
              'The B and N grids in %s contain invalid values.',workspace_file);
    end

    capital_grid = B_grid + N_grid;
    r_grid = workspace.alpha*workspace.Zeta* ...
             capital_grid.^(workspace.alpha-1) ...
             - workspace.delta ...
             - workspace.sigma2.*capital_grid./N_grid;
    w_grid = (1-workspace.alpha)*workspace.Zeta* ...
             capital_grid.^workspace.alpha;
    C_grid = w_grid + r_grid.*B_grid - h_grid;

    finite_C = C_grid(isfinite(C_grid));
    if isempty(finite_C)
        error('d2_C:NoFiniteConsumption', ...
              'No finite aggregate-consumption values were calculated for %s.', ...
              version);
    end

    C_min = min(C_min,min(finite_C));
    C_max = max(C_max,max(finite_C));

    results(i_version).B = B_grid;
    results(i_version).N = N_grid;
    results(i_version).h = h_grid;
    results(i_version).r = r_grid;
    results(i_version).w = w_grid;
    results(i_version).C = C_grid;
end

if C_min == C_max
    C_padding = max(1,abs(C_min))*1e-6;
    C_limits = [C_min-C_padding,C_max+C_padding];
else
    C_limits = [C_min,C_max];
end

fig = figure('Name','Aggregate consumption functions', ...
             'Color','w', ...
             'Position',[100 100 1350 760]);
layout = tiledlayout(fig,2,3, ...
                     'TileSpacing','compact', ...
                     'Padding','compact');

for i_version = 1:numel(versions)
    nexttile(layout,i_version);
    surf(results(i_version).B, ...
         results(i_version).N, ...
         results(i_version).C, ...
         'EdgeColor','none');
    xlabel('$B$','Interpreter','latex');
    ylabel('$N$','Interpreter','latex');
    zlabel('$C(B,N)$','Interpreter','latex');
    title(versions{i_version},'Interpreter','none');
    clim(C_limits);
    view(45,30);
    axis tight
    grid on
    box on
    colorbar
end

nexttile(layout,6);
axis off
text(0.5,0.55,'$C(B,N)=w(B,N)+r(B,N)B-h(B,N)$', ...
     'Interpreter','latex', ...
     'HorizontalAlignment','center', ...
     'FontSize',14);

title(layout,'Aggregate consumption across versions', ...
      'Interpreter','latex');

savefig(fig,fullfile(script_dir,'fig_d2_C.fig'));
exportgraphics(fig,fullfile(script_dir,'fig_d2_C.pdf'), ...
               'ContentType','vector');
