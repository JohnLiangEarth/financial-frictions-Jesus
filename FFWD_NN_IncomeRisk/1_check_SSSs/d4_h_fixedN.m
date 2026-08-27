% Plot h(B,N) against B at a fixed value of N for versions v1--v5.
%
% PLM_finegrid is the converged, relaxed PLM evaluated on the fine
% (BB_grid,NN_grid) grid. Since the PLM is the perceived law of motion
% for B, PLM_finegrid(iB,iN) = h(BB_grid(iB),NN_grid(iN)).

clearvars
close all
clc

script_dir = fileparts(mfilename('fullpath'));
model_dir = fileparts(script_dir);

N_fixed = 2.0; % Common fixed N used for every version.
versions = {'v1','v2','v3','v4','v5'};
colors = lines(numel(versions));
results = struct('version',versions, ...
                 'B',cell(size(versions)), ...
                 'h',cell(size(versions)));

fig = figure('Name','h(B,N) at fixed N','Color','w');
hold on

for i_version = 1:numel(versions)
    version = versions{i_version};
    workspace_file = fullfile(model_dir,version,'z_FinalWorkspace.mat');

    if ~isfile(workspace_file)
        warning('hFixedN:WorkspaceNotFound', ...
                'Skipping %s because its workspace was not found: %s', ...
                version,workspace_file);
        continue
    end

    workspace = load(workspace_file, ...
                     'BB_grid','NN_grid','PLM_finegrid');
    required_fields = {'BB_grid','NN_grid','PLM_finegrid'};
    missing_fields = required_fields(~isfield(workspace,required_fields));
    if ~isempty(missing_fields)
        warning('hFixedN:MissingVariables', ...
                'Skipping %s because %s is missing from its workspace.', ...
                version,strjoin(missing_fields,', '));
        continue
    end

    B_grid_fine = workspace.BB_grid(:);
    N_grid_fine = workspace.NN_grid(:);
    h_grid = workspace.PLM_finegrid;

    expected_size = [numel(B_grid_fine),numel(N_grid_fine)];
    if ~isnumeric(h_grid) || ~isequal(size(h_grid),expected_size)
        warning('hFixedN:InvalidPLMGrid', ...
                ['Skipping %s because PLM_finegrid has size [%s], ', ...
                 'but [%d %d] was expected.'], ...
                version,num2str(size(h_grid)),expected_size(1),expected_size(2));
        continue
    end

    if N_fixed < min(N_grid_fine) || N_fixed > max(N_grid_fine)
        warning('hFixedN:NOutsideGrid', ...
                ['Skipping %s because N_fixed = %.4f is outside its ', ...
                 'N grid [%.4f, %.4f].'], ...
                version,N_fixed,min(N_grid_fine),max(N_grid_fine));
        continue
    end

    % Interpolate in the N dimension for every value of B. Transposing
    % h_grid makes each B slice a separate column for interp1.
    h_fixedN = interp1(N_grid_fine,h_grid.',N_fixed,'linear').';

    results(i_version).B = B_grid_fine;
    results(i_version).h = h_fixedN;

    plot(B_grid_fine,h_fixedN,'LineWidth',1.6, ...
         'Color',colors(i_version,:), ...
         'DisplayName',version);
end

yline(0,'k:','HandleVisibility','off');
xlabel('$B$','Interpreter','latex');
ylabel('$h(B,N_{\mathrm{fixed}})$','Interpreter','latex');
title(sprintf('$h(B,N)$ at $N_{\mathrm{fixed}}=%.2f$',N_fixed), ...
      'Interpreter','latex');
legend('Location','best','Interpreter','none');
grid on
box on
hold off

savefig(fig,fullfile(script_dir,'fig_d4_h_fixedN.fig'));
exportgraphics(fig,fullfile(script_dir,'fig_d4_h_fixedN.pdf'), ...
               'ContentType','vector');
