% Compare dN/dB along h(B,N) = 0 across versions v1--v5.
%
% Each version's z_FinalWorkspace.mat contains PLM0, whose first two
% columns are B and N along the zero locus. The derivative is calculated
% separately on each contiguous finite section so that missing parts of
% the locus are not bridged.

clearvars
close all
clc

script_dir = fileparts(mfilename('fullpath'));
model_dir = fileparts(script_dir);

versions = {'v1','v2','v3','v5'};
colors = lines(numel(versions));
results = struct('version',versions, ...
                 'B',cell(size(versions)), ...
                 'N',cell(size(versions)), ...
                 'dNdB',cell(size(versions)));

fig = figure('Name','dN/dB along h(B,N) = 0', ...
             'Color','w');
hold on

for i_version = 1:numel(versions)
    version = versions{i_version};
    workspace_file = fullfile(model_dir,version,'z_FinalWorkspace.mat');

    if ~isfile(workspace_file)
        error('dNdB_h0:WorkspaceNotFound', ...
              'Workspace not found: %s',workspace_file);
    end

    workspace = load(workspace_file,'PLM0');
    if ~isfield(workspace,'PLM0')
        error('dNdB_h0:PLM0NotFound', ...
              'PLM0 is missing from %s.',workspace_file);
    end
    if ~isnumeric(workspace.PLM0) || size(workspace.PLM0,2) < 2
        error('dNdB_h0:InvalidPLM0', ...
              'PLM0 in %s must be a numeric array with at least two columns.', ...
              workspace_file);
    end

    B_h0 = workspace.PLM0(:,1);
    N_h0 = workspace.PLM0(:,2);
    valid = isfinite(B_h0) & isfinite(N_h0);
    dNdB = NaN(size(B_h0));

    valid_rows = find(valid);
    if numel(valid_rows) < 2
        warning('dNdB_h0:InsufficientPoints', ...
                '%s has fewer than two valid PLM0 points.',version);
    else
        % A new section starts after any invalid row. Calculate the slope
        % within sections only, following v5/d1_PLM_slope.m.
        new_section = [true; diff(valid_rows) > 1];
        section_id = cumsum(new_section);

        for i_section = 1:max(section_id)
            section_rows = valid_rows(section_id == i_section);

            if numel(section_rows) < 2
                warning('dNdB_h0:SinglePointSection', ...
                        ['%s contains a one-point finite section; its ', ...
                         'slope is left as NaN.'],version);
                continue
            end

            B_section = B_h0(section_rows);
            if any(diff(B_section) == 0)
                error('dNdB_h0:RepeatedB', ...
                      '%s contains repeated B values within a finite PLM0 section.', ...
                      version);
            end

            dNdB(section_rows) = gradient(N_h0(section_rows),B_section);
        end
    end

    results(i_version).B = B_h0;
    results(i_version).N = N_h0;
    results(i_version).dNdB = dNdB;

    plot(B_h0,dNdB,'LineWidth',1.6, ...
         'Color',colors(i_version,:), ...
         'DisplayName',version);
end

yline(0,'k:','HandleVisibility','off');
xlabel('$B$','Interpreter','latex');
ylabel('$dN/dB$','Interpreter','latex');
title('$dN/dB$ along $h(B,N)=0$','Interpreter','latex');
legend('Location','best','Interpreter','none');
grid on
box on
hold off

savefig(fig,fullfile(script_dir,'fig_1_dNdB_h0.fig'));
exportgraphics(fig,fullfile(script_dir,'fig_1_dNdB_h0.pdf'), ...
               'ContentType','vector');
