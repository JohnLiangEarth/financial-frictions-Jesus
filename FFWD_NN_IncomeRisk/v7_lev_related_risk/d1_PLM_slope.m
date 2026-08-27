% Calculate the slope of the PLM(B,N)=0 locus as B increases.
%
% Run figures/c2_PLM_NPLM_phase.m first. That script constructs PLM0 as
% an nval_BB-by-2 array:
%   PLM0(:,1) = B
%   PLM0(:,2) = N such that PLM(B,N) = 0.

% Keep the original row locations so that gaps in PLM0 remain identifiable.
PLM0_valid = all(isfinite(PLM0),2);
PLM0_B = PLM0(PLM0_valid,1);
PLM0_N = PLM0(PLM0_valid,2);
PLM0_slope = NaN(size(PLM0_B));

if numel(PLM0_B) < 2
    warning('d1_PLM_slope:InsufficientPoints', ...
        'PLM0 contains fewer than two valid points; its slope cannot be calculated.');
else
    % Calculate dN/dB separately on each contiguous section. This avoids
    % differentiating across B values at which the zero locus is undefined.
    PLM0_valid_rows = find(PLM0_valid);
    new_section = [true; ...
    PLM0_valid_rows(2:end) - PLM0_valid_rows(1:end-1) > 1];
    section_id = cumsum(new_section);

    for i_section = 1:max(section_id)
        in_section = (section_id == i_section);

        if nnz(in_section) >= 2
            PLM0_slope(in_section) = gradient(PLM0_N(in_section), ...
                                              PLM0_B(in_section));
        else
            warning('d1_PLM_slope:SinglePointSection', ...
                'A valid section of PLM0 has only one point; its slope is NaN.');
        end
    end
end

% Collect the zero locus and its slope in one array for inspection or export.
% Columns are B, N on PLM=0, and dN/dB, respectively.
PLM0_slope_results = [PLM0_B, PLM0_N, PLM0_slope];

disp('Slope of the PLM(B,N)=0 locus as B increases:')
disp('             B                  N                 dN/dB')
disp(PLM0_slope_results)

figure;
subplot(2,1,1)
plot(PLM0_B,PLM0_N,'LineWidth',1.5)
xlabel('B','Interpreter','latex')
ylabel('N','Interpreter','latex')
title('$\mathrm{PLM}(B,N)=0$ locus','Interpreter','latex')
grid on

subplot(2,1,2)
plot(PLM0_B,PLM0_slope,'LineWidth',1.5)
yline(0,'k:')
xlabel('B','Interpreter','latex')
ylabel('$dN/dB$','Interpreter','latex')
title('Slope of the zero locus','Interpreter','latex')
grid on
