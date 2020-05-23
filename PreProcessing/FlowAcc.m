function [ DTMF , R , T , AccMatrix ] = FlowAcc( DTM , cellsize )
%FlowAcc Compute flow accomulation for a DTM
%   Inputs:
%   DTM - digital eleveation model [m]
%   cellsize - cellsize [m]
%   Outputs:
%   DTMF - DTM with filled sinks
%   R - the pixel flow direction [rad]. Pixel flow direction is measured counter
%   clockwise from the east-pointing horizontal axis.
%   T - flow matrix
%   AccMatrix - accomulation matrix
%   Reference: Tarboton, "A new method for the determination of flow
%   directions and upslope areas in grid digital elevation models," Water
%   Resources Research, vol. 33, no. 2, pages 309-319, February 1997.
%% Main
DTMF=FillSinks(DTM);
R = demFlow(DTMF,cellsize,cellsize);
T = FlowMatrix( DTMF, R , cellsize , cellsize );
AccMatrix = upslopeArea(DTMF, T);
%% Nested functions
    function [r, s] = facet_flow(e0, e1, e2, d1, d2)
        if nargin < 4
            d1 = 1;
            d2 = 1;
        end
        s1 = (e0 - e1) / d1;           % eqn (1)
        s2 = (e1 - e2) / d2;           % eqn (2)
        r = atan2(s2, s1);             % eqn (3)
        s = hypot(s1, s2);             % eqn (3)
        too_far_south = r < 0;
        r(too_far_south) = 0;
        s(too_far_south) = s1(too_far_south);
        diagonal_angle    = atan2(d2, d1);
        diagonal_distance = hypot(d1, d2);
        too_far_north = r > diagonal_angle;
        r(too_far_north) = diagonal_angle;
        s(too_far_north) = (e0(too_far_north) - e2(too_far_north)) / ...
            diagonal_distance;
    end
    function mask = border_nans(E)
        nan_mask = isnan(E);
        mask = nan_mask & ~imclearborder(nan_mask);
    end
    function DTMF=FillSinks(DTM)
        DTM(isnan(DTM))=-9999;
        DTMF = imfill(DTM,8,'holes');
        DTMF(DTMF<-1000)=NaN;
    end
    function [ R , S ] = demFlow(E, d1, d2)
        Ep = padarray(E, [1 1], 'replicate', 'both');
        [M,N] = size(Ep);
        [i, j] = ndgrid(2:M-1, 2:N-1);
        [R, S] = PixelFlow(Ep, i, j, d1, d2);
    end
    function [R, S] = PixelFlow(E, i, j, d1, d2)
        if nargin < 4
            d1 = 1;
            d2 = 1;
        end
        [M, ~] = size(E);
        highest_value = max(E(:));
        bump = min(1, eps(highest_value));
        E(border_nans(E)) = highest_value + bump;
        e0_idx = (j - 1)*M + i;
        e1_row_offsets = [0 -1 -1  0  0  1  1  0];
        e1_col_offsets = [1  0  0 -1 -1  0  0  1];
        e2_row_offsets = [-1 -1 -1 -1  1  1  1  1];
        e2_col_offsets = [ 1  1 -1 -1 -1 -1  1  1];
        e1_linear_offsets = e1_col_offsets*M + e1_row_offsets;
        e2_linear_offsets = e2_col_offsets*M + e2_row_offsets;
        E0 = E(e0_idx);
        E1 = E(e0_idx + e1_linear_offsets(1));
        E2 = E(e0_idx + e2_linear_offsets(1));
        [R, S] = facet_flow(E0, E1, E2, d1, d2);
        ac = [0  1  1  2  2  3  3  4];
        af = [1 -1  1 -1  1 -1  1 -1];
        positive_S = S > 0;
        R(positive_S) = (af(1) * R(positive_S)) + (ac(1) * pi / 2);  % Equation (6)
        R(~positive_S) = NaN;
        for k = 2:8
            E1 = E(e0_idx + e1_linear_offsets(k));
            E2 = E(e0_idx + e2_linear_offsets(k));
            [Rk, Sk] = facet_flow(E0, E1, E2, d1, d2);
            new_R = (Sk > S) & (Sk > 0);
            S = max(S, Sk);
            R(new_R) = (af(k) * Rk(new_R)) + (ac(k) * pi / 2);
        end
    end
    function T = FlowMatrix(E, R, d1, d2)
        [M, N] = size(R);
        pixel_labels = reshape(1:numel(R), M, N);
        have_neighbor.north.rows = 2:M;
        have_neighbor.north.cols = 1:N;
        have_neighbor.northeast.rows = 2:M;
        have_neighbor.northeast.cols = 1:N-1;
        have_neighbor.east.rows = 1:M;
        have_neighbor.east.cols = 1:N-1;
        have_neighbor.southeast.rows = 1:M-1;
        have_neighbor.southeast.cols = 1:N-1;
        have_neighbor.south.rows = 1:M-1;
        have_neighbor.south.cols = 1:N;
        have_neighbor.southwest.rows = 1:M-1;
        have_neighbor.southwest.cols = 2:N;
        have_neighbor.west.rows = 1:M;
        have_neighbor.west.cols = 2:N;
        have_neighbor.northwest.rows = 2:M;
        have_neighbor.northwest.cols = 2:N;
        offset.north     = -1;
        offset.northeast = M - 1;
        offset.east      = M;
        offset.southeast = M + 1;
        offset.south     = 1;
        offset.southwest = -M + 1;
        offset.west      = -M;
        offset.northwest = -M - 1;
        directions = {'north', 'northeast', 'east', 'southeast', ...
            'south', 'southwest', 'west', 'northwest'};
        ii = [];
        jj = [];
        vv = [];
        for k = 1:numel(directions)
            direction = directions{k};
            i = pixel_labels(have_neighbor.(direction).rows, ...
                have_neighbor.(direction).cols);
            j = i + offset.(direction);
            p = -directional_weight(direction, R(j), d1, d2);
            non_zero_weight = p ~= 0;
            ii = [ii; i(non_zero_weight)]; %#ok<AGROW>
            jj = [jj; j(non_zero_weight)]; %#ok<AGROW>
            vv = [vv; p(non_zero_weight)]; %#ok<AGROW>
        end
        border_mask = border_nans(E);
        border_labels = pixel_labels(border_mask);
        delete_mask = ismember(ii, border_labels) | ...
            ismember(jj, border_labels);
        ii(delete_mask) = [];
        jj(delete_mask) = [];
        vv(delete_mask) = [];
        ii = [ii; pixel_labels(:)];
        jj = [jj; pixel_labels(:)];
        vv = [vv; ones(numel(R), 1)];
        [ip, jp, vp] = plateau_flow_weights(E, R, border_mask);
        ii = [ii; ip];
        jj = [jj; jp];
        vv = [vv; vp];
        T = sparse(ii, jj, double(vv));
    end
    function [ip, jp, vp] = plateau_flow_weights(E, R, border_mask)
        S = size(R);
        E(border_mask) = Inf;
        [nan_list_r, nan_list_c] = find(isnan(R) & ~imregionalmin(E) & ...
            ~border_mask);
        done = false;
        num_nans = numel(nan_list_r);
        ip = zeros(8*num_nans, 1);
        jp = zeros(8*num_nans, 1);
        vp = zeros(8*num_nans, 1);
        total_count = 0;
        while ~done
            done = true;
            delete_from_list = false(numel(nan_list_r), 1);
            rr = zeros(numel(nan_list_r), 1);
            cc = zeros(numel(nan_list_c), 1);
            ww = zeros(8,1);
            zz = zeros(8,1);
            nan_count = 0;
            for k = 1:numel(nan_list_r)
                r = nan_list_r(k);
                c = nan_list_c(k);
                neighbor_count = 0;
                for w = max(1, r-1) : min(S(1), r+1)
                    for z = max(1, c-1) : min(S(2), c+1)
                        if ~isnan(R(w,z)) && (E(r,c) == E(w,z))
                            neighbor_count = neighbor_count + 1;
                            ww(neighbor_count) = w;
                            zz(neighbor_count) = z;
                        end
                    end
                end
                if neighbor_count > 0
                    done = false;
                    nan_count = nan_count + 1;
                    rr(nan_count) = r;
                    cc(nan_count) = c;
                    i = (zz(1:neighbor_count) - 1)*S(1) + ww(1:neighbor_count);
                    j = ((c - 1)*S(1) + r) * ones(neighbor_count, 1);
                    v = -ones(neighbor_count, 1) / neighbor_count;
                    slots = total_count + (1:neighbor_count);
                    ip(slots) = i;
                    jp(slots) = j;
                    vp(slots) = v;
                    total_count = total_count + neighbor_count;
                    delete_from_list(k) = true;
                end
            end
            rr = rr(1:nan_count);
            cc = cc(1:nan_count);
            if ~done
                nan_list_r(delete_from_list) = [];
                nan_list_c(delete_from_list) = [];
                R((cc - 1)*S(1) + rr) = 0;
            end
        end
        ip = ip(1:total_count);
        jp = jp(1:total_count);
        vp = vp(1:total_count);
    end
    function w = directional_weight(direction, R, d1, d2)
        theta_d = atan2(d2, d1);
        angles = [0, theta_d, pi/2, pi-theta_d, pi, -pi+theta_d, -pi/2, -theta_d];
        direction_indices = struct('east', 1, 'northeast', 2, 'north', 3, ...
            'northwest', 4, 'west', 5, 'southwest', 6, ...
            'south', 7, 'southeast', 8);
        inward_direction_index = mod(direction_indices.(direction) + 3, 8) + 1;
        inward_direction_index_plus_1 = mod(inward_direction_index, 8) + 1;
        inward_direction_index_minus_1 = mod(inward_direction_index - 2, 8) + 1;
        inward_angle = angles(inward_direction_index);
        next_inward_angle = angles(inward_direction_index_plus_1);
        prev_inward_angle = angles(inward_direction_index_minus_1);
        interp_table_x = [-angular_difference(inward_angle, prev_inward_angle), ...
            0, ...
            angular_difference(next_inward_angle, inward_angle)];
        interp_table_y = [0 1 0];
        w = interp1(interp_table_x, interp_table_y, angular_difference(R, inward_angle));
        w(isnan(w)) = 0;
    end
    function d = angular_difference(theta1, theta2)
        d = mod(theta1 - theta2 + pi, 2*pi) - pi;
    end
    function A = upslopeArea(E, T)
        rhs = ones(numel(E), 1);
        mask = border_nans(E);
        rhs(mask(:)) = 0;
        A = T \ rhs;
        A = reshape(A, size(E));
    end
end