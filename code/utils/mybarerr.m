function [fo, po] = mybarerr(x, y, yCI, varargin)

% [fo, po] = mybarerr(x,y,yCI,[yi],[colored bars? (default = 0)],[bar "radius"],[x or y bars])

nBoxes  = length(x);        % number of bars/boxes
lusv    = 1.5;              % lighten-up scale value
sucv    = 0.6;             % single-use color value
cvu     = [sucv sucv sucv]; % color value
dcvu    = cvu/lusv;         % dark color value
xy      = 0;                % x or y bars?
inci    = 0;                % two-error intervals?
cf      = 0;                % colored bars?
brm     = [];               % manual bar "radius"

switch nargin
    case 4
        inci    = ~isscalar(varargin{1});
    case 5
        if(~isempty(varargin{1}))
            inci    = ~isscalar(varargin{1});
        end
        if(isscalar(varargin{2}))
            cf      = varargin{2};
        else
            if(~isempty(varargin{2}))
                switch size(varargin{2},1)
                    case 1
                        cvu = varargin{2};
                    otherwise
                        cv  = varargin{2};
                        dcv = cv/lusv;
                end
            end
        end
    case 6
        if(~isempty(varargin{1}))
            inci    = ~isscalar(varargin{1});
        end
        if(isscalar(varargin{2}))
            cf      = varargin{2};
        else
            if(~isempty(varargin{2}))
                switch size(varargin{2},1)
                    case 1
                        cvu = varargin{2};
                    otherwise
                        cv  = varargin{2};
                        dcv = cv/lusv;
                end
            end
        end
        brm     = varargin{3};
    case 7
        if ~isempty(varargin{1})
            inci    = ~isscalar(varargin{1});
        end
        if(isscalar(varargin{2}))
            cf      = varargin{2};
        else
            if(~isempty(varargin{2}))
                switch size(varargin{2},1)
                    case 1
                        cvu = varargin{2};
                    otherwise
                        cv  = varargin{2};
                        dcv = cv/lusv;
                end
            end
        end
        if(~isempty(varargin{3}))
            brm     = varargin{3};
        end
        xy      = varargin{4};
end

% make sure "yCI" rows / columns match "x"
[nCiRows, nCiCols] = size(yCI);
if nCiRows ~= nBoxes
    yCI = yCI';
end

% if there's going to be two intervals drawn (like a mean-ci as well as a
% 95% distribution interval, assign the input argument.
if(inci)
    ymci = varargin{1};
end

% if colored bars are going to be drawn, do some color stuff.
temp = whos;
if(cf)
    switch length(varargin{1})
        case 1
            if(varargin{1}<1)
                pct1 = 1-varargin{1};
                pct2 = 1-pct1;
            else
                pct1 = 0.2;
                pct2 = 0.8;
            end
            cv      = hsv(nBoxes)*pct1 + repmat(pct2,nBoxes,3);
            dcv     = cv/lusv;
        otherwise
            if(isempty(varargin{1}))
                cv = cool(nBoxes);
                cf = 0;
            else
                pct1 = 0.2;
                pct2 = 0.8;
                cv      = repmat(varargin{1}*pct1,nBoxes,1) + ...
                    repmat(pct2,nBoxes,3);
                dcv     = cv/lusv;
            end
    end
elseif(any(strcmp({temp(:).name},'cv')))
    cf = 1;
end

% chose the bar-radius: manual, automatic based on many bars, or
% fixed based on a single bar.
if(length(x)>1 || ~isempty(brm))
    if(isempty(brm))
        br = min([0.4*min(diff(sort(x(:)))), 0.4]);
    else
        br = brm;
    end
else
    br = 0.4;
end

% plotting
if(ishold)
    hflg = 0;
else
    hold on
    hflg = 1;
end
hold on

for i = 1:nBoxes
    if(cf)
        cvu     = cv(i,:);
        dcvu    = dcv(i,:);
    end
    x1 = x(i)-br;
    x2 = x(i)+br;
    
    if(inci)
        try
            h = drawbar(x(i),[yCI(i,1) ymci(i,1) y(i) ymci(i,2) yCI(i,2)],br);
        catch me
            keyboard
        end
        if(cf)
            keyboard
        end
    else
        try
            if ~xy
                Xfill = [x2 x2 x1 x1];
                Yfill = [yCI(i,:) fliplr(yCI(i,:))];
                Xplot = [x1+br/10 x2-br/10];
                Yplot = [y(i) y(i)];
            else
                Yfill = [x2 x2 x1 x1];
                Xfill = [yCI(i,:) fliplr(yCI(i,:))];
                Yplot = [x1+br/10 x2-br/10];
                Xplot = [y(i) y(i)];
            end
            
            fo(i) = fill(Xfill, Yfill,cvu);
            po(i) = plot(Xplot, Yplot, ...
                'LineWidth', 3, 'Color', dcvu);
            set(fo(i), 'EdgeColor', cvu)
        catch me
            
            disp('something went wrong with mybarrerr...')
%             keyboard
        end
    end
end
if(hflg)
    hold off
end