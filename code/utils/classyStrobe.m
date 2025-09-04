classdef classyStrobe < handle
    % hacky class to strobe values to PLEXON ephys.
    % key methods:
    %   addValue - adds a vlaue to the valueList, which will be strobed
    %              once the 'strobe' method is called
    %   strobe - when called strobes all values that are in the valueList.
    %   xxxxx
    %   xxxxx
    
    % IMPORTANT:
    % right now, this only works with Plexon, not any other ephys system.
    
    % 20180606 - lnk & jph wrote it
    
    properties
        valueList       = []; % list of values staged to be strobed upon command (either the 'strobeList', or 'strobeNow' methods)
        armedToStrobe   = []; % boolean.    
        vetoList        = []; % list of values that will not be strobed. The only way to override this is by using methods 'flushVetoList' or 'removeFromVetoList'
        strobedList     = []; % a list of values that were strobed for bookkeeping or analysis
        % define here any strobe codes that are internal to the class.
        % This code should not overlap with any code in the 'codes'
        % file (e.g. in pds.initCodes). This is CRITICAL. 
        internalStrobeCodes = struct('isCell', 32123, 'cellLength', 32124);
        
    end
    
    
    %%
    methods
        
        function self = strb(self)
            % constructor function
            
        end
        
        %%
        function self = addValue(self, value)
            % adds a value to the list of values 'valueList'. Any value on
            % this list will get strobed once method 'strobe' is issued.
            
            
            if isnumeric(value)
                % if 'value' is numeric (either scalar or vecotr) its contents
                % will be added to the list:
                self.valueList      = [self.valueList; value];
            
            elseif iscell(value)
                % however, if it is a cell, we loop through and add the
                % contents of each cell to the list. 
                % Each cell is preceded by a triplette of strobes, for 
                % error-checking: 
                % 1- codes.isCell - to indicate the beginning of a cell
                % 2- codes.cellLength - to indicate length of cell
                % 3- the numerical length of the cell. 
                % Thus, user can check for errors when decoding. 
                
                cellOfValues = value; % just for readability
                % for each cell:
                for ii = 1:numel(cellOfValues)
                    cellLength = numel(cellOfValues{ii});
                    
                    % add the triplette to the list:
                    strobeTriplette = [self.internalStrobeCodes.isCell; ...
                                       self.internalStrobeCodes.cellLength; ...
                                       cellLength]; 
                    self.valueList  = [self.valueList; strobeTriplette]; 
                    
                    % and now add the cell's contents to the valueList:
                    self.valueList      = [self.valueList; value{ii}];
                end
                
            elseif islogical(value)
                % is the value happens to be a logical, convert it into
                % int16 and then add to value list:
                self.valueList      = [self.valueList; double(value)];
            else
                error('THIS IS NOT THE STROBE YOU''RE LOOKING FOR')
            end
            self.armedToStrobe  = true;
                
        end
        
        %%
        function self = strobeList(self)
            % strobes all values that are in the list 'valueList'
            
            % if valueList is empty, get outa here
            if isempty(self.valueList)
                return;
            end
            
            % strobe each of the values stored in valueList:
            nValues = numel(self.valueList);
            for iV = 1:nValues
                value = self.valueList(iV);
                strobe(value);
                % and store for bookeeping:
                self.strobedList(end+1) = value;
            end
            
            % clear the list and un-arm the trigger:
            self.valueList      = [];
            self.armedToStrobe  = false;
            
        end
        
        %%
        function self = strobeNow(self, value)
            % strobes the value supplied as input, independantly of values
            % in 'valueList'
            
            strobe(value);
            % and store for bookeeping:
            self.strobedList(end+1) = value;
        end
        
        %%
        function self = addValueOnce(self, value)
            % adds a value to the 'valiueList' such that it may be strobed
            % once the 'strobeList' method is issued, BUT, it also adds the 
            % value to a 'vetoList' such that it doesn't get strobed again.
            %
            % * In order to override this and
            % strobe the value once more, user will have to remove the
            % value from the vetoList with the method
            % 'removeFromVetoList'.
            
            % if value is in 'vetoList', do nothing. Otherwise, add it to 
            % 'valueList' for subsequent strobing (via 'strobe') but then,
            % critically, add the value to vetoList:
            if any(value == self.vetoList)
%                d disp([num2str(value) ' HAS BEEN VEOTED MWU HAHAHA!'])
            else
                self.addValue(value);
                self.vetoList = [self.vetoList; value];
            end
        end
        
        %% 
        function self = flushVetoList(self)
            % clear all values in the vetoList

            self.vetoList = [];
            
        end
        
        %% 
        function self = flushStrobedList(self)
            % clear all values in the strobedList

            self.strobedList = [];
            
        end
        
        
        %% 
        function self = removeFromVetoList(self, value)
            % if the input 'value' exists in the vetoList, this
            % function removes it from the list such that it may be strobed
            % with the strobeOnce method.
            
            if find(value == self.vetoList)
                ptr     = value == self.vetoList;
                self.vetoList(ptr) = [];
            end
               
            
        end
        
    end
end

%%

function [] = strobe(value, verbose)
% the core function that issues the Datapixx commands that strobe the value
% or list of values to the ephys recording maching.
% *** Currently only tailored to Plexon.

if ~exist('verbose', 'var')
    verbose = false;
end

% strobe it like its hot:
Datapixx('SetDoutValues', value, 32767)    % set word in first 15 bits hex2dec('007fff')
Datapixx('RegWr');
Datapixx('SetDoutValues', 2^16, 65536)   % set STRB to 1 (true) hex2dec('010000')
Datapixx('RegWr');
Datapixx('SetDoutValues', 0, 98303)      % reset strobe and all 15 bits to 0. hex2dec('017fff')
Datapixx('RegWr');

% for debugging, display the values strobed:
if verbose
    disp(['boom - ' num2str(value)])
end
end



