function MList = CreateMoleculeList(numElements, varargin)
%--------------------------------------------------------------------------
% MStruct = CreateMoleculeList(numElements)
% This function creates an array of empty molecule structures. It is used
% to allocate memory for molecule lists.  
%--------------------------------------------------------------------------
% Outputs:
% MList/array of molecule structure: An empty array of molecule structures
%--------------------------------------------------------------------------
% Inputs:
% numElements/integer(1): The number of elements to include in the molecule
%   list
%--------------------------------------------------------------------------
% Variable Inputs:
% 'compact'/boolean (false): A flag which controls whether the molecule
%   list is an array of structures or a structure with array elements
% 'fieldsToLoad'/cell (all fields): Create only the subset of fields provided
%   in this option. 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% October 3, 2012
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;

fieldNames = {'x','y','xc','yc','h','a','w','phi','ax','bg','i','c','density',...
    'frame','length','link','z','zc'};

fieldTypes = {'single','single','single','single','single','single','single',...
    'single','single','single','single','int32','int32','int32','int32',...
    'int32','single','single','single'};

defaultValues = {'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan', 'nan',...
    'nan', 'nan', 'nan', 'nan', 'nan'};

%--------------------------------------------------------------------------
% Default Parameters
%--------------------------------------------------------------------------
compact = true;
fieldsToLoad = fieldNames;

%--------------------------------------------------------------------------
% Parse Required Input
%--------------------------------------------------------------------------
if nargin < 1
    numElements = 0; %Return empty, compact MList
end

%--------------------------------------------------------------------------
% Parse Variable Input
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'compact'
                compact = CheckParameter(parameterValue, 'boolean', parameterName);
            case 'fieldsToLoad'
                fieldsToLoad = CheckParameter(parameterValue, 'cell', parameterName);
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Parse Variable Input
%--------------------------------------------------------------------------
fieldIndsToLoad = find(ismember(fieldNames, fieldsToLoad));

%--------------------------------------------------------------------------
% Create Molecule Structure
%--------------------------------------------------------------------------
if ~compact
    for i=1:length(fieldIndsToLoad)
        MList.(fieldNames{fieldIndsToLoad(i)}) = ...
            eval([fieldTypes{fieldIndsToLoad(i)} '(' defaultValues{fieldIndsToLoad(i)} ')']);
    end
    MList = repmat(MList, [1 numElements]);
else
    for i=1:length(fieldIndsToLoad)
        MList.(fieldNames{fieldIndsToLoad(i)}) = ...
            eval([fieldTypes{fieldIndsToLoad(i)} '(' defaultValues{fieldIndsToLoad(i)} '*ones(' ...
            num2str(numElements) ', 1) )']);
    end
end

