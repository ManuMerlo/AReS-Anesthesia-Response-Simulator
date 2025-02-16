% Hypovolemia indicates the status of fluid loss, and Normovolemia
% represents no fluid loss
classdef VolumeStatus
    enumeration
        HYPOVOLEMIA(struct('sv', 0.74, 'co', 0.868, 'map', 1.04, 'hr', 1.16))
        NORMOVOLEMIA(struct('sv', 1.0, 'co', 1.0, 'map', 1.0, 'hr', 1.0))
    end
end
