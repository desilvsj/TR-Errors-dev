# BAM Tag Legend

- PH:i  Phase indicator
        1 = Phase-1 exact match
        2 = Phase-2 parasail match
        3 = Unmapped after refinement (fallback failed)

- BP:i  Phi (offset into double-consensus used to crop single)

- CL:i  Intended single-consensus length (d)

- MT:i  True matches from Parasail traceback (equals count of identical bases)

- NM:i  Edit distance (mismatches + insertions + deletions)

- CH:i  Chimera heuristic flag
        0 = within threshold, 1 = likely chimera
        Defined as (CL - MT)/CL > 0.05

- RC:i  Original orientation in input BAM
        1 = original read was reverse; 0 = forward
