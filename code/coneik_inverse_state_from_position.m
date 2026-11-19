function [q,qdot,qddot,info] = coneik_inverse_state_from_position(gamma,gamma_dot,gamma_ddot,options,policy,model,candidates,position_info)
%CONEIK_INVERSE_STATE_FROM_POSITION Internal shared post-IK motion evaluator.
% Preconditions: public adapters already validated finite double 2-vectors,
% normalized legacy options/policy and resolved any explicit model. MODEL is
% [] only for the implicit default. CANDIDATES and POSITION_INFO must be the
% unchanged two base branches and metadata for this gamma/model, produced by
% inverse_motion_analysis within the current call (not an external cache).
% This is an implementation helper, not a supported public/cache interface.
% PREVIOUS_Q selects/reconstructs a representative using the legacy rule; it
% is deliberately not treated as an already selected position. No position
% IK or cross-call memoization occurs here. Rank-one motion recovery delegates
% to the checked polynomial interval helper without changing task tolerances.

    [q,qdot,qddot,info]=coneik_inverse_state_resolved(...
        gamma,gamma_dot,gamma_ddot,options,policy,model,candidates,position_info);
end
