"""
Module that allows to perform ground heat exchanger (GHE) sizing to obtain the borehole length
required to cover the heat loads of a building, while keeping the circulating fluid within certain
operating temperature limits. This module only contains the sizing equations, and rely on other
modules or codes to provide heat loads and GHE simulations models (such as QBuildings.jl to provide
a way to compute heat load profiles and GHEModels.jl to provide analytical models to simulate GHE).

Features included as for now are the alternative sizing equation in AlternativeSizing.jl that
is based on g-functions.

Author: Gabriel Dion (dion.gabriel100@gmail.com)
Date: 2025-04
Julia version: 1.11.3
"""

module GHESizings
using Revise

includet("AlternativeSizing.jl")
includet("QAnalysis.jl")

# From AlternativeSizing
export alternative_sizing,
    alternative_sizing_L2,
    alternative_sizing_L3,
    alternative_sizing_L4

# From QAnalysis
export QLoads,
    Q_analysis,
    Q_L3_resample
    building_to_ground_loads

end
