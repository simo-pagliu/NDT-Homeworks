# Current plan:
1) Code all known laws and correlations  
2) Program all functions and test them indipendently with guess values for unkown
3) Investigation on gap thickness and cladding thickness optimum (2 variable optimization)
4) Investigation on plenum height requirments (1 variable optimization, depends on the previous!)
5) Possibly proper numerical optimization by restricting the domain from previous investigations 

# Progress Tracking:
To Fix:
- Proper thermal laws for gap and fuel

To Find data on:
- Fuel roughness --> DONE, code it
- Cladding roughness --> DONE, code it
- Fuel Emissivity --> DONE, write it down properly below
- Cladding Emissivity --> DONE, write it down properly below

To Find data on (not now):
- Cladding hardness (only if we have to consider contact between fuel and cladding)
  
To Do (functions):
- Thermal expansion function: Input: Geometry_Data & T_map --> Output: Geometry_Data
- Pu redistribution: does it needs to be calculated by itself?
- Central void formation function: Input: Geometry_Data & T_map --> Output: Geometry_Data
- Stress computations: ???
  
To Do (code):
For each evaluation of the gap & cladding size we have to Loop for:
- T_fuel_out in gap HTC
- Thermal expansion & Void (hot geom)

To Do (considerations):
- What can we neglect?

To Do ("side quests"):
- Reformat of function arguments by "classifing" free variables
- Code optimization
- Code commenting

To Test (code-wise):
- Power Profile

To Validate (physical result wise):
- Coolant Velociy
- Peak Power
- K coolant
- HTC cladding-coolant
- Contact heat transfer is negligible (hp: if they don't touch in the hot geometry then it's good)
- Coolant thermal resistance
- Cladding thermal resistance
- Gap thermal resistance
- Fuel thermal resistance

Validated:
- Sodium Coolant Density, Viscosity, $C_p$   
[IAEA](https://inis.iaea.org/collection/NCLCollectionStore/_Public/14/776/14776927.pdf)
- Pellet roughness (found ourselves) [NEA](https://www.oecd-nea.org/upload/docs/application/pdf/2019-12/6291-mox.pdf)
- Stainless steel roughness (found ourselves) [Engineering Toolbox](https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html): roughness goes from 6 $\mu m$ to 0.8 $\mu m$
Data on Steel:
[TheWorldMaterial Blog](https://www.theworldmaterial.com/type-304-grade-stainless-steel/)

Possible other refs:
[Sodium Handbook](https://www-pub.iaea.org/MTCD/publications/PDF/CRCP_SOD_003web.pdf)

Done (code-wise):
- Hydraulics parameters
- Axial power profile
- Temperature profile radially
- Temperature profile axially
- All axial temperature profiles

Notes for the report:
- When evaluating K we use a fixed value of Pu concentration as the correlation uses the average concentration and doesn't have any local physical meaning
- Radiative heat transfer in the gap is ~ 14% of the total