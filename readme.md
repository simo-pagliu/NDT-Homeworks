# Current plan:
1) Code all known laws and correlations  
2) Program all functions and test them indipendently with guess values for unkown
3) Investigation on gap thickness and cladding thickness optimum (2 variable optimization)
4) Investigation on plenum height requirments (1 variable optimization, depends on the previous!)
5) Possibly proper numerical optimization by restricting the domain from previous investigations 

# Progress Tracking:
To Fix:
- Temperature axially coolant @ h=0 == t_initial_coolant (to be forced)
- Temperature profile radially in the gap (?) from handouts deltaT should be much larger

To Find data on:
- Fuel roughness
- Cladding roughness
- Fuel Emissivity
- Cladding Emissivity

To Find data on (not now):
- Cladding hardness (only if we have to consider contact between fuel and cladding)
  
To Do:
- All axial temperature profiles
- What can we neglect?
- Thermal Expansion
- Temperature profile in hot geometry (with proper neglectionsz)
- Stress calculations (?)
- Pu redistribution

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

Data we had to get by ourselves:
- Fission Cross section ENDF/B-VI.8
- Energy per Fission

Possible other refs:
[Sodium Handbook](https://www-pub.iaea.org/MTCD/publications/PDF/CRCP_SOD_003web.pdf)

Done (code-wise):
- Hydraulics parameters
- Axial power profile
- Temperature profile radially
- Temperature profile axially (ongoing)