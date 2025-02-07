# DONE
I will update this readme asap below old version.

---


## Current plan:
1) Code all known laws and correlations  
2) Program all functions and test them indipendently with guess values for unkown
3) Investigation on gap thickness and  plenum height requirments (2 variable optimization)

# Progress Tracking:
To Fix:
- 

To Find data on:
- Cladding hardness (only if we have to consider contact between fuel and cladding)
  
To Do (Physical Functions):
- Fission Gas Production
- Pu redistribution
- Stress computations: ???
  
To Do (Iterative Implementation):
- Fission Gas Production
- Void Swelling

To Do (Considerations):
- What can we neglect?

To Do ("side quests"):
- Code optimization
- Code commenting

To Test (code-wise):
- 

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
- Radial profile
- Fuel Restructuring
- Thermal Expansion
- Void Formation
- Void Swelling

Notes for the report:
- When evaluating K we use a fixed value of Pu concentration as the correlation uses the average concentration and doesn't have any local physical meaning
- Radiative heat transfer in the gap is ~ 14% of the total so we preferred not to neglect it
- For the restructuring densities we got data from the handouts so we don't have any other source to confirm that
- For all the non-validated quantities we simply cross-checked reults with handouts and other groups to see if the orders of magnitude where reasonable (which means no external source for those either)
