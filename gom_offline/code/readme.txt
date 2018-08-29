% File Names and Purposes
SIZE.h - Define the domain dimension (x,y,z) (sNx, sNy, Nr), also this is where you can define the subgrid dimensions for
         parrallelization.

PTRACERS_SIZE.h - Define the number of Tracers that will be used in the simulation (part of the pTracer package)

ptracers_apply_forcing.F - Commented out lines so that I can apply the correction in gchem_forcing_sep.f 

do_the_model_io.F - Edited so that the KPP variables calculated by the model are output into the state file when using the offline
package


