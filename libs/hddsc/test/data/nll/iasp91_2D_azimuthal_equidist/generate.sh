mkdir -p model
Vel2Grid run/phaseP.run
Vel2Grid run/phaseS.run
mkdir -p time
Grid2Time run/phaseP.run
Grid2Time run/phaseS.run

