# Generate nll grids for tests. Due to the size of the grids, 
# not all of them can be added to git
cd iasp91_2D_global  ; ./generate.sh ; cd ..
cd iasp91_2D_sdc     ; ./generate.sh ; cd ..
cd iasp91_2D_simple  ; ./generate.sh ; cd ..
cd iasp91_2D_azimuthal_equidist; ./generate.sh ; cd ..
cd iasp91_2D_lambert           ; ./generate.sh ; cd ..
cd iasp91_2D_merc              ; ./generate.sh ; cd ..
cd iasp91_3D_sdc     ; ./generate.sh ; cd ..
cd iasp91_3D_simple  ; ./generate.sh ; cd ..
cd iasp91_3D_azimuthal_equidist; ./generate.sh ; cd ..
cd iasp91_3D_lambert           ; ./generate.sh ; cd ..
cd iasp91_3D_merc              ; ./generate.sh ; cd .. 
