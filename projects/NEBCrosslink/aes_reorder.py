from ase.io import read, write
from ase.geometry import analysis

init = read("initial.geom")
final = read("final.geom")

# match atoms by minimal distance
mapping = init.get_all_distances(final).argmin(axis=1)

reordered = final[mapping]
write("final_reordered.geom", reordered)
