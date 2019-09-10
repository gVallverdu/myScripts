
import os
import shutil
import stat
from pymatgen.io.vasp.inputs import Kpoints

# points k
kx, ky, kz = 4, 4, 4

# job name :
#    * doit finir par .job
#    ou
#    * doit commencer par j
job = "test.job"
basename = "test"

lancer = "#!/bin/bash\n\n"

for U in [0, 2, 4, 6, 8, 10]:
    print("U = %f" % U)

    # make job directory
    dirname = "U_%d" % U
    dirname = os.path.join(os.getcwd(), dirname)
    os.mkdir(dirname)

    # copy files
    shutil.copy("POSCAR", os.path.join(dirname, "POSCAR"))
    shutil.copy("POTCAR", os.path.join(dirname, "POTCAR"))

    kpoints = Kpoints.gamma_automatic((kx, ky, kz))
    kpoints.write_file(os.path.join(dirname, "KPOINTS"))

    jobname = basename + ("_U%d" % U) + ".job"
    shutil.copy(job, os.path.join(dirname, jobname))

    # complete and copy INCAR file
    with open("INCAR", "r") as fincar:
        incar_lines = ""
        for line in fincar:
            if "LDAUU" in line:
                line = "  LDAUU = %3.1f 0.0\n" % U
            incar_lines += line
    with open(os.path.join(dirname, "INCAR"), "w") as fincar:
        fincar.write(incar_lines)

    lancer += "cd %s\n" % dirname
    lancer += "sbatch %s\n\n" % jobname

with open("lancer.sh", "w") as f:
    f.write(lancer)
os.chmod("lancer.sh", 0o755)
