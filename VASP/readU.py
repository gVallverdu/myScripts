
import os
from pymatgen.io.vasp.outputs import Vasprun, Outcar

class UData(object):

    def __init__(self, U=None):
        self.U = U

        self.a = None
        self.b = None
        self.c = None

        self.alpha = None
        self.beta = None
        self.gamma = None

        self.charges = None
        self.ave_charges = None
        self.magnetization = None
        self.ave_magnetization = None

        self.gap = None

        self.energy = None

        self.dirname = None

    def as_dict(self):

        d = {}
        d["a"] = self.a
        d["b"] = self.b
        d["c"] = self.c
        d["alpha"] = self.alpha
        d["beta"] = self.beta
        d["gamma"] = self.gamma
        d["charges"] = [q for at, q in self.charges]
        d["atoms"] = [at for at, q in self.charges]
        d["magnetization"] = [m for at, m in self.magnetization]
        d["average"] = {"magnetization": self.ave_magnetization,
                        "charges": self.ave_charges}
        d["gap"] = self.gap
        d["energy"] = self.energy

        d["dirname"] = self.dirname

        return d

    def write_line(self):
        line = ""
        def write_val(val, fmt):
            if not val:
                val = 0.0
            return fmt % val
        line += write_val(self.U, "%5.1f")
        line += write_val(self.a, "%10.4f")
        line += write_val(self.b, "%10.4f")
        line += write_val(self.c, "%10.4f")
        line += write_val(self.alpha, "%10.4f")
        line += write_val(self.beta, "%10.4f")
        line += write_val(self.gamma, "%10.4f")
        line += write_val(self.energy, "%16.8e")

        for q in self.ave_charges.values():
            line += write_val(q, "%10.4f")

        for m in self.ave_magnetization.values():
            line += write_val(m, "%10.4f")

        line += write_val(self.gap, "%8.3f")

        return line + "\n"

    def head(self):
        head = "#   U"
        head += "a".rjust(10)
        head += "b".rjust(10)
        head += "c".rjust(10)
        head += "alpha".rjust(10)
        head += "beta".rjust(10)
        head += "gamma".rjust(10)
        head += "energy (eV)".rjust(16)
        for el in self.ave_charges:
            head += ("q_%s" % el).rjust(10)
        for el in self.ave_magnetization:
            head += ("m_%s" % el).rjust(10)

        head += "gap (eV)".rjust(10)

        return head + "\n"

    @staticmethod
    def from_calc(U):
        data = UData(U=U)

        dirname = os.path.join(os.getcwd(), "U_%d" % U)
        data.dirname = dirname

        # read lattice
        optxml = os.path.join(dirname, "vasprun.xml")
        optrun = Vasprun(optxml, parse_dos=False, parse_eigen=False, parse_potcar_file=False)

        data.a, data.b, data.c = optrun.final_structure.lattice.abc
        data.alpha, data.beta, data.gamma = optrun.final_structure.lattice.angles

        # read charges
        fcharges = os.path.join(dirname, "Bader/charges.dat")
        ave_charges = {}
        charges = []
        with open(fcharges, "r") as f:
            charge_read = False
            for line in f:
                if not charge_read and "-------" in line:
                    line = f.readline()
                    while "------" not in line:
                        val = line.split()
                        charges.append((val[1], float(val[4])))
                        line = f.readline()
                    charge_read = True
                if "Averages :" in line:
                    f.readline()
                    f.readline()
                    val = f.readline().split()
                    ave_charges[val[0]] = float(val[2])
                    val = f.readline().split()
                    ave_charges[val[0]] = float(val[2])

        data.ave_charges = ave_charges
        data.charges = charges

        # read OUTCAR
        foutcar = os.path.join(dirname, "Bader/OUTCAR")
        outcar = Outcar(foutcar)

        # read energy from OUTCAR
        data.energy = outcar.final_energy

        # read magnetization from OUTCAR
        mag = [(site.specie.symbol, atmag["tot"]) for atmag, site in zip(outcar.magnetization, optrun.final_structure)]
        data.magnetization = mag
        avemag = {el.name: 0. for el in optrun.final_structure.composition.elements}
        for el, atmag in mag:
            avemag[el] += atmag
        for el, count in optrun.final_structure.composition.as_dict().items():
            avemag[el] /= float(count)

        data.ave_magnetization = avemag

        # gap
        dosxml = os.path.join(dirname, "DOS/vasprun.xml")
        dosrun = Vasprun(dosxml, parse_eigen=False, parse_potcar_file=False)
        data.gap = dosrun.complete_dos.get_gap()

        return data

if __name__ == "__main__":

    all_data = [UData.from_calc(U) for U in [0.]]
    # all_data = [UData.from_calc(U) for U in [0., 2., 4., 6., 8., 10.]]

    print(all_data[0].head())

    for data in all_data:
        print(data.write_line())


    # for U in [0., 2., 4., 6., 8., 10.]:
    #
    #     all_data.adata = UData.from_calc(U)
    #
    #     print(3 * "%10.4f" % (data.a, data.b, data.a))
    #     print(3 * "%10.4f" % (data.alpha, data.beta, data.gamma))
    #
    #     print("average charges")
    #     for el, q in data.ave_charges.items():
    #         print("%4s %10.4f" % (el, q))
    #
    #
    #     print("average mu_B:")
    #     for el, atmag in data.ave_magnetization.items():
    #         print("%4s %10.4f" % (el, atmag))
    #
    #     print("gap = %10.4f" % data.gap)
    #
    #     print("energy = %16.6e" % data.energy)
    #
