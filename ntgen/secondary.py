import ViennaRNA


def compute_mfe(sequence, energy_range=10000):
    pred_structure, mfe = ViennaRNA.fold(sequence)
    energies = sorted([r.energy for r in ViennaRNA.subopt(sequence, energy_range)])
    # assert energies[0] == mfe
    delta = .01 * energy_range if len(energies) == 1 else energies[1] - energies[0]
    return pred_structure, mfe, delta
