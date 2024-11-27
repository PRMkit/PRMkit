import re

aa = {
    "A": 71.03711,
    "R": 156.10111,
    "N": 114.04293,
    "D": 115.02694,
    "C": 103.00919,
    "E": 129.04259,
    "Q": 128.05858,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "L": 113.08406,
    "K": 128.09496,
    "M": 131.04049,
    "F": 147.06841,
    "P": 97.05276,
    "S": 87.03203,
    "T": 101.04768,
    "W": 186.07931,
    "Y": 163.06333,
    "V": 99.06841,
    #'X' : 0.00000     #  for unknown amino acids
}  # dictionary of amino acids 'aa' and their monoisotopic mass

OH = 17.00274
H_ = 1.007825


# p=re.compile(r'[{}][\[\+\d\.\]]*'.format(''.join([*aa])))
p = re.compile("([A-Z])(?:[^A-Z\d\+])*(\+\d+\.\d+)?")


def ion_mass(ids, q1_charge):
    ionmass_pair = []
    for name0 in ids:
        seq = p.findall(name0)
        seqmass = [aa[x] + (float(y) if y else 0) for x, y in seq]
        q1, charge = q1_charge
        # check if q1 match
        if q1 is not None:
            ionmass = H_ + sum(seqmass) + OH
            ionmass = ionmass / charge + H_
            # ionmass_pair.append((ionmass,name0,'y'+str(len(seqmass))+'+'*charge))
            if abs(ionmass - q1) > 0.005:
                print(
                    "{}\t{}\t{:.2f}\t{:.2f}\t{:.2f}".format(
                        name0, charge, ionmass, q1, ionmass - q1
                    )
                )
        for i in range(1, len(seqmass) + 1):
            ionmass = H_ + sum(seqmass[:i])
            count_pos = sum(1 for x, _ in seq[:i] if x in "HKR") + 1
            # count_pos=9
            for c in range(1, min(charge, count_pos) + 1):
                # for c in range(1,1+1):
                ionmass_pair.append(
                    ((ionmass + (c - 1) * H_) / c, name0, "b" + str(i) + "+" * c)
                )
            ionmass = OH + sum(seqmass[-i:]) + H_
            count_pos = sum(1 for x, _ in seq[-i:] if x in "HKR") + 1
            # count_pos=9
            for c in range(1, min(charge, count_pos) + 1):
                # for c in range(1,1+1):
                ionmass_pair.append((ionmass / c + H_, name0, "y" + str(i) + "+" * c))
    print(len(ionmass_pair))
    return sorted(ionmass_pair)


if __name__ == "__main__":
    import sys

    name0 = "RRPES[+79.966331]APAESSPSK"
    seq = p.findall(name0)
    print(seq)
    ions = ion_mass([name0], (789.867226, 2))
    for i in ions:
        print(i)
    # count_pos=sum(1 for x in 'DHFGLEGDEESTMLEDSVSPKK' if x in 'HKR')+1
    # print(count_pos)
    # for mm,_,ion in ion_mass(['SHSAGK'],(None,2)):
    #    print(mm,ion)
