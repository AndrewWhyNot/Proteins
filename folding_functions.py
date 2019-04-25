from Bio.PDB.Vector import calc_angle, calc_dihedral
import Bio.PDB.vectors as vect


def chi1Atoms(residue):
    if residue in ['ALA', 'GLY']:
        return []
    res = ['N', 'CA', 'CB']
    if residue.resname in ['ILE', 'VAL']:
        return res + ['CG1']
    if residue.resname in ['CYS']:
        return res + ['SG']
    if residue.resname in ['THR']:
        return res + ['OG1']
    if residue.resname in ['SER']:
        return res + ['OG']
    return res + ['CG']

def chi2Atoms(residue):
    if residue.resname not in ['ARG', 'ASN', 'ASP', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'TRP', 'TYR']:
        return []
    res = ['CA', 'CB']
    if residue.resname in ['ILE']:
        res.append('CG1')
    else:
        res.append('CG')
    if residue.resname in ['ARG', 'GLN', 'GLU', 'LYS', 'PRO']:
        return res + ['CD']
    if residue.resname in ['LEU', 'PHE', 'TRP', 'TYR', 'ILE']:
        return res + ['CD1']
    if residue.resname in ['ASN', 'ASP']:
        return res + ['OD1']
    if residue.resname in ['MET']:
        return res + ['SD']
    return res + ['ND1']

def chi3Atoms(residue):
    if residue.resname not in ['ARG', 'GLN', 'GLU', 'LYS', 'MET']:
        return []
    res = ['CB', 'CG']
    if residue.resname in ['MET']:
        return res + ['SD', 'CE']
    res.append('CD')
    if residue.resname in ['GLN', 'GLU']:
        return res + ['OE1']
    if residue.resname in ['ARG']:
        return res + ['NE']
    return res + ['CE']

def chi4Atoms(residue):
    if residue.resname not in ['ARG', 'LYS']:
        return []
    if residue.resname in ['ARG']:
        return ['CG', 'CD', 'NE', 'CZ']
    else:
        return ['CG', 'CD', 'CE', 'NZ']

def chi5Atoms(residue):
    if residue.resname not in ['ARG']:
        return []
    return ['CD', 'NE', 'CZ', 'NH1']

def phiAtoms(residue):
    return ['C', 'N', 'CA', 'C']

def psiAtoms(residue):
    return ['N', 'CA', 'C', 'N']

def getDihedralAngleAtoms(residues, angleName):
    if angleName == 'chi1':
        return chi1Atoms(residues[0])
    elif angleName == 'chi2':
        return chi2Atoms(residues[0])
    elif angleName == 'chi3':
        return chi3Atoms(residues[0])
    elif angleName == 'chi4':
        return chi4Atoms(residues[0])
    elif angleName == 'chi5':
        return chi5Atoms(residues[0])
    elif angleName == 'phi':
        return phiAtoms(residues)
    elif angleName == 'psi':
        return psiAtoms(residues)
    elif angleName == '':
        return []
    else:
        raise ValueError('No such dihedral angle!')

def getNumberOfChiAngles(resName):
    if resName in ['ALA', 'GLY']:
        return 0
    if resName in ['ARG']:
        return 5
    if resName in ['LYS']:
        return 4
    if resName in ['GLN', 'GLU', 'MET']:
        return 3
    if resName in ['ASN', 'ASP', 'HIS', 'ILE', 'LEU', 'PHE', 'PRO', 'TRP', 'TYR']:
        return 2
    if resName in ['CYS', 'SER', 'THR', 'VAL']:
        return 1
    raise Exception('No such residue name!')
    


def getDihedralAngle(residues, angleName):
    atoms = getDihedralAngleAtoms(residues, angleName)
    if len(atoms) == 0:
        return None
    try:
        if angleName not in ['phi', 'psi']:
            atoms = [residues[0][atomName].get_vector() for atomName in atoms]
        elif angleName == 'phi':
            atoms = [residues[0][atoms[0]].get_vector()] + [residues[1][atoms[i]].get_vector() for i in range(1,4)]
        else:
            atoms = [residues[0][atoms[i]].get_vector() for i in range(3)] + [residues[1][atoms[3]].get_vector()]
        
        angleVal = vect.calc_dihedral(atoms[0], atoms[1], atoms[2], atoms[3]) / 3.1415 * 180
        return vect.calc_dihedral(atoms[0], atoms[1], atoms[2], atoms[3]) / 3.1415 * 180
    except KeyError as e:
        #print('KeyError in', residues[0].resname, ':', e)
        return None

def plotAngleResults(angleValues, angleName, resName, prevNeighbour, nextNeighbour):
    if len(angleValues) == 0:
        return
    hist = plt.hist(angleValues, bins=72, range=(-180, 180))
    plt.title('Значения углов ' + angleName + ' в аминокислоте ' + resName + ' (' + str(len(angleValues)) + ' значений)\nЦепочка: ' + prevNeighbour + '-' + resName + '-' + nextNeighbour)
    plt.show()
    #print('Пики:', getPeakPoints(angleValues))

def aminoacids():
    return ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

def secStructs():
    return ['H', 'B', 'E', 'G', 'I', 'T', 'S', '-']

def getPolarity(resName):
    if resName in ['ALA', 'GLY', 'VAL', 'LEU', 'ILE', 'PRO']:
        return 'UNPOLAR'
    if resName in ['SER', 'THR', 'CYS', 'MET', 'ASN', 'GLN']:
        return 'POLAR_NONCHARGED'
    if resName in ['PHE', 'TYR', 'TRP']:
        return 'AROMATIC'
    if resName in ['ASP', 'GLU']:
        return '-CHARGED'
    if resName in ['LYS', 'ARG', 'HIS']:
        return '+CHARGED'
    raise Exception('No such amino acid!')

