import numpy as np
import subprocess

def run_fixfield(output_fn, b0f, b0d, m, xs, tilt_deg):
    check = False
    for line in open('proton3mev.beam', 'r'):
        if "//number of keywords to be read" in line:
            assert int(line[0]) == 3, "computing closed orbit at wrong place, check proton3mev.beam for keywords, must be 3!"
            check = True
    assert check, "invalid proton3mev.beam, missing number of keywords line"
    print ('*** Running fixfield.exe...')
    bash_command = f'./fixfield.exe {b0f:.15f} {b0d:.15f} {m:.15f} {xs:.15f} {tilt_deg:.15f}'
    with open(output_fn, 'w') as outfile:
        process = subprocess.Popen(bash_command.split(), stdout=outfile)
        output, error = process.communicate()
    assert not error, error
    print ('*** Finished running fixfield.exe, output file written!')

def run_fixfield_da(output_fn, b0f, b0d, m, xs, tilt_deg, x_clo, z_clo, r, u, v):
    check = False
    for line in open('proton3mev.beam', 'r'):
        if "//number of keywords to be read" in line:
            if int(line[0]) != 1:
                print ("unnecessarily checking for closed orbit, check proton3mev.beam for keywords, should be 1!")
            check = True
    if not check:
        print("Did not find 'number of keywords' line in proton3mev.beam!")
    print ('*** Running fixfield-da.exe...')
    dec_matrix = " ".join(map(str, r.real.ravel()))
    bash_command = f'./fixfield-da.exe {b0f:.15f} {b0d:.15f} {m:.15f} {xs:.15f} {tilt_deg:.15f}'
    bash_command += " " + dec_matrix + f" {x_clo:.15f} 0 {z_clo:.15f} 0 {u:.15f} {v:.15f}"
    with open(output_fn, 'w') as outfile:
        process = subprocess.Popen(bash_command.split(), stdout=outfile)
        output, error = process.communicate()
    assert not error, error
    print ('*** Finished running fixfield.exe, output file written!')


def co_transfermatrix(output_fn):
    """Run fixfield.exe and return the computed closed orbit
    and linear (coupled) transfer matrix.

    Args:
       - string output_fn: file name of output file from run_fixfield

    Returns:
       - closed orbit tuple: (x_clo, xprime_clo, z_clo, zprime_clo)
       - transfer matrix np.array
    """

    nl = False

    output = open(output_fn, "r").read()
    for line in str(output).split('\n'):
        if 'errorstop procedure...' in line or '!!! ERROR in moveback_totheboun:' in line or 'Numerical Recipes run-time error...' in line:
            x_clo = np.nan
            xprime_clo = np.nan
            z_clo = np.nan
            zprime_clo = np.nan
            transfer_matrix = np.nan
            print ('*** INFO: no closed orbit found!')
            break
        if 'Closed orbit found:' in line:
            x_clo = float(line.split('x_clo = ')[1].split(' [m]')[0])
            xprime_clo = float(line.split('xprime_clo = ')[1].split(' [deg]')[0])
            nl = True
            continue
        if nl:
            z_clo = float(line.split('z_clo = ')[1].split(' [m]')[0])
            zprime_clo = float(line.split('zprime_clo = ')[1].split(' [deg]')[0])
            #nl = False
            break

    if nl:
        transfer_matrix = str(output).split(
            "--------------------------------- Get matrix first order -------------------------------------\n")[1].split(
            "\nbias")[0]

        transfer_matrix = np.array([(l.split('\t')[:-1]) for l in transfer_matrix.split('\n')], dtype=float)

    return ((x_clo, xprime_clo, z_clo, zprime_clo), transfer_matrix)


def extract_da(output_fn):
    """Run fixfield.exe and return the computed closed orbit
    and linear (coupled) transfer matrix.

    Args:
       - string output_fn: file name of output file from run_fixfield_da

    Returns:
       - 2-tuple of decoupled amplitudes: (u, v)
       - boolean: survival of particle at (u, v)
    """

    nl = False

    output = open(output_fn, "r").read()
    for line in str(output).split('\n'):
        if 'errorstop procedure...' in line or '!!! ERROR in moveback_totheboun:' in line or 'Numerical Recipes run-time error...' in line:
            return (np.nan, np.nan), np.nan
            print ('*** ERROR: no closed orbit found!')
            break
        if 'acceptance in u,v:' in line:
            rem = line.split('au = ')[1]
            u = float(rem.split(',')[0])
            rem = line.split('av = ')[1]
            v = float(rem.split(',')[0])
            nl = True
            continue
        if nl and 'DA' in line:
            survival = 'DA OK' in line
            if not survival:
                assert 'DA NO' in line
            break

    assert nl, f'Something went wrong, didn\'t find expected lines in {output_fn}!'

    return (u, v), survival
