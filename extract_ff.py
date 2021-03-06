import numpy as np
import subprocess

def run_fixfield():
    print ('Running fixfield.exe...')
    # output = subprocess.check_output(["./fixfield.exe"], shell=True)

    # bashCommand = './../fixfield-source/fixfield.exe'
    bashCommand = './../fixfield-source/fixfield.exe -1 1.15 1.31 0 0'
    outfile = open('output.txt', 'w')
    process = subprocess.Popen(bashCommand.split(), stdout=outfile)
    output, error = process.communicate()
    outfile.close()

def co_transfermatrix():
    """Run fixfield.exe and return the computed closed orbit
    and linear (coupled) transfer matrix.

    Returns:
       - closed orbit tuple: (x_clo, xprime_clo, z_clo, zprime_clo)
       - transfer matrix np.array
    """

    nl = False

    file = open("output.txt", "r")
    output = file.read()
    for line in str(output).split('\n'):
        if 'Closed orbit found:' in line:
            x_clo = float(line.split('x_clo = ')[1].split(' [m]')[0])
            xprime_clo = float(line.split('xprime_clo = ')[1].split(' [deg]')[0])
            nl = True
            continue
        if nl:
            z_clo = float(line.split('z_clo = ')[1].split(' [m]')[0])
            zprime_clo = float(line.split('zprime_clo = ')[1].split(' [deg]')[0])
            nl = False
            break

    transfer_matrix = str(output).split(
        "--------------------------------- Get matrix first order -------------------------------------\n")[1].split(
        "\nbias")[0]

    transfer_matrix = np.array([(l.split('\t')[:-1]) for l in transfer_matrix.split('\n')], dtype=float)

    return ((x_clo, xprime_clo, z_clo, zprime_clo), transfer_matrix)
