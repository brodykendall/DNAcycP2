import keras
import pandas as pd
import numpy as np
from numpy import array
from Bio import SeqIO

detrend_int = 0.001641373848542571
detrend_slope = 1.0158132314682007
# Mean and stdev of smoothed C0 for Tiling library:
# (calculated from/on the scale of normalized Cn values)
normal_mean = -0.011196041799376931
normal_std = 0.651684644408004

def dnaOneHot(sequence):
    seq_array = array(list(sequence))
    code = {"A": [0], "C": [1], "G": [2], "T": [3], "N": [4],
            "a": [0], "c": [1], "g": [2], "t": [3], "n": [4]}
    onehot_encoded_seq = []
    for char in seq_array:
        onehot_encoded = np.zeros(5)
        onehot_encoded[code[char]] = 1
        onehot_encoded_seq.append(onehot_encoded[0:4])
    return onehot_encoded_seq

def cycle_fasta(inputfile, folder_path):
    network_final = keras.models.load_model(folder_path)
    genome_file = SeqIO.parse(open(inputfile),'fasta')
    ret = {}
    for fasta in genome_file:
        chrom = fasta.id
        genome_sequence = str(fasta.seq)
        onehot_sequence = dnaOneHot(genome_sequence)
        onehot_sequence = array(onehot_sequence)
        onehot_sequence = onehot_sequence.reshape((onehot_sequence.shape[0],4,1))
        print(f"Sequence length for ID {chrom}: {str(onehot_sequence.shape[0])}")
        print("Predicting cyclizability...")
        fit = []
        fit_reverse = []
        for ind_local in np.array_split(range(25, onehot_sequence.shape[0]-24), 100):
            onehot_sequence_local = []
            for i in ind_local:
                s = onehot_sequence[(i-25):(i+25),]
                onehot_sequence_local.append(s)
            onehot_sequence_local = array(onehot_sequence_local)
            onehot_sequence_local = onehot_sequence_local.reshape((onehot_sequence_local.shape[0],50,4,1))
            onehot_sequence_local_reverse = np.flip(onehot_sequence_local,[1,2])
            fit_local = network_final.predict(onehot_sequence_local, verbose=0)
            fit_local = fit_local.reshape((fit_local.shape[0]))
            fit.append(fit_local)
            fit_local_reverse = network_final.predict(onehot_sequence_local_reverse, verbose=0)
            fit_local_reverse = fit_local_reverse.reshape((fit_local_reverse.shape[0]))
            fit_reverse.append(fit_local_reverse)
        fit = [item for sublist in fit for item in sublist]
        fit = array(fit)
        fit_reverse = [item for sublist in fit_reverse for item in sublist]
        fit_reverse = array(fit_reverse)
        fit = detrend_int + (fit + fit_reverse) * detrend_slope/2
        fit2 = fit * normal_std + normal_mean
        n = fit.shape[0]
        fitall = np.vstack((range(25,25+n),fit,fit2))
        fitall = pd.DataFrame([*zip(*fitall)])
        fitall.columns = ["position","c_score_norm","c_score_unnorm"]
        fitall = fitall.astype({"position": int})
        ret["cycle_"+chrom] = fitall
    
    return ret

def cycle(sequences, folder_path):
    network_final = keras.models.load_model(folder_path)
    X = []
    all50 = True
    print("Reading sequences...")
    for sequence_nt in sequences:
        if len(sequence_nt) != 50:
            all50=False
        X.append(dnaOneHot(sequence_nt))

    if all50:
        print("Predicting cyclizability...")
        X = array(X)
        X = X.reshape((X.shape[0],50,4,1))
        X_reverse = np.flip(X,[1,2])

        model_pred = network_final.predict(X)
        model_pred_reverse = network_final.predict(X_reverse)

        model_pred = detrend_int + (model_pred + model_pred_reverse) * detrend_slope/2
        output_cycle = model_pred.flatten()
        output_cycle2 = np.array([item * normal_std + normal_mean for item in output_cycle])
    else:
        print("Not all sequences are length 50, predicting every subsequence...")
        output_cycle = []
        lenX = len(X)
        for j, onehot_loop in enumerate(X):
            l = len(onehot_loop)
            onehot_loop = array(onehot_loop)
            onehot_loop = onehot_loop.reshape((l,4,1))
            onehot_loops = []
            for i in range(l-49):
                onehot_loops.append(onehot_loop[i:i+50])
            onehot_loops = array(onehot_loops)
            onehot_loops_reverse = np.flip(onehot_loops,[1,2])
            if l > 1000:
                # Provide status bar for long sequences:
                cycle_local = network_final.predict(onehot_loops)
                cycle_local_reverse = network_final.predict(onehot_loops_reverse)
            else:
                # No status bar for short sequences (verbose=0):
                cycle_local = network_final.predict(onehot_loops, verbose=0)
                cycle_local_reverse = network_final.predict(onehot_loops_reverse, verbose=0)
            cycle_local = detrend_int + (cycle_local + cycle_local_reverse) * detrend_slope/2
            cycle_local = cycle_local.reshape(cycle_local.shape[0])
            output_cycle.append(cycle_local)
            if j%10==9:
                print(f"Completed {j+1} out of {lenX} total sequences")
        output_cycle2 = [item * normal_std + normal_mean for item in output_cycle]

    ret = {"norm": output_cycle,
           "unnorm": output_cycle2}

    return ret
