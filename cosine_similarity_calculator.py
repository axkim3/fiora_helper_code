import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.spatial.distance import cosine
from collections import defaultdict
from matplotlib.ticker import FuncFormatter

def parse_pluskal_mgf_file(mgf_filename):
    spectra = []
    with open(mgf_filename, 'r') as mgf_file:
        current_spectrum = None
        for line in mgf_file:
            line = line.strip()
            if line.startswith('BEGIN IONS'):
                current_spectrum = {
                    'title': None,
                    'adduct': None,
                    'collision_energy': None,
                    'peaks': []
                }
            elif line.startswith('NAME='):
                current_spectrum['title'] = line.split('=', 1)[1]
            elif line.startswith('ADDUCT='):
                current_spectrum['adduct'] = line.split('=', 1)[1]
            elif line.startswith('Collision energy='):
                current_spectrum['collision_energy'] = line.split('=', 1)[1]
            elif line.startswith('END IONS'):
                spectra.append(current_spectrum)
                current_spectrum = None
            elif current_spectrum is not None and line:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        if float(parts[0]) and float(parts[1]):
                            mz, intensity = float(parts[0]), float(parts[1])
                            current_spectrum['peaks'].append((mz, intensity))
                    except ValueError:
                        pass
    return spectra
def parse_fiora_mgf_file(mgf_filename):
    spectra = []
    with open(mgf_filename, 'r') as mgf_file:
        current_spectrum = None
        for line in mgf_file:
            line = line.strip()
            if line.startswith('BEGIN IONS'):
                current_spectrum = {
                    'title': None,
                    'adduct': None,
                    'collision_energy': None,
                    'peaks': []
                }
            elif line.startswith('TITLE='):
                current_spectrum['title'] = line.split('=', 1)[1]
            elif line.startswith('PRECURSORTYPE='):
                current_spectrum['adduct'] = line.split('=', 1)[1]
            elif line.startswith('COLLISIONENERGY='):
                current_spectrum['collision_energy'] = line.split('=', 1)[1]
            elif line.startswith('END IONS'):
                spectra.append(current_spectrum)
                current_spectrum = None
            elif current_spectrum is not None and line:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        if float(parts[0]) and float(parts[1]):
                            mz, intensity = float(parts[0]), float(parts[1])
                            current_spectrum['peaks'].append((mz, intensity))
                    except ValueError:
                        pass
    return spectra

def cosine_similarity(spectrum1, spectrum2, mz_tolerance = 0.02):
    # Create dictionaries to store m/z and intensity values
    spec1_dict = dict(spectrum1)
    spec2_dict = dict(spectrum2)
    
    # Get all unique m/z values
    all_mz = sorted(set(list(spec1_dict.keys()) + list(spec2_dict.keys())))
    
    # Initialize intensity vectors
    intensity1 = []
    intensity2 = []
    
    for mz in all_mz:
        # Find matching peaks within tolerance
        intensity1.append(max([spec1_dict.get(m, 0) for m in all_mz if abs(m - mz) <= mz_tolerance], default=0))
        intensity2.append(max([spec2_dict.get(m, 0) for m in all_mz if abs(m - mz) <= mz_tolerance], default=0))
    
    # Convert to numpy arrays
    intensity1 = np.array(intensity1)
    intensity2 = np.array(intensity2)
    
    # Calculate cosine similarity
    if np.sum(intensity1) == 0 or np.sum(intensity2) == 0:
        return 0.0  # Return 0 if either spectrum is empty
    
    similarity = 1000 * (1 - cosine(intensity1, intensity2))
    
    return similarity

def compare_spectra(experimental_spectra, fiora_spectra):
    similarities = []
    for exp_spec in experimental_spectra:
        for fiora_spec in fiora_spectra:
            if (exp_spec['title'] == fiora_spec['title'] and 
                exp_spec['adduct'] == fiora_spec['adduct'] and 
                exp_spec['collision_energy'] == fiora_spec['collision_energy']):
                similarity = cosine_similarity(exp_spec['peaks'], fiora_spec['peaks'])
                b = True
                for i in similarities:
                    if i['title'] == exp_spec['title']:
                        if similarity < i['similarity']:
                            b = False
                        else:
                            similarities.remove(i)
                if b:
                    similarities.append({
                        'title': exp_spec['title'],
                        'adduct': exp_spec['adduct'],
                        'collision_energy': exp_spec['collision_energy'],
                        'similarity': similarity,
                        'exp_peaks': exp_spec['peaks'],
                        'fiora_peaks': fiora_spec['peaks']
                    })
    return similarities

def plot_spectrum(ax, mz, intensity, color, label, direction = 1):
    """Plot a single MS/MS spectrum with intensity direction."""
    intensity = np.array(intensity) * direction
    ax.stem(mz, intensity, linefmt=color, markerfmt=" ", basefmt=" ", label=label)

    for m, i in zip(mz, intensity):
        if abs(i)  > 10:
            ax.text(m, i, f'{m:.2f}', va='bottom' if i > 0 else 'top', ha='center', fontsize=8)


def abs_formatter(x, pos):
    return f'{abs(int(x))}'

def graph_spectra(spectra):
    title = spectra['title']
    exp_mz = [i[0] for i in spectra['exp_peaks']]
    exp_intensity = [i[1] for i in spectra['exp_peaks']]
    fiora_mz = [i[0] for i in spectra['fiora_peaks']]
    m = 0
    for i in spectra['fiora_peaks']:
        m = max(i[1], m)
    fiora_intensity = [100 * i[1] / m for i in spectra['fiora_peaks']]
    
    fig, ax = plt.subplots()

    # Plot experimental spectrum upwards
    plot_spectrum(ax, exp_mz, exp_intensity, 'b', 'Experimental', direction=1)
    
    # Plot theoretical spectrum downwards
    plot_spectrum(ax, fiora_mz, fiora_intensity, 'r', 'Fiora', direction=-1)

    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity (%)")
    ax.legend()
    ax.set_title(f"Fiora Results for {spectra['title']} \nCosine Similarity: {spectra['similarity']:.4f}")
    
    # Set zero line
    ax.axhline(0, color='black', linewidth=0.5)
    ax.yaxis.set_major_formatter(FuncFormatter(abs_formatter))
    plt.show()



if __name__ == "__main__":
    experimental_mgf = 'Pluskal/MS2/pos/20240411_mcebio_library_pos_all_lib_MS2.mgf'  # Replace with your experimental .mgf file path
    fiora_mgf = 'Pluskal/MS2/pos/Pluskal-mcebio-fioraoutput.mgf'  # Replace with your Fiora-generated .mgf file path
    experimental_spectra = parse_pluskal_mgf_file(experimental_mgf)
    fiora_spectra = parse_fiora_mgf_file(fiora_mgf)
    similarities = compare_spectra(experimental_spectra, fiora_spectra)

    for result in similarities:
        #if result['title'] == 'STL568828':
            #graph_spectra(result)
        print(f"TITLE: {result['title']}, Adduct: {result['adduct']}, Collision Energy: {result['collision_energy']}, Cosine Similarity: {result['similarity']:.4f}")
    cos_similarities = [i['similarity'] for i in similarities]
    plt.hist(cos_similarities, bins=30, color='skyblue', edgecolor='black')
    # Adding labels and title
    plt.xlabel('Cosine Similarity Values')
    plt.ylabel('Frequency')
    plt.title('Histogram of Cosine Similarities')
    plt.show(block=True)
