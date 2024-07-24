import numpy as np
from scipy.spatial.distance import cosine
from collections import defaultdict

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
    """
    Calculate the cosine similarity between two MS/MS spectra.
    
    :param spectrum1: List of (m/z, intensity) tuples for the first spectrum
    :param spectrum2: List of (m/z, intensity) tuples for the second spectrum
    :param mz_tolerance: Mass-to-charge ratio tolerance for matching peaks
    :return: Cosine similarity score between 0 and 1
    """
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
    
    similarity = 1 - cosine(intensity1, intensity2)
    
    return similarity

def compare_spectra(experimental_spectra, fiora_spectra):
    similarities = []
    for exp_spec in experimental_spectra:
        for fiora_spec in fiora_spectra:
            if (exp_spec['title'] == fiora_spec['title'] and 
                exp_spec['adduct'] == fiora_spec['adduct'] and 
                exp_spec['collision_energy'] == fiora_spec['collision_energy']):
                #print(f"TITLE: {exp_spec['title']}, Adduct: {exp_spec['adduct']}, Collision Energy: {exp_spec['collision_energy']}")
                #print(i for i in exp_spec['peaks'])
                #print(i for i in fiora_spec['peaks'])
                similarity = cosine_similarity(exp_spec['peaks'], fiora_spec['peaks'])
                similarities.append({
                    'title': exp_spec['title'],
                    'adduct': exp_spec['adduct'],
                    'collision_energy': exp_spec['collision_energy'],
                    'similarity': similarity
                })
    return similarities

if __name__ == "__main__":
    experimental_mgf = 'Pluskal/20231031_nihnp_library_pos_all_lib_MS2.mgf'  # Replace with your experimental .mgf file path
    fiora_mgf = 'Pluskal-60-pos-Orbitrap-100-OUTPUT.mgf'  # Replace with your Fiora-generated .mgf file path
    experimental_spectra = parse_pluskal_mgf_file(experimental_mgf)
    fiora_spectra = parse_fiora_mgf_file(fiora_mgf)
    similarities = compare_spectra(experimental_spectra, fiora_spectra)
    
    for result in similarities:
        print(f"TITLE: {result['title']}, Adduct: {result['adduct']}, Collision Energy: {result['collision_energy']}, Cosine Similarity: {result['similarity']:.4f}")
