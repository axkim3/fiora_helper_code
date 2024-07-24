import csv

def parse_mgf_file(mgf_filename, csv_filename):
    # Open the input .mgf file and output .csv file
    with open(mgf_filename, 'r') as mgf_file, open(csv_filename, 'w', newline='') as csv_file:
        # Create a CSV writer object
        csv_writer = csv.writer(csv_file)
        # Write the header row to the CSV file
        csv_writer.writerow(['Name', 'SMILES', 'Precursor_type', 'CE', 'Instrument_type'])

        # Initialize variables to hold the extracted data
        name = smiles = precursor_type = collision_energy = instrument_type = None

        # Iterate through each line in the .mgf file
        for line in mgf_file:
            line = line.strip()
            if line.startswith('NAME='):
                name = line.split('=', 1)[1]
            elif line.startswith('SMILES='):
                smiles = line.split('=', 1)[1]
            elif line.startswith('ADDUCT='):
                precursor_type = line.split('=', 1)[1]
            elif line.startswith('Collision energy='):
                collision_energy = line.split('=', 1)[1]
            elif line.startswith('INSTRUMENT_TYPE='):
                instrument_type = line.split('=', 1)[1]
            elif line.startswith('END IONS'):
                # Write the extracted data to the CSV file when 'END IONS' is encountered
                if name or smiles or precursor_type or collision_energy or instrument_type:
                    csv_writer.writerow([
                        name if name else '',
                        smiles if smiles else '',
                        precursor_type if precursor_type else '',
                        collision_energy if collision_energy else '',
                        instrument_type if instrument_type else ''
                    ])
                # Reset the variables for the next spectrum
                name = smiles = precursor_type = collision_energy = instrument_type = None

if __name__ == "__main__":
    mgf_filename = '20231031_nihnp_library_pos_all_lib_MS2.mgf'  # Replace with your .mgf file path
    csv_filename = 'input2.csv'  # Replace with your desired .csv file path
    parse_mgf_file(mgf_filename, csv_filename)
    print(f"Data successfully extracted to {csv_filename}")
