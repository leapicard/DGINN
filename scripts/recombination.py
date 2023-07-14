import sys, os, json
from Bio import AlignIO


if __name__ == "__main__":


    # Load the configuration file
    with open(sys.argv[1], 'r') as config_in:
        config_dict = json.load(config_in)

    parameters = config_dict["parameters"]
    data = config_dict["data"]

    if parameters["recombination"]:
        alignment_file_path = data['aln']

        lines = AlignIO.read(open(alignment_file_path, 'r'), 'fasta')
        cutting_point = 2

        sequence_sets = [[] for _ in range(cutting_point)]

        for line in lines:
            sequence_length = len(line.seq)
            for i in range(cutting_point):
                n = min((i + 1) * (sequence_length // cutting_point), sequence_length)
                sequence_sets[i].append('>' + line.id + '\n')
                sequence_sets[i].append(format(line.seq[i * (sequence_length // cutting_point):n]) + '\n')

        base_directory = data.get('o', '')  
        recomb_files = []
        for i in range(cutting_point):
            file_path = os.path.abspath(os.path.join(base_directory, f'recomb_file_{i}.fasta'))  # Generate file paths dynamically
            recomb_files.append(file_path)
            with open(file_path, 'w') as out:
                out.writelines(sequence_sets[i])
        
        data['recomb_files'] = recomb_files

        # Update the configuration dictionary
        config_dict["data"] = data
        config_dict["parameters"] = parameters

        filename = data["o"]+"recomb_"+parameters["infile"].split("data/")[1].split(".")[0]
        
        recomb = open(filename+"_files_list.txt", "w")

        for file in data["recomb_files"]:
            if len(data["recomb_files"]): 
                recomb.write(file+"\n")

        with open(sys.argv[1], 'w') as config_out:
            json.dump(config_dict, config_out, indent=" ")
