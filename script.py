import pandas as pd
import numpy as np
import sys
import os
import matplotlib.pyplot as plt



def parse():
    with open('Text.txt') as f:
        contents = f.readlines()
    freq = {}
    for i, line in enumerate(contents):
        line = line.rstrip().replace(',', '.')
        if i == 0:
            key = line
            freq[key] = []
        else:
            freq[key].append(float(line))
    freq_df = pd.DataFrame().from_dict(freq)
    # print(freq_df.head())

    with open('Text2.txt') as f:
        contents = f.readlines()

    hla = {}
    for i, line in enumerate(contents):
        line = line.rstrip()
        if i == 0:
            key = line
            hla[key] = []
        else:
            hla[key].append(line)
    hla_df = pd.DataFrame().from_dict(hla)
    # print(hla_df.head())
    db = pd.concat([hla_df, freq_df], axis=1)
    db.to_csv('database.csv', index=False)


def main():
    args = sys.argv[1:]
    if not args or len(args) > 3:
        raise AttributeError
    elif len(args) == 1:
        input_file = args[0]
        rank_threshold, score_threshold = 1, 0.9
    elif len(args) == 3:
        input_file, rank_threshold, score_threshold = args[0], float(
            args[1]), float(args[2])

    file_name, file_type = input_file.split('_')[0], input_file.split('_')[1]
    if file_type == 'mhcI.csv':
        file_diag = 'mhcI'
        df = pd.read_csv(input_file)
        output = pd.DataFrame([
            df['allele'], df['core'],
            df['peptide'], df['score'], df['rank']]
        ).T
    elif file_type == 'mhcII.csv':
        file_diag = 'mhcII'
        df = pd.read_csv(input_file)
        output = pd.DataFrame([
            df['allele'], df['core_peptide'],
            df['peptide'], df['score'], df['rank']]
        ).T
    else:
        raise TypeError('Wrong input format')

    output.allele = output.allele.str.split('-', expand=True)[1]
    if file_type == 'mhcII.csv':
        split_df = output.allele.str.split('/')
        output.allele = split_df.apply(lambda x: x[-1] if len(x) > 0 else None)
        print('II')
        


    output_path = 'output'
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    # parse()
    db = pd.read_csv('database.csv')
    
    output = output[
        (output['rank'].astype(float) <= rank_threshold) &
        (output['score'].astype(float) >= score_threshold)
    ]
    output.reset_index(drop=True, inplace=True)
    print(output)

    freq = {'Frequency': []}
    for allele in output.allele:
        fr = np.nan
        search = db[db.HLA.str[:len(allele)] == allele]

        if search.size != 0:
            fr = search.Frequency.sum()
        else:
            fr = 0
        freq['Frequency'].append(fr)
    freq_df = pd.DataFrame().from_dict(freq)
    output = pd.concat([output, freq_df], axis=1)
    print(output.head())
    output.to_csv(f'{output_path}/Freq_{file_name}_{file_type}', index=False)
    count_freq = ((output['Frequency'] >= 0.1)&(output['rank'] > 0)).sum()
    count = (output['rank'] > 0).sum()

    def diagram():
        common = count_freq
        rare = count - common
        labels = [common, rare]
        sizes = [common, rare]
        plt.pie(sizes, labels = labels, textprops={'fontsize': 25})
        plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        # Display the circle diagram
        diagram_path = 'diagrams/' + file_name + '_' + file_diag + '.png'
        plt.savefig(diagram_path)



    if count != 0:
        diagram()
    else:
        print('There is not any optimal alleles with your parametrs')


if __name__ == '__main__':
    main()
