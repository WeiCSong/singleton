import argparse
import pandas as pd
from bgen_reader import open_bgen
from tqdm import tqdm

def bgen_to_sparse_matrix(input_path, output_path):

    bgen = open_bgen(input_path, verbose=False)


    variant_df = bgen.read_variants()
    total_variants = len(variant_df)


    records = []


    for variant, genotype in tqdm(zip(variant_df.itertuples(), bgen.iter_genotype()),
                                   total=total_variants,
                                   desc="Processing variants"):

        if genotype.phased is False and genotype.probability is None:
            hardcalls = genotype.values  # shape: (n_samples,)

            het_indices = (hardcalls == 1).nonzero()[0]
            if len(het_indices) != 1:
                raise ValueError(f"Variant {variant.id} is not singleton (has {len(het_indices)} hets).")

            sample_index = het_indices[0]
            records.append((variant.id, sample_index))


    df = pd.DataFrame(records, columns=["variant_id", "sample_index"])
    df.to_csv(output_path, sep="\t", index=False)
    print(f"\n? Saved sparse matrix to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Convert singleton-only BGEN to sparse matrix format")
    parser.add_argument("--input", required=True, help="Input .bgen file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    bgen_to_sparse_matrix(args.input, args.output)

if __name__ == "__main__":
    main()
